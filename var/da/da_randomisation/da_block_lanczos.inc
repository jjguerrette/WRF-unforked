subroutine da_block_lanczos(grid, mz, config_flags,        &
                           cv_size, xbx, be, iv, cv,   &
                           re, y, xhat, nsamp, j_cost,  &
                           qhat, eignval, eignvec, neign, hessian_info)

   !-------------------------------------------------------------------------
   ! Purpose:         Main Block Lanczos routine, separated into 2 stages, and 2 type of inner iteration
   !
   ! Here 
   !    cv   is updated in outer-loop.
   !    xhat is the control variable in inner-loop.
   !
   ! Called from da_solve
   !
   ! History: 09/01/2016  Creation (JJ Guerrette)
   !
   ! References: Halko et al. (2011); Bousserez and Henze (2016)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   integer, intent(in)               :: mz(:)
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: cv_size                  ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                      ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                       ! background error structure.
   type (iv_type), intent(inout)     :: iv                       ! ob. increment vector.
   real,   intent(inout)             :: cv(1:cv_size)            ! control variable (local).
   type (y_type), intent(inout)      :: re                       ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                        ! y = H(x_inc) structure.
   real,   intent(out)               :: xhat(1:cv_size)          ! control variable (local).
   integer, intent(in)               :: nsamp                     ! Block size

   type (j_type), optional, intent(out) :: j_cost                ! cost function
   real, optional, intent(out)      :: qhat(cv_size,0:nmodes_global)  ! Hessian basis vectors
   real, optional, intent(out)      :: eignvec(nmodes_global, nmodes_global)
   real, optional, intent(out)      :: eignval(nmodes_global)
   integer, optional, intent(out)   :: neign                     ! Number of non-truncated eigenpairs
   type(hessian_type), optional, intent(inout) :: hessian_info

#if defined(LAPACK)
   real, allocatable :: qhat_(:,:)  ! Hessian basis vectors
   real                         :: LAM(1:nmodes_global)
   real, allocatable            :: yhat(:,:)  ! yhat = A * omega
   real, allocatable            :: omega(:)   ! Current inner iteration hessian basis vectors
   real, allocatable            :: qB(:,:), qM(:,:)
   real, allocatable            :: MMAT(:,:,:), BMAT(:,:,:), TMAT(:,:), RMAT(:,:)

 
   integer                      :: iens, jens, i, j, nmodes, info
   integer                      :: iiens, jjens, block_start, block_end, iter

   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real                         :: ghat(1:cv_size)        ! cv copy.

   real*8                       :: temp_l
   real*8                       :: temp_e(nmodes_global)
   real                         :: beta
   real                         :: j0_total
   real                         :: hess_trace_approx
   character*10                 :: cproc, cens, cinner
   character(len=filename_len)  :: filename
   integer                      :: mat_unit
   character(len=5) :: cnsamp=""
   logical          :: swapped

   if (trace_use) call da_trace_entry("da_block_lanczos")

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
   je_start   = be % cv % size_jb + 1
   je_end     = be % cv % size_jb + be % cv % size_je
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
   jl_start   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + 1
   jl_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + be % cv % size_jl

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

   xhat = 0.D0
   write(cproc,fmt='(i4.4)') myproc
   write(cens,fmt='(i4.4)') ensmember
   write(cinner,fmt='(i4.4)') rand_inner_it

   if ( rand_stage .eq. 1 ) then
   !-------------------------------------------------------------------------
   ! [2.0] Stage 1: Calculate the ensemember gradient vector realization
   !-------------------------------------------------------------------------
      write(unit=stdout,fmt='(A)') 'Stage 1 of Block Lanczos'
      write(unit=stdout,fmt=*) ' '

      if ( ensmember.gt.0 .and. ensmember.le.nsamp) then
         ! Calculate J and grad J (ghat)
         call da_calculate_j(rand_outer, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                             be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)

         j0_total = j_cost%total
         if (j0_total == 0.0) return

         beta = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))

         write(unit=stdout,fmt='("Starting outer iteration : ",i3)') rand_outer
         write(unit=stdout,fmt=11) j0_total, beta
11       format('Starting cost function: ' ,1PD15.8,', Gradient: ',1PD15.8)

         ghat = - ghat

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
         ! Write ghat realization to file
         filename = trim(randvecdir)//'ghat.e'//trim(adjustl(cens))// &
                        '.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, ghat, grid, mz )
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_ghat')
#endif
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nsamp, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if

   else if ( rand_stage .eq. 2 ) then
   !-------------------------------------------------------------------------
   ! [3.0] Stage 2: Orthonormalize all gradient vector realizations,
   !    store in qhat vector files for Block Lanczos Initialization
   !-------------------------------------------------------------------------
      write(unit=stdout,fmt='(A)') 'Stage 2 of Block Lanczos'
      write(unit=stdout,fmt=*) ' '

      allocate(qhat_(1:cv_size,0:nsamp))
      qhat_ = 0.D0

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      do iens = 1, nsamp
         write(cens,fmt='(i4.4)') iens
         write(unit=stdout,fmt='("Reading ghat vectors; iens,iter = ",i3)') iens
         filename = trim(randvecdir)//'ghat.e'//trim(adjustl(cens))// &
                        '.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, qhat_(1:cv_size,iens), grid, mz )
      end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_ghat')
#endif

!      if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end) = 0.D0
!      if ( be%cv%size_je.gt.0) omega(je_start:je_end) = 0.D0
!      if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end) = 0.D0

!      if (rotate_omega .and. spectral_precon .and. rand_outer.gt.1) &
!         call da_rotate_omega(grid, cv_size, mz, jp_start, jp_end, &
!                                 rand_outer-1, omega, omega)

      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                      qhat_,cv_size,nsamp)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      do iens = 1, nsamp
         write(cens,fmt='(i4.4)') iens
         write(unit=stdout,fmt='("Writing qhat ensemble member ",i3," of ",i3)') iens,nsamp
         filename = trim(randvecdir)//'qhat.e'//trim(adjustl(cens))// &
                        '.p'//trim(adjustl(cproc))// &
                        '.iter'//trim(adjustl(cinner))
         call da_cv_io ('w', filename, cv_size, qhat_(1:cv_size,iens), grid, mz )
      end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_qhat')
#endif

      deallocate(qhat_)

   else if ( rand_stage .eq. 3 ) then
   !-------------------------------------------------------------------------
   ! [4.0] Stage 3: Multiply Hessian by current inner iteration ensemember 
   !                qhat vector
   !-------------------------------------------------------------------------
      write(unit=stdout,fmt='(A)') 'Stage 3 of Block Lanczos'
      write(unit=stdout,fmt=*) ' '

      if ( ensmember.gt.0 .and. ensmember.le.nsamp) then
         allocate(omega(1:cv_size))
         write(unit=stdout,fmt='("Reading qhat ensemble member ",i3," of ",i3)') ensmember,nsamp
         filename = trim(randvecdir)//'qhat.e'//trim(adjustl(cens))// &
                        '.p'//trim(adjustl(cproc))// &
                        '.iter'//trim(adjustl(cinner))

         call da_cv_io ('r', filename, cv_size, omega(1:cv_size), grid, mz )

         allocate(yhat(1:cv_size,1))
         yhat = 0.D0

         call da_calculate_gradj(rand_outer, rand_inner_it, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size),y,yhat(1:cv_size,1),grid,config_flags)

         hess_trace_approx = &
            da_dot_cv(cv_size, omega(:), yhat(:,1), grid, mz, jp_start, jp_end)
         write(unit=stdout,fmt='(A,E27.16E3)') 'HESS_TRACE: ',hess_trace_approx

         deallocate(omega)

         ! Write yhat to file
         filename = trim(randvecdir)//'yhat.e'//trim(adjustl(cens))// &
                           '.p'//trim(adjustl(cproc))// &
                           '.iter'//trim(adjustl(cinner))

         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1), grid, mz )

         deallocate(yhat)
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nsamp, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if
   else if ( rand_stage .eq. 4 ) then
   !-------------------------------------------------------------------------
   ! [5.0] Stage 4: Gather qhat/yhat from previous/current inner iterations
   !                Carry out Block Lanczos algorithm for eigendecomposition
   !                  and increment 
   !-------------------------------------------------------------------------
      write(unit=stdout,fmt='(A)') 'Stage 4 of Block Lanczos'
      write(unit=stdout,fmt=*) ' '

      nmodes = nmodes_global

      allocate(qhat_(1:cv_size,0:nmodes+nsamp))
      allocate(MMAT(nsamp,nsamp,rand_inner_it))
      allocate(BMAT(nsamp,nsamp,rand_inner_it))
      allocate(TMAT(nmodes,nmodes))

      qhat_=0.D0

      TMAT = 0.D0
      MMAT = 0.D0

      do iter = 1, rand_inner_it
         block_start = nsamp*(iter-1) + 1
         block_end = block_start + nsamp - 1

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

         write(cinner,fmt='(i4.4)') iter
         write(unit=stdout,fmt='("Reading qhat vectors; iter = ",i3)') iter
         do iens = 1, nsamp
            write(cens,fmt='(i4.4)') iens
            filename = trim(randvecdir)//'qhat.e'//trim(adjustl(cens))// &
                           '.p'//trim(adjustl(cproc))// &
                           '.iter'//trim(adjustl(cinner))
            call da_cv_io ('r', filename, cv_size, qhat_(1:cv_size,block_start+iens-1), grid, mz )
         end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_qhat')
#endif

         ! Populate TMAT with MMAT and BMAT blocks
         if (iter .lt. rand_inner_it) then
            !Read MMAT(:,:,iter) from file
!            write(cinner,fmt='(i4.4)') iter
            filename = trim(randvecdir)//'MMAT.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_mat_io ('r', filename, nsamp, nsamp, MMAT(:,:,iter) )
         else
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
            !Calculate MMAT(:,:,rand_inner_it) from latest yhat and qhat vectors
!            write(cinner,fmt='(i4.4)') iter

            write(unit=stdout,fmt='("Reading yhat vectors; iter = ",i3)') iter
!            allocate(yhat(1:cv_size,block_start:block_end))
            allocate(yhat(1:cv_size,1:nsamp))

            do iens = 1, nsamp
               write(cens,fmt='(i4.4)') iens
               filename = trim(randvecdir)//'yhat.e'//trim(adjustl(cens))// &
                              '.p'//trim(adjustl(cproc))// &
                              '.iter'//trim(adjustl(cinner))
!               call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,block_start+iens-1), grid, mz )
               call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,iens), grid, mz )
            end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_yhat')
#endif

!#ifdef DM_PARALLEL
!            call mpi_barrier(comm, ierr)
!#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
            write(unit=stdout,fmt='("Generating current MMAT")')
            !MMAT = Q^T * Y
            do iens = 1,nsamp
               iiens = block_start+iens-1
               do jens = 1, nsamp
!                  jjens = block_start+jens-1
                  MMAT(iens,jens,iter) = da_dot_cv(cv_size, qhat_(:,iiens), &
                                                    yhat(:,jens), grid, mz, jp_start, jp_end)
!                                                    yhat(:,jjens), grid, mz, jp_start, jp_end)
               end do
            end do
            MMAT(:,:,iter) = (MMAT(:,:,iter) + TRANSPOSE(MMAT(:,:,iter))) / 2.D0

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_MMAT')
#endif

            !Write MMAT(:,:,iter) to file
            filename = trim(randvecdir)//'MMAT.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_mat_io ('w', filename, nsamp, nsamp, MMAT(:,:,iter) )

            if (rand_inner_it .eq. max_rand_inner) deallocate(yhat)
         end if

         TMAT(block_start:block_end,block_start:block_end) = MMAT(:,:,iter)


         if (iter .gt. 1) then
            !Read BMAT(:,:,iter) from file
            write(cinner,fmt='(i4.4)') iter-1
            filename = trim(randvecdir)//'BMAT.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_mat_io ('r', filename, nsamp, nsamp, BMAT(:,:,iter-1) )

            TMAT(block_start-nsamp:block_end-nsamp,block_start:block_end) = TRANSPOSE(BMAT(:,:,iter-1))
            TMAT(block_start:block_end,block_start-nsamp:block_end-nsamp) = BMAT(:,:,iter-1)
         end if

      end do

!      !Write TMAT to file (not necessary, only for debugging)
!      write(cinner,fmt='(i4.4)') rand_inner_it
!      filename = trim(randvecdir)//'TMAT.iter'//trim(adjustl(cinner))// &
!                     '.p'//trim(adjustl(cproc))
!      call da_mat_io ('w', filename, nmodes, nmodes, TMAT, .true.)

      ! Construct eigendecomposition of TMAT
      if (rand_inner_it.eq.max_rand_inner .or. var4d_inc_out) then

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

!BEGIN ENCAPUSLATE THIS EIGENDECOMP --> Same as second half of da_rsvd56
       ! Determine eigenvalues and eigenvectors of the symmetric TMAT matrix
       !-----------------------------------------------------------------
         write(unit=stdout,fmt='("Decomposing TMAT")')
         info = 0
         call SYEV(TMAT, LAM, JOBZ = 'V' , INFO = info) !F95
#ifdef DM_PARALLEL
         call mpi_barrier(comm, ierr)
#endif
         if (info /=0) then
            write(stdout,*) 'Error in Block Lanczos: SYEV returned ',info
            CALL wrf_error_fatal('Canceling Block Lanczos')
         end if

         do iens = 1, nmodes
            do jens = 1, nmodes
               eignvec(iens,jens) = TMAT(iens,jens)
            end do
         end do
         deallocate(TMAT)

!         temp_e = 1.D0
!         if (all(LAM .le. 0.D0)) then
!            write(unit=stdout,fmt='(A)') &
!               'WARNING: Preconditioned Hessian: all(eigs <= 0), setting them positive.'
!            LAM = - LAM
!         else if ( sum(sign(temp_e,LAM)) .lt. 0.D0 ) then
!            write(unit=stdout,fmt='(A,F10.1,A)') &
!               'WARNING: Preconditioned Hessian: More negative than positive eigenvalues,', sum(sign(temp_e,LAM)), ', setting LAM=-LAM'
!            LAM = - LAM
!         end if

         !Sort eigenvalues and eigenvectors in descending order
         do j = nmodes-1, 1, -1
           swapped = .false.
           do i = 1, j
             if (LAM(i) .lt. LAM(i+1)) then
               temp_l = LAM(i)
               LAM(i) = LAM(i+1)
               LAM(i+1) = temp_l

               temp_e = eignvec(:,i)
               eignvec(:,i) = eignvec(:,i+1)
               eignvec(:,i+1) = temp_e

               swapped = .true.
             end if
           end do
           if (.NOT. swapped) exit
         end do

         neign = nmodes

         !Note: neign should be less than nobs or rank of Hessian update
         ! (not a concern for large systems, but is for single/low # obs systems)

!        ! Optional: Adjust neign according to how many eigenvalues are positive (very small impact)
!         do i = 1, nmodes
!            if (i.gt.1 .and. LAM(i) .lt. 0.D0) then
!               neign = i-1
!               exit
!            end if
!         end do

!         !Preconditioned observation Hessian should be positive semi-definite
!         do i = 1, nmodes
!            LAM(i) = max(LAM(i),0.D0)
!         end do
         eignval = LAM
!END ENCAPUSLATE THIS EIGENDECOMP --> Same as second half of da_rsvd56

         ! Add 1.0 to get eigenvalues of full preconditioned Hessian
         eignval = eignval + 1.D0

         ! Print current eigenpairs
         ! ----------------------------------------------
         write(cnsamp,fmt='(I0)') nmodes
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'Eigenvalues(r): ',eignval(1:nmodes)

!         do i = 1, nmodes
!            write(unit=stdout,fmt='(A,'//trim(cnsamp)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(i,1:nmodes)
!         end do

         ! Read gradient
         filename = trim(randvecdir)//'ghat.e0001.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, ghat, grid, mz )

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_EIGENDECOMP')
#endif


         if ( var4d_inc_out &
              .and. present(hessian_info) &
              .and. allocated(hessian_info%ghat) ) then
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
            hessian_info%ghat = ghat

            hessian_info%nmodes_max = nmodes
            hessian_info%qhat = 0.D0
            hessian_info%qhat = qhat_

            allocate(hessian_info%eigs(1)%eignval(1:nmodes))
            allocate(hessian_info%eigs(1)%eignvec(1:nmodes,1:nmodes))

            !neign should be less than nobs or rank of Hessian update 
            ! (only a concern for small systems)
            hessian_info%eigs(1)%neign = neign
            hessian_info%eigs(1)%eignval = eignval
            hessian_info%eigs(1)%eignvec = eignvec
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: increment_diagnostics1')
#endif
         end if

         if (rand_inner_it.eq.max_rand_inner) then
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
            !! Calculate increment
            write(unit=stdout,fmt='(A)') ' '
            write(unit=stdout,fmt='(A,I5,A)') 'Using ', nmodes, ' eigenpairs'
            write(unit=stdout,fmt='(A)') ' '

!            ! Adaptive increment
!            if (minval(eignval(1:neign)) .lt. 2.D0) then
!               write(unit=stdout,fmt='(A)') 'Carrying out RIOT w/ Block Lanczos, LRU'
!               do i = 1, nmodes
!                  LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
!               end do
!            else
               write(unit=stdout,fmt='(A)') 'Calculating Block Lanczos Increment'
               do i = 1, nmodes
                  LAM(i) =  1.D0 / eignval(i)
               end do
!            end if

            call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                   nmodes,nmodes,neign, &
                                   LAM,eignvec,qhat_(:,0:nmodes), &
                                   ghat,xhat)

!            if (minval(eignval(1:neign)) .lt. 2.D0) xhat = xhat + ghat

            if (spectral_precon .and. rand_outer.gt.1) &
               call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                       rand_outer-1, xhat, xhat, 'TL' )

            if (spectral_precon) &
               call da_hessian_io ('w',grid, cv_size,rand_outer,nmodes,nmodes,&
                                  eignvec,eignval,qhat_(:,0:nmodes),mz,ghat)
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_INCREMENT')
#endif
         end if
      end if

      !Perhaps put a test here for convergence (similar to Lanczos/CG)
      ! Increment max_rand_inner if test fails
      !if (INCREMENT_FAILS_CONVERGENCE_TEST) max_rand_inner = max_rand_inner + 1

      !if (present(qhat)) qhat = qhat_(:,0:nmodes)
      qhat = qhat_(:,0:nmodes)

      if (rand_inner_it .lt. max_rand_inner) then
         write(unit=stdout,fmt= &
         '("Generate RMAT basis residual matrix")')
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
         !This section pepares for next inner iteration
         iter = rand_inner_it
         block_start = nsamp*(iter-1) + 1
         block_end = block_start + nsamp - 1

         ! Construct new residual matrix, RMAT
         allocate(RMAT(1:cv_size,0:nsamp))
         RMAT = 0.D0

         !qM = Q * MMAT
         allocate(qM(1:cv_size,nsamp))
         qM = 0.D0
         do iens = 1,nsamp
            iiens = block_start+iens-1
            do jens = 1, nsamp
               qM(:,jens) = qM(:,jens) + qhat_(:,iiens) * MMAT(iens,jens,iter)
            end do
         end do
         do iens = 1,nsamp
!            iiens = block_start+iens-1
!            RMAT(:,iens) = yhat(:,iiens) - qM(:,iens)
            RMAT(:,iens) = yhat(:,iens) - qM(:,iens)
         end do
         deallocate(qM)
         deallocate(MMAT)
         deallocate(yhat)

         !qB = Q * B^T
         if (iter .gt. 1) then
            allocate(qB(1:cv_size,nsamp))
            qB = 0.D0
            do iens = 1,nsamp
               iiens = block_start+iens-1 - nsamp
               do jens = 1, nsamp
                  qB(:,jens) = qB(:,jens) + qhat_(:,iiens) * BMAT(jens,iens,iter-1)
               end do
            end do
            do iens = 1,nsamp
               RMAT(:,iens) = RMAT(:,iens) - qB(:,iens)
            end do
            deallocate(qB)
         end if

         block_start = nmodes + 1
         block_end = block_start + nsamp - 1

         write(unit=stdout,fmt=&
         '("Generate qhat and BMAT for next inner iteraion")')
         !Decompose RMAT into QR = qhat_(:,block_start:block_end) * BMAT(:,:,iter)
         call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         Q=RMAT,cv_size=cv_size,nvec=nsamp, &
                         R=BMAT(:,:,iter))
         qhat_(:,block_start:block_end) = RMAT(:,1:nsamp)
         deallocate(RMAT)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_RMAT')
#endif

         !Write BMAT(:,:,iter) to file
         write(cinner,fmt='(i4.4)') iter
         filename = trim(randvecdir)//'BMAT.iter'//trim(adjustl(cinner))// &
                        '.p'//trim(adjustl(cproc))
         call da_mat_io ('w', filename, nsamp, nsamp, BMAT(:,:,iter) )

         deallocate(BMAT)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
         !Orthonormalize most recent qhat vectors w.r.t. previous vectors and write to disk
         call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         qhat_,cv_size,nmodes+nsamp,nmodes+1)
         write(cinner,fmt='(i4.4)') rand_inner_it+1
         do iens = 1, nsamp
            write(cens,fmt='(i4.4)') iens
            write(unit=stdout,fmt='("Writing qhat vector for next iteration; iens,iter = ",i3,",",i3)') iens, rand_inner_it+1
            filename = trim(randvecdir)//'qhat.e'//trim(adjustl(cens))// &
                           '.p'//trim(adjustl(cproc))// &
                           '.iter'//trim(adjustl(cinner))
            call da_cv_io ('w', filename, cv_size, qhat_(1:cv_size,block_start+iens-1), grid, mz )
         end do
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_cv_qhat')
#endif

      end if !end prep for next iteration

      deallocate(qhat_)

   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: For Block Lanczos, rand_stage should be 0, 1, 2, 3, or 4, not ',rand_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_block_lanczos")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_Lanczos')
#endif


#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for block Lanczos"/))
#endif

end subroutine da_block_lanczos


