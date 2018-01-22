subroutine da_block_lanczos(grid, mz, config_flags,        &
                           cv_size, xbx, be, iv, cv,   &
                           re, y, xhat, nens, j_cost,  &
                           qhat_out, eignval, eignvec, neign, hessian_info)

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
   integer, intent(in)               :: nens                     ! Block size

   type (j_type), optional, intent(out) :: j_cost                ! cost function
   real, optional, alocatable, intent(out) :: qhat(:,:)  ! Hessian basis vectors
   real,      optional, intent(out) :: eignvec(nens*rand_inner, nens*rand_inner)
   real,      optional, intent(out) :: eignval(nens*rand_inner)
   integer,   optional, intent(out) :: neign                         ! Number of non-truncated eigenpairs
   type(hessian_type), optional, intent(inout) :: hessian_info

#if defined(LAPACK)
   real                         :: LAM(1:nens*rand_inner)
   real, allocatable            :: yhat(:,:)  ! yhat = A * omega
   real, allocatable            :: omega(:)   ! Current inner iteration hessian basis vectors
   real, allocatable            :: qB(:,:), qM(:,:)

 
   integer                      :: iens, jens, i, j, nmodes, info
   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real                         :: ghat(1:cv_size)        ! cv copy.

   real*8                       :: temp_l
   real*8                       :: temp_e(nens*rand_inner)
   real*8                       :: temp_mat(nens,nens)
   real                         :: beta
   real                         :: j0_total
   character*10                 :: cproc, cens, cinner
   character(len=filename_len)  :: filename
   integer                      :: mat_unit
   character(len=3) :: cnens=""
   logical          :: swapped

   if (trace_use) call da_trace_entry("da_block_lanczos")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using randomized SVD method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
   je_start   = be % cv % size_jb + 1
   je_end     = be % cv % size_jb + be % cv % size_je
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
   jl_start   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + 1
   jl_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + be % cv % size_jl

   xhat = 0.D0
   write(cproc,fmt='(i4.4)') myproc

   if ( rand_stage .eq. 1) then
      ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = L^T * H^T * R^-1 * H * L 
      write(cens,fmt='(i4.4)') ensmember
      write(cinner,fmt='(i4.4)') rand_inner

      if ( ensmember.gt.0 .and. ensmember.le.nens) then
         allocate(omega(1:cv_size))

         if ( rand_inner .eq. 1 ) then 
!!!!!!!!!!!!!!
         ! The random vectors used in perturbed initial conditions must be orthogonal!
         ! Make sure this happens in da_main.inc
!Also remove irrellevant portions of omega for adjoint perturbation
!      if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,:) = 0.D0
!      if ( be%cv%size_je.gt.0) omega(je_start:je_end,:) = 0.D0
!      if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,:) = 0.D0


!!!!!!!!!!!!!!

            ! Calculate J and grad J (ghat)
!            if (ensmember .eq. 1) use_randomblock = .false.
            call da_calculate_j(rand_outer, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                                be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)
!            if (ensmember .eq. 1) use_randomblock = .true.

            j0_total = j_cost%total
            if (j0_total == 0.0) return

!            beta = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
!
!            write(unit=stdout,fmt='("Starting outer iteration : ",i3)') rand_outer
!            write(unit=stdout,fmt=11) j0_total, beta
!11          format('Starting cost function: ' ,1PD15.8,', Gradient: ',1PD15.8)

            ghat = - ghat

            ! Write ghat to file
            filename = 'qhat.e'//trim(adjustl(cens))// &
                        '.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_cv_io ('w', filename, cv_size, ghat, grid, mz)

!            if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end) = 0.D0
!            if ( be%cv%size_je.gt.0) omega(je_start:je_end) = 0.D0
!            if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end) = 0.D0

!            if (rotate_omega .and. spectral_precon .and. rand_outer.gt.1) &
!               call da_rotate_omega(grid, cv_size, mz, jp_start, jp_end, &
!                                       rand_outer-1, omega, omega)

            omega = ghat            
         else
            write(unit=stdout,fmt='("Reading qhat ensemble member ",i3," of ",i3)') ensmember,nens
            filename = 'qhat.e'//trim(adjustl(cens))// &
                           '.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size), grid, mz)
         end if

         allocate(yhat(1:cv_size,1))
         yhat = 0.D0

         call da_calculate_gradj(rand_outer, 1, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size),y,yhat(1:cv_size,1),grid,config_flags)

         deallocate(omega)

         ! Write yhat to file
         filename = 'yhat.e'//trim(adjustl(cens))// &
                           '.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1), grid, mz)

         deallocate(yhat)

      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nens, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if
   else if ( rand_stage .eq. 2 ) then
      nmodes = nens*rand_inner

      allocate(qhat(1:cv_size,0:nmodes+nens))
      allocate(MMAT(nens,nens,rand_inner))
      allocate(BMAT(nens,nens,rand_inner))
      allocate(TMAT(nmodes,nmodes))

      qhat=0.D0

      TMAT = 0.D0
      MMAT = 0.D0

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      ! Read gradient
      filename = 'yhat.e0001.iter0001.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat, grid, mz)

      do iter = 1, rand_inner
         vstart = nens*(iter-1) + 1
         vend = vstart + nens - 1

         write(cinner,fmt='(i4.4)') iter
         do iens = 1, nens
            write(cens,fmt='(i4.4)') iens
            write(unit=stdout,fmt='("Reading qhat vectors; iens,iter = ",i3,",",i3)') iens, iter
            filename = 'qhat.e'//trim(adjustl(cens))// &
                           '.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, qhat(1:cv_size,vstart+iens-1), grid, mz)
         end do

         if (iter .gt. 1) then
            !Read BMAT(:,:,iter) from file
            write(cinner,fmt='(i4.4)') iter-1
            filename = 'BMAT.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_mat_io ('r', filename, nens, nens, BMAT(:,:,iter-1))

            TMAT(vstart-nens:vend-nens,vstart:vend) = BMAT(:,:,iter-1)
            TMAT(vstart:vend,vstart-nens:vend-nens) = TRANSPOSE(BMAT(:,:,iter-1))
         end if

         if (iter .lt. rand_inner)
            !Read MMAT(:,:,iter) from file
            filename = 'MMAT.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_mat_io ('r', filename, nens, nens, MMAT(:,:,iter))
         else
            allocate(yhat(1:cv_size,vstart:vend))
            do iens = 1, nens
               write(cens,fmt='(i4.4)') iens
               write(unit=stdout,fmt='("Reading yhat vectors; iens,iter = ",i3,",",i3)') iens, iter
               filename = 'yhat.e'//trim(adjustl(cens))// &
                              '.iter'//trim(adjustl(cinner))// &
                              '.p'//trim(adjustl(cproc))
               call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,vstart+iens-1), grid, mz)
            end do

            !MMAT = Q^T * Y
            do iens = 1,nens
               iiens = vstart+iens-1
               do jens = 1, nens
                  jjens = vstart+jens-1
                  MMAT(iens,jens,iter) = da_dot_cv(cv_size, qhat(:,iiens), &
                                                    yhat(:,jjens), grid, mz, jp_start, jp_end)
               end do
            end do
            MMAT(:,:,iter) = (MMAT(:,:,iter) + TRANSPOSE(MMAT(:,:,iter))) / 2.D0

            deallocate(yhat)

         end if

         TMAT(vstart:vend,vstart:vend) = MMAT(:,:,iter)
      end do

      !Write MMAT(:,:,iter) to file
      filename = 'MMAT.iter'//trim(adjustl(cinner))// &
                     '.p'//trim(adjustl(cproc))
      call da_mat_io ('w', filename, nens, nens, MMAT(:,:,rand_inner))


      !Everything is in place to approximate Hessian and increment.  Do it now.

      ! Construct eigendecomposition of TMAT

!      if ( var4d_inc_out &
!           .and. present(hessian_info) &
!           .and. allocated(hessian_info%ghat) ) then
!
!         hessian_info%ghat = ghat
!
!         hessian_info%nmodes_max = nens
!         do iens = 0, nens
!            hessian_info%qhat(:,iens) = qhat(:,iens)
!         end do
!
!         do iens = nens, 1, -inc_out_interval
!            allocate(hessian_info%eigs(iens)%eignval(1:iens))
!            allocate(hessian_info%eigs(iens)%eignvec(1:iens,1:iens))
!
!!!            call da_rsvd56(grid, mz, jp_start, jp_end, &
!!!                          omega(:,1:iens), yhat(:,1:iens), qhat(:,0:iens), cv_size, &
!!!                          hessian_info%eigs(iens)%eignval, hessian_info%eigs(iens)%eignvec, &
!!!                          iens, hessian_info%eigs(iens)%neign)
!!!
!            hessian_info%eigs(iens)%eignval = hessian_info%eigs(iens)%eignval + 1.D0
!         end do
!      end if

!!!      call da_rsvd56(grid, mz, jp_start, jp_end, &
!!!                     omega, yhat, qhat, cv_size, &
!!!                     eignval, eignvec, nens, neign)


      if (rand_inner .eq. max_rand_inner) then

!ENCAPUSLATE THIS EIGENDECOMP --> Same as second half of rsvd56
    ! Determine eigenvalues and eigenvectors of the symmetric TMAT matrix
    !-----------------------------------------------------------------
      info = 0
      call SYEV(TMAT, LAM, JOBZ = 'V' , INFO = info) !F95

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

         temp_e = 1.D0
         if (all(LAM .le. 0.D0)) then
            if(.not.present(print_messages)) write(unit=stdout,fmt='(A)') &
               'WARNING: Preconditioned Hessian: all(eigs <= 0), setting them positive.'
            LAM = - LAM
         else if ( sum(sign(temp_e,LAM)) .lt. 0.D0 ) then
            if(.not.present(print_messages)) write(unit=stdout,fmt='(A,F10.1,A)') &
               'WARNING: Preconditioned Hessian: More negative than positive eigenvalues,', sum(sign(temp_e,LAM)), ', setting LAM=-LAM'
            LAM = - LAM
         end if

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

         !Preconditioned observation Hessian should be positive semi-definite
         do iens = 1, nmodes
            LAM(iens) = max(LAM(iens),0.D0)
         end do
         eignval = LAM

         ! Add 1.0 to get eigenvalues of full preconditioned Hessian
         eignval = eignval + 1.D0

         ! Print current eigenpairs
         ! ----------------------------------------------
         write(cnens,fmt='(I0)') nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'Eigenvalues(r): ',eignval(1:nens)

         do iens = 1, nens
            write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nens)
         end do

         write(unit=stdout,fmt='(A)') ' '
         write(unit=stdout,fmt='(A,I5,A)') 'Using ', nmodes, ' eigenpairs'
         write(unit=stdout,fmt='(A)') ' '


         !! Calculate increment
         do i = 1, nmodes
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do

         write(unit=stdout,fmt='(A)') 'Carrying out RIOT w/ Block Lanczos, LRU'

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                nmodes,nmodes,nmodes, &
                                LAM,eignvec,qhat, &
                                ghat,xhat)

         xhat = xhat + ghat

         if (spectral_precon .and. rand_outer.gt.1) &
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                    rand_outer-1, xhat, xhat, 'TL' )

         if (riot_precon.ge.1 .or. spectral_precon) &
            call da_hessian_io ('w',grid, cv_size,rand_outer,nmodes,nmodes,&
                               eignvec,eignval,qhat,mz,ghat)

         
      end if

      !Perhaps put a test here for convergence (similar to Lanczos/CG)
      ! Increment max_rand_inner if test fails

      !if (INCREMENT_FAILS_CONVERGENCE_TEST) max_rand_inner = max_rand_inner + 1

      if (rand_inner .lt. max_rand_inner) then
         !This sections pepares for next inner iteration

         iter = rand_inner
         vstart = nens*(iter-1) + 1
         vend = vstart + nens - 1

         ! Construct new residual matrix, RMAT
         allocate(RMAT(1:cv_size,nens))

         !qM = Q * MMAT
         allocate(qM(1:cv_size,nmodes)
         qM = 0.D0
         do iens = 1,nens
            iiens = vstart+iens-1
            do jens = 1, nens
               qM(:,jens) = qM(:,jens) + qhat(:,iiens) * MMAT(iens,jens,iter)
            end do
         end do
         do iens = 1,nens
            iiens = vstart+iens-1
            RMAT(:,iens) = yhat(:,iens) - qM(:,iens)
         end do
         deallocate(qM)
         deallocate(MMAT)


         !qB = Q * B^T
         if (iter .gt. 1) then
            allocate(qB(1:cv_size,nmodes)
            qB = 0.D0
            do iens = 1,nens
               iiens = vstart+iens-1
               do jens = 1, nens
                  qB(:,jens) = qB(:,jens) + qhat(:,iiens) * BMAT(jens,iens,iter-1)
               end do
            end do
            do iens = 1,nens
               RMAT(:,iens) = RMAT(:,iens) - qB(:,iens)
            end do
            deallocate(qB)
         end if

         vstart = nens*(rand_inner) + 1
         vend = vstart + nens - 1

!!!!!!!!!
         !Decompose RMAT into QR = qhat(:,vstart:vend) * BMAT(:,:,iter)
          !NEED QR decomposition that produces next sections of qhat and BMAT
!!!!!!!!!

         deallocate(RMAT)

         !Write BMAT(:,:,iter) to file
         write(cinner,fmt='(i4.4)') iter
         filename = 'BMAT.iter'//trim(adjustl(cinner))// &
                        '.p'//trim(adjustl(cproc))
         call da_mat_io ('w', filename, nens, nens, BMAT(:,:,iter))
         deallocate(BMAT)

         !Orthonormalize most recent qhat vectors w.r.t. previous vectors
         call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         qhat,cv_size,nmodes+nens,nmodes+1)


         write(cinner,fmt='(i4.4)') rand_inner+1
         do iens = 1, nens
            write(cens,fmt='(i4.4)') iens
            write(unit=stdout,fmt='("Writing qhat vector for next iteration; iens,iter = ",i3,",",i3)') iens, rand_inner+1

            filename = 'qhat.e'//trim(adjustl(cens))// &
                           '.iter'//trim(adjustl(cinner))// &
                           '.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, qhat(1:cv_size,vstart+iens-1), grid, mz)
         end do
      end if !end prep for next iteration

      deallocate(qhat)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: Block_Lanczos')
#endif

   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: For RSVD5.6, rand_stage should be 0, 1, or 2, not ',rand_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_block_lanczos")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for block Lanczos"/))
#endif

end subroutine da_block_lanczos


