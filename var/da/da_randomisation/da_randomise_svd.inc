subroutine da_randomise_svd(grid, mz, config_flags,        &
                           cv_size, xbx, be, iv, cv,   &
                           re, y, xhat, nsamp, j_cost,  &
                           qhat, eignval, eignvec, neign, hessian_info)

   !-------------------------------------------------------------------------
   ! Purpose:         Main randomized SVD routine, separated into 2 stages 
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
   integer, intent(in)               :: nsamp                     ! Total number of eigenpairs
   type (j_type), optional, intent(out) :: j_cost                ! cost function
   real, optional, intent(out)       :: qhat(1:cv_size, 0:nsamp)  ! orthonormal basis set of A = Q * B * Q^T
   real,      optional, intent(out) :: eignvec(nsamp, nsamp)
   real,      optional, intent(out) :: eignval(nsamp)
   integer,   optional, intent(out) :: neign                         ! Number of non-truncated eigenpairs
   type(hessian_type), optional, intent(inout) :: hessian_info

#if defined(LAPACK)
   real                         :: LAM(1:nsamp)
   real                         :: hess_trace(1:nsamp)
   real, allocatable            :: yhat(:,:)  ! yhat = A * omega
   real, allocatable            :: omega(:,:) ! Gaussian draw
  
   integer                      :: iens, jens, i, j, itype, np
   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real                         :: ghat(1:cv_size)        ! cv copy.

   real                         :: beta
   real                         :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=4) :: cnsamp=""
   logical          :: swapped

   if (trace_use) call da_trace_entry("da_randomise_svd")

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

   if ( rand_stage .eq. 1 .and. ensmember .gt. 0) then
      ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = L^T * H^T * R^-1 * H * L 
      write(cens,fmt='(i4.4)') ensmember

      if ( ensmember .eq. nsamp+1 ) then 
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

         ! Write ghat to file
         filename = trim(randvecdir)//'ghat.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, ghat, grid, mz)

      else if (ensmember .gt. 0 .and. ensmember .le. nsamp) then
         !Generate OMEGA
         allocate(omega(1:cv_size,1))
         if (read_omega) then
            ! Read in omega vector for this ensemble member
            write(unit=stdout,fmt='("Reading omega ensemble member ",i4," of ",i4)') ensmember,nsamp
            filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size,1), grid, mz)
         else
            call da_set_randomcv(cv_size, omega(:,1))
         end if
         if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,1) = 0.D0
         if ( be%cv%size_je.gt.0) omega(je_start:je_end,1) = 0.D0
         if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,1) = 0.D0

         if (rotate_omega .and. spectral_precon .and. rand_outer.gt.1) then
            call da_rotate_omega(grid, cv_size, mz, jp_start, jp_end, &
                                    rand_outer-1, omega(:,1), omega(:,1))

!            !Included in da_rotate_omega to save hessian I/O time
!            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
!                                    rand_outer-1, omega(:,1), omega(:,1), 'IT' )
         end if

         if (.not. read_omega .or. (rotate_omega .and. spectral_precon .and. rand_outer.gt.1)) then
            !Write omega to file
            filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1), grid, mz)
         end if

         ! Multiply A * OMEGA
         allocate(yhat(1:cv_size,1))
         yhat = 0.D0
         call da_calculate_gradj(rand_outer, 1, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size,1),y,yhat(1:cv_size,1),grid,config_flags)

         deallocate(omega)

         ! Write yhat to file
         filename = trim(randvecdir)//'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1), grid, mz)

         deallocate(yhat)
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nsamp+1, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if
   else if ( rand_stage .eq. 2 ) then
      allocate(yhat(1:cv_size,nsamp))
      allocate(omega(1:cv_size,nsamp))


#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      ! Read gradient
      filename = trim(randvecdir)//'ghat.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat, grid, mz)

      ! Read in yhat and omega vectors from all ensembles on each processor
      do iens = 1, nsamp
         write(unit=stdout,fmt='("Reading yhat and omega ensemble member : ",i4)') iens
         write(cens,fmt='(i4.4)') iens
         filename = trim(randvecdir)//'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,iens), grid, mz)

         filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, omega(1:cv_size,iens), grid, mz)
      end do
      if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,:) = 0.D0
      if ( be%cv%size_je.gt.0) omega(je_start:je_end,:) = 0.D0
      if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,:) = 0.D0

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: READ_ghat_yhat_omega')
#endif

      qhat(:,1:nsamp) = yhat(:,1:nsamp)
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nsamp)

!REMOVE EVENTUALLY or add option since very useful for debugging
!      ! Write qhat vectors to disk
!      do iens = 1, nsamp
!         write(unit=stdout,fmt='("Writing qhat vectors to disk: ",i4)') iens
!         write(cens,fmt='(i4.4)') iens
!         filename = trim(trim(randvecdir))//'qhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', filename, cv_size, qhat(1:cv_size,iens), grid, mz)
!      end do
!REMOVE EVENTUALLY or add option since very useful for debugging
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega, yhat, qhat, cv_size, &
                     eignval, eignvec, nsamp, hess_trace_approx = hess_trace)

      ! Add 1.0 to get eigenvalues of full preconditioned Hessian
      eignval = eignval + 1.D0

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: rsvd56')
#endif

      if ( var4d_inc_out &
           .and. present(hessian_info) &
           .and. allocated(hessian_info%ghat) ) then
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif
         write(cnsamp,fmt='(I0)') nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'HESS_TRACE: ',hess_trace(1:nsamp)

         hessian_info%ghat = ghat

         hessian_info%nmodes_max = nsamp
         do iens = 0, nsamp
            hessian_info%qhat(:,iens) = qhat(:,iens)
         end do

         do iens = nsamp, 1, -inc_out_interval
            allocate(hessian_info%eigs(iens)%eignval(1:iens))
            allocate(hessian_info%eigs(iens)%eignvec(1:iens,1:iens))

            hessian_info%eigs(iens)%neign = iens
            if (iens.eq.nsamp) then
               hessian_info%eigs(iens)%eignval = eignval
               hessian_info%eigs(iens)%eignvec = eignvec
            else
               call da_rsvd56(grid, mz, jp_start, jp_end, &
                             omega(:,1:iens), yhat(:,1:iens), qhat(:,0:iens), cv_size, &
                             hessian_info%eigs(iens)%eignval, hessian_info%eigs(iens)%eignvec, &
                             iens)
               hessian_info%eigs(iens)%eignval = hessian_info%eigs(iens)%eignval + 1.D0
            end if
         end do
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: increment_diagnostics2')
#endif
      end if


      ! Print current eigenpairs
      ! ----------------------------------------------
      write(cnsamp,fmt='(I0)') nsamp
      write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'Eigenvalues(r): ',eignval(1:nsamp)

!      do iens = 1, nsamp
!         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nsamp)
!      end do

      neign = nsamp
      write(unit=stdout,fmt='(A)') ' '
      write(unit=stdout,fmt='(A,I5,A,I5,A)') 'Using ', neign, ' of ', nsamp, ' eigenpairs'
      write(unit=stdout,fmt='(A)') ' '

      np = neign
      if (svd_p .gt. 0 .and. neign .gt. svd_p) np = neign - svd_p
      if ((minval(eignval(1:np)) .lt. 2.D0 .and. adapt_svd.eq.1) .or. adapt_svd.gt.1) then
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
         itype = 2

         write(unit=stdout,fmt='(A)') 'Carrying out RIOT-56 amat_mul, LRU'
      else
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
         itype = 1

         write(unit=stdout,fmt='(A)') 'Carrying out RIOT-56 amat_mul, LRA'
      end if

      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             nsamp,nsamp,neign, &
                             LAM,eignvec,qhat, &
                             ghat,xhat)

      if (itype .eq. 2) then !LRU
         xhat = xhat + ghat
      end if

      if (spectral_precon .and. rand_outer.gt.1) &
         call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 rand_outer-1, xhat, xhat, 'TL' )

      if (spectral_precon) &
         call da_hessian_io ('w',grid, cv_size,rand_outer,neign,neign,&
                            eignvec(1:neign,1:neign),eignval(1:neign),qhat(:,0:neign),mz,ghat)

      deallocate(yhat)
      deallocate(omega)
   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: For RSVD5.6, rand_stage should be 0, 1, or 2, not ',rand_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_randomise_svd")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd


