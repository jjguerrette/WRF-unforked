subroutine da_randomise_svd(grid, mz, config_flags,        &
                           cv_size, xbx, be, iv, cv,   &
                           re, y, xhat, nens, j_cost,  &
                           qhat, eignval, eignvec, neign)

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
   integer, intent(in)               :: cv_size                ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                    ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                     ! background error structure.
   type (iv_type), intent(inout)     :: iv                     ! ob. increment vector.
   real*8, intent(inout)             :: cv(1:cv_size)          ! control variable (local).
   type (y_type), intent(inout)      :: re                     ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                      ! y = H(x_inc) structure.
   real*8, intent(out)               :: xhat(1:cv_size)        ! control variable (local).
   integer, intent(in)               :: nens                   ! Total number of eigenpairs
   type (j_type), optional, intent(out) :: j_cost                        ! cost function

   real*8, optional, intent(out)     :: qhat(1:cv_size, 0:nens) ! orthonormal basis set of A = Q * B * Q^T
   real*8, optional, intent(out)     :: eignvec(nens, nens)
   real*8, optional, intent(out)     :: eignval(nens)
   integer, optional, intent(out)    :: neign                         ! Number of non-truncated eigenpairs


#if defined(LAPACK)

   real(wp), dimension(1:nens) :: LAM
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
  
   integer                      :: iens, jens, i, j, itype, np
   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real*8                       :: ghat(1:cv_size)        ! cv copy.
   real(wp)                     :: F1(nens, nens) ! intermediate matrix
   real(wp)                     :: F2(nens, nens) ! intermediate matrix
   real(wp)                     :: U(nens, nens)
   real(wp)                     :: VT(nens, nens)
   real(wp)                     :: SIG(nens)
   real(wp), allocatable        :: TEMP1(:,:), TEMP2(:,:)


   real(wp)                     :: KMAT(nens, nens)  ! B = Q^T * A * Q

   real*8                       :: beta
   real*8                       :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=3) :: cnens=""
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

   if ( svd_stage .eq. 1 .and. ensmember .gt. 0) then
      ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = U^T * H^T * R^-1 * H * U 
      allocate(yhat(1:cv_size,1))
      yhat = 0.D0

      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember

      if ( ensmember .eq. nens+1 ) then 
         ! Calculate J and grad J (ghat)
         use_randomsvd = .false.
         call da_calculate_j(svd_outer, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                             be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)
         use_randomsvd = .true.

         j0_total = j_cost%total
         if (j0_total == 0.0) return

         beta = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))

         write(unit=stdout,fmt='("Starting outer iteration : ",i3)') svd_outer
         write(unit=stdout,fmt=11) j0_total, beta
11       format('Starting cost function: ' ,1PD15.8,', Gradient: ',1PD15.8)

         ghat = - ghat

         ! Write ghat to file
         filename = 'ghat.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, ghat, grid, mz)

         if (any(riot_precon.eq.(/1,2,3,4/)) .and. svd_outer.gt.1 .and. sum(ntmax(1:svd_outer-1)).gt.nens) &
            call da_precon_omega(grid, cv_size, mz, jp_start, jp_end, nens, ghat)

      else if (ensmember .gt. 0 .and. ensmember .le. nens) then
         allocate(omega(1:cv_size,1))
         if (read_omega) then
            ! Read in omega vector from this ensemble on each processor
            write(unit=stdout,fmt='("Reading omega ensemble member ",i3," of ",i3)') iens,nens
            filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size,1), grid, mz)
         else
            ! call da_set_randomcv(cv_size, omega(:,1)) << LT
            call da_gen_omega(omega, cv_size, 1, grid, mz)
            if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,1) = 0.D0
            if ( be%cv%size_je.gt.0) omega(je_start:je_end,1) = 0.D0
            if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,1) = 0.D0
         end if

         if (rotate_omega .and. spectral_precon .and. svd_outer.gt.1) then
            call da_rotate_omega(grid, cv_size, mz, jp_start, jp_end, &
                                    svd_outer-1, omega(:,1), omega(:,1))

!            !Included in da_rotate_omega to save I/O time
!            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
!                                    svd_outer-1, omega(:,1), omega(:,1), 'IT' )
         end if

         call da_calculate_gradj(svd_outer, 1, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size,1),y,yhat(1:cv_size,1),grid,config_flags)

         !Write omega to file
         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1), grid, mz)
         deallocate(omega)

         ! Write yhat to file
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1), grid, mz)
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nens+1, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')

      end if

      deallocate(yhat)

   else if ( svd_stage .eq. 2 ) then
      allocate(yhat(1:cv_size,nens))
      allocate(omega(1:cv_size,nens))

      write(cproc,fmt='(i4.4)') myproc

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      ! Read gradient
      filename = 'ghat.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat, grid, mz)

      ! Read in yhat and omega vectors from all ensembles on each processor
      do iens = 1, nens
         write(unit=stdout,fmt='("Reading yhat and omega ensemble member : ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,iens), grid, mz)

         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, omega(1:cv_size,iens), grid, mz)
      end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: READ_ghat_yhat_omega')
#endif

      if (prepend_rsvd_basis(svd_outer) .eq. 1) then
         qhat(:,2:nens)  = yhat(:,1:nens-1) 
         qhat(:,1) = ghat
      else
         qhat(:,1:nens) = yhat(:,1:nens)
      end if

      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nens)

!REMOVE EVENTUALLY or add option
!      ! Write qhat vectors to disk
!      do iens = 1, nens
!         write(unit=stdout,fmt='("Writing qhat vectors to disk: ",i3)') iens
!         write(cens,fmt='(i4.4)') iens
!         filename = 'qhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', filename, cv_size, qhat(1:cv_size,iens), grid, mz)
!      end do
!REMOVE EVENTUALLY or add option

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega, yhat, qhat, cv_size, &
                     eignval, eignvec, nens, neign)

      ! Add 1.0 to get eigenvalues of full preconditioned Hessian
      eignval = eignval + 1.D0

      ! Print current eigenpairs
      ! ----------------------------------------------
      write(cnens,fmt='(I0)') nens
      write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'Eigenvalues(r): ',eignval(1:nens)-1.D0

      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nens)
      end do

      write(unit=stdout,fmt='(A)') ' '
      write(unit=stdout,fmt='(A,I5,A,I5,A)') 'Using ', neign, ' of ', nens, ' eigenpairs'
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
                             nens,nens,neign, &
                             LAM,eignvec,qhat, &
                             ghat,xhat)

      if (itype .eq. 2) then !LRU
         xhat = xhat + ghat
      end if

      if (spectral_precon .and. svd_outer.gt.1) &
         call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 svd_outer-1, xhat, xhat, 'TL' )

      if (riot_precon.ge.1 .or. spectral_precon) &
         call da_hessian_io ('w',grid, cv_size,svd_outer,nens,neign,&
                            eignvec,eignval,qhat,mz,ghat)

      if (riot_precon.gt.0 .and. sum(ntmax(1:svd_outer)).lt.ntmax(svd_outer+1)) then
         if ( riot_precon.eq.1 ) riot_precon=14
         if ( riot_precon.eq.3 ) riot_precon=15
      end if

      if ( any(riot_precon.eq.(/12,13,14,15/)) ) &
         call da_precon_omega(grid, cv_size, mz, jp_start, jp_end, ntmax(svd_outer+1), ghat)

      deallocate(yhat)
      deallocate(omega)
   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: For RSVD5.6, svd_stage should be 0, 1, or 2, not ',svd_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_randomise_svd")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd


