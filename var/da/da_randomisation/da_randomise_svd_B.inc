subroutine da_randomise_svd_B(grid, mz, config_flags,            &
                           cv_size, xbx, be, iv, cv, &
                           re, y, xhat, nsamp, j_cost, &
                           qhat, eignval, eignvec, neign)

   !-------------------------------------------------------------------------
   ! Purpose:         Main randomized SVD routine, separated into 2 stages (B only)
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
   real, intent(inout)               :: cv(1:cv_size)            ! control variable (local).
   type (y_type), intent(inout)      :: re                       ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                        ! y = H(x_inc) structure.
   real, intent(out)                 :: xhat(1:cv_size)          ! control variable (local).
   integer, intent(in)               :: nsamp                     ! Total number of eigenpairs
   type (j_type), optional, intent(out) :: j_cost                ! cost function
   real, optional, intent(out)       :: qhat(1:cv_size, 0:nsamp)  ! orthonormal basis set of A = Q * B * Q^T
   real*8, optional, intent(out)     :: eignvec(nsamp, nsamp)
   real*8, optional, intent(out)     :: eignval(nsamp)
   integer, optional, intent(out)    :: neign                         ! Number of non-truncated eigenpairs

#if defined(LAPACK)

   real(wp), dimension(1:nsamp) :: LAM, LAMinv
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
   real(wp), allocatable        :: nuhat(:,:) ! Eigenvectors of A
   real(wp)                     :: dot_cv
  
   integer                      :: iens, jens, i, j, ccount, p
   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real                         :: ghat(1:cv_size)        ! cv copy.
   real(wp)                     :: F1(nsamp, nsamp) ! intermediate matrix
   real(wp)                     :: F2(nsamp, nsamp) ! intermediate matrix
   real(wp)                     :: U(nsamp, nsamp)
   real(wp)                     :: VT(nsamp, nsamp)
   real(wp)                     :: SIG(nsamp)
   real(wp), allocatable        :: TEMP1(:,:), TEMP2(:,:)


   real(wp)                     :: KMAT(nsamp, nsamp)  ! B = Q^T * A * Q
   integer                      :: piv(nsamp)
   integer                      :: info

   real                         :: gdot, beta, temp_l
   real                         :: temp_e(nsamp)
   real                         :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=4) :: cnsamp=""
   logical          :: swapped
  
   if (trace_use) call da_trace_entry("da_randomise_svd_B")

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

   if ( rand_stage .eq. 1 .and. ensmember .gt. 0) then
      ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = U^T * H^T * R^-1 * H * U 
      allocate(yhat(1:cv_size,1))
      yhat = 0.D0

      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember

      if ( ensmember .eq. nsamp+1 ) then 
         use_randomblock = .false.
         call da_calculate_j(rand_outer, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                             be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)
         use_randomblock = .true.

         j0_total = j_cost%total
         if (j0_total == 0.0) return
         yhat(1:cv_size,1) = - ghat

         ! Write ghat to file
         filename = 'ghat.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1))

      else if (ensmember .gt. 0 .and. ensmember .le. nsamp) then
         allocate(omega(1:cv_size,1))
         if (read_omega) then
            ! Read in omega vector from this ensemble on each processor
            write(unit=stdout,fmt='("Reading omega ensemble member : ",i4)') iens
            filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size,1))
         else
            call da_set_randomcv(cv_size, omega(:,1))
            if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,1) = 0.D0
            if ( be%cv%size_je.gt.0) omega(je_start:je_end,1) = 0.D0
            if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,1) = 0.D0
         end if

#if (WRF_CHEM == 1)
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                omega(1:be%cv%size_jb,1), grid%vv, grid%vp, grid%vch)
#else
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                omega(1:be%cv%size_jb,1), grid%vv, grid%vp)
#endif
         yhat(1:cv_size,1) = 0.D0
#if (WRF_CHEM == 1)
         call da_transform_vtox_adj(grid, be%cv%size_jb, xbx, be, grid%ep, grid%vp, grid%vv, &
                                yhat(1:be%cv%size_jb,1), grid%vch)
#else
         call da_transform_vtox_adj(grid, be%cv%size_jb, xbx, be, grid%ep, grid%vp, grid%vv, &
                                yhat(1:be%cv%size_jb,1))
#endif


!         call da_calculate_gradj(rand_outer, 1, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
!                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size,1),y,yhat(1:cv_size,1),grid,config_flags)

         !Write omega to file
         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1))
         deallocate(omega)

         ! Write yhat to file
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1))
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nsamp+1, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')

      end if

      deallocate(yhat)

   else if ( rand_stage .eq. 2 ) then
      allocate(yhat(1:cv_size,nsamp))
      allocate(omega(1:cv_size,nsamp))
      qhat(:,0) = 0.D0

      write(cproc,fmt='(i4.4)') myproc

      ! Read in yhat and omega vectors from all ensembles on each processor
      do iens = 1, nsamp
         write(unit=stdout,fmt='("Reading yhat and omega ensemble member : ",i4)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,iens))

         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, omega(1:cv_size,iens))
      end do

!      ! Read gradient
!      filename = 'ghat.p'//trim(adjustl(cproc))
!      call da_cv_io ('r', filename, cv_size, ghat)

      ! Generate qhat basis vectors by Modified Gram-Schmidt orthonormalization
      do iens = 1, nsamp
         qhat(:,iens) = yhat(:,iens)
         if (iens .gt. 1) then
            do i = 1,2 !Second iteration improves orthogonalization for numerical precision (Bjorck, 1994)
            do jens = 1, iens-1
               gdot = da_dot_cv(cv_size, qhat(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
               qhat(:,iens) = qhat(:,iens) - gdot * qhat(:,jens)
            end do
            end do
         end if

         beta   = SQRT(da_dot_cv(cv_size, qhat(:,iens), qhat(:,iens), grid, mz, jp_start, jp_end))
         qhat(:,iens) = qhat(:,iens) / beta
      end do

      ! Write qhat vectors to disk
      do iens = 1, nsamp
         write(unit=stdout,fmt='("Writing qhat vectors to disk: ",i4)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'qhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, qhat(1:cv_size,iens))
      end do

      ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
      do iens = 1, nsamp
         do jens = 1, nsamp
            F1(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), omega(:,jens), grid, mz, jp_start, jp_end)
            F2(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), yhat(:,jens), grid, mz, jp_start, jp_end)
         end do
      end do

    ! Solve system D * x = b, &
    !  where B * F1^T = F2^T, F1^T = Q^T * OMEGA, F2^T = Q^T * Y
    !  D = F1
    !  b = F2(:,iens)
    !  x = KMAT(iens,:)
    !-----------------------------------------------------------
!      write(cnsamp,fmt='(I0)') nsamp
!      do iens = 1, nsamp
!         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'F1: ',F1(iens,1:nsamp)
!      end do
!      do iens = 1, nsamp
!         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'F2: ',F2(iens,1:nsamp)
!      end do

      info = 0
      p = -1
      if (p .eq. 0) then
         call GESV(F1, F2, piv, info) !F95

         ! GESV returns B^T in F2
         do iens = 1, nsamp
            do jens = 1, nsamp
               KMAT(iens,jens) = F2(jens,iens)
            end do
         end do
      else
         if (p .lt. 0) p = 0
         call GESVD(F1(1:nsamp,1:nsamp-p), &
                    SIG(1:nsamp-p), &
                    U = U(1:nsamp,1:nsamp-p), &
                    VT = VT(1:nsamp-p,1:nsamp-p), &
                    INFO = info)

         allocate(TEMP1(1:nsamp-p,1:nsamp-p))
         allocate(TEMP2(1:nsamp-p,1:nsamp-p))
         do iens = 1, nsamp-p
            TEMP2(iens,:) = (1.D0 / SIG(iens)) * VT(iens,:)
            do jens = 1, nsamp-p
               TEMP1(jens,iens) = sum(F2(1:nsamp,iens) * U(1:nsamp,jens))
            end do
         end do


         do iens = 1, nsamp-p
            do jens = 1, nsamp-p
                  KMAT(iens,jens) = sum(TEMP1(:,iens) * TEMP2(:,jens))
            end do
         end do
         deallocate(TEMP1)
         deallocate(TEMP2)
      end if
      if (info /=0) then
         write(stdout,*) 'Error in randomised SVD: GESV returned ',info
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if
      ! Force symmetry in B
      do iens = 1, nsamp
         do jens = min(iens+1,nsamp), nsamp
            KMAT(iens,jens) = (KMAT(iens,jens) + KMAT(jens,iens)) / 2.D0
            KMAT(jens,iens) = KMAT(iens,jens)
         end do
      end do

      write(cnsamp,fmt='(I0)') nsamp
      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'E27.16E3)') 'KMAT: ',KMAT(iens,1:nsamp)
      end do

    ! Determine eigenvalues and eigenvectors of the symmetric B matrix
    !-----------------------------------------------------------------
!      ndgewrk           = MAX(4*nsamp-1,1)
      info              = 0
      call SYEV(KMAT, LAM, JOBZ = 'V' , INFO = info) !F95

      if (info /=0) then
         write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if

      do iens = 1, nsamp
         do jens = 1, nsamp 
            eignvec(iens,jens) = KMAT(iens,jens)
         end do
      end do
      
      !First sort eigenvalues and eigenvectors in descending order
      do j = nsamp-1, 1, -1
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

      neign = nsamp
      do iens = 1, nsamp
         eignval(iens) = LAM(iens) + 1
         if (neign .eq. nsamp .and. LAM(iens) .lt. 0.D0) then
            if (iens .gt. 1) then
               neign = iens - 1
            else
               write(unit=stdout,fmt='(A)') &
                  'FATAL ERROR: There need to be some positive eigenvalues of KMAT'
               CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
            end if
         end if
      end do

      ! Print current eigenpairs
      ! ----------------------------------------------
      write(cnsamp,fmt='(I0)') nsamp
      write(unit=stdout,fmt='(A,'//trim(cnsamp)//'D25.16)') 'Eigenvalues(r): ',LAM(1:nsamp)
      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,'//trim(cnsamp)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nsamp)
      end do

      write(unit=stdout,fmt='(A)') ' '
      write(unit=stdout,fmt='(A,I5,A,I5,A)') 'Showing ', neign, ' of ', nsamp, ' eigenpairs'
      write(unit=stdout,fmt='(A)') ' '

      deallocate(yhat)
      deallocate(omega)
   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: rand_stage should be 0, 1, or 2, not ',rand_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_randomise_svd_B")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd_B


