subroutine da_randomise_svd_B11(grid, mz, config_flags,            &
                           cv_size, xbx, be, iv, cv, &
                           re, y, xhat, j_cost, &
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
   integer, intent(in)               :: cv_size                       ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                           ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                            ! background error structure.
   type (iv_type), intent(inout)     :: iv                            ! ob. increment vector.
   real, intent(inout)               :: cv(1:cv_size)                 ! control variable (local).
   type (y_type), intent(inout)      :: re                            ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                             ! y = H(x_inc) structure.
   real, intent(out)                 :: xhat(1:cv_size)               ! control variable (local).
   type (j_type), optional, intent(out) :: j_cost                        ! cost function

   real, optional, intent(out)       :: qhat(1:cv_size, 0:ensdim_svd) ! orthonormal basis set of A = Q * B * Q^T
   real*8, optional, intent(out)     :: eignvec(ensdim_svd, ensdim_svd)
   real*8, optional, intent(out)     :: eignval(ensdim_svd)
   integer, optional, intent(out)    :: neign                         ! Number of non-truncated eigenpairs

#if defined(LAPACK)

   real(wp), dimension(1:ensdim_svd) :: LAM, LAMinv
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
   real(wp), allocatable        :: nuhat(:,:) ! Eigenvectors of A
   real(wp)                     :: dot_cv
  
   integer                      :: iens, jens, nens, i, j, ccount, p
   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   real                         :: ghat(1:cv_size)        ! cv copy.
   real(wp)                     :: F1(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: F2(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: U(ensdim_svd, ensdim_svd)
   real(wp)                     :: VT(ensdim_svd, ensdim_svd)
   real(wp)                     :: SIG(ensdim_svd)
   real(wp), allocatable        :: TEMP1(:,:), TEMP2(:,:)

   integer                      :: svd_unit
   real(wp)                     :: KMAT(ensdim_svd, ensdim_svd)  ! B = Q^T * A * Q
   real(wp)                     :: KMAT_antidiag(ensdim_svd,2)
   real(wp)                     :: KMAT_maxerr(ensdim_svd)
   integer                      :: KMAT_maxerri(ensdim_svd)

   integer                      :: piv(ensdim_svd)
   integer                      :: info

   real                         :: gdot, beta, temp_l
   real                         :: temp_e(ensdim_svd)
   real                         :: j0_total, erri
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=5) :: cnens=""
   logical          :: swapped
  
   if (trace_use) call da_trace_entry("da_randomise_svd_B11")

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

   nens = ensdim_svd

   xhat = 0.D0

   !============================================================================================
   !SVD STAGE1
   ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = U^T * H^T * R^-1 * H * U 
   allocate(yhat(1:cv_size,nens))

   write(cproc,fmt='(i4.4)') myproc

   allocate(omega(1:cv_size,nens))

   do iens = 1,nens
      write(cens,fmt='(i4.4)') iens

      yhat(1:cv_size,iens) = 0.D0

!      if (read_omega) then
!         ! Read in omega vector from this ensemble on each processor
!         write(unit=stdout,fmt='("Reading omega ensemble member : ",I5)') iens
!         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('r', filename, cv_size, omega(1:cv_size,iens))
!      else
         write(unit=stdout,fmt='("Multiplying UU^T * omega for ens member",I5)') iens
         call da_set_randomcv(cv_size, omega(:,iens))
         if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,1) = 0.D0
         if ( be%cv%size_je.gt.0) omega(je_start:je_end,1) = 0.D0
         if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,1) = 0.D0
!      end if

#if (WRF_CHEM == 1)
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                             omega(1:be%cv%size_jb,iens), grid%vv, grid%vp, grid%vch)
#else
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                             omega(1:be%cv%size_jb,iens), grid%vv, grid%vp)
#endif
#if (WRF_CHEM == 1)
      call da_transform_vtox_adj(grid, be%cv%size_jb, xbx, be, grid%ep, grid%vp, grid%vv, &
                             yhat(1:be%cv%size_jb,iens), grid%vch)
#else
      call da_transform_vtox_adj(grid, be%cv%size_jb, xbx, be, grid%ep, grid%vp, grid%vv, &
                             yhat(1:be%cv%size_jb,iens))
#endif


!      !Write omega to file
!      filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!      call da_cv_io ('w', filename, cv_size, omega(1:cv_size,iens))

!      ! Write yhat to file
!      filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!      call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,iens))
   end do

   !============================================================================================
   !!SVD STAGE2
   qhat(:,0) = 0.D0

   write(cproc,fmt='(i4.4)') myproc


   ! Generate qhat basis vectors by Modified Gram-Schmidt orthonormalization
   write(unit=stdout,fmt='("Generating qhat")')
   do iens = 1, nens
      write(unit=stdout,fmt='("iens = ",I5)') iens
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

!   ! Write qhat vectors to disk
!   do iens = 1, nens
!      write(unit=stdout,fmt='("Writing qhat vectors to disk: ",I5)') iens
!      write(cens,fmt='(i4.4)') iens
!      filename = 'qhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!      call da_cv_io ('w', filename, cv_size, qhat(1:cv_size,iens))
!   end do

   write(unit=stdout,fmt='("Generating F1 and F2")')
   ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
   do iens = 1, nens
      write(unit=stdout,fmt='("iens = ",I5)') iens
      do jens = 1, nens
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
   write(cnens,fmt='(I0)') nens
!      do iens = 1, nens
!         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F1: ',F1(iens,1:nens)
!      end do
!      do iens = 1, nens
!         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F2: ',F2(iens,1:nens)
!      end do

   info = 0
   p = -1
   if (p .eq. 0) then
      write(unit=stdout,fmt='("Running GESV")')
      call GESV(F1, F2, piv, info) !F95

      ! GESV returns B^T in F2
      do iens = 1, nens
         do jens = 1, nens
            KMAT(iens,jens) = F2(jens,iens)
         end do
      end do
   else
      if (p .lt. 0) p = 0
      write(unit=stdout,fmt='("Running GESVD")')
      call GESVD(F1(1:nens,1:nens-p), &
                 SIG(1:nens-p), &
                 U = U(1:nens,1:nens-p), &
                 VT = VT(1:nens-p,1:nens-p), &
                 INFO = info)

      allocate(TEMP1(1:nens-p,1:nens-p))
      allocate(TEMP2(1:nens-p,1:nens-p))
      do iens = 1, nens-p
         TEMP2(iens,:) = (1.D0 / SIG(iens)) * VT(iens,:)
         do jens = 1, nens-p
            TEMP1(jens,iens) = sum(F2(1:nens,iens) * U(1:nens,jens))
         end do
      end do


      do iens = 1, nens-p
         do jens = 1, nens-p
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

   write(unit=stdout,fmt='("Printing KMAT symmetry")')
   do iens = 1, nens
      KMAT_antidiag(iens,1)=KMAT(nens-iens+1,iens)
      KMAT_antidiag(iens,2)=KMAT(iens,nens-iens+1)
      KMAT_maxerr(iens) = 0.D0
!!      do jens = min(iens+1,nens), nens
      do jens = 1, max(1,iens-1)
         erri = abs(2.D0 * (KMAT(iens,jens) - KMAT(jens,iens)) / &
                                (KMAT(iens,jens) + KMAT(jens,iens)))
         if ( erri .gt. KMAT_maxerr(iens) ) then
            KMAT_maxerr(iens) = erri
            KMAT_maxerri(iens) = jens
         end if
      end do
   end do
   call da_get_unit (svd_unit)
   filename = 'KMAT_symm.txt'
   open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
   do iens = 1, nens
      write(unit=svd_unit,fmt='(I5,3D25.16,I6)') iens, KMAT_antidiag(iens,1:2), KMAT_maxerr(iens), KMAT_maxerri(iens)
   end do
   call da_free_unit (svd_unit)

   ! Force symmetry in B
   do iens = 1, nens
      do jens = min(iens+1,nens), nens
         KMAT(iens,jens) = (KMAT(iens,jens) + KMAT(jens,iens)) / 2.D0
         KMAT(jens,iens) = KMAT(iens,jens)
      end do
   end do

   ! Determine eigenvalues and eigenvectors of the symmetric B matrix
   !-----------------------------------------------------------------

   write(unit=stdout,fmt='("Running SYEV")')
   info              = 0
   call SYEV(KMAT, LAM, JOBZ = 'V' , INFO = info) !F95

   if (info /=0) then
      write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if

   do iens = 1, nens
      do jens = 1, nens 
         eignvec(iens,jens) = KMAT(iens,jens)
      end do
   end do

   write(unit=stdout,fmt='("Sorting Eigenmodes")')  
   !First sort eigenvalues and eigenvectors in descending order
   do j = nens-1, 1, -1
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

   neign = nens
   do iens = 1, nens
      eignval(iens) = LAM(iens)
!      if (neign .eq. nens .and. LAM(iens) .lt. 0.D0) then
!         if (iens .gt. 1) then
!            neign = iens - 1
!         else
!            write(unit=stdout,fmt='(A)') &
!               'FATAL ERROR: There need to be some positive eigenvalues of KMAT'
!            CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
!         end if
!      end if
   end do

   ! Print current eigenpairs
   ! ----------------------------------------------
   write(unit=stdout,fmt='("Printing Eigenpairs")')
   if (nens .le. 100) then
      write(cnens,fmt='(I0)') nens
      write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'Eigenvalues(r): ',LAM(1:nens)
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nens)
      end do
   else
      if(rootproc) then
         call da_get_unit (svd_unit)
         filename = 'eigenvalues.txt'
         open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
         write(unit=svd_unit,fmt='(D25.16)') LAM(1:nens)
         call da_free_unit (svd_unit)

         do iens = 1, nens
            call da_get_unit (svd_unit)
            write(cens,fmt='(i4.4)') iens
            filename = 'eigenvector.e'//trim(adjustl(cens))//'.txt'
            open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
            write(unit=svd_unit,fmt='(D25.16)') eignvec(1:nens,iens)
            call da_free_unit (svd_unit)
         end do
      end if
#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
#endif
   end if


   write(unit=stdout,fmt='(A)') ' '
   write(unit=stdout,fmt='(A,I5,A)') 'Showing ', nens, ' eigenpairs'
   write(unit=stdout,fmt='(A)') ' '

   deallocate(yhat)
   deallocate(omega)

   if (trace_use) call da_trace_exit("da_randomise_svd_B11")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd_B11


