subroutine da_randomise_svd(grid, config_flags,            &
                           cv_size, xbx, be, iv, cv, &
                           re, y, xhat, j_cost) !!, eignvec, LAM, neign)

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

!   real, optional, intent(out)       :: qhat(1:cv_size, 1:ensdim_svd) ! orthonormal basis set.
!   real*8, optional, intent(out)     :: eignvec(ensdim_svd, ensdim_svd)
!   real*8, optional, intent(out)     :: nuhat(1:cv_size, 1:ensdim_svd)
!   real*8, optional, intent(out)     :: LAM(1:ensdim_svd)
!   real*8, optional, intent(out)     :: LAMi(1:ensdim_svd)
!   real*8, optional, intent(out)     :: eignvec(ensdim_svd, ensdim_svd)
!!   real*8, optional, intent(out)     :: nuhat(1:cv_size, 1:ensdim_svd)
!   real*8, optional, intent(out)     :: LAM(1:ensdim_svd)
!   real*8, optional, intent(out)     :: LAMi(1:ensdim_svd)
!   integer, optional, intent(out)    :: neign                         ! Number of non-truncated eigenpairs

#if defined(LAPACK)

   real(wp)                     :: eignvec(ensdim_svd, ensdim_svd)
   real(wp)                     :: LAM(1:ensdim_svd), LAMi(1:ensdim_svd)
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: qhat(:,:)  ! orthonormal basis set of A = Q * B * Q^T = qhat * BMAT * qhat^T
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
   real(wp), allocatable        :: nuhat(:,:) ! Eigenvectors of A
   real(wp)                     :: dot_cv

   complex*16, allocatable      :: nuhat_z(:,:) ! Eigenvectors of A
   complex*16, allocatable      :: qhat_z(:,:)
   complex*16, allocatable      :: xhat_z(:)   ! complex control variable (local).
   complex*16, allocatable      :: ghat_z(:)   ! complex control variable (local).
   complex*16                   :: eignvec_z(ensdim_svd, ensdim_svd)
   complex*16                   :: LAM_z(1:ensdim_svd)
   complex*16                   :: dot_cv_z

   integer                      :: neign                         ! Number of non-truncated eigenpairs
   
   integer                      :: iens, jens, nens, i, j, ccount
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   integer                      :: mz(11)
#else
   integer                      :: mz(7)
#endif
#else
#ifdef CLOUD_CV
   integer                      :: mz(13)
#else
   integer                      :: mz(9)
#endif
#endif
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   real                         :: ghat(1:cv_size)        ! cv copy.
   real(wp)                     :: F1(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: F2(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: BMAT(ensdim_svd, ensdim_svd)  ! B = Q^T * A * Q
   integer                      :: piv(ensdim_svd)
!   real                         :: dgewrk(4*nensdim_svd)
!   real                         :: dummymat(1, ensdim_svd)
!   integer                      :: ndgewrk, info
   integer                      :: info

   real                         :: gdot, beta
   real                         :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=3) :: cnens=""
  
   if (trace_use) call da_trace_entry("da_randomise_svd")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using randomized SVD method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,  be%v6%mz, be%v7%mz,  be%v8%mz, be%v9%mz, be%alpha%mz, be % ne /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,be%alpha%mz, be % ne /)
#endif

#else

#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz,  be%v6%mz, be%v7%mz,  be%v8%mz, be%v9%mz, be%alpha%mz, be%ne, be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#else
   mz =  (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne, be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#endif
#endif

   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   nens = ensdim_svd

   xhat = 0.D0

   if ( svd_stage .eq. 1 .and. ensmember .gt. 0) then
      ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = U^T * H^T * R^-1 * H * U 
      allocate(yhat(1:cv_size,1))

      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember

      if ( ensmember .eq. nens+1 ) then 
         call da_calculate_j(svd_outer, 0, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                             be%cv%size_jl, xbx, be, iv, xhat, cv, re, y, j_cost, ghat, grid, config_flags)
         j0_total = j_cost%total
         if (j0_total == 0.0) return
         yhat(1:cv_size,1) = - ghat

         ! Write ghat to file
         filename = 'ghat.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1))

      else if (ensmember .gt. 0 .and. ensmember .le. nens) then
         allocate(omega(1:cv_size,1))
         if (read_omega) then
            ! Read in omega vector from this ensemble on each processor
            write(unit=stdout,fmt='("Reading omega ensemble member : ",i3)') iens
            filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size,1))
         else
            call da_set_randomcv(cv_size, omega(:,1))
         end if

         call da_calculate_gradj(svd_outer, 1, cv_size,be%cv%size_jb,be%cv%size_je,be%cv%size_jp, &
                                 be%cv%size_jl, xbx,be,iv,omega(1:cv_size,1),y,yhat(1:cv_size,1),grid,config_flags)

         !Write omega to file
         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1))
         deallocate(omega)

         ! Write yhat to file
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, yhat(1:cv_size,1))
      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nens+1, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')

      end if

      deallocate(yhat)

   else if ( svd_stage .eq. 2 ) then
      allocate(yhat(1:cv_size,nens))
      allocate(omega(1:cv_size,nens))
      allocate(qhat(1:cv_size,0:nens)) !! zero first for amat_mul later
      qhat(:,0) = 0.D0

      write(cproc,fmt='(i4.4)') myproc

      ! Read in yhat and omega vectors from all ensembles on each processor
      do iens = 1, nens
         write(unit=stdout,fmt='("Reading yhat and omega ensemble member : ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'yhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, yhat(1:cv_size,iens))

         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, omega(1:cv_size,iens))
      end do

      ! Read gradient
      filename = 'ghat.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat)

      ! Generate qhat basis vectors
      do iens = 1, nens
         qhat(:,iens) = yhat(:,iens)
         if (iens .gt. 1) then
            do jens = iens-1, 1, -1
               gdot = da_dot_cv(cv_size, qhat(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
               qhat(:,iens) = qhat(:,iens) - gdot * qhat(:,jens)
            end do
         end if

         beta   = SQRT(da_dot_cv(cv_size, qhat(:,iens), qhat(:,iens), grid, mz, jp_start, jp_end))
         qhat(:,iens) = qhat(:,iens) / beta
      end do

      ! Write qhat vectors to disk
      do iens = 1, nens
         write(unit=stdout,fmt='("Writing qhat vectors to disk: ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'qhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, qhat(1:cv_size,iens))
      end do

      ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
      do iens = 1, nens
         do jens = 1, nens
            F1(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), omega(:,jens), grid, mz, jp_start, jp_end)
            F2(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), yhat(:,jens), grid, mz, jp_start, jp_end)
         end do
      end do

    ! Solve system D * x = b, &
    !  where B * F1^T = F2^T, F1^T = Q^T * OMEGA, F2^T = Q^T * Y
    !  D = F1
    !  b = F2(:,iens)
    !  x = BMAT(iens,:)
    !-----------------------------------------------------------
      write(cnens,fmt='(I0)') ensdim_svd
      do iens = 1, ensdim_svd
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'F1: ',F1(iens,1:ensdim_svd)
      end do
      do iens = 1, ensdim_svd
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'F2: ',F2(iens,1:ensdim_svd)
      end do

      info = 0
!      call GESV(nens, nens, F1, nens, piv, F2, nens, info) !F77
      call GESV(F1, F2, piv, info) !F95

      if (info /=0) write(stdout,*) 'Error in randomised SVD: GESV returned ',info

      do iens = 1, nens
         do jens = 1, nens
            if (svd_symm_type .eq. 1) then
               BMAT(iens,jens) = F2(jens,iens)  !BMAT is not necessarily symmetric, and complex eigs can results
            else
               BMAT(iens,jens) = (F2(jens,iens) + F2(iens,jens)) / 2
            end if
         end do
      end do

      do iens = 1, ensdim_svd
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'BMAT: ',BMAT(iens,1:ensdim_svd)
      end do

    ! Determine eigenvalues and eigenvectors of the B matrix
    !---------------------------------------------------------
!      ndgewrk           = MAX(4*nens-1,1)
      info              = 0
      if (svd_symm_type .eq. 1) then
!      call GEEV('N', 'V', nens, BMAT, nens, LAM, LAMi, dummymat, 1, eignvec, nens, &
!              dgewrk(1:ndgewrk), ndgewrk, info) !F77
         call GEEV(BMAT, LAM, LAMi, VR = eignvec, INFO = info) !F95
         if (info /=0) write(stdout,*) 'Error in randomised SVD: GEEV returned ',info
      else
         call SYEV(BMAT, LAM, JOBZ = 'V' , INFO = info) !F95
         do iens = 1, nens
            do jens = 1, nens 
               eignvec(iens,jens) = BMAT(iens,jens)
            end do
         end do
         LAMi = 0.D0
         if (info /=0) write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
      end if


      neign = nens

      ! Print current eigenpairs
      ! ----------------------------------------------
      write(cnens,fmt='(I0)') ensdim_svd
      write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'Eigenvalues(r): ',LAM(1:ensdim_svd)
      write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'Eigenvalues(i): ',LAMi(1:ensdim_svd)
      do iens = 1, ensdim_svd
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'EIGENVECTOR: ',eignvec(iens,1:ensdim_svd)
      end do

      ! Note: da_amat_mul does not properly handle imaginary eigenpairs,
      !       nor does it perform proper truncation when neign < nens.
      ! Although symmetric real matrices should not have complex eigenvalues or eigenvectors


      if (any(LAMi .gt. 0) .and. svd_amat_type .le. 1) then
         ccount = 0
         do jens = 1, max(neign,nens)
            LAM_z(jens) = (1.D0, 0.D0) / ( dcmplx( LAM(jens) + 1.D0, LAMi(jens) ) )
            if (LAMi(jens) .gt. 0. .or. LAMi(jens) .lt. 0. ) then
              ccount = ccount + 1
              write(unit=stdout,fmt='(A,I5)') 'ccount = ', ccount
            end if
            do i = 1, max(neign,nens)
               if ((LAMi(jens) .gt. 0. .or. LAMi(jens) .lt. 0.)) then
                  if (mod(ccount,2) .eq. 1) then
                     eignvec_z(i,jens  ) = dcmplx(eignvec(i,jens),   eignvec(i,jens+1))
                     eignvec_z(i,jens+1) = dcmplx(eignvec(i,jens), - eignvec(i,jens+1))
                  end if
               else
                  eignvec_z(i,jens) = dcmplx(eignvec(i,jens), 0.D0)
               end if
            end do
         end do

         allocate(qhat_z(1:cv_size,1:nens))
         allocate(ghat_z(1:cv_size))
         allocate(xhat_z(1:cv_size))

         qhat_z = dcmplx (qhat, 0.D0)
         ghat_z = dcmplx (ghat, 0.D0)
         xhat_z = dcmplx(0.D0, 0.D0)

! LAST CHANGE INCASE OF ERROR
!         do i = 1, cv_size
!            xhat_z(i) = dcmplx(0.D0, 0.D0)
!            ghat_z(i) = dcmplx( ghat(i) , 0.D0 )
!            do jens = 1, nens
!               qhat_z(i,jens) = dcmplx( qhat(i,jens) , 0.D0 )
!            end do
!         end do

         if (svd_amat_type .eq. 1) then
            ! Generate eigenvectors of A
            write(unit=stdout,fmt='(A)') 'Carrying out truncated complex amat_mul'
            allocate(nuhat_z(cv_size, nens))

            call GEMM(qhat_z, eignvec_z, nuhat_z)

            ! Alternatively, if BLAS isn't available
!            do jens = 1, nens
!               do i = 1, cv_size
!                  nuhat_z(i,jens) = da_dot_z(neign, qhat_z(i,1:neign), eignvec_z(1:neign,jens))
!               end do
!            end do

            ! For a symmetric matrix, P = V*D*V^-1 = V*D*V^*, where ^* represents the conjugate transpose
            do i = 1, neign
               dot_cv_z = da_dot_cv_z(cv_size, conjg(nuhat_z(:,i)), ghat_z, grid, mz, jp_start, jp_end)
               write(unit=stdout,fmt='("i, dot_cv = ",I3,D25.16,SP,D25.16,"i")') i, real(dot_cv_z,8),real(imag(dot_cv_z),8)
               write(unit=stdout,fmt='("i, LAM_z = ",I3,D25.16,SP,D25.16,"i")') i, real(LAM_z(i),8),real(imag(LAM_z(i)),8)

               xhat_z = xhat_z + LAM_z(i) * nuhat_z(:,i) * dot_cv_z
            end do

            deallocate(nuhat_z)

         else if (svd_amat_type .eq. 0) then
            ! Copied directly from da_amat_mul, seems to work for complex conjugate eigen pairs
            write(unit=stdout,fmt='(A)') 'Carrying out full complex amat_mul'
            do j = 1, neign
               dot_cv_z = da_dot_cv_z(cv_size, qhat_z(:,j), ghat_z, grid, mz, jp_start, jp_end)
               write(unit=stdout,fmt='("j, dot_cv = ",I3,D25.16,SP,D25.16,"i")') j, real(dot_cv_z,8),real(imag(dot_cv_z),8)
               write(unit=stdout,fmt='("j, LAM_z = ",I3,D25.16,SP,D25.16,"i")') j, real(LAM_z(j),8),real(imag(LAM_z(j)),8)

               do i = 1, neign
                  xhat_z = xhat_z + qhat_z(:,i) * &
                         SUM(eignvec_z(i,1:neign)*LAM_z(1:neign)*eignvec_z(j,1:neign)) * dot_cv_z
               end do
            end do 
         end if

         ! Since xhat is in real space, the complex parts must cancel.
         do i = 1, cv_size
            xhat(i) = REAL( xhat_z(i) ,8)
            if ( REAL(IMAG( xhat_z(i) ),8) .gt. 1.D-10 ) then
               write(unit=stdout,fmt='(A,2I2)') &
                'ERROR: There should be no imaginary component to XHAT',i,cv_size
            end if
         end do

         deallocate(qhat_z)
         deallocate(ghat_z)
         deallocate(xhat_z)

      else
         ! Eigenvalues of the inverse Hessian (including all background terms)
         ccount = 0
         do i = 1, nens
            LAM(i) = 1.D0 / ( LAM(i) + 1.D0 )
            if ( LAMi(i) .lt. 0. .or. LAMi(i) .gt. 0. ) then
               if (ccount .eq. 0) neign = i - 1
               ccount = ccount + 1
            end if
         end do

         if (neign .eq. nens) then
            write(unit=stdout,fmt='(A)') 'Carrying out full real amat_mul'
            call da_amat_mul(be,grid,cv_size,nens,neign,LAM,eignvec,qhat,ghat,xhat)
         else
            !  Truncate to neign modes Here we use full eigenvec matrix to generate nuhat=qhat*eigenvec,
            !  then truncate by only using neign of the nuhat vectors 

            write(unit=stdout,fmt='(A)') 'Carrying out truncated amat_mul'

            ! (1) Generate eigenvectors of A
            allocate(nuhat(cv_size, neign))
            nuhat = 0.D0

            call GEMM(qhat(1:cv_size,1:nens), eignvec(1:nens,1:neign), nuhat)

            ! Alternatively, if BLAS isn't available
!            do jens = 1, neign
!               do i = 1, cv_size
!                  nuhat(i,jens) = da_dot(neign, qhat(i,1:neign), eignvec(1:neign,jens))
!               end do
!            end do

            ! Write nuhat vectors to disk
            do i = 1, neign
               write(unit=stdout,fmt='("Writing nuhat vectors to disk: ",i3)') i
               write(cens,fmt='(i4.4)') i
               filename = 'nuhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
               call da_cv_io ('w', filename, cv_size, nuhat(1:cv_size,i))
            end do


            ! (2) Construct low rank estimate of A^-1 * ghat with leading neign modes of A
            do i = 1, neign
               dot_cv = da_dot_cv(cv_size, nuhat(:,i), ghat, grid, mz, jp_start, jp_end)
               write(unit=stdout,fmt='("i, dot_cv = ",I3,D25.16)') i,dot_cv

               xhat = xhat + LAM(i) * nuhat(:,i) * dot_cv
            end do

            deallocate(nuhat)
         end if
      end if

      write(unit=stdout,fmt='(A,I5)') 'svd_amat_type = ', svd_amat_type
      write(unit=stdout,fmt='(A,I5,A,I5)') 'Number of eigenpairs used: ', neign, ' of ', nens


      deallocate(yhat)
      deallocate(omega)
      deallocate(qhat)

   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: svd_stage should be 0, 1, or 2, not ',svd_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_randomise_svd")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd


