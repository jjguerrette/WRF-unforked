subroutine da_randomise_svd(grid, config_flags,            &
                           cv_size, xbx, be, iv, cv, &
                           re, y, xhat, j_cost, &
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
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   integer          :: mz(13)
#else
   integer          :: mz(7)
#endif
#else
#ifdef CLOUD_CV
   integer          :: mz(15)
#else
   integer          :: mz(9)
#endif
#endif
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   real                         :: ghat(1:cv_size)        ! cv copy.
   real(wp)                     :: F1(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: F2(ensdim_svd, ensdim_svd) ! intermediate matrix
   real(wp)                     :: U(ensdim_svd, ensdim_svd)
   real(wp)                     :: VT(ensdim_svd, ensdim_svd)
   real(wp)                     :: SIG(ensdim_svd)
   real(wp), allocatable        :: TEMP1(:,:), TEMP2(:,:)


   real(wp)                     :: BMAT(ensdim_svd, ensdim_svd)  ! B = Q^T * A * Q
   integer                      :: piv(ensdim_svd)
   integer                      :: info

   real                         :: gdot, beta, temp_l
   real                         :: temp_e(ensdim_svd)
   real                         :: j0_total
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
#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
           be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /)
#endif
#else
#ifdef CLOUD_CV
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
           be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#else
   mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne, &
           be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
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
            do jens = 1, iens-1
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
      write(cnens,fmt='(I0)') nens
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'F1: ',F1(iens,1:nens)
      end do
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'F2: ',F2(iens,1:nens)
      end do

      info = 0
      p = svd_p
      if (p .eq. 0) then
         call GESV(F1, F2, piv, info) !F95

         ! GESV returns B^T in F2
         do iens = 1, nens
            do jens = 1, nens
               BMAT(iens,jens) = F2(jens,iens)
            end do
         end do
      else
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
                  BMAT(iens,jens) = sum(TEMP1(:,iens) * TEMP2(:,jens))
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
      do iens = 1, nens
         do jens = iens+1, nens
               BMAT(iens,jens) = (BMAT(iens,jens) + BMAT(jens,iens)) / 2.D0
               BMAT(jens,iens) = BMAT(iens,jens)
         end do
      end do


      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'BMAT: ',BMAT(iens,1:nens)
      end do

    ! Determine eigenvalues and eigenvectors of the symmetric B matrix
    !-----------------------------------------------------------------
!      ndgewrk           = MAX(4*nens-1,1)
      info              = 0
      call SYEV(BMAT, LAM, JOBZ = 'V' , INFO = info) !F95

      if (info /=0) then
         write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if

      do iens = 1, nens
         do jens = 1, nens 
            eignvec(iens,jens) = BMAT(iens,jens)
         end do
      end do
      
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
         eignval(iens) = LAM(iens) + 1
         if (neign .eq. nens .and. LAM(iens) .lt. 0.D0) then
            if (iens .gt. 1) then
               neign = iens - 1
            else
               write(unit=stdout,fmt='(A)') &
                  'FATAL ERROR: There need to be some positive eigenvalues of BMAT'
               CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
            end if
         end if
      end do

      ! Print current eigenpairs
      ! ----------------------------------------------
      write(cnens,fmt='(I0)') nens
      write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'Eigenvalues(r): ',LAM(1:nens)
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'D25.16)') 'EIGENVECTOR_MATRIX: ',eignvec(iens,1:nens)
      end do

      write(unit=stdout,fmt='(A)') ' '
      write(unit=stdout,fmt='(A,I5,A,I5,A)') 'Using ', neign, ' of ', nens, ' eigenpairs'
      write(unit=stdout,fmt='(A)') ' '


      ! Note: da_amat_mul does not properly handle imaginary eigenpairs,
      !       nor does it perform proper truncation when neign < nens.
      ! Although symmetric real matrices should not have complex eigenvalues or eigenvectors

      ! Eigenvalues of the inverse Hessian (including all background terms)

      if (neign .eq. nens .and. .not. adapt_svd) then
         do i = 1, nens
            LAMinv(i) = 1.D0 / ( LAM(i) + 1.D0 )
         end do
         write(unit=stdout,fmt='(A)') 'Carrying out full real amat_mul'
         call da_amat_mul(be,grid,cv_size,nens,neign,LAMinv,eignvec,qhat,ghat,xhat)
      else
         !  Truncate to neign modes Here we use full eigenvec matrix to generate nuhat=qhat*eigenvec,
         !  then truncate by only using neign of the nuhat vectors 

         if (minval(LAM(1:neign)) .lt. 1.D0 .and. adapt_svd) then
            ! (1) Generate eigenvectors of A
            allocate(nuhat(1:cv_size, 1:neign))
            nuhat = 0.D0

            do jens = 1, neign
               do i = 1, cv_size
                  nuhat(i,jens) = da_dot(nens, qhat(i,1:nens), eignvec(1:nens,jens))
               end do
            end do

            ! Write nuhat vectors to disk
            do i = 1, neign
               write(unit=stdout,fmt='("Writing nuhat vectors to disk: ",i3)') i
               write(cens,fmt='(i4.4)') i
               filename = 'nuhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
               call da_cv_io ('w', filename, cv_size, nuhat(1:cv_size,i))
            end do

            xhat = ghat

            ! (2) Construct low rank estimate of A^-1 * ghat with leading neign modes of A
            do i = 1, neign
               dot_cv = da_dot_cv(cv_size, nuhat(:,i), ghat, grid, mz, jp_start, jp_end)
               write(unit=stdout,fmt='("i, dot_cv = ",I3,D25.16)') i,dot_cv

!!               if (minval(LAM(1:neign)) .lt. 1.D0 .and. adapt_svd) then
                  xhat = xhat - LAM(i) / (LAM(i) + 1.D0) * nuhat(:,i) * dot_cv
!!               else
!!                  xhat = xhat + 1.D0 / (LAM(i) + 1.D0) * nuhat(:,i) * dot_cv
!!               end if
            end do

            deallocate(nuhat)
         else
            do i = 1, nens
               LAMinv(i) = 1.D0 / ( LAM(i) + 1.D0 )
            end do
            call da_amat_mul_trunc(be,grid,cv_size,nens,nens,neign,LAMinv,eignvec,qhat,ghat,xhat)
         end if

         write(unit=stdout,fmt='(A)') 'Carrying out truncated amat_mul'

      end if

      deallocate(yhat)
      deallocate(omega)
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


