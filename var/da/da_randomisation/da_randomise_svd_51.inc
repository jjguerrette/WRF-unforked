subroutine da_randomise_svd_51(grid, mz, config_flags,     &
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
   integer, intent(in)               :: cv_size                  ! Total cv size
   type (xbx_type),intent(inout)     :: xbx                      ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                       ! background error structure.
   type (iv_type), intent(inout)     :: iv                       ! ob. increment vector.
   real*8, intent(inout)             :: cv(1:cv_size)            ! control variable (local).
   type (y_type), intent(inout)      :: re                       ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                        ! y = H(x_inc) structure.
   real*8, intent(out)               :: xhat(1:cv_size)          ! control variable (local).
   integer, intent(in)               :: nens                     ! Total number of eigenpairs
   type (j_type), optional, intent(out) :: j_cost                ! cost function

   real*8, optional, intent(out)       :: qhat(1:cv_size, 0:nens) ! orthonormal basis set of A = Q * B * Q^T
   real*8, optional, intent(out)     :: eignvec(nens, nens)
   real*8, optional, intent(out)     :: eignval(nens)
   integer, optional, intent(out)    :: neign                         ! Number of non-truncated eigenpairs


#if defined(LAPACK)

   real(wp)                     :: SIG(1:nens)
   real(wp)                     :: LAM(1:nens)

   real(wp), allocatable   :: KMAT(:,:)  ! B = Q ^T * A  or  B^T = A^T * Q, where A = R^-1/2 * H * U
   real(wp), allocatable   :: KMAT_g(:,:)  ! B = Q ^T * A  or  B^T = A^T * Q, where A = R^-1/2 * H * U

   real(wp), allocatable   :: omega(:,:) ! Gaussian draw
   real(wp), allocatable   :: nuhat_T(:,:) ! Eigenvectors of A
   real(wp), allocatable   :: nuhat_T_g(:,:) ! Eigenvectors of A

   real(wp), allocatable   :: yhat_obs(:)  ! yhat = A * omega
   real(wp), allocatable   :: YHAT(:,:)    ! yhat = A * omega

   real(wp), allocatable   :: qhat_obs(:)  ! orthonormal basis set of A approx= Q * B
   real(wp), allocatable   :: qhat_ens(:,:)    ! orthonormal basis set of A approx= Q * B

   real(wp), allocatable       :: dumcv(:)
   real(wp), allocatable       :: dumcv_g(:)


   type (y_type)                :: jo_grad_y    ! Grad_y(jo)
   real*8                       :: grad_jo(cv_size)

   integer                      :: nall(num_ob_indexes)
   integer, allocatable         :: nobs(:)
   
   integer                      :: iens, jens, i, j, ccount, itype
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
   real                         :: ghat(1:cv_size)        ! cv copy.
   integer                      :: info

   real*8                       :: gdot, beta
   real*8                       :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=3) :: cnens=""

   if (trace_use) call da_trace_entry("da_randomise_svd_51")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using randomized SVD method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------

   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   xhat = 0.D0

   if ( svd_stage .eq. 1 .and. ensmember .gt. 0) then
      
      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember

      if ( ensmember .eq. nens+1 ) then 
         ! Calculate grad J for xhat = 0
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

         if (riot_precon.eq.1 .and. svd_outer.gt.1 .and. sum(ntmax(1:svd_outer-1)).gt.nens) &
            call da_precon_omega(grid, cv_size, mz, jp_start, jp_end, nens, ghat)

      else if (ensmember .gt. 0 .and. ensmember .le. nens) then
         ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = R^-1/2 * H * U 
         allocate(omega(1:cv_size,1))

!         omega(:,1) = 0.D0

         if (read_omega) then
            ! Read in omega vector from this ensemble on each processor
            write(unit=stdout,fmt='("Reading omega ensemble member : ",i3)') ensmember
            filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('r', filename, cv_size, omega(1:cv_size,1), grid, mz)
         else
            write(unit=stdout,fmt='("Generating omega ensemble member : ",i3)') ensmember
            call da_set_randomcv(cv_size, omega(:,1))
         end if

         call da_transform_vtoy(cv_size, be, grid%ep, omega(1:cv_size,1), iv, grid%vp, &
                 grid%vv, grid%vp6, grid%vv6, xbx, y, grid, config_flags &
#if (WRF_CHEM == 1)
                 , grid%vch &
#endif
                 )

         call da_transform_ytoyhat (iv, y, yhat_obs, nall)

         if (allocated(yhat_obs)) then
            if (size(yhat_obs) .ne. sum(nall) .or. sum(nall).le.0) then
               write(unit=stdout,fmt='(A)') 'FATAL ERROR: Incorrect size yhat_obs'
               CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
            end if
         end if

         allocate(nobs(1))
         nobs(1) = sum(nall)
         filename = 'yhat_obs.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_yhat_io('w', filename, nall, nobs(1), yhat_obs)
         deallocate(nobs)

         if (allocated(yhat_obs)) deallocate(yhat_obs)

!         !Write omega to file (only necessary if same random omega is desired in each outer loop)
!         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1), grid, mz)
         deallocate(omega)

      else
         write(unit=stdout,fmt='(A,I4)') 'FATAL ERROR: ensmember should be between 1 and nens+1, not ',ensmember
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')

      end if

   else if ( svd_stage .eq. 2 ) then

      allocate(nobs(nens))
      write(cproc,fmt='(i4.4)') myproc

      ! Read in yhat vector from all ensembles on each processor
      do iens = 1, nens
         write(unit=stdout,fmt='("Reading yhat ensemble member : ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'yhat_obs.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))

         call da_yhat_io ('r', filename, nall, nobs(iens), yhat_obs)

         if (iens .eq. 1) then
            if (nobs(iens) .gt. 0) then
               allocate(YHAT(1:nobs(1),nens))
            else
!               deallocate(yhat_obs)
               exit
            end if
         else
            if (nobs(iens) .ne. nobs(iens-1)) then
               write(unit=stdout,fmt='(A,2I6)') &
                  'FATAL ERROR: nobs does not match across ensembles',nobs(iens),nobs(iens-1)
               CALL wrf_error_fatal('Canceling RANDOMIZED SVD')

               !May need to add a test for nall also...
            end if
         end if

         YHAT(:,iens) = yhat_obs

         deallocate(yhat_obs)

      end do

      call mpi_barrier(comm, ierr)

!      if (nobs(1) .gt. 0) then
         allocate(qhat_ens(1:nobs(1),nens))
!      end if

      ! Generate qhat basis vectors by Modified Gram-Schmidt orthonormalization
      do iens = 1, nens
         if (nobs(1) .gt. 0) qhat_ens(:,iens) = YHAT(:,iens)
         if (iens .gt. 1) then
            do i = 1,2 !Second iteration improves orthogonalization for numerical precision (Bjorck, 1994)
            do jens = 1, iens-1
               gdot = da_dot_obs(nobs(1), qhat_ens(:,iens), qhat_ens(:,jens))
               
               if (nobs(1) .gt. 0) qhat_ens(:,iens) = qhat_ens(:,iens) - gdot * qhat_ens(:,jens)
            end do
            end do
         end if

         beta   = SQRT(da_dot_obs(nobs(1), qhat_ens(:,iens), qhat_ens(:,iens)))
         if (nobs(1) .gt. 0) qhat_ens(:,iens) = qhat_ens(:,iens) / beta
      end do

      allocate(qhat_obs(1:nobs(1)))
      ! Write qhat vectors to disk
      do iens = 1, nens
         write(unit=stdout,fmt='("Writing qhat vectors to disk: ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'qhat_obs.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))

         qhat_obs = qhat_ens(:,iens)
         call da_yhat_io('w', filename, nall, nobs(1), qhat_obs)
      end do
      deallocate(qhat_obs)

      if (allocated(YHAT)) deallocate(YHAT)
      if (allocated(qhat_ens)) deallocate(qhat_ens)
      deallocate(nobs)
   else if ( svd_stage .eq. 3 ) then

      call da_allocate_y(iv, jo_grad_y)
#if (WRF_CHEM == 1)
      call da_allocate_y_chem(iv, jo_grad_y)
#endif

      allocate(nobs(1))
      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember
      filename = 'qhat_obs.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))

      call da_yhat_io('r', filename, nall, nobs(1), qhat_obs)

      call da_transform_ytoyhat_adj (iv, qhat_obs, nall, jo_grad_y)

      deallocate(nobs)
      
      if (allocated(qhat_obs)) deallocate(qhat_obs)

      call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
              grid%vp, grid%vv, grid%vp6, grid%vv6, &
              xbx, jo_grad_y, grid, config_flags, .true. &
#if (WRF_CHEM == 1)
              , grid%vch &
#endif
              )
      grad_jo = - grad_jo    !! Compensate for sign in calculation of grad_v (Jo)

      ! Write KMAT row to file
      filename = 'bhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
      call da_cv_io ('w', filename, cv_size, grad_jo, grid, mz)

      call da_deallocate_y(jo_grad_y)

   else if ( svd_stage .eq. 4 ) then
      ! Eigen decomposition of A (Algorithm 5.1 of Halko)

      allocate(KMAT(nens,1:cv_size))
      allocate(nuhat_T(nens,1:cv_size))
      nuhat_T = 0.D0

      write(cproc,fmt='(i4.4)') myproc

      ! Read in KMAT for all rows on each processor, KMAT = (A^T * Q)^T
      do iens = 1, nens
         write(unit=stdout,fmt='("Reading KMAT row : ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'bhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, KMAT(iens,1:cv_size), grid, mz)
      end do

      ! Read gradient
      filename = 'ghat.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat, grid, mz)

#ifdef DM_PARALLEL
      !Perform global SVD of KMAT matrix
      cv_size_domain = wrf_dm_sum_integer(cv_size)
      if (.not.rootproc) cv_size_domain = 1

!      if (rootproc) then
         allocate(KMAT_g(nens,1:cv_size_domain))
         allocate(nuhat_T_g(nens,1:cv_size_domain))
!      else
!         allocate(KMAT_g(1:nens,1))
!         allocate(nuhat_T_g(1:nens,1))
!      end if
      allocate(dumcv(1:cv_size))
      allocate(dumcv_g(1:cv_size_domain))
      do iens = 1, nens
         dumcv = KMAT(iens,1:cv_size)

!REMOVE EVENTUALLY or add option
!         write(cens,fmt='(i4.4)') iens
!         filename = 'dumcv.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', filename, cv_size, dumcv, grid, mz)
!REMOVE EVENTUALLY or add option

         call da_cv_to_global(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)

         if(rootproc) KMAT_g(iens,:) = dumcv_g

         call mpi_barrier(comm, ierr)
      end do
      deallocate(dumcv)
      deallocate(dumcv_g)

      write(unit=stdout,fmt='(A)') 'Calculate SVD of B = Q^T * R^-1/2 * H * U'
      if (rootproc) then
         !SVD of global KMAT on monitor process
         call GESVD(KMAT_g, SIG, VT = nuhat_T_g, INFO = info)
      else
         SIG = 0.D0
      end if
      deallocate(KMAT_g)

      call mpi_barrier(comm, ierr)

      call wrf_dm_bcast_real(SIG,nens)

      call wrf_dm_bcast_integer(info,1)

      allocate(dumcv(1:cv_size))
      allocate(dumcv_g(1:cv_size_domain))
      do iens = 1, nens
         if(rootproc) then
            dumcv_g = nuhat_T_g(iens,:)
         else
            dumcv_g = 0.D0
         end if
         dumcv = 0.D0

         call da_global_to_cv(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)
         nuhat_T(iens,:) = dumcv

         call mpi_barrier(comm, ierr)
      end do
      deallocate(dumcv)
      deallocate(dumcv_g)
      deallocate(nuhat_T_g)
#else
      !SVD of global KMAT
      call GESVD(KMAT, SIG, VT = nuhat_T, INFO = info)
#endif

      if (info /=0) then
         write(stdout,*) 'Error in randomised SVD: GESVD returned ',info
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if

      !Eigenvalues of A^T * A = U^T H^T R^-1 H U are the square of singular values of A
      do i = 1, nens
         LAM(i) = SIG(i) * SIG(i)
      end do

      neign = nens

      write(cnens,fmt='(I0)') neign
      write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'Eigenvalues(r): ',LAM(1:nens)

!REMOVE EVENTUALLY or add option
!      ! Write nuhat vectors to disk
!      do i = 1, neign
!         write(unit=stdout,fmt='("Writing nuhat vectors to disk: ",i3)') i
!         write(cens,fmt='(i4.4)') i
!         filename = 'nuhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', filename, cv_size, nuhat_T(i,1:cv_size), grid, mz)
!      end do
!REMOVE EVENTUALLY or add option

      ! Store output vectors
      eignvec(:,:) = 0.D0
      do i = 1, neign
         eignvec(i,i) = 1.D0
         eignval(i) = LAM(i) + 1.D0
      end do
      qhat(1:cv_size,1:nens) = TRANSPOSE(nuhat_T)

      ! (2) Construct low rank estimate of A^-1 * ghat with leading neign modes of A
      if (minval(eignval(1:neign)) .lt. 2.D0 .and. adapt_svd) then
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
         itype = 2

         write(unit=stdout,fmt='(A)') 'Carrying out RIOT-51 LRU'
      else
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
         itype = 1

         write(unit=stdout,fmt='(A)') 'Carrying out RIOT-51 LRA'
      end if

      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             nens,nens,neign, &
                             LAM,eignvec,qhat, &
                             ghat,xhat,.true.)

      if (itype .eq. 2) then !LRU
         xhat = xhat + ghat
      end if

      if (riot_precon.ge.1) &
         call da_hessian_io ('w',grid, cv_size,svd_outer,nens,neign,&
                            eignvec,eignval,qhat,mz,ghat)

      if (riot_precon.eq.1 .and. sum(ntmax(1:svd_outer)).lt.ntmax(svd_outer+1)) riot_precon=4

      if (any(riot_precon.eq.(/2,3,4/))) &
         call da_precon_omega(grid, cv_size, mz, jp_start, jp_end, ntmax(svd_outer+1), ghat)

      deallocate(KMAT)
      deallocate(nuhat_T)

   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: For RSVD5.1, svd_stage should be 0, 1, 2, 3, or 4, not ',svd_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_randomise_svd_51")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_randomise_svd_51


