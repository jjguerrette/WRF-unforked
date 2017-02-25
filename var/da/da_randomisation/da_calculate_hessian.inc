subroutine da_calculate_hessian(grid, config_flags, &
                                cv_size, xbx, be, iv, cv, &
                                re, y, xhat, j_cost, &
                                qhat, eignval, eignvec, neign)

   !-------------------------------------------------------------------------
   ! Purpose: Direct Hessian calculation, separated into 2 stages 
   !          + Comparisons of alternative Hessian decomps and increments
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

   real(wp)                     :: SIG(1:ensdim_svd)
   real(wp)                     :: LAM(1:ensdim_svd)

   real(wp), allocatable   :: AMAT(:,:)  ! B = Q ^T * A  or  B^T = A^T * Q, where A = R^-1/2 * H * U
   real(wp), allocatable   :: AMAT_g(:,:)  ! B = Q ^T * A  or  B^T = A^T * Q, where A = R^-1/2 * H * U

   real(wp), allocatable   :: omega(:,:) ! Gaussian draw
   real(wp), allocatable   :: nuhat(:,:) ! Eigenvectors of A
   real(wp), allocatable   :: nuhat_g(:,:) ! Eigenvectors of A

   real(wp), allocatable       :: dumcv(:)
   real(wp), allocatable       :: dumcv_g(:)

   type (y_type)                :: jo_grad_y    ! Grad_y(jo)
   real                         :: grad_jo(cv_size)

   integer                      :: nall(num_ob_indexes)
   integer, allocatable         :: nobs(:)
   
   integer                      :: iens, jens, nyobs, i, j, ccount
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
   integer                      :: info

   real*8                       :: beta
   real*8                       :: j0_total
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename
   character(len=3) :: cnyobs=""
  
   if (trace_use) call da_trace_entry("da_calculate_hessian")

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

   nyobs = ensdim_svd

   xhat = 0.D0

   if ( svd_stage .eq. 1 .and. ensmember .gt. 0) then
      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i4.4)') ensmember


      if ( ensmember .eq. nyobs+1 ) then
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
         call da_cv_io ('w', filename, cv_size, ghat)

      else if (ensmember .gt. 0 .and. ensmember .le. nyobs) then
      
         call da_allocate_y(iv, jo_grad_y)
#if (WRF_CHEM == 1)
         call da_allocate_y_chem(iv, jo_grad_y)
#endif
         call da_force_grad_hess (iv, y, nall, ensmember, jo_grad_y)

         call da_transform_vtoy_adj(cv_size, be, grid%ep, grad_jo, iv, &
                 grid%vp, grid%vv, grid%vp6, grid%vv6, &
                 xbx, jo_grad_y, grid, config_flags, .true. &
#if (WRF_CHEM == 1)
                 , grid%vch &
#endif
                 )

   ! Should not be necessary since forcing is positive
   !      grad_jo = - grad_jo    !! Compensate for sign in calculation of grad_v (Jo)

         ! Write AMAT row to file
         filename = 'ahat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, grad_jo)

         call da_deallocate_y(jo_grad_y)
      end if
   else if ( svd_stage .eq. 2 ) then
      ! Eigen decomposition of A

      allocate(AMAT(1:cv_size,nyobs))
      allocate(nuhat(1:cv_size,nyobs))
      nuhat = 0.D0

      write(cproc,fmt='(i4.4)') myproc

      ! Read in AMAT for all rows on each processor, AMAT = (A^T * I)
      do iens = 1, nyobs
         write(unit=stdout,fmt='("Reading AMAT row : ",i3)') iens
         write(cens,fmt='(i4.4)') iens
         filename = 'ahat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', filename, cv_size, AMAT(1:cv_size,iens))
      end do

      ! Read gradient
      filename = 'ghat.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, cv_size, ghat)

#ifdef DM_PARALLEL
      !Perform global SVD of AMAT matrix
      cv_size_domain = wrf_dm_sum_integer(cv_size)
      if (rootproc) then
         allocate(AMAT_g(1:cv_size_domain,nyobs))
         allocate(nuhat_g(1:cv_size_domain,nyobs))
      else
         allocate(AMAT_g(1,1:nyobs))
         allocate(nuhat_g(1,1:nyobs))
      end if
      allocate(dumcv(1:cv_size))
      allocate(dumcv_g(1:cv_size_domain))
      do iens = 1, nyobs
!         dumcv = AMAT(iens,1:cv_size)
         dumcv = AMAT(1:cv_size,iens)


         write(cens,fmt='(i4.4)') iens
         filename = 'dumcv.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, dumcv)

         call da_cv_to_global(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)

!         if(rootproc) AMAT_g(iens,:) = dumcv_g
         if(rootproc) AMAT_g(:,iens) = dumcv_g

         call mpi_barrier(comm, ierr)
      end do
      deallocate(dumcv)
      deallocate(dumcv_g)

      write(unit=stdout,fmt='(A)') 'Calculate SVD of B = Q^T * R^-1/2 * H * U'
      if (rootproc) then
         !SVD of global AMAT on monitor process
         call GESVD(AMAT_g, SIG, U = nuhat_g, INFO = info)
!         call GESVD(AMAT_g, SIG, VT = nuhat_g, INFO = info)

      else
         SIG = 0.D0
      end if
      deallocate(AMAT_g)

      call mpi_barrier(comm, ierr)

      call wrf_dm_bcast_real(SIG,nyobs)

      call wrf_dm_bcast_integer(info,1)

      allocate(dumcv(1:cv_size))
      allocate(dumcv_g(1:cv_size_domain))
      do iens = 1, nyobs
         if(rootproc) then
            dumcv_g = nuhat_g(:,iens)
         else
            dumcv_g = 0.D0
         end if
         dumcv = 0.D0

         call da_global_to_cv(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)
         nuhat(:,iens) = dumcv

         call mpi_barrier(comm, ierr)
      end do
      deallocate(dumcv)
      deallocate(dumcv_g)
      deallocate(nuhat_g)
#else
      !SVD of global AMAT
      call GESVD(AMAT, SIG, U = nuhat, INFO = info)
!      call GESVD(AMAT, SIG, VT = nuhat, INFO = info)
#endif

      if (info /=0) then
         write(stdout,*) 'Error in randomised SVD: GESVD returned ',info
         CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
      end if

      !Eigenvalues of A^T * A = U^T H^T R^-1 H U are the square of singular values of A
      do i = 1, nyobs
         LAM(i) = SIG(i) * SIG(i)
      end do

      neign = nyobs

      ! Write nuhat vectors to disk
      do i = 1, neign
         write(unit=stdout,fmt='("Writing nuhat vectors to disk: ",i3)') i
         write(cens,fmt='(i4.4)') i
         filename = 'nuhat.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', filename, cv_size, nuhat(1:cv_size,i))
      end do

      write(cnyobs,fmt='(I0)') neign
      write(unit=stdout,fmt='(A,'//trim(cnyobs)//'E25.16E3)') 'Eigenvalues(r): ',LAM(1:ensdim_svd)

      ! Store output vectors
      eignvec(:,:) = 0.D0
      do i = 1, neign
         eignvec(i,i) = 1.D0
         eignval(i) = LAM(i) + 1.D0
      end do
      qhat(:,0) = 0.D0
      qhat(1:cv_size,1:nyobs) = nuhat

      ! (2) Construct low rank estimate of A^-1 * ghat with leading neign modes of A
      if (minval(eignval(1:neign)) .lt. 2.D0 .and. adapt_svd) then
         do i = 1, neign
            LAM(i) = -(eignval(i) - 1.D0) / eignval(i)
         end do
         xhat = ghat

         write(unit=stdout,fmt='(A)') 'Carrying out Full-SVD LRU'
      else
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do

         write(unit=stdout,fmt='(A)') 'Carrying out Full-SVD LRA'
      end if

      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             neign,neign,neign,&
                             LAM,eignvec,qhat,&
                             ghat,xhat)

      deallocate(AMAT)
      deallocate(nuhat)

      if (nens_compare .gt. 0) then
         call da_compare_decomp_methods(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                          neign, eignval, eignvec, qhat, &
                                          ghat, xhat)
      end if

   else
      write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: svd_stage should be 0, 1, or 2, not ',svd_stage
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   if (trace_use) call da_trace_exit("da_calculate_hessian")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for full Hessian Calculation"/))
#endif

end subroutine da_calculate_hessian


