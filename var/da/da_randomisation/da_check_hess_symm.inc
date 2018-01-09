subroutine da_check_hess_symm( grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                      config_flags, iv, re, y, &
                                      rankA, eignval_hess, eignvec_hess, nu_hess, AMAT)

   !-------------------------------------------------------------------------
   ! Purpose:         Check Hessian Symmetry
   !
   ! Here 
   !    cv   is updated in outer-loop.
   !    xhat is the control variable in inner-loop.
   !
   ! Called from da_calculate_hessian.inc
   !
   ! History: 06/27/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   type (xbx_type),intent(inout) :: xbx     ! Header & non-gridded vars.
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: rankA
   real*8,         intent(in) :: eignvec_hess(rankA, rankA)
   real*8,         intent(in) :: eignval_hess(rankA)
   real*8,         intent(in) :: nu_hess(cv_size, 0:rankA) ! Ritz vectors

   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)
   real*8,         intent(in) :: AMAT(cv_size, 1:rankA) ! Half Hessian


   type(grid_config_rec_type), intent(inout) :: config_flags
   type (be_type), intent(in)    :: be      ! background error structure.
   type (iv_type), intent(inout)     :: iv                       ! ob. increment vector.
   type (y_type), intent(inout)      :: re                       ! residual (o-a) structure.
   type (y_type), intent(inout)      :: y                        ! y = H(x_inc) structure.


#if defined(LAPACK)
   real(wp), allocatable   :: omega(:) ! Gaussian draw
   real(wp), allocatable   :: omega_g(:) ! Eigenvectors of A

   real(wp), allocatable   :: yhat_obs(:)  ! yhat = A * omega
   real(wp), allocatable   :: yhat_obs_g(:)  ! yhat = A * omega

   real(wp), allocatable       :: dumcv(:)
   real(wp), allocatable       :: dumcv_g(:)
   real(wp)                    :: ACOL(rankA), temp

   integer                      :: nall(num_ob_indexes)
   integer, allocatable         :: nobs(:)
   
   integer                      :: iens, jens, i, j, imax(nens_compare), svd_unit, tempi
   logical, allocatable         :: mask(:)

   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   if (trace_use) call da_trace_entry("da_check_hess_symm")
   write(unit=stdout,fmt=*) ' '
   write(unit=stdout,fmt='(A)') '**Checking Hessian symmetry using TL model and left-half decomposition**'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------

   ! Generate vectors in Y = A*OMEGA (Halko, Algorithm 4.1), where A = R^-1/2 * H * U 
   allocate(omega(1:cv_size))

#ifdef DM_PARALLEL
   cv_size_domain = wrf_dm_sum_integer(cv_size)

   if (.not.rootproc) cv_size_domain = 1
#else
   cv_size_domain = cv_size
#endif

   allocate(mask(1:cv_size_domain))
   mask(:) = .true.

   imax(:) = 0
   do iens = 1, nens_compare
      omega(:) = 0.D0
#ifdef DM_PARALLEL
      allocate(dumcv(1:cv_size))
      allocate(dumcv_g(1:cv_size_domain))
      dumcv_g = 0.D0
      dumcv = nu_hess(1:cv_size,iens)
      call da_cv_to_global(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)

      allocate(omega_g(1:cv_size_domain))
      omega_g = 0.D0
      if(rootproc) then
         imax(iens:iens) = maxloc(abs(dumcv_g),mask)
!         imax(iens:iens) = maxloc(dumcv_g,mask)
!         imax(iens:iens) = minloc(dumcv_g,mask)

         mask(imax(iens)) = .false.

         omega_g(imax(iens)) = 1.D0
         tempi = imax(iens)
      end if
      call mpi_barrier(comm, ierr)

      call wrf_dm_bcast_integer(tempi,1)
      imax(iens) = tempi

      dumcv = 0.D0
      call da_global_to_cv(cv_size, cv_size_domain, dumcv, grid, mz, omega_g)
      omega = dumcv
      call mpi_barrier(comm, ierr)

      deallocate(dumcv)
      deallocate(dumcv_g)
      deallocate(omega_g)
#else
      imax(iens:iens) = maxloc(nu_hess(1:cv_size,iens),mask)
      mask(imax(iens)) = .false.

      omega(imax(iens)) = 1.D0
#endif

      call da_transform_vtoy(cv_size, be, grid%ep, omega, iv, grid%vp, &
              grid%vv, grid%vp6, grid%vv6, xbx, y, grid, config_flags &
#if (WRF_CHEM == 1)
              , grid%vch &
#endif
              )

      !da_transform_ytoyhat needs non-chem observations before it and da_check_hess_symm are generalizable
      call da_transform_ytoyhat (iv, y, yhat_obs, nall)

      if (allocated(yhat_obs)) then
         if (size(yhat_obs) .ne. sum(nall) .or. sum(nall).le.0) then
            write(unit=stdout,fmt='(A)') 'FATAL ERROR: Incorrect size yhat_obs'
            CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
         end if
      end if

      write(cproc,fmt='(i4.4)') myproc
      write(cens,fmt='(i8.8)') imax(iens)

      allocate(nobs(1))
      nobs(1) = sum(nall)
      filename = 'right_hess_sqrt_column.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
      call da_yhat_io('w', filename, nall, nobs(1), yhat_obs)
      deallocate(nobs)

      if (allocated(yhat_obs)) deallocate(yhat_obs)

   !         !Write omega to file (only necessary if same random omega is desired in each outer loop)
   !         filename = 'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
   !         call da_cv_io ('w', filename, cv_size, omega(1:cv_size), grid, mz)

#ifdef DM_PARALLEL
      if (rootproc) then
         call da_get_unit (svd_unit)
         filename = 'left_hess_sqrt_row.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
         open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
      end if

      do jens = 1, rankA
         allocate(dumcv(1:cv_size))
         allocate(dumcv_g(1:cv_size_domain))
         dumcv_g = 0.D0
         dumcv = AMAT(1:cv_size,jens)
         call da_cv_to_global(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)

         if(rootproc) then
            temp = dumcv_g(imax(iens))
         else
            temp = 0.D0
         end if
         call mpi_barrier(comm, ierr)
         call wrf_dm_bcast_real(temp,1)

         ACOL(jens) = temp

         deallocate(dumcv)
         deallocate(dumcv_g)
#else
         ACOL(jens) = AMAT(imax(iens),jens)
#endif

         if (rootproc) write(unit=svd_unit,fmt='(E25.16E3,I12)') ACOL(jens), imax(iens)
      end do

!!      write(cens,fmt='(I0)') rankA
!!      write(unit=svd_unit,fmt='(A,I8,A'//trim(cens)//'E25.16E3)') 'AMAT(:,', imax(iens), ') = ',ACOL(1:rankA)

      if (rootproc) then
         close(unit=svd_unit)
         call da_free_unit (svd_unit)
      end if

      !NOW DO SOMETHING WITH ACOL...compare to yhat_obs somehow (condense yhat_obs to 241 instead of 1196...not obvious which ones should be removed)

      ! OR make sure to output i,j,k coordinate of imax(iens) for comparison to AMAT in MATLAB

      ! MAYBE JUST REMOVE ZEROS (or < tolerance) FROM BOTH SINCE THEY ARE IN THE SAME ORDER...according to model time for acft obs

   end do

   deallocate(omega)
   deallocate(mask)

   if (trace_use) call da_trace_exit("da_check_hess_symm")

#else
   call da_error(__FILE__,__LINE__, &
       (/"Must compile with $LAPACK option for randomised svd"/))
#endif

end subroutine da_check_hess_symm


