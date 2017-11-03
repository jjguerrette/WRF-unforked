subroutine da_output_increment_by_modes(grid,config_flags,mz, &
                                        it,cv_size,xbx,be, &
                                        hessian_info)

   !-------------------------------------------------------------------------
   ! Purpose:   Output increment for all inner loop truncations using Hessian
   !            spectra and qhat basis vectors.  This subroutine is for
   !            diagnostic purposes and should not be used when computational
   !            expense need be minimized.
   !
   ! Called from da_main
   !
   ! History: 10/31/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout) :: grid
   type (grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: mz(:)

   integer,         intent(in)   :: it    ! outer-loop index
   integer,        intent(in)    :: cv_size ! Size of cv array.
   type(xbx_type), intent(in)    :: xbx  ! For header & non-grid arrays.
   type(be_type),  intent(in)    :: be   ! background errors.
   type(hessian_type), intent(in) :: hessian_info

!   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.
!   integer                      :: jl_start, jl_end       ! Start/end indices of Jl.

   integer :: i, imode, nmodes, neign, itype, iend, icount
   real*8, allocatable :: LAM(:)
   real*8, allocatable :: eignval(:)
   real*8, allocatable :: eignvec(:,:)
   real*8, allocatable :: qhat(:,:)
   real*8 :: ghat(1:cv_size)
   real*8 :: xhat(1:cv_size)
   real*8 :: xhat_REF(1:cv_size)
   real*8 :: dvnorm_(1:cv_size)
   real*8 :: EUCNORM_xPa(2)
!   real*8 :: minlam(1:1)
   type(hessian_eig_type) :: hess_REF

   character(len=10) :: PREFIX
   character(len=256) :: filnam
   character(len=256) :: timestr
   character(len=256) :: dummy_name
  

   if (trace_use) call da_trace_entry("da_output_increment_by_modes")

!   je_start   = be % cv % size_jb + 1
!   je_end     = be % cv % size_jb + be % cv % size_je
   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp
!   jl_start   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + 1
!   jl_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp + be % cv % size_jl

   nmodes = hessian_info%nmodes_max

   ghat = hessian_info%ghat

   allocate(qhat(cv_size,0:nmodes))
   allocate(LAM(1:nmodes))

   dummy_name = trim(config_flags%auxhist18_outname)
   call nl_set_io_form_auxhist18 ( grid%id, 2 )
   config_flags%io_form_auxhist18 = 2
   config_flags%frames_per_auxhist18 = 1

   iend = 2
   if (use_lanczos) iend = 1

   do imode = 1, nmodes
      qhat(:,imode) = hessian_info%qhat(:,imode)
   end do

   icount = 0

   do imode = nmodes, 1, -inc_out_interval
      allocate(eignval(1:imode))
      allocate(eignvec(1:imode,1:imode))

      neign = hessian_info%eigs(imode)%neign
      eignval = hessian_info%eigs(imode)%eignval
      eignvec = hessian_info%eigs(imode)%eignvec

      EUCNORM_xPa = 0.D0
   
      do itype = 1, iend
         icount = icount+ 1
!!!      if ((minval(eignval(1:np)) .lt. 2.D0 .and. adapt_svd.eq.1) .or. adapt_svd.gt.1) then
         if (itype .eq. 2) then
            do i = 1, neign
               LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
            end do

            PREFIX = 'LRU'

         else
            do i = 1, neign
               LAM(i) = 1.D0 / eignval(i)
            end do

            PREFIX = 'LRA'
            if (use_lanczos) PREFIX = 'CGLZ'

         end if
         write(unit=filnam,fmt='(2a,i2.2,a,i4.4)') trim(PREFIX),'inc_it',it,'_nmode',imode

         xhat = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                imode,imode,neign, &
                                LAM(1:imode),eignvec,qhat(:,0:imode), &
                                ghat,xhat)

         if (itype .eq. 2) then !LRU
            xhat = xhat + ghat
         end if

         call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                 xhat(1:be%cv%size_jb),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                                 , grid%vch &
#endif
                                 )

         call domain_clock_get( grid, current_timestr=timestr )
         config_flags%auxhist18_outname = trim(filnam)
         call da_transfer_xatowrftl(grid, config_flags, 'DUMM', timestr, .true.)

         !-----------------------------------------
         ! Calculate ||dx - dx_REF||^2_{Pa_REF^-1}
         !-----------------------------------------
         !Define reference case (alternatively, could read in from file)
         if (icount .eq. 1) then
            xhat_REF = xhat

            allocate(hess_REF%eignval(1:nmodes))
            allocate(hess_REF%eignvec(1:nmodes,1:nmodes))
            hess_REF%neign = neign
            hess_REF%eignval = eignval
            hess_REF%eignvec = eignvec
         end if

         !diffvec
         xhat = xhat - xhat_REF

         do i = 1, hess_REF%neign
            LAM(i) = hess_REF%eignval(i) - 1.D0
         end do

         !Multiply Pav_REF^-1 * (dv - dv_REF) in two steps
         !(1) (L^T H^T R^-1 H L)_REF * diffvec
         dvnorm_ = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                                nmodes,nmodes,hess_REF%neign, &
                                LAM,hess_REF%eignvec,qhat,&
                                xhat,dvnorm_)

         !(2) (I + L^T H^T R^-1 H L) * diffvec
         dvnorm_ = dvnorm_ + xhat

         !Multiply (dv - dv_REF)^T * Pav_REF^-1 * (dv - dv_REF) = ||dx - dx_REF||^2_{Pa_REF^-1}
         EUCNORM_xPa(itype) = da_dot_cv(cv_size, xhat, dvnorm_, grid, mz, jp_start, jp_end)

      end do

!      minlam = minval(eignval)

      write(unit=stdout,fmt='(A,2I6,3E27.16E3)') &
         'imode,neign,min(lambda),EUCNORM_xPa[LRA,LRU] = ',&
         imode,neign,minval(eignval),EUCNORM_xPa

      deallocate(eignval)
      deallocate(eignvec)

   end do

   deallocate(LAM)
   deallocate(hess_REF%eignval)
   deallocate(hess_REF%eignvec)
   deallocate(qhat)

   config_flags%auxhist18_outname = trim(dummy_name)
   call nl_set_io_form_auxhist18 ( grid%id, 0 )
   config_flags%io_form_auxhist18 = 0

   if (trace_use) call da_trace_exit("da_output_increment_by_modes")

end subroutine da_output_increment_by_modes
