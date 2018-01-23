subroutine da_evaluate_j(grid, config_flags, cv_size, cv_size_jb, &
                         cv_s, cv_e, eta, jout, &
                         xbx, be, iv, xhat, cv, ob, &
                         it, mz, iter)

   !---------------------------------------------------------------------------
   ! Purpose: Relinearize cost function around Lb * (cv + eta*dv)
   !          cv, previous accumulated cv increments
   !          xhat, current cv increment
   !          eta, scalar to multiply preconditioned cv members
   !          +Applies to MET and CHEM
   !          +Only calculates Jo and Jb, no others
   !
   !          WARNING: this J evaluation will cause errors for var4d_lbc
   !
   ! Called from da_dgn
   !---------------------------------------------------------------------------

   implicit none

   type (domain),       intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)                :: cv_size    ! Total cv size.
   integer, intent(in)                :: cv_size_jb ! Jb cv size.
   type (xbx_type),intent(inout)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)         :: be     ! background error structure.
   type (iv_type), intent(inout)      :: iv     ! innovation vector (o-b).
   real, intent(in)                   :: xhat(1:cv_size) ! control variables.
   real, intent(in)                   :: cv(1:cv_size)   ! control variables.
   type (y_type) , intent(inout)      :: ob          ! Observation structure.
   integer, intent(in)                :: it, iter 
   integer, intent(in)                :: mz(:)
!   type (xb_type)      :: shuffleb

   ! Need to allow ncv_type to be +1 when var4d_lbc
#if (WRF_CHEM == 1)
!   type (xch_type)      :: shuffleb_ch
   integer, parameter    :: ncv_type = 2
#else
   integer, parameter    :: ncv_type = 1
#endif
   integer, intent(in), dimension(1:ncv_type) :: cv_s, cv_e
   real, intent(in), dimension(1:ncv_type)    :: eta
   real, intent(out), dimension(1:ncv_type)   :: jout
   integer          :: num_qcstat_conv(2,num_ob_indexes,num_ob_vars,npres_print+1)
   real             :: jo_partial                   ! jo for this processor
   real             :: joc_partial                  ! jo for cross-validation for this processor
   real             :: cv_xhat_jb(cv_size_jb)
   real             :: xhat_temp(1:cv_size) ! initial control variables.
   type (j_type)    :: j
   type (y_type)    :: jo_grad_y 
   type (y_type)    :: re     ! residual vector (o-a).
   type (y_type)    :: y      ! y = H(x_inc) = 0.; just a placeholder here
   integer          :: i
   integer, dimension(1:ncv_type) :: ncv
   integer          :: rc
   character(len=filename_len) :: fg_filename
   logical :: original_disable
   integer :: jl_start, jl_end
   character*256 :: timestr, timestr1
   type(x_type) :: shuffle

   if (trace_use) call da_trace_entry("da_evaluate_j")

   !-------------------------------------------------------------------------
   ! [0.0] initialization:
   !-------------------------------------------------------------------------

   if (it .eq. 1) then
      fg_filename = 'fg'
   else
      write(fg_filename,fmt='("wrfvar_output_",I2.2)') it-1
   end if

   xhat_temp = xhat

   !Separable line search when chem is uncoupled, need to add coupling options for DA
   do i = 1, ncv_type
      ncv(i) = cv_e(i) - cv_s(i) + 1
      if ( ncv(i) .gt. 0 ) xhat_temp(cv_s(i):cv_e(i)) = eta(i) * xhat(cv_s(i):cv_e(i))
   end do

   ! Store xb, which get_innov_vector modifies below
!   shuffleb = grid%xb
!#if (WRF_CHEM == 1)
!   shuffleb_ch = grid%xbch
!#endif

   evalj = .true.

   !-------------------------------------------------------------------------
   ! [1.0] calculate jo:
   !-------------------------------------------------------------------------

   call da_transform_vtox (grid,cv_size_jb,xbx,be,grid%ep,&
                            xhat_temp(1:cv_size_jb),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                            , grid%vch &
#endif
                            )

   call da_transform_xatowrftl_lbc_wrapper( &
           it, cv_size, be, grid%ep, &
           xhat, iv, grid%vp6, grid%vv6, xbx, &
           grid, config_flags &
!#if (WRF_CHEM == 1)
!          , grid%vch6 &
!#endif
           )

   ! Final transfers to analysis increment --> analysis --> and model_grid
   call da_transform_xtoxa (grid)

   call da_transfer_xatowrf (grid, config_flags)

!!None of wrftoxb transfers should be needed since xb is unused in the update to the model
!! then overwritten in da_get_innov_vector
!#if (WRF_CHEM == 1)
!   if (cv_options_chem .gt. 0 .and. chem_opt(1) .gt. 0) &
!      call da_transfer_wrftoxb_chem(grid)
!   if (cv_options .gt. 0) &
!#endif
!      call da_transfer_wrftoxb(xbx, grid, config_flags)

   call da_transfer_headtomodel (grid)

   ! Turn off model boundary reading as we already provide a new one.
   if ( checkpoint_interval.eq.0 .or. var4d_lbc .or. .not.boundary_io) &
      call da_model_lbc_off

   if ( checkpoint_interval .gt. 0 ) write_checkpoints = 3
   original_disable = disable_traj
   disable_traj = .false.
   call da_nl_model(-1)
   disable_traj = original_disable
   if ( checkpoint_interval .gt. 0 ) write_checkpoints = 0

   ! Calculate new cost function
#if (WRF_CHEM == 1)
   if ( use_chemobs .and. ncv(2).gt.0 ) call da_retrieve_chem_hx( model_grid, iv )
#endif

   num_qcstat_conv=0
   call da_get_innov_vector (it+1+100*iter, num_qcstat_conv, ob, iv, grid , config_flags)

   call da_allocate_y (iv, re)
   call da_allocate_y (iv, y)
   call da_allocate_y (iv, jo_grad_y)
#if (WRF_CHEM == 1)
   call da_allocate_y_chem (iv, re)
   call da_allocate_y_chem (iv, y)
   call da_allocate_y_chem (iv, jo_grad_y)
#endif

!! The remainder of jo could be carried out by da_calculate_j (consider replacing)

   call da_calculate_residual(iv, y, re)

   call da_jo_and_grady(iv, re, jo_partial, j % jo, jo_grad_y &
#if (WRF_CHEM == 1)
                        , joc_partial &
#endif
                        )

   call da_deallocate_y (re)
   call da_deallocate_y (y)
   call da_deallocate_y (jo_grad_y)

#ifdef DM_PARALLEL
   j % jo % total = wrf_dm_sum_real(jo_partial)
#else 
   j % jo % total = jo_partial
#endif

   jout = 0.D0

   jout(1) = j % jo % total
#if (WRF_CHEM == 1)
   if ( ncv(2) .gt. 0 ) then
      jo_partial = &
                   j%jo%chem_acft &
                   + j%jo%chem_surf 
      jout(2) = wrf_dm_sum_real(jo_partial)
   end if
#endif

   write(unit=stdout,fmt='("DGN: jo = ",2F11.3)') jout

   !-------------------------------------------------------------------------
   ! [2.0] calculate jb:
   !-------------------------------------------------------------------------
   j % jb = 0.0
   if (cv_size_jb > 0) then
      cv_xhat_jb(1:cv_size_jb) = cv(1:cv_size_jb) + xhat_temp(1:cv_size_jb)
      j % jb = jb_factor * 0.5 * da_dot_cv(cv_size_jb,  cv_xhat_jb, cv_xhat_jb, grid, mz)

      write(unit=stdout,fmt='("DGN: jb = ",F11.3)') j%jb

      jout(1) = jout(1) + j % jb

      write(unit=stdout,fmt='("DGN: jtotal = ",2F11.3)') jout

#if (WRF_CHEM == 1)
      if ( ncv(1) .le. 0 ) then
         jout(2) = jout(1)
         jout(1) = 0.D0
      else if ( ncv(2) .gt. 0 ) then
         jout(2) = jout(2) + jb_factor * 0.5 * &
              da_dot_cv(ncv(2),  cv_xhat_jb(cv_s(2):cv_e(2)), cv_xhat_jb(cv_s(2):cv_e(2)), grid, mz)
         jout(1) = jout(1) - jout(2)
      end if
#endif

      write(unit=stdout,fmt='("DGN: jtotal = ",2F11.3)') jout

   end if

   !-------------------------------------------------------------------------
   ! [3.0] restore grid to appropriate state:
   !-------------------------------------------------------------------------

   evalj = .false.

!   !FROM MEMORY, not dependable since xtraj may be re-initialized
!   call domain_clock_get (grid, start_timestr=timestr)
!   call domain_clock_set( grid, current_timestr=timestr )
!   call da_read_basicstates ( xbx, grid, config_flags, timestr ) 

   write(unit=message(1),fmt='(A)') 'Resetting state variables to values from before da_evaluate_j.'
   write(unit=message(2),fmt='(3A)') '--> reading ', trim(fg_filename),' as INPUT stream'
   call wrf_message(trim(message(1)))
   call wrf_message(trim(message(2)))

   ! Reread fg from file since xtraj is overwritten by da_nl_model
   call da_med_initialdata_input (grid, config_flags, trim(fg_filename))

   ! Reset xb, which get_innov_vector modifies above
   call da_transfer_wrftoxb(xbx, grid, config_flags)

   ! Needed for chem IC's when they are added as CVs
#if (WRF_CHEM == 1)
   call da_transfer_wrftoxb_chem(grid)
#endif


   !Reset boundary values in grid, which were overwritten in da_transfer_xatowrf
   if (it.eq.1 .or. .not.var4d_lbc) then
      if (checkpoint_interval .eq. 0 .or. &
          (checkpoint_interval*60.lt.interval_seconds &
           .and. run_hours*3600.eq.interval_seconds)) then
         CALL WRFU_AlarmRingerOn( grid%alarms( BOUNDARY_ALARM ), rc=rc )
         if (it .eq. 1) then
            fg_filename = "wrfbdy_d<domain>"
         else
            write(fg_filename,fmt='("wrfvar_bdyout_",I2.2)') it-1
         end if
         write(unit=message(1),fmt='(A)') 'Resetting boundary variables to values from before da_evaluate_j.'
         write(unit=message(2),fmt='(3A)') '--> reading ', trim(fg_filename),' as BOUNDARY stream'
         call wrf_message(trim(message(1)))
         call wrf_message(trim(message(2)))

         config_flags%bdy_inname = trim(fg_filename)
         call med_latbound_in ( grid, config_flags )
         call close_dataset ( grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )
         config_flags%bdy_inname = "wrfbdy_d<domain>"
      end if
!!   else
!!      if (checkpoint_interval .eq. 0) then
!!         !NEED SOLUTION for var4d_lbc after first iteration (I/O of ana02_$it?)
!!      else
!!         !This solution might be more complicated
!      end if
   end if

   if (trace_use) call da_trace_exit("da_evaluate_j")

end subroutine da_evaluate_j

