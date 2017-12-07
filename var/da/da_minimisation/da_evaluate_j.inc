subroutine da_evaluate_j(grid, config_flags, cv_size, cv_size_jb, &
                         cv_s, cv_e, eta, jout, &
                         xbx, be, iv, xhat, cv, ob, &
                         it, mz)

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
   integer, intent(in)                :: it 
   integer, intent(in)                :: mz(:)
#if (WRF_CHEM == 1)
   integer, parameter    :: ncv_type = 2
#else
   integer, parameter    :: ncv_type = 1
#endif
   integer, intent(in), dimension(1:ncv_type) :: cv_s, cv_e
   real, intent(in), dimension(1:ncv_type)    :: eta
   real, intent(out), dimension(1:ncv_type)   :: jout
!#if (WRF_CHEM == 1)
!   integer, intent(in)                :: cv_s_chem, cv_e_chem
!   real, intent(in)                   :: eta_chem
!   real, intent(out)                  :: jout_chem
!#endif

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

   character(len=filename_len) :: filename

   if (trace_use) call da_trace_entry("da_evaluate_j")

   evalj = .true.

   !-------------------------------------------------------------------------
   ! [0.0] initialization:
   !-------------------------------------------------------------------------

   call da_allocate_y (iv, re)
   call da_allocate_y (iv, y)
   call da_allocate_y (iv, jo_grad_y)
#if (WRF_CHEM == 1)
   call da_allocate_y_chem (iv, re)
   call da_allocate_y_chem (iv, y)
   call da_allocate_y_chem (iv, jo_grad_y)
#endif

   xhat_temp = xhat

   !Separable line search when chem is uncoupled, need to add coupling options for DA
   do i = 1, ncv_type
      ncv(i) = cv_e(i) - cv_s(i) + 1
      if ( ncv(i) .gt. 0 ) xhat_temp(cv_s(i):cv_e(i)) = eta(i) * xhat(cv_s(i):cv_e(i))
   end do

   !-------------------------------------------------------------------------
   ! [1.0] calculate jo:
   !-------------------------------------------------------------------------

   call da_transform_vtox (grid,cv_size_jb,xbx,be,grid%ep,&
                            xhat_temp(1:cv_size_jb),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                            , grid%vch &
#endif
                            )

   call da_transform_xtoxa (grid)

   call da_transfer_xatoanalysis (it, xbx, grid, config_flags)

   call da_transfer_headtomodel (grid)

   ! Turn off model boundary reading as we already provide a new one.
   if ( checkpoint_interval.eq.0 .or. var4d_lbc .or. .not.boundary_io) &
      call da_model_lbc_off

   if ( checkpoint_interval .gt. 0 ) write_checkpoints = 3
   call da_nl_model(-1)
   if ( checkpoint_interval .gt. 0 ) write_checkpoints = 0

   ! Calculate new cost function
#if (WRF_CHEM == 1)
   if ( use_chemobs .and. ncv(2).gt.0 ) call da_retrieve_chem_hx( model_grid, iv )
#endif
   call da_get_innov_vector (2, num_qcstat_conv, ob, iv, grid , config_flags)

!! The remainder of jo could be carried out by da_calculate_j (consider replacing)

   call da_calculate_residual(iv, y, re)

   call da_jo_and_grady(iv, re, jo_partial, j % jo, jo_grad_y &
#if (WRF_CHEM == 1)
                        , joc_partial &
#endif
                        )

   if (test_dm_exact) then
      ! jo_partial has been already summed at lower level
      j % jo % total = jo_partial
   else
      j % jo % total = wrf_dm_sum_real(jo_partial)
   end if

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

   !-------------------------------------------------------------------------
   ! [2.0] calculate jb:
   !-------------------------------------------------------------------------
   j % jb = 0.0
   if (cv_size_jb > 0) then
      cv_xhat_jb(1:cv_size_jb) = cv(1:cv_size_jb) + xhat_temp(1:cv_size_jb)
      j % jb = jb_factor * 0.5 * da_dot_cv(cv_size_jb,  cv_xhat_jb, cv_xhat_jb, grid, mz)

      jout(1) = jout(1) + j % jb

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
   end if


   call da_deallocate_y (re)
   call da_deallocate_y (y)
   call da_deallocate_y (jo_grad_y)

   !-------------------------------------------------------------------------
   ! [3.0] restore grid to appropriate state:
   !-------------------------------------------------------------------------
   if (it .eq. 1) then
#if (WRF_CHEM == 1)
      if ( ncv(2) .gt. 0 ) then
         filename = 'fgc'
      else
#endif
         filename = 'fg'
#if (WRF_CHEM == 1)
      end if
#endif
   else
      write(filename,fmt='("wrfvar_output_",I2.2)') it-1
   end if

   write(unit=message(1),fmt='(A)') 'Resetting state variables to values from previous iteration.'
   write(unit=message(2),fmt='(3A)') '--> reading ', filename,' as INPUT stream'
   call wrf_message(trim(message(1)))
   call wrf_message(trim(message(2)))

   call da_med_initialdata_input (grid, config_flags, trim(filename))

!setup_firstguess would probably be more general, but also has extra initialization code
!   call da_setup_firstguess(xbx, grid, config_flags, .false.)

   call da_transfer_wrftoxb(xbx, grid, config_flags)
#if (WRF_CHEM == 1)
   call da_transfer_wrftoxb_chem(grid)
#endif


   !NEED TO RESTORE BOUNDARY with file READ or recall from memory, 
   ! plus another firstguess setup (see da_solve.inc)

   evalj = .false.

   if (trace_use) call da_trace_exit("da_evaluate_j")

end subroutine da_evaluate_j

