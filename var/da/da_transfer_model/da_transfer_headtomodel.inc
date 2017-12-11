subroutine da_transfer_headtomodel ( grid )

   !-----------------------------------------------------------------------------------------
   ! Purpose:        Copy head_grid state to model_grid for initializing model with updated state
   !
   ! Called from:    da_solve, da_evaluate_j
   !
   ! History: 12/10/2017  Creation (JJ Guerrette )
   !
   !-----------------------------------------------------------------------------------------

   implicit none

   type(domain),   intent(inout)     :: grid

#ifdef VAR4D

   integer                           :: i

   if (trace_use) call da_trace_entry("da_transfer_headtomodel")

#if (WRF_CHEM == 1)
   if (.not. use_nonchemobs .and. cv_options.eq.0) then
      ! Reset MET variables in model_grid to their initial values
      call da_init_model_input
   else
#endif
      ! Set MET variables in model_grid to their posterior values
      ! This functionality is not available yet, since grid%XX is not updated for MET.
      call kj_swap (grid%u_2, model_grid%u_2, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      call kj_swap (grid%v_2, model_grid%v_2, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      call kj_swap (grid%w_2, model_grid%w_2, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      call kj_swap (grid%t_2, model_grid%t_2, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      call kj_swap (grid%ph_2, model_grid%ph_2, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      call kj_swap (grid%p, model_grid%p, &
                    grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)

      model_grid%mu_2 = grid%mu_2
      model_grid%t2 = grid%t2
      model_grid%th2 = grid%th2
      model_grid%q2 = grid%q2
      model_grid%u10 = grid%u10
      model_grid%v10 = grid%v10
      model_grid%tsk = grid%tsk
      model_grid%psfc = grid%psfc

      do i = PARAM_FIRST_SCALAR, num_moist
         call kj_swap (grid%moist(:,:,:,i), model_grid%moist(:,:,:,i), &
                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      enddo
#if (WRF_CHEM == 1)
   end if

   if (cv_options_chem .gt. 0 .and. chem_opt(1) .gt. 0) then

      ! Transferring chem from grid to model_grid
      ! --> Gives same result as da_init_model_input, unless chem is included
      !      as a control variable, in which case the code below is correct.
      do i = PARAM_FIRST_SCALAR, num_chem
         call kj_swap (grid%chem(:,:,:,i), model_grid%chem(:,:,:,i), &
                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, grid%xp%kms, grid%xp%kme)
      enddo

      ! Transfer posterior scaling factors to model_grid
      do i = PARAM_FIRST_SCALAR, num_scaleant
         call kj_swap (grid%scaleant(:,:,:,i), model_grid%scaleant(:,:,:,i), &
                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, 1, num_ant_steps)
      enddo
      do i = PARAM_FIRST_SCALAR, num_scalebb
         call kj_swap (grid%scalebb(:,:,:,i), model_grid%scalebb(:,:,:,i), &
                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, 1, num_bb_steps)
      enddo

   end if

   if (cv_options.gt.0) then
#endif

   ! Update boundary condition of model
      model_grid%u_bxs = grid%u_bxs
      model_grid%u_bxe = grid%u_bxe
      model_grid%u_bys = grid%u_bys
      model_grid%u_bye = grid%u_bye
      model_grid%v_bxs = grid%v_bxs
      model_grid%v_bxe = grid%v_bxe
      model_grid%v_bys = grid%v_bys
      model_grid%v_bye = grid%v_bye
      model_grid%w_bxs = grid%w_bxs
      model_grid%w_bxe = grid%w_bxe
      model_grid%w_bys = grid%w_bys
      model_grid%w_bye = grid%w_bye
      model_grid%t_bxs = grid%t_bxs
      model_grid%t_bxe = grid%t_bxe
      model_grid%t_bys = grid%t_bys
      model_grid%t_bye = grid%t_bye
      model_grid%mu_bxs = grid%mu_bxs
      model_grid%mu_bxe = grid%mu_bxe
      model_grid%mu_bys = grid%mu_bys
      model_grid%mu_bye = grid%mu_bye
      model_grid%ph_bxs = grid%ph_bxs
      model_grid%ph_bxe = grid%ph_bxe
      model_grid%ph_bys = grid%ph_bys
      model_grid%ph_bye = grid%ph_bye
      model_grid%moist_bxs = grid%moist_bxs
      model_grid%moist_bxe = grid%moist_bxe
      model_grid%moist_bys = grid%moist_bys
      model_grid%moist_bye = grid%moist_bye

      model_grid%u_btxs = grid%u_btxs
      model_grid%u_btxe = grid%u_btxe
      model_grid%u_btys = grid%u_btys
      model_grid%u_btye = grid%u_btye
      model_grid%v_btxs = grid%v_btxs
      model_grid%v_btxe = grid%v_btxe
      model_grid%v_btys = grid%v_btys
      model_grid%v_btye = grid%v_btye
      model_grid%w_btxs = grid%w_btxs
      model_grid%w_btxe = grid%w_btxe
      model_grid%w_btys = grid%w_btys
      model_grid%w_btye = grid%w_btye
      model_grid%t_btxs = grid%t_btxs
      model_grid%t_btxe = grid%t_btxe
      model_grid%t_btys = grid%t_btys
      model_grid%t_btye = grid%t_btye
      model_grid%mu_btxs = grid%mu_btxs
      model_grid%mu_btxe = grid%mu_btxe
      model_grid%mu_btys = grid%mu_btys
      model_grid%mu_btye = grid%mu_btye
      model_grid%ph_btxs = grid%ph_btxs
      model_grid%ph_btxe = grid%ph_btxe
      model_grid%ph_btys = grid%ph_btys
      model_grid%ph_btye = grid%ph_btye
      model_grid%moist_btxs = grid%moist_btxs
      model_grid%moist_btxe = grid%moist_btxe
      model_grid%moist_btys = grid%moist_btys
      model_grid%moist_btye = grid%moist_btye
#if (WRF_CHEM == 1)
   end if

   if (cv_options_chem .gt. 0 .and. chem_opt(1) .gt. 0) then
      model_grid%chem_bxs = grid%chem_bxs
      model_grid%chem_bxe = grid%chem_bxe
      model_grid%chem_bys = grid%chem_bys
      model_grid%chem_bye = grid%chem_bye

      model_grid%chem_btxs = grid%chem_btxs
      model_grid%chem_btxe = grid%chem_btxe
      model_grid%chem_btys = grid%chem_btys
      model_grid%chem_btye = grid%chem_btye
   end if
#endif

   if (trace_use) call da_trace_exit("da_transfer_headtomodel")

#endif

end subroutine da_transfer_headtomodel

