subroutine da_nl_model (it)

   !-----------------------------------------------------------------------
   ! Purpose: Acrivate nonlinear model
   !    Author: Xin Zhang,  MMM/ESSL/NCAR,  Date: 04/28/2010
   !-----------------------------------------------------------------------

   implicit none
   
   integer, intent(in) :: it ! outer looper counter.
   integer :: time_step, alarmid
   character(len=256) :: rst_outname
#if (WRF_CHEM == 1)
   integer :: ioform8, ioform16, io_form
#endif

   if (trace_use_dull) call da_trace_entry("da_nl_model")

#if (WRF_CHEM != 1)
   if ( it == 1 ) then
#else
   if ( it == 1 .or. init_osse_chem ) then
#endif
#ifdef _OPENMP
      call setfeenv()
#endif

#ifdef DM_PARALLEL
      ! We need the communicator generated in WRFDA
      local_communicator = comm
#endif

      ! Set up WRF model.
      call wrf_init ( .false. )
      model_config_rec%var4d_run = .true.
      config_flags%var4d_run = .true.
   endif

   ! 4dvar init
!   if (.not. evalj) call init_4dvar
   call init_4dvar

   if (it .gt. 0) then
      nl_called = 1 
   else
      nl_called = nl_called + 1
   end if


#if (WRF_CHEM == 1)
   if ( use_chemobs .and. cv_options_chem .gt. 0 ) then
      call da_copy_chem_to_model( head_grid, it )

      !Adds scaleant, scalebb, or other variables initialized in WRFDA to fg file
      if ( it == 1 .and. init_chem_scale ) then
!         call da_init_model_output()
         call da_add_var_to_firstguess()
      end if
   end if
#endif

   if (write_checkpoints .eq. 2) return

   ! Turn off  history output
   call nl_set_history_begin_y (head_grid%id, 9999)
   call nl_set_history_end_y (head_grid%id, 9999)
   call nl_set_history_interval (head_grid%id, 0)

   ! Set up basic states output
   if ( multi_inc == 1 ) then
     call nl_set_auxhist6_interval_s ( head_grid%id, var4d_bin )
   else
     call nl_get_time_step ( head_grid%id, time_step )
     call nl_set_auxhist6_interval_s ( head_grid%id, time_step ) 
   endif
   call nl_set_io_form_auxhist6 ( head_grid%id, 2 ) 
   call nl_set_frames_per_auxhist6 ( head_grid%id, 1 ) 

   !Turn off xtraj writing
!   if (write_checkpoints.gt.0 .or. evalj) call nl_set_auxhist6_interval_s ( head_grid%id, 0 )
   if (write_checkpoints.gt.0) call nl_set_auxhist6_interval_s ( head_grid%id, 0 )

   ! Set up basic states writing for offline obs. operators
!   if ((any(write_checkpoints.eq.(/1,3/)) .or. evalj) &
   if (any(write_checkpoints.eq.(/1,3/)) &
#if (WRF_CHEM == 1)
        .and. use_offlineobs &
#endif
        ) then
      model_config_rec%auxhist16_outname = trim(xtrajprefix)//"_d<domain>_<date>"
      call nl_set_auxhist16_interval_s (head_grid%id, var4d_bin )
      call nl_get_io_form_auxhist16 ( head_grid%id, ioform16 )
      if (ioform16 .le. 0) call nl_set_io_form_auxhist16 ( head_grid%id, io_form_xtraj )
      call nl_set_frames_per_auxhist16 ( head_grid%id, 1 )
   end if

   ! Set up checkpoint file writing (restart files) for offline and online obs. operators
   if ( write_checkpoints .eq. 1 ) then
      call nl_get_restart_interval_h ( head_grid%id, original_restart_interval)
#if (WRF_CHEM == 1)
      if( .not. calc_hx_only) then
#endif
         !Turn on Output for standard checkpoint files
         call nl_set_restart_interval_m ( head_grid%id, checkpoint_interval)
         call nl_set_restart_interval_h ( head_grid%id, 0)
         call nl_get_io_form_restart ( head_grid%id, io_form )
         if (io_form .le. 0) call nl_set_io_form_restart ( head_grid%id, io_form_xtraj )
         rst_outname = "wrf_checkpoint_d<domain>_<date>"
         call nl_set_rst_outname ( head_grid%id, rst_outname)
!Should not be needed, equal to default
!         call nl_set_restart_begin_m( head_grid%id, 0 )
#if (WRF_CHEM == 1)
      endif
#endif
      ! input_inname used to initialize first checkpoints in TL/AD models
      if (it .eq. 1) then
         model_config_rec%input_inname = "fg"
      else if (it .gt. 1) then
         write(model_config_rec%input_inname,fmt='("wrfvar_output_",I2.2)') it-1
      end if
#if (WRF_CHEM == 1)
   else
      ! Turn off gocart_bg field reading
!      call nl_get_io_form_auxinput8(head_grid%id, ioform8)
      if ( model_config_rec%start_hour(head_grid%id) .ne. 0 ) then
         call nl_set_io_form_auxinput8(head_grid%id, 0)
      endif
#endif
   endif


   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   call Setup_Timekeeping ( head_grid )

   ! Set the physics options to what they are in TL and AD
   !call nl_set_mp_physics (head_grid%id, 0)
   !call nl_set_ra_lw_physics (head_grid%id, 0)
   !call nl_set_ra_sw_physics (head_grid%id, 0)
   !call nl_set_sf_sfclay_physics (head_grid%id, 0)
   !call nl_set_bl_pbl_physics (head_grid%id, 0)
   !call nl_set_cu_physics (head_grid%id, 0)
   !call nl_set_mp_zero_out (head_grid%id, 0)
   !call nl_set_sf_surface_physics (head_grid%id, 1)
   !call nl_set_ifsnow (head_grid%id, 1)
   !call nl_set_icloud (head_grid%id, 0)
   !call nl_set_isfflx (head_grid%id, 0)

   ! WRF model time-stepping.  Calls integrate().
!   ! Automatic restarting after each cycle_interval.  Useful for long assimilation windows.
!   ! Not fully functional yet. -JJGuerrette
!   if ( cycle_interval .gt. 0 ) then 
!      call wrf_run_cycle
!   else
      call wrf_run
!   end if

   ! restore the physics options to the originals
   !call nl_set_mp_physics (head_grid%id, original_mp_physics)
   !call nl_set_ra_lw_physics (head_grid%id, original_ra_lw_physics)
   !call nl_set_ra_sw_physics (head_grid%id, original_ra_sw_physics)
   !call nl_set_sf_sfclay_physics (head_grid%id, original_sf_sfclay_physics)
   !call nl_set_bl_pbl_physics (head_grid%id, original_bl_pbl_physics)
   !call nl_set_cu_physics (head_grid%id, original_cu_physics)
   !call nl_set_mp_zero_out (head_grid%id, original_mp_zero_out)
   !call nl_set_sf_surface_physics (head_grid%id, original_sf_surface_physics)
   !call nl_set_ifsnow (head_grid%id, original_ifsnow)
   !call nl_set_icloud (head_grid%id, original_icloud)
   !call nl_set_isfflx (head_grid%id, original_isfflx)

   ! Turn off basic states output
   call nl_set_io_form_auxhist6 ( head_grid%id, 0 ) 
   call nl_set_auxhist6_interval_s ( head_grid%id, 0 ) 

!   if ( (any(write_checkpoints.eq.(/1,3/)) .or. evalj) &
   if (any(write_checkpoints.eq.(/1,3/)) &
#if (WRF_CHEM == 1)
        .and. use_offlineobs &
#endif
        ) then
      model_config_rec%auxhist16_outname = "auxhist16_d<domain>_<date>"
      call nl_set_auxhist16_interval_s (head_grid%id, 0 )
      call nl_set_io_form_auxhist16 ( head_grid%id, ioform16 )
   end if

   if ( write_checkpoints .eq. 1 ) then
      call nl_set_restart_interval_h ( head_grid%id, original_restart_interval)
      call nl_set_restart_interval_m ( head_grid%id, 0)
      call nl_set_io_form_restart ( head_grid%id, 0 )
      rst_outname = "wrfrst_d<domain>_<date>"
      call nl_set_rst_outname ( head_grid%id, rst_outname)
      call nl_set_restart_begin_h( head_grid%id, original_restart_interval)
   end if

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   call Setup_Timekeeping ( head_grid )

   ! Turn off boundary input from now.
   if ( checkpoint_interval.eq.0 .or. var4d_lbc ) &
      CALL nl_set_io_form_boundary( head_grid%id, 0 )

   if (trace_use_dull) call da_trace_exit("da_nl_model")

end subroutine da_nl_model
