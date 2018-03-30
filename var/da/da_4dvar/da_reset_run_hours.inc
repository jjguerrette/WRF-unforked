subroutine da_reset_run_hours ()

   !-----------------------------------------------------------------------
   ! Purpose: Reset run hours
   !    Author: JJ Guerrette,  CU Boulder,  Date: 07/25/2016
   !-----------------------------------------------------------------------

   implicit none
   
   integer :: alarmid

   if (trace_use_dull) call da_trace_entry("da_reset_run_hours")

   if ( adtl_run_hours.gt.0 &
        .and. run_hours.gt.adtl_run_hours &
        .and. run_days.eq.0 ) then

      model_config_rec%run_hours = run_hours

      model_config_rec%end_year(head_grid%id) = end_year0
      model_config_rec%end_month(head_grid%id) = end_month0
      model_config_rec%end_day(head_grid%id) = end_day0
      model_config_rec%end_hour(head_grid%id) = end_hour0

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

   end if

   if (trace_use_dull) call da_trace_exit("da_reset_run_hours")

end subroutine da_reset_run_hours
