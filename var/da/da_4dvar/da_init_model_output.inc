subroutine da_init_model_output(filename)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   character(*), intent(in), optional    :: filename               ! I/O filename

   if (trace_use_dull) call da_trace_entry("da_init_model_output")

   if ( present(filename) ) then
      config_flags%input_outname = trim(filename)
   else
      config_flags%input_outname = "fg"
   end if
   config_flags%write_input = .true.

   call med_filter_out (head_grid, config_flags)

   config_flags%input_outname = "wrfinput_d<domain>"
   config_flags%write_input = .false.

   if (trace_use_dull) call da_trace_exit("da_init_model_output")

end subroutine da_init_model_output
