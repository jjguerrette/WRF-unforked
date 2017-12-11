subroutine da_init_model_input

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   !-----------------------------------------------------------------------

   implicit none

   character(len=256) :: original_input

   if (trace_use_dull) call da_trace_entry("da_init_model_input")

   head_grid%itimestep = 0

   original_input = trim(config_flags%input_inname)

   config_flags%input_inname = "fg"

   call med_initialdata_input (head_grid, config_flags)

   config_flags%input_inname = trim(original_input)

   if (trace_use_dull) call da_trace_exit("da_init_model_input")

end subroutine da_init_model_input
