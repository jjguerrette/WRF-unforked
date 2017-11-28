subroutine da_read_basicstates ( xbx, grid, config_flags, timestr, fgat )

   !-------------------------------------------------------------------------
   ! Purpose: Read basic state at time = timestr
   !
   ! History: 10/01/2010  Creation (Xin Zhang )
   !
   !-------------------------------------------------------------------------

   implicit none

   character(len=256),         intent(in)    ::   timestr
   type(domain),               intent(inout) ::   grid
   type(xbx_type),             intent(inout) ::   xbx
   type(grid_config_rec_type), intent(inout) ::   config_flags
   character(len=120), intent(in), optional  ::   fgat

   integer :: io_form

#ifdef VAR4D
   if ( grid%var4d ) then
      if ( grid%trajectory_io ) then
         if ( checkpoint_interval .le. 0 ) then
            call input_nl_xtraj (timestr)
!            call da_swap_xtraj (grid)
         else
!            if ( xtrajswitch .gt. 0) then
               call input_nl_xtraj_checkpt (timestr)
!               call da_swap_xtraj (grid)
!            else
!
!               call nl_get_io_form_auxinput16 ( grid%id, io_form )
!               if (io_form .le. 0) then
!                  call nl_set_io_form_auxinput16 ( grid%id, io_form_xtraj )
!                  config_flags%io_form_auxinput16 = io_form_xtraj
!               else
!                  config_flags%io_form_auxinput16 = io_form
!               end if
!
!               config_flags%frames_per_auxinput16 = 1
!
!               config_flags%auxinput16_inname = trim(xtrajprefix)//"_d<domain>_<date>"
!
!               call med_auxinput_in ( grid, auxinput16_only, config_flags )
!
!               call nl_set_io_form_auxinput16 ( grid%id, io_form )
!            end if

            call da_swap_xtraj (grid)
         end if
      else
         config_flags%auxinput16_inname = "auxhist6_d<domain>_<date>"
         call med_auxinput_in ( grid, auxinput16_only, config_flags )
      endif
      call da_transfer_wrftoxb(xbx, grid, config_flags)
   endif
#endif
! Following codes are for FGAT
   if ( .not. grid%var4d .and. grid%num_fgat_time > 1 .and. present(fgat) ) then
      call da_med_initialdata_input (grid, config_flags, trim(fgat))
      call da_transfer_wrftoxb(xbx, grid, config_flags)
   end if
   return

end subroutine da_read_basicstates

