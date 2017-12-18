subroutine da_transform_xatowrftl_lbc_wrapper(it, cv_size, be, ep, cv, iv, vp6, vv6, xbx, &
                              grid, config_flags &
!#if (WRF_CHEM == 1)
!                        ,vch6 &
!#endif
                        )

   !----------------------------------------------------------------------------------------
   ! Purpose:  Does transform of control variable xa to model TL space for boundary updating
   !----------------------------------------------------------------------------------------

   implicit none

   integer,                    intent(in)    :: it

   integer,                    intent(in)    :: cv_size ! Size of cv array.
   type(be_type),              intent(in)    :: be     ! background error structure.
   type(ep_type),              intent(in)    :: ep     ! Ensemble perturbation structure.
   real,                       intent(in)    :: cv(1:cv_size)     ! control variables.
   type(iv_type),              intent(inout) :: iv     ! innovation vector (o-b).
   type(vp_type),              intent(inout) :: vp6     ! Grdipt/level CV for 6h.
   type(vp_type),              intent(inout) :: vv6     ! Grdipt/EOF CV for 6h.
!#if (WRF_CHEM == 1)
!   type(xch_type), optional, intent(out) :: vch6   ! CHEM CV
!#endif
   type(xbx_type),             intent(inout) :: xbx    ! For header & non-grid arrays.
   type(domain),               intent(inout) :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags

   type(x_type) :: shuffle
   integer :: nobwin, jl_start, jl_end
   character(len=filename_len) :: fg_filename

   character(len=256) :: timestr, timestr1

   if (trace_use) call da_trace_entry("da_transform_xatowrftl_lbc_wrapper")

   call domain_clock_get( grid, start_timestr=timestr )
   call da_transfer_xatowrftl(grid, config_flags, 'tl01', timestr)

   if ( var4d_lbc ) then 
      call domain_clock_get (grid, stop_timestr=timestr1)
      call domain_clock_set( grid, current_timestr=timestr1 )
      grid%u_2 = u6_2 ; grid%v_2 = v6_2; grid%t_2 = t6_2;
      grid%w_2 = w6_2 ; grid%mu_2 = mu6_2 ; grid%ph_2 =ph6_2
      grid%moist = moist6; grid%p = p6; grid%psfc = psfc6
      call da_transfer_wrftoxb(xbx, grid, config_flags)

      shuffle = grid%xa
      jl_start    = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + 1
      jl_end      = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_jl
      grid%xa  = grid%x6a
      call da_transform_vtox(grid, be%cv%size_jl, xbx, be, ep, cv(jl_start:jl_end), vv6, vp6)
      grid%xa  = shuffle

      call da_transfer_xatowrftl_lbc(grid, config_flags, 'tl01')

      call domain_clock_get( grid, start_timestr=timestr1 )
      call domain_clock_set( grid, current_timestr=timestr1 )

      if (nl_called .eq. 1) then
         call da_read_basicstates ( xbx, grid, config_flags, timestr1 )
      else
         if (it .eq. 1) then
            fg_filename = 'fg'
         else
            write(fg_filename,fmt='("wrfvar_output_",I2.2)') it-1
         end if

         write(unit=message(1),fmt='(3A)') '--> reading ', trim(fg_filename),' as INPUT stream'
         call wrf_message(trim(message(1)))

         ! Reread fg from file since xtraj may have been overwritten by da_nl_model
         call da_med_initialdata_input (grid, config_flags, trim(fg_filename))

         ! Reset xb, which get_innov_vector modifies above
         call da_transfer_wrftoxb(xbx, grid, config_flags)

         ! Needed for chem IC's when they are added as CVs
#if (WRF_CHEM == 1)
         call da_transfer_wrftoxb_chem(grid)
#endif
      end if
   else
      call da_transfer_xatowrftl_lbc(grid, config_flags, 'tl01')
   endif


   if (trace_use) call da_trace_exit("da_transform_xatowrftl_lbc_wrapper")

end subroutine da_transform_xatowrftl_lbc_wrapper


