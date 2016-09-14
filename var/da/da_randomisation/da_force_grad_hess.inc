subroutine da_force_grad_hess(iv, y, nall, hessmember, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Does part of the obs gradient operation   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y           ! Residual vector (O-A).
   integer,           intent(out)   :: nall (num_ob_indexes)
   integer,           intent(in)    :: hessmember
   type (y_type),     intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: istart, iend, itype

   if (trace_use) call da_trace_entry("da_force_grad_hess")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------

   nall(:) = 0


   do itype = 1, num_ob_indexes
      if ( iv%info(itype)%nlocal .gt. 0 ) then
         select case (itype)
!            case (sound)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_sound    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (sonde_sfc)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_sonde_sfc(iv, y, jo_grad_y, hessmember, nall(itype))
!            case (synop)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_synop    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (geoamv)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_geoamv   (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (polaramv)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_polaramv (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (airep)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_airep    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (pilot)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_pilot    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (profiler)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_profiler (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (satem)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_satem    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (metar)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_metar    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (ships)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_ships    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (buoy)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_buoy     (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (gpspw)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_gpspw    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (gpsref)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_gpsref   (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (ssmi_tb)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_ssmi_tb  (iv, y, jo_grad_y, hessmember, nall(itype)) 
!            case (ssmi_rv)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_ssmi_rv  (iv, y, jo_grad_y, hessmember, nall(itype)) 
!            case (ssmt1)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_ssmt1    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (ssmt2)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_ssmt2    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (pseudo)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_pseudo   (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (bogus)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_bogus    (iv, y, jo_grad_y, hessmember, nall(itype))  
!            case (qscat)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_qscat    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (radar)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_radar    (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (mtgirs)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_mtgirs   (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (tamdar)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_tamdar   (iv, y, jo_grad_y, hessmember, nall(itype))
!            case (tamdar_sfc)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_tamdar_sfc(iv, y, jo_grad_y, hessmember, nall(itype))
!            case (rain)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_rain     (iv, y, jo_grad_y, hessmember, nall(itype))
!
!            case (airsr)
!               if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_airsr    (iv, y, jo_grad_y, hessmember, nall(itype))

#if (WRF_CHEM == 1)
            case (chem_surf)
               if (sum(nall) .lt. hessmember) &
                  call da_force_grad_chem_surf (iv, y, jo_grad_y, hessmember, nall(itype))
            case (chem_acft)
               if (sum(nall) .lt. hessmember) &
                  call da_force_grad_chem_acft (iv, y, jo_grad_y, hessmember, nall(itype))
#endif
         end select

!#if defined(CRTM) || defined(RTTOV)
!      else if ( itype .eq. radiance .and. iv%num_inst .gt. 0 ) then
!              if (sum(nall) .lt. hessmember) &
!                  call da_force_grad_rad      (iv, y, jo_grad_y, hessmember, nall(itype))
!#endif

      end if
   end do

   if (trace_use) call da_trace_exit("da_force_grad_hess")

end subroutine da_force_grad_hess


