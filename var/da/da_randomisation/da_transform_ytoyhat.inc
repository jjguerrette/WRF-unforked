subroutine da_transform_ytoyhat(iv, re, yhat_all, nall)

   !-------------------------------------------------------------------------
   ! Purpose: Does part of the obs gradient operation   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: re          ! Residual vector (O-A).
   real, allocatable, intent(out)   :: yhat_all(:)
   integer,           intent(out)   :: nall (num_ob_indexes)

   type (vec_type), dimension(num_ob_indexes) :: yhat
   integer :: istart, iend

   if (trace_use) call da_trace_entry("da_transform_ytoyhat")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------

   do itype = 1, num_ob_indexes
      nall(itype) = 0

      if ( iv%info(itype)%nlocal .gt. 0 ) then
         select case (itype)
!            case (sound)
!               call da_ytoyhat_sound    (iv, re, yhat(itype)%vec)
!            case (sonde_sfc)
!               call da_ytoyhat_sonde_sfc(iv, re, yhat(itype)%vec)
!            case (synop)
!               call da_ytoyhat_synop    (iv, re, yhat(itype)%vec)
!            case (geoamv)
!               call da_ytoyhat_geoamv   (iv, re, yhat(itype)%vec)
!            case (polaramv)
!               call da_ytoyhat_polaramv (iv, re, yhat(itype)%vec)
!            case (airep)
!               call da_ytoyhat_airep    (iv, re, yhat(itype)%vec)
!            case (pilot)
!               call da_ytoyhat_pilot    (iv, re, yhat(itype)%vec)
!            case (profiler)
!               call da_ytoyhat_profiler (iv, re, yhat(itype)%vec)
!            case (satem)
!               call da_ytoyhat_satem    (iv, re, yhat(itype)%vec)
!            case (metar)
!               call da_ytoyhat_metar    (iv, re, yhat(itype)%vec)
!            case (ships)
!               call da_ytoyhat_ships    (iv, re, yhat(itype)%vec)
!            case (buoy)
!               call da_ytoyhat_buoy     (iv, re, yhat(itype)%vec)
!            case (gpspw)
!               call da_ytoyhat_gpspw    (iv, re, yhat(itype)%vec)
!            case (gpsref)
!               call da_ytoyhat_gpsref   (iv, re, yhat(itype)%vec)
!            case (ssmi_tb)
!               call da_ytoyhat_ssmi_tb  (iv, re, yhat(itype)%vec) 
!            case (ssmi_rv)
!               call da_ytoyhat_ssmi_rv  (iv, re, yhat(itype)%vec) 
!            case (ssmt1)
!               call da_ytoyhat_ssmt1    (iv, re, yhat(itype)%vec)
!            case (ssmt2)
!               call da_ytoyhat_ssmt2    (iv, re, yhat(itype)%vec)
!            case (pseudo)
!               call da_ytoyhat_pseudo   (iv, re, yhat(itype)%vec)
!            case (bogus)
!               call da_ytoyhat_bogus    (iv, re, yhat(itype)%vec)  
!            case (qscat)
!               call da_ytoyhat_qscat    (iv, re, yhat(itype)%vec)
!            case (radar)
!               call da_ytoyhat_radar    (iv, re, yhat(itype)%vec)
!            case (mtgirs)
!               call da_ytoyhat_mtgirs   (iv, re, yhat(itype)%vec)
!            case (tamdar)
!               call da_ytoyhat_tamdar   (iv, re, yhat(itype)%vec)
!            case (tamdar_sfc)
!               call da_ytoyhat_tamdar_sfc(iv, re, yhat(itype)%vec)
!            case (rain)
!               call da_ytoyhat_rain     (iv, re, yhat(itype)%vec)
!
!            case (airsr)
!               call da_ytoyhat_airsr    (iv, re, yhat(itype)%vec)

#if (WRF_CHEM == 1)
            case (chem_surf)
               call da_ytoyhat_chem_surf (iv, re, yhat(itype)%vec)
            case (chem_acft)
               call da_ytoyhat_chem_acft (iv, re, yhat(itype)%vec)
#endif
         end select

!#if defined(CRTM) || defined(RTTOV)
!      else if ( itype .eq. radiance .and. iv%num_inst .gt. 0 ) then
!              call da_ytoyhat_rad      (iv, re, yhat(itype)%vec)
!#endif


      end if

      if (allocated(yhat(itype)%vec) then
         nall(itype) = size(yhat(itype)%vec)
      end if
   end do

   if (sum(nall) .gt. 0) then
      allocate( yhat_all(sum(nall)) )

      istart = 1
      do itype = 1, num_ob_indexes
         if ( nall(itype) .gt. 0 ) then
            iend = istart + nall(itype) - 1

            yhat_all(istart:iend) = yhat(itype)%vec

            istart = iend + 1
         end if
      end do
   end if

   if (trace_use) call da_trace_exit("da_transform_ytoyhat")

end subroutine da_transform_ytoyhat


