subroutine da_transform_ytoyhat(iv, y, yhat_obs, nall)

   !-------------------------------------------------------------------------
   ! Purpose: Does part of the obs gradient operation   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector (O-B).
   type (y_type),     intent(inout) :: y           ! Residual vector (O-A).
   real, allocatable, intent(out)   :: yhat_obs(:)
   integer,           intent(out)   :: nall (num_ob_indexes)

   type (yhat_type) :: yhat(num_ob_indexes)
   integer :: istart, iend, itype

   if (trace_use) call da_trace_entry("da_transform_ytoyhat")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------

   do itype = 1, num_ob_indexes
      nall(itype) = 0

write(unit=stdout,fmt='(A,2I6)') &
   'DEBUG_A01',itype,nall(itype)

      if ( iv%info(itype)%nlocal .gt. 0 ) then
         select case (itype)
!            case (sound)
!               call da_ytoyhat_sound    (iv, y, yhat(itype)%obs)
!            case (sonde_sfc)
!               call da_ytoyhat_sonde_sfc(iv, y, yhat(itype)%obs)
!            case (synop)
!               call da_ytoyhat_synop    (iv, y, yhat(itype)%obs)
!            case (geoamv)
!               call da_ytoyhat_geoamv   (iv, y, yhat(itype)%obs)
!            case (polaramv)
!               call da_ytoyhat_polaramv (iv, y, yhat(itype)%obs)
!            case (airep)
!               call da_ytoyhat_airep    (iv, y, yhat(itype)%obs)
!            case (pilot)
!               call da_ytoyhat_pilot    (iv, y, yhat(itype)%obs)
!            case (profiler)
!               call da_ytoyhat_profiler (iv, y, yhat(itype)%obs)
!            case (satem)
!               call da_ytoyhat_satem    (iv, y, yhat(itype)%obs)
!            case (metar)
!               call da_ytoyhat_metar    (iv, y, yhat(itype)%obs)
!            case (ships)
!               call da_ytoyhat_ships    (iv, y, yhat(itype)%obs)
!            case (buoy)
!               call da_ytoyhat_buoy     (iv, y, yhat(itype)%obs)
!            case (gpspw)
!               call da_ytoyhat_gpspw    (iv, y, yhat(itype)%obs)
!            case (gpsref)
!               call da_ytoyhat_gpsref   (iv, y, yhat(itype)%obs)
!            case (ssmi_tb)
!               call da_ytoyhat_ssmi_tb  (iv, y, yhat(itype)%obs) 
!            case (ssmi_rv)
!               call da_ytoyhat_ssmi_rv  (iv, y, yhat(itype)%obs) 
!            case (ssmt1)
!               call da_ytoyhat_ssmt1    (iv, y, yhat(itype)%obs)
!            case (ssmt2)
!               call da_ytoyhat_ssmt2    (iv, y, yhat(itype)%obs)
!            case (pseudo)
!               call da_ytoyhat_pseudo   (iv, y, yhat(itype)%obs)
!            case (bogus)
!               call da_ytoyhat_bogus    (iv, y, yhat(itype)%obs)  
!            case (qscat)
!               call da_ytoyhat_qscat    (iv, y, yhat(itype)%obs)
!            case (radar)
!               call da_ytoyhat_radar    (iv, y, yhat(itype)%obs)
!            case (mtgirs)
!               call da_ytoyhat_mtgirs   (iv, y, yhat(itype)%obs)
!            case (tamdar)
!               call da_ytoyhat_tamdar   (iv, y, yhat(itype)%obs)
!            case (tamdar_sfc)
!               call da_ytoyhat_tamdar_sfc(iv, y, yhat(itype)%obs)
!            case (rain)
!               call da_ytoyhat_rain     (iv, y, yhat(itype)%obs)
!
!            case (airsr)
!               call da_ytoyhat_airsr    (iv, y, yhat(itype)%obs)

#if (WRF_CHEM == 1)
            case (chem_surf)
               call da_ytoyhat_chem_surf (iv, y, yhat(itype)%obs)
            case (chem_acft)
               call da_ytoyhat_chem_acft (iv, y, yhat(itype)%obs)
#endif
         end select

!#if defined(CRTM) || defined(RTTOV)
!      else if ( itype .eq. radiance .and. iv%num_inst .gt. 0 ) then
!              call da_ytoyhat_rad      (iv, y, yhat(itype)%obs)
!#endif

         if (allocated(yhat(itype)%obs)) then
            nall(itype) = size(yhat(itype)%obs)
         end if
      end if


write(unit=stdout,fmt='(A,2I6)') &
   'DEBUG_A02',itype,nall(itype)

   end do

   if (sum(nall) .gt. 0) then
      allocate( yhat_obs(sum(nall)) )

      istart = 1
      do itype = 1, num_ob_indexes

write(unit=stdout,fmt='(A,2I6)') &
   'DEBUG_A03',itype,nall(itype)


         if ( nall(itype) .gt. 0 .and. allocated(yhat(itype)%obs) ) then
            iend = istart + nall(itype) - 1

            yhat_obs(istart:iend) = yhat(itype)%obs(1:nall(itype))

            istart = iend + 1
         end if

write(unit=stdout,fmt='(A,2I6)') &
   'DEBUG_A04',itype,nall(itype)

      end do
   end if

   if (trace_use) call da_trace_exit("da_transform_ytoyhat")

end subroutine da_transform_ytoyhat


