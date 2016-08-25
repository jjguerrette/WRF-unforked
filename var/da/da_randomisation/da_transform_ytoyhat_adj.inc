subroutine da_transform_ytoyhat_adj(iv, qhat_obs, nall, jo_grad_y)

   !-------------------------------------------------------------------------
   ! Purpose: Does part of the obs gradient operation   
   !-------------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector (O-B).
   real, allocatable, intent(in)    :: qhat_obs(:)
   integer,           intent(in)    :: nall(num_ob_indexes)
   type (y_type),     intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   type (yhat_type) :: qhat(num_ob_indexes)
   integer :: istart, iend, itype

   if (trace_use) call da_trace_entry("da_transform_ytoyhat_adj")

   !-------------------------------------------------------------------------
   ! [1.0] Compute components of Grad_y(Jo):
   !-------------------------------------------------------------------------

   if (sum(nall) .gt. 0 .and. allocated(qhat_obs) ) then

      istart = 1
      do itype = 1, num_ob_indexes
         if ( nall(itype) .gt. 0 ) then
            iend = istart + nall(itype) - 1
            allocate( qhat(itype)%obs(nall(itype)) )
            qhat(itype)%obs(1:nall(itype)) = qhat_obs(istart:iend)
   
            istart = iend + 1
         end if
      end do

      do itype = 1, num_ob_indexes
         if ( iv%info(itype)%nlocal .gt. 0 ) then
            select case (itype)
   !            case (sound)
   !               call da_ytoyhat_sound_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (sonde_sfc)
   !               call da_ytoyhat_sonde_sfc_adj (iv, qhat(itype)%obs, jo_grad_y)
   !            case (synop)
   !               call da_ytoyhat_synop_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (geoamv)
   !               call da_ytoyhat_geoamv_adj    (iv, qhat(itype)%obs, jo_grad_y)
   !            case (polaramv)
   !               call da_ytoyhat_polaramv_adj  (iv, qhat(itype)%obs, jo_grad_y)
   !            case (airep)
   !               call da_ytoyhat_airep_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (pilot)
   !               call da_ytoyhat_pilot_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (profiler)
   !               call da_ytoyhat_profiler_adj  (iv, qhat(itype)%obs, jo_grad_y)
   !            case (satem)
   !               call da_ytoyhat_satem_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (metar)
   !               call da_ytoyhat_metar_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (ships)
   !               call da_ytoyhat_ships_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (buoy)
   !               call da_ytoyhat_buoy_adj      (iv, qhat(itype)%obs, jo_grad_y)
   !            case (gpspw)
   !               call da_ytoyhat_gpspw_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (gpsref)
   !               call da_ytoyhat_gpsref_adj    (iv, qhat(itype)%obs, jo_grad_y)
   !            case (ssmi_tb)
   !               call da_ytoyhat_ssmi_tb_adj   (iv, qhat(itype)%obs, jo_grad_y)
   !            case (ssmi_rv)
   !               call da_ytoyhat_ssmi_rv_adj   (iv, qhat(itype)%obs, jo_grad_y)
   !            case (ssmt1)
   !               call da_ytoyhat_ssmt1_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (ssmt2)
   !               call da_ytoyhat_ssmt2_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (pseudo)
   !               call da_ytoyhat_pseudo_adj    (iv, qhat(itype)%obs, jo_grad_y)
   !            case (bogus)
   !               call da_ytoyhat_bogus_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (qscat)
   !               call da_ytoyhat_qscat_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (radar)
   !               call da_ytoyhat_radar_adj     (iv, qhat(itype)%obs, jo_grad_y)
   !            case (mtgirs)
   !               call da_ytoyhat_mtgirs_adj    (iv, qhat(itype)%obs, jo_grad_y)
   !            case (tamdar)
   !               call da_ytoyhat_tamdar_adj    (iv, qhat(itype)%obs, jo_grad_y)
   !            case (tamdar_sfc)
   !               call da_ytoyhat_tamdar_sfc_adj (iv, qhat(itype)%obs, jo_grad_y)
   !            case (rain)
   !               call da_ytoyhat_rain_adj      (iv, qhat(itype)%obs, jo_grad_y)
   !
   !!#if defined(CRTM) || defined(RTTOV)
   !!   if  (iv%num_inst               > 0) call da_ytoyhat_rad_adj      (iv, jo_grad_y)
   !!#endif
   !            case (airsr)
   !               call da_ytoyhat_airsr_adj     (iv, qhat(itype)%obs, jo_grad_y)
   
   #if (WRF_CHEM == 1)
               case (chem_surf)
                  call da_ytoyhat_chem_surf_adj  (iv, qhat(itype)%obs, jo_grad_y)
               case (chem_acft)
                  call da_ytoyhat_chem_acft_adj  (iv, qhat(itype)%obs, jo_grad_y)
   #endif
            end select
   
   !#if defined(CRTM) || defined(RTTOV)
   !      else if ( itype .eq. radiance .and. iv%num_inst .gt. 0 ) then
   !              call da_ytoyhat_rad_adj        (iv, qhat(itype)%obs, jo_grad_y)
   !#endif
   
         end if
      end do
   
   end if

   if (trace_use) call da_trace_exit("da_transform_ytoyhat_adj")

end subroutine da_transform_ytoyhat_adj


