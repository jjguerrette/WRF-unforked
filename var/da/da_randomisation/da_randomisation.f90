module da_randomisation

   !---------------------------------------------------------------------------
   ! Purpose: Collection of routines associated with randomisation. 
   !---------------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use module_domain, only : domain !&
   use module_state_description, only : &
#if (WRF_CHEM == 1)
      num_scaleant, num_scalebb, &
#endif
      PARAM_FIRST_SCALAR

   use da_control, only : svd_stage, ensmember, ensdim_svd, svd_outer, &
       myproc, filename_len, test_dm_exact, rootproc, cv_size_domain, &
       stdout, trace_use, adapt_svd, &
       num_ob_indexes, read_omega, svd_p, ierr, comm, &
       use_randomsvd, nens_compare, tsvd_compare, cglz_compare, rsvd_compare, &
       riot_precon, ntmax, &
       hess_comp_x, &
       hess_comp_Pav_LRU, hess_comp_Pav_LRA, &
       hess_comp_Pa_LRU, hess_comp_Pa_LRA, &
       hess_comp_PaPa_LRU, hess_comp_PaPa_LRA, &
       scale_hcheck, &
#if (WRF_CHEM == 1)
       chem_surf, chem_acft, &
       num_ant_steps, num_bb_steps, &
#endif
       sound, mtgirs, sonde_sfc, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
       satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep,tamdar, tamdar_sfc, rain, &
       bogus, buoy, qscat, pseudo, radiance, &
       its, ite, jts, jte

   use da_minimisation, only: da_transform_vtoy, da_transform_vtoy_adj, &
       da_calculate_grady, da_calculate_j, da_calculate_gradj, &
       da_amat_mul
   use da_vtox_transforms, only : da_transform_vtox, da_transform_vtox_adj
   use da_define_structures, only : iv_type, y_type, j_type, be_type, xbx_type, &
#if defined(LAPACK)
       yhat_type, &
#endif
#if (WRF_CHEM == 1)
      da_allocate_y_chem, &
#endif
      da_allocate_y, da_deallocate_y
   use da_par_util, only : da_cv_to_global, da_global_to_cv, da_cv_to_vv, da_vv_to_cv
   use da_reporting, only : da_message, da_warning, da_error, message
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tools, only: da_set_randomcv
   use da_tracing, only : da_trace_entry, da_trace_exit,da_trace
!#ifdef VAR4D
!#endif
   use da_lapack, only : dsteqr
#if defined(LAPACK)
!   use mkl95_precision, only: WP => DP
!   use mkl95_lapack, only: gesv, syev, gesvd
   use f95_precision, only: WP => DP
   use lapack95, only: gesv, syev, gesvd
!   use blas95, only: gemm

!   use da_airep, only : da_ytoyhat_airep, da_ytoyhat_airep_adj
!   use da_airsr , only : da_ytoyhat_airsr, da_ytoyhat_airsr_adj
!   use da_bogus, only : da_ytoyhat_bogus, da_ytoyhat_bogus_adj
!   use da_buoy , only : da_ytoyhat_buoy, da_ytoyhat_buoy_adj
!   use da_geoamv, only : da_ytoyhat_geoamv, da_ytoyhat_geoamv_adj
!   use da_gpspw, only : da_ytoyhat_gpspw, da_ytoyhat_gpspw_adj
!   use da_gpsref, only : da_ytoyhat_gpsref, da_ytoyhat_gpsref_adj
!   use da_metar, only : da_ytoyhat_metar, da_ytoyhat_metar_adj
!   use da_pilot, only : da_ytoyhat_pilot, da_ytoyhat_pilot_adj
!   use da_polaramv, only : da_ytoyhat_polaramv, da_ytoyhat_polaramv_adj
!   use da_profiler, only : da_ytoyhat_profiler, da_ytoyhat_profiler_adj
!   use da_pseudo, only : da_ytoyhat_pseudo, da_ytoyhat_pseudo_adj
!   use da_qscat, only : da_ytoyhat_qscat, da_ytoyhat_qscat_adj
!   use da_mtgirs, only : da_ytoyhat_mtgirs, da_ytoyhat_mtgirs_adj
!   use da_tamdar, only : da_ytoyhat_tamdar, da_ytoyhat_tamdar_adj, &
!      da_ytoyhat_tamdar_sfc, da_ytoyhat_tamdar_sfc_adj
!   use da_radiance, only : da_ytoyhat_rad, da_ytoyhat_rad_adj
!   use da_radar, only :  da_ytoyhat_radar, da_ytoyhat_radar_adj
!   use da_rain, only :  da_ytoyhat_rain, da_ytoyhat_rain_adj
!   use da_satem, only : da_ytoyhat_satem, da_ytoyhat_satem_adj
!   use da_ships, only : da_ytoyhat_ships, da_ytoyhat_ships_adj
!   use da_sound, only : da_ytoyhat_sound, da_ytoyhat_sound_adj
!   use da_ssmi, only : da_ytoyhat_ssmi_tb, da_ytoyhat_ssmi_tb_adj, &
!      da_ytoyhat_ssmi_rv, da_ytoyhat_ssmi_rv_adj, &
!      da_ytoyhat_ssmt1, da_ytoyhat_ssmt1_adj, &
!      da_ytoyhat_ssmt2, da_ytoyhat_ssmt2_adj
!   use da_synop, only : da_ytoyhat_synop, da_ytoyhat_synop_adj
#if (WRF_CHEM == 1)
   use da_chem, only: da_ytoyhat_chem_surf, da_ytoyhat_chem_surf_adj, &
      da_ytoyhat_chem_acft, da_ytoyhat_chem_acft_adj, &
      da_force_grad_chem_acft, da_force_grad_chem_surf

#endif

#endif

   implicit none

#ifdef DM_PARALLEL
    include 'mpif.h'
#endif

   private :: da_dot, da_dot_cv, da_dot_z, da_dot_cv_z, da_dot_obs

contains

#include "da_dot.inc"
#include "da_dot_cv.inc"     
#include "da_dot_z.inc"
#include "da_dot_cv_z.inc"
#include "da_dot_obs.inc"
#include "da_cv_io.inc"
#include "da_cv_io_int.inc"
#include "da_yhat_io.inc"
#include "da_hessian_io.inc"
#include "da_precon_omega.inc"
#if defined(LAPACK)
#include "da_transform_ytoyhat.inc"
#include "da_transform_ytoyhat_adj.inc"
#endif
#include "da_amat_mul_trunc.inc"
#include "da_randomise_svd.inc"
#include "da_randomise_svd_51.inc"
#include "da_force_grad_hess.inc"
#include "da_gram_schmidt.inc"
#include "da_rsvd56.inc"

#include "da_calculate_hessian.inc"
#include "da_compare_decomp_methods.inc"
#include "da_evaluate_decomp.inc"
#include "da_evaluate_decomp_init.inc"
#include "da_write_norms.inc"
#include "da_evaluate_increment.inc"
#include "da_evaluate_hessian.inc"
#include "da_gen_omega.inc"

#include "da_randomise_svd_B.inc"
#include "da_randomise_svd_B11.inc"



end module da_randomisation
