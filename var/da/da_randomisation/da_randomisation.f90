module da_randomisation

   !---------------------------------------------------------------------------
   ! Purpose: Collection of routines associated with randomisation. 
   !---------------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use module_domain, only : domain !&
   use module_timing
   use module_state_description, only : &
#if (WRF_CHEM == 1)
      num_scaleant, num_scalebb, &
#endif
      PARAM_FIRST_SCALAR

   use da_control, only : ensmember, &
       rand_outer, rand_stage, rand_inner_it, max_rand_inner, &
       myproc, filename_len, rootproc, &
       var4d_lbc, stdout, trace_use, adapt_svd, &
       read_omega, svd_p, ierr, comm, &
       use_global_cv_io, nmodes_global, var4d_inc_out, inc_out_interval, &
       spectral_precon, rotate_omega, spectral_trunc, &
       enforce_posdef, &
#if (WRF_CHEM == 1)
       chem_surf, chem_acft, &
       num_ant_steps, num_bb_steps, &
#endif
       sound, mtgirs, sonde_sfc, synop, profiler, gpsref, gpspw, polaramv, geoamv, ships, metar, &
       satem, radar, ssmi_rv, ssmi_tb, ssmt1, ssmt2, airsr, pilot, airep,tamdar, tamdar_sfc, rain, &
       bogus, buoy, qscat, pseudo, radiance, &
       its, ite, jts, jte
   use da_minimisation, only: da_transform_vtoy, da_transform_vtoy_adj, &
       da_calculate_grady, da_calculate_j, da_calculate_gradj
   use da_linear_ops, only:  da_gram_schmidt, da_amat_mul_trunc, &
       da_cv_io, da_cv_io_global, da_spectral_precon, da_hessian_io, &
       da_dot_cv, da_dot, da_mat_io
   use da_vtox_transforms, only : da_transform_vtox, da_transform_vtox_adj
   use da_define_structures, only : iv_type, y_type, j_type, be_type, xbx_type, &
       hessian_type
   use da_par_util, only : da_cv_to_global, da_global_to_cv, da_cv_to_vv, da_vv_to_cv
   use da_reporting, only : da_message, da_warning, da_error, message
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tools, only: da_set_randomcv
   use da_tracing, only : da_trace_entry, da_trace_exit,da_trace
   use da_lapack, only : dsteqr !WRFDA-specific version of dsteqr (not MKL)

#if defined(LAPACK)
   use f95_precision, only: WP => DP
   use lapack95, only: gesv, syev, gesvd
#endif

   implicit none

#ifdef DM_PARALLEL
    include 'mpif.h'
#endif

!   private :: da_dot_obs

   character(len=14) :: randvecdir

contains

#include "da_gen_omega.inc"
#include "da_rotate_omega.inc"
#include "da_randomise_svd.inc"
#include "da_rsvd56.inc"
#include "da_block_lanczos.inc"

#if (WRF_CHEM == 1)
#include "da_evaluate_increment.inc"
#include "da_evaluate_hessian.inc"
#endif

end module da_randomisation
