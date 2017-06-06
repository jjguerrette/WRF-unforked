module da_chem_tools

!---------------------------------------------------------------------------
! Purpose: Collection of routines useful for chem inversions
!---------------------------------------------------------------------------

   use module_dm, only : wrf_dm_sum_reals, wrf_dm_sum_real, &
         wrf_dm_sum_integer
   use da_reporting, only : da_message,da_warning,message
#if (WRF_CHEM == 1)
   use module_domain, only : domain
   use module_state_description, only : &
         PARAM_FIRST_SCALAR, num_moist, num_chem, &
         num_scaleant, num_scalebb, &
         num_chem_acft

   use da_control, only : eta_emiss, trace_use, max_ext_its, stdout, &
         checkpoint_interval, write_checkpoints, num_qcstat_conv, &
         test_dm_exact, jb_factor, rootproc, ierr, comm, cv_size_domain, &
         num_ant_steps, num_bb_steps, use_nonchemobs, &
         osse_chem_rel, osse_chem_abs, &
         myproc, filename_len, read_omega, &
         num_ant_steps, num_bb_steps, &
         missing_r, &
         num_platform, num_ts, &
         its, ite, jts, jte, kts, kte, &
         ims, ime, jms, jme, kms, kme, &
         ids, ide, jds, jde, kds, kde 

!         pinterp_option, num_acft_obs, num_surf_obs, 
!         num_ts, num_platform, &&
!         ims, ime, jms, jme, kms, kme
   use da_chem, only: da_retrieve_chem_hx
   use da_define_structures, only : &
         iv_type, y_type, j_type, be_type, xbx_type, &
         da_allocate_y, da_allocate_y_chem, da_deallocate_y, &
         da_gauss_noise, da_random_seed
!   use da_obs, only: da_fill_obs_structures_chem
   use da_obs_io, only: da_read_obs_chem_multiplat
#endif
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tools, only : da_set_randomcv
   use da_tracing, only: da_trace_entry, da_trace_exit
   use da_vtox_transforms, only : da_transform_vtox
#if (WRF_CHEM == 1)
   use da_transfer_model, only: da_transfer_xatowrf_temp
#endif
   use da_minimisation, only:   da_get_innov_vector, &
                             da_calculate_residual, da_jo_and_grady
   use da_randomisation, only:   da_cv_io
   use da_par_util, only : da_cv_to_global
   use da_wrf_interfaces, only : wrf_dm_bcast_real
!!WRFPLUS domain type
!   use module_wrf_top, only: domain
!WRFPLUS model grid is usable only after da_nl_model(1) has been called
   use da_4dvar, only: &
#if (WRF_CHEM == 1)
      da_init_model_input, &
#endif
      model_grid, kj_swap, da_nl_model
!WRFDA domain type
   use module_configure, only : grid_config_rec_type

   implicit none

#ifdef DM_PARALLEL
    include 'mpif.h'
#endif
#if (WRF_CHEM == 1)
   private :: da_dot, da_dot_cv
#endif

contains

#if (WRF_CHEM == 1)

#include "da_hdgn.inc"
#include "da_dgn.inc"
#include "da_evaluate_j.inc"
#include "da_dot.inc"
#include "da_dot_cv.inc"
#include "da_setup_osse_chem.inc"


#endif

end module da_chem_tools
