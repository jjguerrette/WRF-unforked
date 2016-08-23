module da_randomisation

   !---------------------------------------------------------------------------
   ! Purpose: Collection of routines associated with randomisation. 
   !---------------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use module_domain, only : domain !&
   use module_state_description, only : PARAM_FIRST_SCALAR      
   use da_control, only : svd_stage, ensmember, ensdim_svd, svd_outer, &
       myproc, filename_len, test_dm_exact, rootproc, cv_size_domain, &
       stdout, trace_use, svd_amat_type, svd_symm_type
   use da_minimisation, only: da_transform_vtoy, da_transform_vtoy_adj, &
       da_calculate_grady, da_calculate_j, da_calculate_gradj, &
       da_amat_mul
   use da_define_structures, only : iv_type, y_type, j_type, be_type, xbx_type, &
#if (WRF_CHEM == 1)
      da_allocate_y_chem, &
#endif
      da_allocate_y, da_deallocate_y
   use da_par_util, only : da_cv_to_global
   use da_reporting, only : da_message, da_warning, da_error, message
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tools, only: da_set_randomcv
   use da_tracing, only : da_trace_entry, da_trace_exit,da_trace
!#ifdef VAR4D
!#endif
#if defined(LAPACK)
!   use mkl95_precision, only: WP => DP
!   use mkl95_lapack, only: gesv, geev
   use f95_precision, only: WP => DP
   use lapack95, only: gesv, geev, syev
   use blas95, only: gemm
#endif

   implicit none

#ifdef DM_PARALLEL
    include 'mpif.h'
#endif

   private :: da_dot, da_dot_cv, da_dot_z, da_dot_cv_z

contains

#include "da_dot.inc"
#include "da_dot_cv.inc"     
#include "da_dot_z.inc"
#include "da_dot_cv_z.inc"
#include "da_randomise_svd.inc"
#include "da_cv_io.inc"

end module da_randomisation
