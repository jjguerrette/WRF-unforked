module da_linear_ops

   !---------------------------------------------------------------------------
   ! Purpose: Collection of routines associated with randomisation. 
   !---------------------------------------------------------------------------

   use module_configure, only : grid_config_rec_type
   use module_dm, only : wrf_dm_sum_real, wrf_dm_sum_integer
   use module_domain, only : domain, domain_clock_get
   use module_timing

   use da_control, only : &
#if defined(LAPACK)
       use_randomblock, &
#endif
       myproc, filename_len, test_dm_exact, rootproc, cv_size_domain, &
       stdout, trace_use, ierr, comm, &
       use_lanczos, use_global_cv_io, ntmax, inc_out_interval, &
       read_hess_REF, nmodes_hess_REF, &
       spectral_precon
   use da_define_structures, only : iv_type, y_type, j_type, be_type, xbx_type, &
       hessian_type, hessian_eig_type

   use da_par_util, only : da_cv_to_global, da_global_to_cv, da_cv_to_vv, da_vv_to_cv
   use da_vtox_transforms, only : da_transform_vtox
   use da_transfer_model, only : da_transfer_xatowrftl
   use da_reporting, only : da_message, da_warning, da_error, message
   use da_tools_serial, only : da_get_unit,da_free_unit
   use da_tools, only: da_set_randomcv
   use da_tracing, only : da_trace_entry, da_trace_exit,da_trace

   implicit none

#ifdef DM_PARALLEL
    include 'mpif.h'
#endif

!   private :: da_dot, da_dot_cv

contains

#include "da_dot.inc"
#include "da_dot_cv.inc"     

#include "da_cv_io.inc"
#include "da_cv_io_global.inc"
#include "da_hessian_io.inc"
#include "da_hessian_io_global.inc"

#include "da_gram_schmidt.inc"
#include "da_amat_mul_trunc.inc"
#include "da_spectral_precon.inc"
#include "da_output_increments.inc"

end module da_linear_ops
