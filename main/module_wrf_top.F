!WRF:DRIVER_LAYER:TOP
!

!TBH:  $$$  move this to ../frame?  

MODULE module_wrf_top
!<DESCRIPTION>
! This module defines top-level wrf_init(), wrf_adtl_check(), wrf_run(), and wrf_finalize() 
! routines.  
!</DESCRIPTION>

   USE module_machine
   USE module_domain
   USE module_integrate
   USE module_driver_constants
   USE module_configure
   USE module_check_a_mundo
   USE module_streams

   USE module_timing
   USE module_wrf_error
   USE module_adtl_grid_utilities, ONLY : allocate_grid, deallocate_grid, copy_grid_to_s, copy_grid_to_b, &
                                          restore_grid, inner_dot_g, add_grid, inner_dot, init_domain_size, &
                                          copy_s_to_g_adjtest, copy_g_to_b_adjtest, inner_dot_g_adjtest, &
                                          copy_g_to_a_adjtest, inner_dot_a_b_adjtest, zero_out_ad, zero_out_tl, &
                                          get_forc_locations, allocate_locations, gen_scenario_matrix, &
                                          spot_force_ad, spot_force_tl, spot_force_nl, &
                                          extract_ad_der, extract_tl_der, extract_nl_den, extract_nl_num
   USE mediation_pertmod_io, ONLY : xtraj_pointer, xtraj_head, xtraj_tail, xtraj_io_initialize, adtl_initialize
   USE module_nesting
   USE module_io_domain, ONLY: close_dataset, wrf_inquire_opened
   USE module_io, ONLY: io_form_for_dataset
#if (EM_CORE==1)
   USE module_state_description, ONLY : SURFDRAGSCHEME, ACMPBLSCHEME, &
                                        SFCLAYSCHEME, PXSFCSCHEME, &
                                        SLABSCHEME, PXLSMSCHEME, &
                                        PARAM_FIRST_SCALAR, num_moist, &
#if (WRF_CHEM == 1)
                                        num_chem, num_scaleant, num_scalebb, &
#endif
                                        num_tracer
#endif

#ifdef DM_PARALLEL
   USE module_dm, ONLY : wrf_dm_initialize,wrf_get_hostid,domain_active_this_task,mpi_comm_allcompute,wrf_dm_sum_real, &
                         local_communicator
#else
   USE module_dm, ONLY : domain_active_this_task
#endif

   USE module_date_time, ONLY : current_date, start_date
   USE module_io_domain, ONLY : open_w_dataset, output_boundary
   USE module_cpl, ONLY : coupler_on, cpl_finalize, cpl_defdomain
   USE esmf_timemod, ONLY : ESMF_TimeInc

   IMPLICIT NONE

#include <wrf_io_flags.h>

   REAL    :: time

   INTEGER :: loop , &
              levels_to_process

   TYPE (domain) , POINTER :: keep_grid, grid_ptr, null_domain, grid
   TYPE (domain) , pointer :: parent_grid, new_nest
   LOGICAL                                :: a_nest_was_opened
   TYPE (grid_config_rec_type), SAVE :: config_flags
   INTEGER        :: kid, nestid
   INTEGER                 :: number_at_same_level
   INTEGER                 :: time_step_begin_restart

   INTEGER :: max_dom , domain_id , fid , oid , idum1 , idum2 , ierr
   INTEGER :: debug_level, alarmid
   INTEGER :: xtraj_boundary_form=2, input_boundary_form

   LOGICAL :: input_from_file

#ifdef DM_PARALLEL
   INTEGER                 :: nbytes
   INTEGER, PARAMETER      :: configbuflen = 4* CONFIG_BUF_LEN
   INTEGER                 :: configbuf( configbuflen )
   LOGICAL , EXTERNAL      :: wrf_dm_on_monitor
#endif

   CHARACTER (LEN=256)     :: rstname
   CHARACTER (LEN=80)      :: message
   LOGICAL                 :: gradient_out = .FALSE.
   CHARACTER (LEN=256) , PRIVATE :: a_message

   INTERFACE 
     SUBROUTINE Setup_Timekeeping( grid )
      USE module_domain
      TYPE(domain), POINTER :: grid
     END SUBROUTINE Setup_Timekeeping

! #if (EM_CORE == 1)
     SUBROUTINE wrf_dfi_write_initialized_state( )
     END SUBROUTINE wrf_dfi_write_initialized_state
 
     SUBROUTINE wrf_dfi_startfwd_init( )
     END SUBROUTINE wrf_dfi_startfwd_init
     
     SUBROUTINE wrf_dfi_startbck_init( )
     END SUBROUTINE wrf_dfi_startbck_init
     
     SUBROUTINE wrf_dfi_bck_init( )
     END SUBROUTINE wrf_dfi_bck_init
     
     SUBROUTINE wrf_dfi_fwd_init( )
     END SUBROUTINE wrf_dfi_fwd_init
     
     SUBROUTINE wrf_dfi_fst_init( )
     END SUBROUTINE wrf_dfi_fst_init
     
     SUBROUTINE wrf_dfi_array_reset ( )
     END SUBROUTINE wrf_dfi_array_reset
! #endif

     SUBROUTINE med_nest_initial ( parent , grid , config_flags )
       USE module_domain
       USE module_configure
       TYPE (domain), POINTER ::  grid , parent
       TYPE (grid_config_rec_type) config_flags
     END SUBROUTINE med_nest_initial

   END INTERFACE


CONTAINS


   SUBROUTINE wrf_init( no_init1 )
!<DESCRIPTION>
!     WRF initialization routine.
!</DESCRIPTION>
#ifdef _OPENMP
     use omp_lib
#endif
#ifdef _ACCEL
     use accel_lib
#endif
     LOGICAL, OPTIONAL, INTENT(IN) :: no_init1
     INTEGER i, myproc, nproc, hostid, loccomm, ierr, buddcounter, mydevice, save_comm
     INTEGER, ALLOCATABLE :: hostids(:), budds(:)
     CHARACTER*512 hostname
     CHARACTER*512 mminlu_loc
#ifdef _ACCEL
     INTEGER :: it, nt, in, devnum
#endif
#if defined(DM_PARALLEL) && !defined(STUBMPI) && ( defined(RUN_ON_GPU) || defined(_ACCEL))
     include "mpif.h"
#endif
#include "version_decl"


!<DESCRIPTION>
! Program_name, a global variable defined in frame/module_domain.F, is
! set, then a routine <a href=init_modules.html>init_modules</a> is
! called. This calls all the init programs that are provided by the
! modules that are linked into WRF.  These include initialization of
! external I/O packages.   Also, some key initializations for
! distributed-memory parallelism occur here if DM_PARALLEL is specified
! in the compile: setting up I/O quilt processes to act as I/O servers
! and dividing up MPI communicators among those as well as initializing
! external communication packages such as RSL or RSL_LITE.
!
!</DESCRIPTION>

   program_name = "WRF " // TRIM(release_version) // " MODEL"

   ! Initialize WRF modules:  
   ! Phase 1 returns after MPI_INIT() (if it is called)
   CALL init_modules(1)
   IF ( .NOT. PRESENT( no_init1 ) ) THEN
     ! Initialize utilities (time manager, etc.)
#ifdef NO_LEAP_CALENDAR
     CALL WRFU_Initialize( defaultCalKind=WRFU_CAL_NOLEAP )
#else
     CALL WRFU_Initialize( defaultCalKind=WRFU_CAL_GREGORIAN )
#endif
   ENDIF
   ! Phase 2 resumes after MPI_INIT() (if it is called)
   CALL init_modules(2)

!<DESCRIPTION>
! The wrf namelist.input file is read and stored in the USE associated
! structure model_config_rec, defined in frame/module_configure.F, by the
! call to <a href=initial_config.html>initial_config</a>.  On distributed
! memory parallel runs this is done only on one processor, and then
! broadcast as a buffer.  For distributed-memory, the broadcast of the
! configuration information is accomplished by first putting the
! configuration information into a buffer (<a
! href=get_config_as_buffer.html>get_config_as_buffer</a>), broadcasting
! the buffer, then setting the configuration information (<a
! href=set_config_as_buffer.html>set_config_as_buffer</a>).
!
!</DESCRIPTION>

#ifdef DM_PARALLEL
   CALL wrf_get_dm_communicator( save_comm )
   CALL wrf_set_dm_communicator( mpi_comm_allcompute )
   IF ( wrf_dm_on_monitor() ) THEN
     CALL initial_config
   ENDIF
   CALL get_config_as_buffer( configbuf, configbuflen, nbytes )
   CALL wrf_dm_bcast_bytes( configbuf, nbytes )
   CALL set_config_as_buffer( configbuf, configbuflen )
   CALL wrf_dm_initialize
   CALL wrf_set_dm_communicator( save_comm )
#else
   CALL initial_config
#endif

   CALL set_derived_rconfigs
   CALL check_nml_consistency
   CALL set_physics_rconfigs

#ifdef _ACCEL
   buddcounter = 1
   mydevice = 0
# if defined(DM_PARALLEL) && !defined(STUBMPI) 
   CALL wrf_get_myproc( myproc )
   CALL wrf_get_nproc( nproc )
   CALL wrf_get_hostid ( hostid )
   CALL wrf_get_dm_communicator ( loccomm )

   ALLOCATE( hostids(nproc) )
   ALLOCATE( budds(nproc) )
   CALL mpi_allgather( hostid, 1, MPI_INTEGER, hostids, 1, MPI_INTEGER, loccomm, ierr )
   if ( ierr .NE. 0 ) print * ,'error in mpi_allgather ',ierr
   budds = -1
   buddcounter = 0
   ! mark the ones i am on the same node with
   DO i = 1, nproc
      IF ( hostid .EQ. hostids(i) ) THEN
         budds(i) = buddcounter
         buddcounter = buddcounter + 1
      ENDIF
   ENDDO
   mydevice = budds(myproc+1)
   DEALLOCATE( hostids )
   DEALLOCATE( budds )
# endif
   in = acc_get_num_devices(acc_device_nvidia)
   if (in .le. 0) print *, 'error:  No GPUS present: ',in
# ifdef _OPENMP
   !$OMP PARALLEL SHARED(mydevice,in) PRIVATE(it,nt,devnum)
   it = omp_get_thread_num()
   nt = omp_get_num_threads()
   devnum = mod(mod(mydevice*nt,in) + it, in)
# ifdef _ACCEL_PROF
   print *, "Process, Thread, Device: ",mydevice, it, devnum
# endif
   call acc_set_device_num(devnum, acc_device_nvidia)

   !$OMP END PARALLEL
# else
   it = 0
   nt = 1
   devnum = mod(mod(mydevice*nt,in) + it, in)
#  ifdef _ACCEL_PROF
   print *, "Process, Thread, Device: ",mydevice, it, devnum
#  endif
   call acc_set_device_num(devnum, acc_device_nvidia)
# endif
#endif

#ifdef RUN_ON_GPU
   CALL wrf_get_myproc( myproc )
   CALL wrf_get_nproc( nproc )
# ifdef DM_PARALLEL
   CALL wrf_get_hostid ( hostid ) 
   CALL wrf_get_dm_communicator ( loccomm )
   ALLOCATE( hostids(nproc) )
   ALLOCATE( budds(nproc) )
   CALL mpi_allgather( hostid, 1, MPI_INTEGER, hostids, 1, MPI_INTEGER, loccomm, ierr )
   IF ( ierr .NE. 0 ) THEN
      write(a_message,*)__FILE__,__LINE__,'error in mpi_allgather ',ierr
      CALL wrf_message ( a_message ) 
   END IF
   budds = -1
   buddcounter = 0 
   ! mark the ones i am on the same node with
   DO i = 1, nproc 
      IF ( hostid .EQ. hostids(i) ) THEN
         budds(i) = buddcounter 
         buddcounter = buddcounter + 1
      ENDIF
   ENDDO
   mydevice = budds(myproc+1)
   DEALLOCATE( hostids )
   DEALLOCATE( budds )
# else
   mydevice = 0
# endif
   CALL wsm5_gpu_init( myproc, nproc, mydevice )
#endif

!<DESCRIPTION>
! Among the configuration variables read from the namelist is
! debug_level. This is retrieved using nl_get_debug_level (Registry
! generated and defined in frame/module_configure.F).  The value is then
! used to set the debug-print information level for use by <a
! href=wrf_debug.html>wrf_debug</a> throughout the code. Debug_level
! of zero (the default) causes no information to be printed when the
! model runs. The higher the number (up to 1000) the more information is
! printed.
! 
!</DESCRIPTION>

   CALL nl_get_debug_level ( 1, debug_level )
   CALL set_wrf_debug_level ( debug_level )

   ! allocated and configure the mother domain

   NULLIFY( null_domain )

!<DESCRIPTION>
! RSL is required for WRF nesting options.
! The non-MPI build that allows nesting is only supported on machines
! with the -DSTUBMPI option.  Check to see if the WRF model is being asked 
! for a for a multi-domain run (max_dom > 1, from the namelist).  If so,
! then we check to make sure that we are under the parallel
! run option or we are on an acceptable machine.
!</DESCRIPTION>

   CALL nl_get_max_dom( 1, max_dom )
   IF ( max_dom > 1 ) THEN
#if ( ! defined(DM_PARALLEL)  &&   ! defined(STUBMPI) )
   CALL wrf_error_fatal( &
     'nesting requires either an MPI build or use of the -DSTUBMPI option' ) 
#endif
   END IF

!<DESCRIPTION>
! The top-most domain in the simulation is then allocated and configured
! by calling <a href=alloc_and_configure_domain.html>alloc_and_configure_domain</a>.
! Here, in the case of this root domain, the routine is passed the
! globally accessible pointer to TYPE(domain), head_grid, defined in
! frame/module_domain.F.  The parent is null and the child index is given
! as negative, signifying none.  Afterwards, because the call to
! alloc_and_configure_domain may modify the model's configuration data
! stored in model_config_rec, the configuration information is again
! repacked into a buffer, broadcast, and unpacked on each task (for
! DM_PARALLEL compiles). The call to <a
! href=setup_timekeeping.html>setup_timekeeping</a> for head_grid relies
! on this configuration information, and it must occur after the second
! broadcast of the configuration information.
! 
!</DESCRIPTION>
   CALL       wrf_message ( program_name )
   CALL       wrf_debug ( 100 , 'wrf: calling alloc_and_configure_domain ' )
   CALL alloc_and_configure_domain ( domain_id  = 1 ,                  &
                               active_this_task = domain_active_this_task(1), &
                                     grid       = head_grid ,          &
                                     parent     = null_domain ,        &
                                     kid        = -1                   )

   CALL       wrf_debug ( 100 , 'wrf: calling model_to_grid_config_rec ' )
   CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
   CALL       wrf_debug ( 100 , 'wrf: calling set_scalar_indices_from_config ' )
   CALL set_scalar_indices_from_config ( head_grid%id , idum1, idum2 )
   CALL       wrf_debug ( 100 , 'wrf: calling init_wrfio' )
   CALL init_wrfio

#ifdef DM_PARALLEL
   CALL wrf_get_dm_communicator( save_comm )
   CALL wrf_set_dm_communicator( mpi_comm_allcompute )
   CALL get_config_as_buffer( configbuf, configbuflen, nbytes )
   CALL wrf_dm_bcast_bytes( configbuf, nbytes )
   CALL set_config_as_buffer( configbuf, configbuflen )
   CALL wrf_set_dm_communicator( save_comm )
#endif

! #if (EM_CORE == 1)
   ! In case we are doing digital filter initialization, set dfi_stage = DFI_SETUP 
   !   to indicate in Setup_Timekeeping that we want forecast start and
   !   end times at this point 
   IF ( head_grid%dfi_opt .NE. DFI_NODFI ) head_grid%dfi_stage = DFI_SETUP
! #endif

   CALL Setup_Timekeeping (head_grid)

!<DESCRIPTION>
! The head grid is initialized with read-in data through the call to <a
! href=med_initialdata_input.html>med_initialdata_input</a>, which is
! passed the pointer head_grid and a locally declared configuration data
! structure, config_flags, that is set by a call to <a
! href=model_to_grid_config_rec.html>model_to_grid_config_rec</a>.  It is
! also necessary that the indices into the 4d tracer arrays such as
! moisture be set with a call to <a
! href=set_scalar_indices_from_config.html>set_scalar_indices_from_config</a>
! prior to the call to initialize the domain.  Both of these calls are
! told which domain they are setting up for by passing in the integer id
! of the head domain as <tt>head_grid%id</tt>, which is 1 for the
! top-most domain.
! 
! In the case that write_restart_at_0h is set to true in the namelist,
! the model simply generates a restart file using the just read-in data
! and then shuts down. This is used for ensemble breeding, and is not
! typically enabled.
! 
!</DESCRIPTION>

 IF ( domain_active_this_task(1) ) THEN
   CALL med_initialdata_input( head_grid , config_flags )

   IF ( config_flags%write_restart_at_0h ) THEN
      CALL med_restart_out ( head_grid, config_flags )
#ifndef AUTODOC_BUILD
! prevent this from showing up before the call to integrate in the autogenerated call tree
      CALL wrf_debug ( 0 , ' 0 h restart only wrf: SUCCESS COMPLETE WRF' )
! TBH:  $$$ Unscramble this later...  
! TBH:  $$$ Need to add state to avoid calling wrf_finalize() twice when ESMF 
! TBH:  $$$ library is used.  Maybe just set clock stop_time=start_time and 
! TBH:  $$$ do not call wrf_finalize here...  
      CALL wrf_finalize( )
#endif
   END IF
  ENDIF  ! domain_active_this_task

   ! set default values for subtimes
   head_grid%start_subtime = domain_get_start_time ( head_grid )
   head_grid%stop_subtime = domain_get_stop_time ( head_grid )

 IF ( domain_active_this_task(1) ) THEN
   !  For EM (but not DA), if this is a DFI run, we can allocate some space.  We are
   !  not allowing anyting tricky for nested DFI.  If there are any nested domains,
   !  they all need to start at the same time.  Otherwise, why even do the DFI?  If
   !  the domains do not all start at the same time, then there will be inconsistencies,
   !  which is what DFI is supposed to address.

#if (EM_CORE == 1)
   IF ( head_grid%dfi_opt .NE. DFI_NODFI ) THEN
      CALL alloc_doms_for_dfi ( head_grid )
   END IF
#endif

   IF (coupler_on) CALL cpl_defdomain( head_grid ) 
  ENDIF  ! domain_active_this_task

   END SUBROUTINE wrf_init

   SUBROUTINE wrf_adtl_check( )

   IF ( config_flags%scenario_type .EQ. 0 ) THEN
      CALL wrf_adtl_check_sum
   ELSE
      CALL wrf_adtl_check_spot
   END IF

   END SUBROUTINE wrf_adtl_check


   SUBROUTINE wrf_adtl_check_sum( )
!<DESCRIPTION>
!     WRF adjoint and tangent linear code check routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  
! 
!</DESCRIPTION>

   REAL :: alpha_m, factor, val_l, val_a, save_l, val_n, coef
   INTEGER :: nt, ij, time_step, rc
   CHARACTER*256 :: timestr

   ! Return immediately if not dyn_em_check
   IF ( config_flags%dyn_opt .NE. dyn_em_check ) RETURN

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   ! Force to read the lateral boundary condition here.
   CALL med_before_solve_io ( head_grid, config_flags )

   ! Close boundary file, as it will be read again from start later
   CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )

   CALL init_domain_size ( head_grid, config_flags )

   ! Release linked list for trajectory
   call xtraj_io_initialize

IF ( config_flags%check_TL .or. config_flags%check_AD ) THEN

   ! Save the initial condition and boundary condition, x
   CALL allocate_grid ( )

   CALL copy_grid_to_s ( head_grid , &
                         head_grid%i_start(1), head_grid%i_end(1),                &
                         head_grid%j_start(1), head_grid%j_end(1)                 )

   CALL       wrf_message ( "wrf: calling nonlinear integrate" )
   model_config_rec%dyn_opt = dyn_em

   ! Set up basic states output
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, time_step )
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )
   IF ( config_flags%check_TL ) THEN
      IF ( config_flags%mp_physics .NE. 0 .and. config_flags%mp_physics .NE. 99 .and. &
           config_flags%mp_physics .NE. 98 ) CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
      IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. config_flags%bl_pbl_physics .GT. 0 ) CALL nl_set_bl_pbl_physics (head_grid%id, 98)
      IF ( config_flags%cu_physics .GT. 0 ) THEN
         CALL nl_set_cu_physics (head_grid%id, 98)
         head_grid%cudt = 0
      ENDIF
      CALL nl_set_ra_lw_physics (head_grid%id, 0)
      CALL nl_set_ra_sw_physics (head_grid%id, 0)
      IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
      IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. config_flags%sf_surface_physics.NE.PXLSMSCHEME) CALL nl_set_sf_surface_physics (head_grid%id, 0)

      ! Force to turn off boundary input as we can only perturb the initial boundary condition.
      CALL nl_set_io_form_boundary( head_grid%id, 0 )
   ENDIF

   CALL wrf_run

   ! Turn off basic states output
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )
   if ( .not. head_grid%trajectory_io ) CALL nl_set_inputout_interval_s ( head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL       wrf_message ( "wrf: back from nonlinear integrate" )

   IF ( config_flags%check_TL ) THEN

   wrf_err_message = '=============== Tangent Linear Check ==================='
   CALL wrf_message(TRIM(wrf_err_message))
   
   wrf_err_message = 'check==== U === V === W == PH === T == MU == MOIST == TRACER ==='
   CALL wrf_message(TRIM(wrf_err_message))

   WRITE(wrf_err_message,'(A,8(1x,L5))') 'check',head_grid%check_u, head_grid%check_v,  &
                               head_grid%check_w, head_grid%check_ph, &
                               head_grid%check_t, head_grid%check_mu, &
                               head_grid%check_moist, head_grid%check_tracer
#if (WRF_CHEM == 1)
   wrf_err_message = 'check==== CHEM == SCALEANT == SCALEBB  ==='
   CALL wrf_message(TRIM(wrf_err_message))

   WRITE(wrf_err_message,'(A,3(1x,L5))') 'check', head_grid%check_chem, &
                               head_grid%check_scaleant, head_grid%check_scalebb
#endif
   CALL wrf_message(TRIM(wrf_err_message))

   ! Save the f(x)
   CALL copy_grid_to_b ( head_grid )

   ! Restore the x and assign the \delta x
   CALL restore_grid ( head_grid )
   CALL copy_s_to_g_adjtest ( head_grid, 1.0 )

   ! Set up basic states reading
   model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxinput16_interval_s (head_grid%id, time_step )
   CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

   CALL wrf_run_tl

   ! Turn off auxinput16 reading
   CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   save_l = 0.0
   ! Calculate the inner dot.
   !$OMP PARALLEL DO    &
   !$OMP DEFAULT (SHARED) PRIVATE ( ij ) &           
   !$OMP REDUCTION (+:save_l) 
   DO ij = 1 , head_grid%num_tiles
     CALL inner_dot_g ( head_grid , save_l,                                        &
                        head_grid%i_start(ij), head_grid%i_end(ij),                &
                        head_grid%j_start(ij), head_grid%j_end(ij)                 )
   END DO
   !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
   save_l = wrf_dm_sum_real ( save_l )
#endif

   alpha_m = 1.
   tangentLinearCheck : DO nt = 1 , 11

      alpha_m = 0.1 * alpha_m
      factor = 1.0 + alpha_m

      CALL add_grid ( head_grid,  factor )

      CALL       wrf_message ( "wrf: calling nonlinear integrate" )
      model_config_rec%dyn_opt = dyn_em
      CALL domain_clock_get( head_grid, start_timestr=timestr )
      CALL domain_clock_set( head_grid, current_timestr=timestr )

      ! Force to turn off history output in this case
      CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )
      ! Force to turn off boundary input as it was perturbated in add_grid.
      CALL nl_set_io_form_boundary( head_grid%id, 0 )

      head_grid%dtbc = 0
      head_grid%itimestep = 0
      CALL wrf_run
      CALL       wrf_message ( "wrf: back from nonlinear integrate" )

      val_n = 0.0
      !$OMP PARALLEL DO    &
      !$OMP PRIVATE ( ij ) &
      !$OMP REDUCTION (+:val_n)
      DO ij = 1 , head_grid%num_tiles
        CALL inner_dot ( head_grid, val_n,                                          &
                         head_grid%i_start(ij), head_grid%i_end(ij),                &
                         head_grid%j_start(ij), head_grid%j_end(ij)                 )
      END DO
      !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
      val_n = wrf_dm_sum_real ( val_n )
#endif
      val_l=save_l*alpha_m**2
      coef=val_n/val_l
      WRITE(wrf_err_message, FMT='(A,E9.4,A,E23.14,A,E14.7,A,E14.7)') &
           'tl_check: alpha_m=',alpha_m,'  coef=',coef, &
           '  val_n=',val_n,'  val_l=',val_l
      CALL wrf_message(TRIM(wrf_err_message))

   ENDDO  tangentLinearCheck 

   END IF !end of Tangent linear check

   IF ( config_flags%check_AD ) THEN

   wrf_err_message = '==================== Adjoint check ====================='
   CALL wrf_message(TRIM(wrf_err_message))

   wrf_err_message = 'check==== U === V === W == PH === T == MU == MOIST == TRACER ==='
   CALL wrf_message(TRIM(wrf_err_message))

   WRITE(wrf_err_message,'(A,8(1x,L5))') 'check',head_grid%check_u, head_grid%check_v,  &
                               head_grid%check_w, head_grid%check_ph, &
                               head_grid%check_t, head_grid%check_mu, &
                               head_grid%check_moist, head_grid%check_tracer
#if (WRF_CHEM == 1)
   wrf_err_message = 'check==== CHEM == SCALEANT == SCALEBB  ==='
   CALL wrf_message(TRIM(wrf_err_message))

   WRITE(wrf_err_message,'(A,3(1x,L5))') 'check', head_grid%check_chem, &
                               head_grid%check_scaleant, head_grid%check_scalebb
#endif

   CALL wrf_message(TRIM(wrf_err_message))

   CALL restore_grid ( head_grid )
   CALL copy_s_to_g_adjtest ( head_grid, 0.01 )

   CALL copy_g_to_b_adjtest ( head_grid )

   ! Set up basic states reading
   model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxinput16_interval_s (head_grid%id, time_step )
   CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

   IF ( config_flags%mp_physics .NE. 0 .and. config_flags%mp_physics .NE. 99 .and. &
        config_flags%mp_physics .NE. 98 ) CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
   IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. config_flags%bl_pbl_physics .GT. 0 ) CALL nl_set_bl_pbl_physics (head_grid%id, 98)
   IF ( config_flags%cu_physics .GT. 0 ) THEN
      CALL nl_set_cu_physics (head_grid%id, 98)
      head_grid%cudt = 0
   ENDIF
   CALL nl_set_ra_lw_physics (head_grid%id, 0)
   CALL nl_set_ra_sw_physics (head_grid%id, 0)
   IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
   IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. config_flags%sf_surface_physics.NE.PXLSMSCHEME) CALL nl_set_sf_surface_physics (head_grid%id, 0)

   CALL wrf_run_tl

   ! Turn off auxinput16 reading
   CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   val_l = 0.0
   !$OMP PARALLEL DO    &
   !$OMP PRIVATE ( ij ) &
   !$OMP REDUCTION (+:val_l)
   DO ij = 1 , head_grid%num_tiles
     CALL inner_dot_g_adjtest ( head_grid , val_l,                                         &
                        head_grid%i_start(ij), head_grid%i_end(ij),                &
                        head_grid%j_start(ij), head_grid%j_end(ij)                 )
   END DO
   !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
   val_l =  wrf_dm_sum_real ( val_l )
#endif

!  CALL restore_grid ( head_grid )
   !$OMP PARALLEL DO    &
   !$OMP DEFAULT (SHARED) PRIVATE ( ij )
   DO ij = 1 , head_grid%num_tiles
      CALL copy_g_to_a_adjtest ( head_grid,                                           &
                        head_grid%i_start(ij), head_grid%i_end(ij),                &
                        head_grid%j_start(ij), head_grid%j_end(ij)                 )
   END DO
   !$OMP END PARALLEL DO

   ! Set the file names and interval for reading basic states.
   model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxinput16_interval_s (head_grid%id, time_step )
   CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

   CALL wrf_run_ad

   ! Turn off auxinput16 reading
   CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   val_a = 0.0
   !$OMP PARALLEL DO    &
   !$OMP PRIVATE ( ij ) &
   !$OMP REDUCTION (+:val_a)
   DO ij = 1 , head_grid%num_tiles
     CALL inner_dot_a_b_adjtest ( head_grid, val_a,                                  &
                          head_grid%i_start(ij), head_grid%i_end(ij),                &
                          head_grid%j_start(ij), head_grid%j_end(ij)                 )
   END DO
   !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
   val_a =  wrf_dm_sum_real ( val_a )
#endif

   WRITE(wrf_err_message, FMT='(A,E23.14)') 'ad_check: VAL_TL: ', val_l
   CALL wrf_message(TRIM(wrf_err_message))
   WRITE(wrf_err_message, FMT='(A,E23.14)') 'ad_check: VAL_AD: ', val_a
   CALL wrf_message(TRIM(wrf_err_message))

   END IF !ADJ test

   ! Deallocate all working space
   CALL deallocate_grid()

ENDIF

   ! Release linked list for trajectory
   call xtraj_io_initialize

   ! WRF model clean-up.  This calls MPI_FINALIZE() for DM parallel runs.
   CALL wrf_finalize
   STOP
   END SUBROUTINE wrf_adtl_check_sum

   SUBROUTINE wrf_adtl_check_spot( )
!<DESCRIPTION>
!     WRF adjoint and tangent linear code check routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  
! 
!</DESCRIPTION>

#ifdef DM_PARALLEL
#ifndef STUBMPI
      INCLUDE 'mpif.h'
#endif
#endif

   REAL :: alpha_m, val_l, val_a, save_l, val_n, coef, nl_num, nl_den, dummy_der
   REAL, DIMENSION(:,:), ALLOCATABLE :: factors
   REAL, DIMENSION(:,:,:), ALLOCATABLE :: ad_derivative, tl_derivative
   REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: nl_derivative
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: locations_f, locations_i, scenario_matrix
   INTEGER :: ni, nf, nsc, nvar, ij, time_step, rc, &
              nfinal, ninitial, firatio, iter, check_type, psign, &
              iloc_f, jloc_f, kloc_f, iloc_i, jloc_i, kloc_i, &
              vars_count, nnumer, ndenom
   INTEGER :: original_restart_interval, original_checkpt_int, ialarm, interval
   CHARACTER*10, DIMENSION(:), ALLOCATABLE :: check_name !large for expectation of many future chem variables
   CHARACTER*256 :: timestr, rst_outname
   CHARACTER*5000 :: wrf_msg

   ! Return immediately if not dyn_em_check
   IF ( config_flags%dyn_opt .NE. dyn_em_check ) RETURN

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   ! Force to read the lateral boundary condition here.
   CALL med_before_solve_io ( head_grid, config_flags )

   ! Close boundary file, as it will be read again from start later
   CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                                            "DATASET=BOUNDARY" )

   CALL init_domain_size ( head_grid, config_flags )

   ! Release linked list for trajectory
   call xtraj_io_initialize

IF ( config_flags%check_TL .or. config_flags%check_AD .or. config_flags%check_NL ) THEN

   ! Save the initial condition and boundary condition, x
   CALL allocate_grid ( )

#if (WRF_CHEM == 1)
   !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
   ! which is necessary to do 4DVAR -JJG 9-23-2014
   head_grid%scaleant = 1.0D0
   head_grid%scalebb = 1.0D0
#endif

   !$OMP PARALLEL DO    &
   !$OMP DEFAULT (SHARED) PRIVATE ( ij ) &           
   DO ij = 1 , head_grid%num_tiles
      CALL copy_grid_to_s ( head_grid , &
                            head_grid%i_start(ij), head_grid%i_end(ij),                &
                            head_grid%j_start(ij), head_grid%j_end(ij)                 )
   ENDDO

   CALL       wrf_message ( "wrf: calling nonlinear integrate" )
   !model_config_rec%dyn_opt = dyn_em

   IF (head_grid%checkpoint_interval .LE. 0) THEN
      IF ( head_grid%check_AD .OR. head_grid%check_TL) THEN
         ! Set up basic states output
         CALL nl_get_time_step ( head_grid%id, time_step )
         CALL nl_set_auxhist6_interval_s ( head_grid%id, time_step )
         CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
         CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )
      ENDIF
      ! Force to turn off trajectory boundary update
      xtraj_boundary_form=0
   ELSE
      call nl_get_restart_interval_m ( head_grid%id, original_restart_interval)
      call nl_set_restart_interval_m ( head_grid%id, head_grid%checkpoint_interval)
      call nl_set_io_form_restart ( head_grid%id, 2 )
      rst_outname = "wrf_checkpoint_d<domain>_<date>"
      call nl_set_rst_outname ( head_grid%id, rst_outname )
      call nl_set_restart_begin_h( head_grid%id, 0)
      IF (.not. head_grid%boundary_io ) xtraj_boundary_form=0
      IF ( model_config_rec%interval_seconds .lt. 60*head_grid%checkpoint_interval .or. &
           mod( model_config_rec%interval_seconds,60*head_grid%checkpoint_interval ) .ne. 0) THEN
         ! Force to turn off boundary input for trajectory calculation
         xtraj_boundary_form = 0
         CALL wrf_message ("WARNING: interval_seconds must be divisible by &
              &60*checkpoint_interval to use boundary updating in a TL or AD &
              &simulation. Turning off boundary updating.")
      ENDIF
   ENDIF
   ! Set trajectory boundary form
   CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )


   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   IF (head_grid%checkpoint_interval .GT. 0) THEN
      DO ialarm = first_auxinput, last_auxinput
         CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
         IF ( interval .gt. 60*head_grid%checkpoint_interval ) THEN
            IF ( mod(interval, 60*head_grid%checkpoint_interval) .ne. 0 ) THEN
               write(wrf_msg,FMT='(A,I2)')&
                 'WARNING: when checkpoint_interval < auxinputXX_interval,&
                 &auxinput must be divisible by checkpoint.  &
                 &Disabling auxinput',&
                 ialarm-first_input
               CALL wrf_message (TRIM(wrf_msg))
               CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
            ENDIF
         ELSEIF ( mod(60*head_grid%checkpoint_interval, interval) .ne. 0 ) THEN
               write(wrf_msg,FMT='(A,I2)')&
                 'WARNING: when checkpoint_interval > auxinputXX_interval,&
                 &checkpoint must be divisible by auxinput.  &
                 &Disabling auxinput',&
                 ialarm-first_input
               CALL wrf_message (TRIM(wrf_msg))
               CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
         ENDIF
      ENDDO
   ENDIF

   IF ( config_flags%mp_physics .NE. 0 .AND. &
        config_flags%mp_physics .NE. 99 .AND. &
        config_flags%mp_physics .NE. 98 ) &
      CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
   IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. &
        config_flags%bl_pbl_physics .GT. 0 ) &
      CALL nl_set_bl_pbl_physics (head_grid%id, 98)
   IF ( config_flags%cu_physics .GT. 0 ) THEN
      CALL nl_set_cu_physics (head_grid%id, 98)
      head_grid%cudt = 0
   ENDIF
   CALL nl_set_ra_lw_physics (head_grid%id, 0)
   CALL nl_set_ra_sw_physics (head_grid%id, 0)
   IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. &
        config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) &
      CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
   IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. &
        config_flags%sf_surface_physics.NE.PXLSMSCHEME) &
      CALL nl_set_sf_surface_physics (head_grid%id, 0)

   !CALL wrf_run
   CALL wrf_run_cycle


   IF (head_grid%checkpoint_interval .EQ. 0) THEN
      ! Turn off basic states output
      CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
      CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )
   ELSE
     call nl_set_restart_interval_m ( head_grid%id, original_restart_interval)
     call nl_set_io_form_restart ( head_grid%id, 0 )
     rst_outname = "wrfrst_d<domain>_<date>"
     call nl_set_rst_outname ( head_grid%id, rst_outname )
     call nl_set_restart_begin_m ( head_grid%id, original_restart_interval)
   ENDIF
   IF ( xtraj_boundary_form .gt. 0 ) THEN
      ! Close boundary file, as it will be read again from start later
      CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                           "DATASET=BOUNDARY" )
   ENDIF

   if ( .not. head_grid%trajectory_io ) CALL nl_set_inputout_interval_s ( head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL       wrf_message ( "wrf: back from nonlinear integrate" )

   CALL copy_grid_to_b ( head_grid )

   wrf_msg = '============== Adjoint - Tangent Linear Check ==============='
   CALL wrf_message(TRIM(wrf_msg))

   ! Get the locations of adjoint and tangent linear forcings
   CALL allocate_locations( 'locations_f', locations_f, nfinal, ierr)
   IF( ierr .GT. 0) THEN
      CALL wrf_error_fatal( 'adtl_check: error opening locations_f for reading' )
   ENDIF
   CALL allocate_locations( 'locations_i', locations_i, ninitial, ierr)
   IF( ierr .GT. 0) THEN
      CALL wrf_error_fatal(&
                   'adtl_check: error opening locations_i for reading' )
   ENDIF

   ierr = 0
   CALL get_forc_locations( 'locations_f', locations_f, nfinal, ierr)
   CALL get_forc_locations( 'locations_i', locations_i, ninitial, ierr)



   vars_count = 6 
   DO nsc = 1, num_moist
      vars_count = vars_count + 1
   ENDDO
#if (WRF_CHEM == 1)
   DO nsc = 1, num_chem
      vars_count = vars_count + 1
   ENDDO
   DO nsc = 1, num_scaleant
      vars_count = vars_count + 1
   ENDDO
   DO nsc = 1, num_scalebb
      vars_count = vars_count + 1
   ENDDO
#endif
   DO nsc = 1, num_tracer
      vars_count = vars_count + 1
   ENDDO

   ALLOCATE(factors(vars_count,2))
   ALLOCATE(scenario_matrix(vars_count,vars_count))
   scenario_matrix = 0
   CALL gen_scenario_matrix(config_flags, scenario_matrix, vars_count, &
                            model_config_rec%numer_vars, model_config_rec%denom_vars)

   ALLOCATE(check_name(vars_count))
   WRITE(check_name(1),FMT='(A)') 'U'
   WRITE(check_name(2),FMT='(A)') 'V'
   WRITE(check_name(3),FMT='(A)') 'W'
   WRITE(check_name(4),FMT='(A)') 'T'
   WRITE(check_name(5),FMT='(A)') 'PH'
   WRITE(check_name(6),FMT='(A)') 'MU'
   DO nsc = 1, num_moist
      IF( nsc .LT. PARAM_FIRST_SCALAR ) THEN
         WRITE(check_name(6+nsc),FMT='(A)')'DUMMY'
      ELSE
         WRITE(check_name(6+nsc),FMT='(A,I2.1)')'MOIST',nsc-PARAM_FIRST_SCALAR+1
      ENDIF
   ENDDO
   DO nsc = 1, num_tracer
      IF( nsc .LT. PARAM_FIRST_SCALAR ) THEN
         WRITE(check_name(6+nsc+num_moist),FMT='(A)')'DUMMY' 
      ELSE
         WRITE(check_name(6+nsc+num_moist),FMT='(A,I2.1)')'TRACER',nsc-PARAM_FIRST_SCALAR+1
      ENDIF
   ENDDO 
#if (WRF_CHEM == 1)
   DO nsc = 1, num_chem
      IF( nsc .LT. PARAM_FIRST_SCALAR ) THEN
         WRITE(check_name(6+nsc+num_moist+num_tracer),FMT='(A)')'DUMMY' 
      ELSE
         WRITE(check_name(6+nsc+num_moist+num_tracer),FMT='(A,I2.1)')'CHEM',nsc-PARAM_FIRST_SCALAR+1
      ENDIF
   ENDDO  
   DO nsc = 1, num_scaleant
      IF( nsc .LT. PARAM_FIRST_SCALAR ) THEN
         WRITE(check_name(6+nsc+num_moist+num_tracer+num_chem),FMT='(A)')'DUMMY' 
      ELSE
         WRITE(check_name(6+nsc+num_moist+num_tracer+num_chem),FMT='(A,I2.1)')'SCALEANT',nsc-PARAM_FIRST_SCALAR+1
      ENDIF
   ENDDO  
   DO nsc = 1, num_scalebb
      IF( nsc .LT. PARAM_FIRST_SCALAR ) THEN
         WRITE(check_name(6+nsc+num_moist+num_tracer+num_chem+num_scaleant),FMT='(A)')'DUMMY' 
      ELSE
         WRITE(check_name(6+nsc+num_moist+num_tracer+num_chem+num_scaleant),FMT='(A,I2.1)')'SCALEBB',nsc-PARAM_FIRST_SCALAR+1
      ENDIF
   ENDDO 
#endif

   CALL wrf_message("AVAILABLE CONTROL VARIABLES FOR VALIDATION")
   DO nvar = 1, vars_count
      WRITE(wrf_msg, FMT='(I2,A,A)') nvar,'  ',check_name(nvar)
      CALL wrf_message(TRIM(wrf_msg)) 
   ENDDO

   CALL wrf_message("")
   CALL wrf_message("SELECTED CONTROL VARIABLES")
   DO nnumer = 0, vars_count      
      WRITE(wrf_msg, FMT='(I2,A)') nnumer, ' ,  '
      DO ndenom = 1, vars_count
         IF ( nnumer .EQ. 0 ) THEN
            WRITE(wrf_msg, FMT='(A,I2,A)') TRIM(wrf_msg),&
                                ndenom, ' ,  '
         ELSE
            WRITE(wrf_msg, FMT='(A,I2,A)') TRIM(wrf_msg),&
                                scenario_matrix(nnumer,ndenom), ' ,  '
         ENDIF
      ENDDO
      CALL wrf_message(TRIM(wrf_msg))
   ENDDO


   IF ( .NOT. head_grid%check_TLNLonly ) THEN
      firatio = ninitial/nfinal

      ALLOCATE(ad_derivative(ninitial,vars_count,vars_count))   
      ad_derivative = 0.0D0

      ALLOCATE(tl_derivative(ninitial,vars_count,vars_count))
      ALLOCATE(nl_derivative(ninitial,3,vars_count,vars_count))
      tl_derivative = 0.0D0
      nl_derivative = 0.0D0

      !Find all the adjoint sensitivities
      IF ( head_grid%check_AD ) THEN
      numer_vars_Loop1: DO nnumer = 1, vars_count
   
         IF( ALL(scenario_matrix(nnumer,:) .EQ. 0) ) CYCLE numer_vars_Loop1
         factors(:,1) = 0.0D0
         factors(nnumer,1) = 1.0D0
      
         iter = 0
      
         DO nf = 1, nfinal
            iloc_f = locations_f(nf,1)
            jloc_f = locations_f(nf,2)
            kloc_f = locations_f(nf,3)
   
            CALL zero_out_tl( head_grid )
            CALL zero_out_ad( head_grid )
            !$OMP PARALLEL DO    &
            !$OMP PRIVATE ( ij )
            DO ij = 1 , head_grid%num_tiles
               IF( head_grid%i_start(ij) .le. iloc_f .AND. &
                   head_grid%i_end(ij)  .ge. iloc_f .AND. &
                   head_grid%j_start(ij) .le. jloc_f .AND. &
                   head_grid%j_end(ij)  .ge. jloc_f ) THEN
                   CALL spot_force_ad( head_grid, iloc_f, jloc_f , kloc_f, factors(:,1), vars_count ) 
               ENDIF
            ENDDO
            !$OMP END PARALLEL DO
   
            !gradient_out=.TRUE.
            model_config_rec%dyn_opt = dyn_em_check
   
            IF ( wrf_dm_on_monitor() ) THEN
               CALL start_timing
            END IF
   
            CALL wrf_run_ad_checkpt
   
            IF ( wrf_dm_on_monitor() ) THEN
               CALL end_timing('Checkpointed AD')
            END IF
            !gradient_out=.FALSE.
   
            DO ni = 1,firatio
               iter = iter + 1
               iloc_i = locations_i(iter,1)
               jloc_i = locations_i(iter,2)
               kloc_i = locations_i(iter,3)
      
               denom_vars_Loop1: DO ndenom = 1,vars_count
                  IF ( scenario_matrix(nnumer,ndenom) .EQ. 0) CYCLE denom_vars_Loop1
                  factors(:,2)=0.0D0
                  factors(ndenom,2) = 1.0D0
                  WRITE(wrf_msg, FMT='(2(A,I3),5(A))') 'nf = ',nf, &
                               ',  ni = ',ni, ',  check_type = d[', &
                               TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']'
                  CALL wrf_message(TRIM(wrf_msg))
     
                  CALL wrf_message("Extracting the AD sensitivity")
                  dummy_der = 0.D0
                  !Store the AD derivative
                  !$OMP PARALLEL DO    &
                  !$OMP PRIVATE ( ij ) &
                  !$OMP REDUCTION (+: dummy_der) 
                  DO ij = 1 , head_grid%num_tiles
                     IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                        head_grid%i_end(ij)  .ge. iloc_i .AND. &
                        head_grid%j_start(ij) .le. jloc_i .AND. &
                        head_grid%j_end(ij)  .ge. jloc_i ) THEN
                        CALL extract_ad_der( head_grid, dummy_der, iloc_i, jloc_i, kloc_i, head_grid%scale_hr, factors (:,2), vars_count)
                     ENDIF
                  END DO
                  !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                  CALL mpi_barrier(local_communicator, ierr)
                  ad_derivative(iter,nnumer,ndenom) = wrf_dm_sum_real ( dummy_der )
#endif
   
               ENDDO denom_vars_Loop1
            ENDDO      
         ENDDO
      ENDDO numer_vars_Loop1
      ENDIF
   
   
      !Find all the nonlinear and TL sensitivities
      iter = 0
      DO nf = 1, nfinal
         iloc_f = locations_f(nf,1)
         jloc_f = locations_f(nf,2)
         kloc_f = locations_f(nf,3)  
         DO ni = 1,firatio
            iter = iter + 1
            iloc_i = locations_i(iter,1)
            jloc_i = locations_i(iter,2)
            kloc_i = locations_i(iter,3)
   
            denom_vars_Loop2: DO ndenom = 1,vars_count
               IF ( ALL(scenario_matrix(:,ndenom) .EQ. 0) ) CYCLE denom_vars_Loop2
               factors(:,2)=0.0D0
               factors(ndenom,2) = 1.0D0
     
               ! Do Finite Difference test of nonlinear model
               IF ( head_grid%check_NL ) THEN
                  nonlinearLoop1: DO psign=-1,1,2
                     CALL wrf_message ( "wrf: calling nonlinear integrate" ) 
   
                     original_checkpt_int = model_config_rec%checkpoint_interval
                     model_config_rec%checkpoint_interval = 0
                     head_grid%dtbc = 0
                     head_grid%itimestep = 0
   
                     CALL domain_clock_get( head_grid, start_timestr=timestr )
                     CALL domain_clock_set( head_grid, current_timestr=timestr )
   
                     IF ( head_grid%checkpoint_interval .eq. 0 ) THEN !This line used for debugging only
                        CALL restore_grid ( head_grid )
                     ELSE
                        !Read initial conditions
                        model_config_rec%restart = .false.
                        head_grid%restart = .false.
                        CALL nl_set_restart ( head_grid%id, .false. )
                        CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
                        CALL med_initialdata_input( head_grid , config_flags )
#if (WRF_CHEM == 1)
                        !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
                        ! which is necessary to do 4DVAR -JJG 9-23-2014
                        head_grid%scaleant = 1.0D0
                        head_grid%scalebb = 1.0D0
#endif
                     ENDIF
   
!!         IF ( head_grid%check_scaleant) THEN !This line used for debugging only
                     !$OMP PARALLEL DO    &
                     !$OMP PRIVATE ( ij )
                     DO ij = 1 , head_grid%num_tiles
                        IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                            head_grid%i_end(ij)  .ge. iloc_i .AND. &
                            head_grid%j_start(ij) .le. jloc_i .AND. &
                            head_grid%j_end(ij)  .ge. jloc_i ) THEN 
                            CALL spot_force_nl( head_grid, iloc_i, jloc_i, kloc_i, head_grid%scale_hr, &
                                                factors(:,2), vars_count, &
                                                REAL(psign,8) * REAL(config_flags%nl_pert, 8))
                        ENDIF
                     ENDDO
                     !$OMP END PARALLEL DO
!!         ENDIF
   
                     ! Set trajectory boundary form
                     CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )
   
                     model_config_rec%io_form_restart = 0
                     model_config_rec%dyn_opt = dyn_em_check
   
                     !Setup clock and I/O alarms
                     IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
                       IF ( head_grid%domain_clock_created ) THEN
                         CALL WRFU_ClockDestroy( head_grid%domain_clock )
                         head_grid%domain_clock_created = .FALSE.
                       ENDIF
                     ENDIF
                     IF ( ASSOCIATED( head_grid%alarms ) .AND. &
                          ASSOCIATED( head_grid%alarms_created ) ) THEN
                       DO alarmid = 1, MAX_WRF_ALARMS
                         IF ( head_grid%alarms_created( alarmid ) ) THEN
                           CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
                           head_grid%alarms_created( alarmid ) = .FALSE.
                         ENDIF
                       ENDDO
                     ENDIF
                     CALL Setup_Timekeeping ( head_grid )
                     head_grid%start_subtime = domain_get_start_time ( head_grid )
                     head_grid%stop_subtime = domain_get_stop_time ( head_grid )
   
                     ! Force to turn off history output in this case
                     CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )
   
                     CALL wrf_run
   
                     model_config_rec%checkpoint_interval = original_checkpt_int
                     IF ( xtraj_boundary_form .gt. 0 ) THEN
                        ! Close boundary file, as it will be read again from start later
                        CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                                             "DATASET=BOUNDARY" )
                     ENDIF
   
                     CALL wrf_message( "wrf: back from nonlinear integrate" ) 

                     CALL wrf_message("Extracting the FD sensitivities")
                     nl_den = 0.0D0
                     !$OMP PARALLEL DO    &
                     !$OMP PRIVATE ( ij ) &
                     !$OMP REDUCTION (+: nl_den) 
                     DO ij = 1 , head_grid%num_tiles
                        IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                           head_grid%i_end(ij)  .ge. iloc_i .AND. &
                           head_grid%j_start(ij) .le. jloc_i .AND. &
                           head_grid%j_end(ij)  .ge. jloc_i ) THEN
                           CALL extract_nl_den( nl_den, &
                                                iloc_i, jloc_i, kloc_i, 1, &
                                                factors(:,2), vars_count, &
                                                REAL(psign,8)*REAL(config_flags%nl_pert, 8))
                        ENDIF
                     END DO
                     !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                     CALL mpi_barrier(local_communicator, ierr)
                     nl_den = wrf_dm_sum_real ( nl_den )
#endif
                     WRITE(iloc_i,FMT='(3(I6),D27.18)') iloc_i,jloc_i,kloc_i,nl_den

                     numer_vars_Loop2: DO nnumer = 1, vars_count
                        nl_num = 0.0D0
 
                        IF( scenario_matrix(nnumer,ndenom) .EQ. 0 ) CYCLE numer_vars_Loop2
                        factors(:,1) = 0.0D0
                        factors(nnumer,1) = 1.0D0
                        WRITE(wrf_msg, FMT='(2(A,I3),5(A))')'nf = ',nf, &
                                     ',  ni = ',ni, ',  check_type = d[', &
                                     TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']'
                        CALL wrf_message(TRIM(wrf_msg))
   
                        !Store the NL derivative
                        !$OMP PARALLEL DO    &
                        !$OMP PRIVATE ( ij ) &
                        !$OMP REDUCTION (+: nl_num) 
                        DO ij = 1 , head_grid%num_tiles
                           IF( head_grid%i_start(ij) .le. iloc_f .AND. &
                              head_grid%i_end(ij)  .ge. iloc_f .AND. &
                              head_grid%j_start(ij) .le. jloc_f .AND. &
                              head_grid%j_end(ij)  .ge. jloc_f ) THEN
                              CALL extract_nl_num( head_grid, nl_num, &
                                                   iloc_f, jloc_f, kloc_f, &
                                                   factors(:,1), vars_count)
                           ENDIF
                        END DO
                        !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                        CALL mpi_barrier(local_communicator, ierr)
                        nl_num = wrf_dm_sum_real ( nl_num )
#endif
                        nl_derivative(iter,2+psign,nnumer,ndenom) = nl_num / nl_den
                     ENDDO numer_vars_Loop2
                  ENDDO nonlinearLoop1
   !!!               DO nnumer = 1, vars_count
   !!!                  nl_derivative(iter,2,nnumer,ndenom) = (nl_derivative(iter,1,nnumer,ndenom) + nl_derivative(iter,3,nnumer,ndenom)) * 0.5
   !!!               ENDDO
               ENDIF
   
               IF ( head_grid%check_TL ) THEN
                  CALL restore_grid ( head_grid )
   
                  !Set the single forcing value for each tangent linear run
                  CALL zero_out_tl( head_grid )
                  CALL zero_out_ad( head_grid )
      
                  !$OMP PARALLEL DO    &
                  !$OMP PRIVATE ( ij )
                  DO ij = 1 , head_grid%num_tiles
                     IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                        head_grid%i_end(ij)  .ge. iloc_i .AND. &
                        head_grid%j_start(ij) .le. jloc_i .AND. &
                        head_grid%j_end(ij)  .ge. jloc_i ) THEN
                        CALL spot_force_tl( head_grid, iloc_i, jloc_i, kloc_i, head_grid%scale_hr, factors(:,2), vars_count) 
                     ENDIF
                  ENDDO
                  !$OMP END PARALLEL DO
   
                  !gradient_out=.TRUE.
                  model_config_rec%dyn_opt = dyn_em_check
                  IF ( wrf_dm_on_monitor() ) THEN
                     CALL start_timing
                  END IF
   
                  CALL wrf_run_tl_checkpt
   
                  IF ( wrf_dm_on_monitor() ) THEN
                     CALL end_timing('Checkpointed TL')
                  END IF
                  !gradient_out=.FALSE.
     
                  numer_vars_Loop3: DO nnumer = 1, vars_count
                     IF( scenario_matrix(nnumer,ndenom) .EQ. 0 ) CYCLE numer_vars_Loop3
                     factors(:,1) = 0.0D0
                     factors(nnumer,1) = 1.0D0
                     WRITE(wrf_msg, FMT='(2(A,I3),5(A))') 'nf = ',nf, &
                                  ',  ni = ',ni, ',  check_type = d[', &
                                  TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']'
                     CALL wrf_message(TRIM(wrf_msg))
       
                     CALL wrf_message("Extracting the TL sensitivity")
                     dummy_der = 0.D0
                     !Store the TL derivative
                     !$OMP PARALLEL DO    &
                     !$OMP PRIVATE ( ij ) &
                     !$OMP REDUCTION (+: dummy_der)
                     DO ij = 1 , head_grid%num_tiles
                        IF( head_grid%i_start(ij) .le. iloc_f .AND. &
                           head_grid%i_end(ij)  .ge. iloc_f .AND. &
                           head_grid%j_start(ij) .le. jloc_f .AND. &
                           head_grid%j_end(ij)  .ge. jloc_f ) THEN
                           CALL extract_tl_der( head_grid, dummy_der, iloc_f, jloc_f, kloc_f, factors(:,1), vars_count)
                        ENDIF
                     END DO
                     !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                     CALL mpi_barrier(local_communicator, ierr)
                     tl_derivative(iter,nnumer,ndenom) = wrf_dm_sum_real ( dummy_der )
#endif
   
                  ENDDO numer_vars_Loop3
               ENDIF
            ENDDO denom_vars_Loop2
         ENDDO      
      ENDDO
   ELSE
      firatio = nfinal/ninitial

      ALLOCATE(tl_derivative(nfinal,vars_count,vars_count))
      ALLOCATE(nl_derivative(nfinal,3,vars_count,vars_count))
      tl_derivative = 0.0D0
      nl_derivative = 0.0D0

      !Find all the nonlinear and TL sensitivities
      DO ni = 1,ninitial
         iloc_i = locations_i(ni,1)
         jloc_i = locations_i(ni,2)
         kloc_i = locations_i(ni,3)

         head_grid%ni_tlnl = ni

         denom_vars_Loop3: DO ndenom = 1,vars_count
            IF ( ALL(scenario_matrix(:,ndenom) .EQ. 0) ) CYCLE denom_vars_Loop3
            factors(:,2)=0.0D0
            factors(ndenom,2) = 1.0D0
  
            ! Do Finite Difference test of nonlinear model
            IF ( head_grid%check_NL ) THEN
               nonlinearLoop2: DO psign=-1,1,2
                  CALL wrf_message ( "wrf: calling nonlinear integrate" ) 

                  original_checkpt_int = model_config_rec%checkpoint_interval
                  model_config_rec%checkpoint_interval = 0
                  head_grid%dtbc = 0
                  head_grid%itimestep = 0

                  CALL domain_clock_get( head_grid, start_timestr=timestr )
                  CALL domain_clock_set( head_grid, current_timestr=timestr )

                  IF ( head_grid%checkpoint_interval .eq. 0 ) THEN !This line used for debugging only
                     CALL restore_grid ( head_grid )
                  ELSE
                     !Read initial conditions
                     model_config_rec%restart = .false.
                     head_grid%restart = .false.
                     CALL nl_set_restart ( head_grid%id, .false. )
                     CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
                     CALL med_initialdata_input( head_grid , config_flags )
#if (WRF_CHEM == 1)
                     !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
                     ! which is necessary to do 4DVAR -JJG 9-23-2014
                     head_grid%scaleant = 1.0D0
                     head_grid%scalebb = 1.0D0
#endif
                  ENDIF

!!      IF ( head_grid%check_scalebb) THEN !This line used for debugging only
                  !$OMP PARALLEL DO    &
                  !$OMP PRIVATE ( ij )
                  DO ij = 1 , head_grid%num_tiles
                     IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                         head_grid%i_end(ij)  .ge. iloc_i .AND. &
                         head_grid%j_start(ij) .le. jloc_i .AND. &
                         head_grid%j_end(ij)  .ge. jloc_i ) THEN 
                         CALL spot_force_nl( head_grid, iloc_i, jloc_i, kloc_i, head_grid%scale_hr, &
                                             factors(:,2), vars_count, &
                                             REAL(psign,8) * REAL(config_flags%nl_pert, 8))
                     ENDIF
                  ENDDO
                  !$OMP END PARALLEL DO
!!      ENDIF

                  ! Set trajectory boundary form
                  CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )

                  model_config_rec%io_form_restart = 0
                  model_config_rec%dyn_opt = dyn_em_check

                  !Setup clock and I/O alarms
                  IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
                    IF ( head_grid%domain_clock_created ) THEN
                      CALL WRFU_ClockDestroy( head_grid%domain_clock )
                      head_grid%domain_clock_created = .FALSE.
                    ENDIF
                  ENDIF
                  IF ( ASSOCIATED( head_grid%alarms ) .AND. &
                       ASSOCIATED( head_grid%alarms_created ) ) THEN
                    DO alarmid = 1, MAX_WRF_ALARMS
                      IF ( head_grid%alarms_created( alarmid ) ) THEN
                        CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
                        head_grid%alarms_created( alarmid ) = .FALSE.
                      ENDIF
                    ENDDO
                  ENDIF
                  CALL Setup_Timekeeping ( head_grid )
                  head_grid%start_subtime = domain_get_start_time ( head_grid )
                  head_grid%stop_subtime = domain_get_stop_time ( head_grid )

                  ! Force to turn off history output in this case
                  CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

                  head_grid%nl_perturbed = .true.
                  CALL wrf_run
                  head_grid%nl_perturbed = .false.

                  model_config_rec%checkpoint_interval = original_checkpt_int
                  IF ( xtraj_boundary_form .gt. 0 ) THEN
                     ! Close boundary file, as it will be read again from start later
                     CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                                          "DATASET=BOUNDARY" )
                  ENDIF

                  CALL wrf_message( "wrf: back from nonlinear integrate" ) 

                  CALL wrf_message("Extracting the FD sensitivities")  
                  nl_den = 0.0D0
                  !$OMP PARALLEL DO    &
                  !$OMP PRIVATE ( ij ) &
                  !$OMP REDUCTION (+: nl_den) 
                  DO ij = 1 , head_grid%num_tiles
                     IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                        head_grid%i_end(ij)  .ge. iloc_i .AND. &
                        head_grid%j_start(ij) .le. jloc_i .AND. &
                        head_grid%j_end(ij)  .ge. jloc_i ) THEN
                        CALL extract_nl_den( nl_den, &
                                             iloc_i, jloc_i, kloc_i, 1, &
                                             factors(:,2), vars_count, &
                                             REAL(psign,8)*REAL(config_flags%nl_pert, 8))

                        WRITE(iloc_i*100+jloc_i,FMT='(3(I6),D27.18)') iloc_i,jloc_i,kloc_i,nl_den
                     ENDIF
                  END DO
                  !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                  CALL mpi_barrier(local_communicator, ierr)
                  nl_den = wrf_dm_sum_real ( nl_den )
#endif

                  numer_vars_Loop4: DO nnumer = 1, vars_count
                     IF( scenario_matrix(nnumer,ndenom) .EQ. 0 ) CYCLE numer_vars_Loop4
                     factors(:,1) = 0.0D0
                     factors(nnumer,1) = 1.0D0

                     iter=(ni-1)*firatio
                     DO nf = 1, firatio
                        nl_num = 0.0D0
                        iter=iter+1
                        iloc_f = locations_f(iter,1)
                        jloc_f = locations_f(iter,2)
                        kloc_f = locations_f(iter,3)
                        !Store the NL derivative
                        !$OMP PARALLEL DO    &
                        !$OMP PRIVATE ( ij ) &
                        !$OMP REDUCTION (+: nl_num) 
                        DO ij = 1 , head_grid%num_tiles
                           IF( head_grid%i_start(ij) .le. iloc_f .AND. &
                              head_grid%i_end(ij)  .ge. iloc_f .AND. &
                              head_grid%j_start(ij) .le. jloc_f .AND. &
                              head_grid%j_end(ij)  .ge. jloc_f ) THEN
                              CALL extract_nl_num( head_grid, nl_num, &
                                                   iloc_f, jloc_f, kloc_f, &
                                                   factors(:,1), vars_count)
                           ENDIF
                        END DO
                        !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                        CALL mpi_barrier(local_communicator, ierr)
                        nl_num = wrf_dm_sum_real ( nl_num )
#endif  
                        nl_derivative(iter,2+psign,nnumer,ndenom) = nl_num / nl_den

                        WRITE(wrf_msg, FMT='(2(A,I3),5(A))')'nf = ',nf, &
                                     ',  ni = ',ni, ',  check_type = d[', &
                                     TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']'
                        CALL wrf_message(TRIM(wrf_msg))
                     ENDDO
                  ENDDO numer_vars_Loop4
               ENDDO nonlinearLoop2
            ENDIF

            IF ( head_grid%check_TL ) THEN
               CALL restore_grid ( head_grid )

               !Set the single forcing value for each tangent linear run
               CALL zero_out_tl( head_grid )
               CALL zero_out_ad( head_grid )
   
               !$OMP PARALLEL DO    &
               !$OMP PRIVATE ( ij )
               DO ij = 1 , head_grid%num_tiles
                  IF( head_grid%i_start(ij) .le. iloc_i .AND. &
                     head_grid%i_end(ij)  .ge. iloc_i .AND. &
                     head_grid%j_start(ij) .le. jloc_i .AND. &
                     head_grid%j_end(ij)  .ge. jloc_i ) THEN
                     CALL spot_force_tl( head_grid, iloc_i, jloc_i, kloc_i, head_grid%scale_hr, factors(:,2), vars_count) 
                  ENDIF
               ENDDO
               !$OMP END PARALLEL DO

               !gradient_out=.TRUE.
               model_config_rec%dyn_opt = dyn_em_check
               IF ( wrf_dm_on_monitor() ) THEN
                  CALL start_timing
               END IF

               head_grid%print_NL = .false.
               CALL wrf_run_tl_checkpt
               head_grid%print_NL = .true.
  
               IF ( wrf_dm_on_monitor() ) THEN
                  CALL end_timing('Checkpointed TL')
               END IF
               !gradient_out=.FALSE.
  
               numer_vars_Loop5: DO nnumer = 1, vars_count
                  IF( scenario_matrix(nnumer,ndenom) .EQ. 0 ) CYCLE numer_vars_Loop5
                  factors(:,1) = 0.0D0
                  factors(nnumer,1) = 1.0D0
    
                  CALL wrf_message("Extracting the TL sensitivity")
                  dummy_der = 0.D0
                  iter=(ni-1)*firatio
                  DO nf = 1, firatio
                     iter=iter+1
                     iloc_f = locations_f(iter,1)
                     jloc_f = locations_f(iter,2)
                     kloc_f = locations_f(iter,3) 
                     !Store the TL derivative
                     !$OMP PARALLEL DO    &
                     !$OMP PRIVATE ( ij ) &
                     !$OMP REDUCTION (+: dummy_der)
                     DO ij = 1 , head_grid%num_tiles
                        IF( head_grid%i_start(ij) .le. iloc_f .AND. &
                           head_grid%i_end(ij)  .ge. iloc_f .AND. &
                           head_grid%j_start(ij) .le. jloc_f .AND. &
                           head_grid%j_end(ij)  .ge. jloc_f ) THEN
                           CALL extract_tl_der( head_grid, dummy_der, iloc_f, jloc_f, kloc_f, factors(:,1), vars_count)
                        ENDIF
                     END DO
                     !$OMP END PARALLEL DO
#ifdef DM_PARALLEL
                     CALL mpi_barrier(local_communicator, ierr)
                     tl_derivative(iter,nnumer,ndenom) = wrf_dm_sum_real ( dummy_der )
#endif
                     WRITE(wrf_msg, FMT='(2(A,I3),5(A))') 'nf = ',nf, &
                                  ',  ni = ',ni, ',  check_type = d[', &
                                  TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']'
                     CALL wrf_message(TRIM(wrf_msg))

                  ENDDO
               ENDDO numer_vars_Loop5
            ENDIF
         ENDDO denom_vars_Loop3
      ENDDO
   ENDIF

   ! Print out the sensitivities
   CALL wrf_message(&
        "====================== Results =======================")
   CALL wrf_message(&
        "Open fort.90456 for full results")

   IF ( .NOT. head_grid%check_TLNLonly ) THEN
      numer_vars_Loop6: DO nnumer = 1, vars_count
         IF( ALL(scenario_matrix(nnumer,:) .EQ. 0) ) CYCLE numer_vars_Loop6
         iter = 0
         DO nf = 1, nfinal
            DO ni = 1, firatio
               iter = iter+1
               iloc_i = locations_i(iter,1)
               jloc_i = locations_i(iter,2)
               kloc_i = locations_i(iter,3)
               iloc_f = locations_f(nf,1)
               jloc_f = locations_f(nf,2)
               kloc_f = locations_f(nf,3)
               WRITE(90456, FMT='(I6,6(1x,I8))')&
                            iter,&
                            iloc_f,&
                            jloc_f,&
                            kloc_f,&
                            iloc_i,&
                            jloc_i,&
                            kloc_i
            ENDDO
         ENDDO

         DO ndenom = 1, vars_count
            IF( scenario_matrix(nnumer,ndenom) .EQ. 0) CYCLE
            WRITE(90456, FMT='(5(A))') '        ------check_type = d[',&
                              TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']------'

            IF ( head_grid%check_NL ) THEN
               WRITE(90456, FMT='(A,2(A3),4(1x,A17))') '      ','i','f','TL','AD','NL_BD','NL_FD'

               iter = 0
               DO nf = 1, nfinal
                  DO ni = 1, firatio
                     iter = iter+1
                     WRITE(90456, FMT='(I6,2(I5),4(D25.16))')&
                                  iter,ni,nf,&
                                  tl_derivative(iter,nnumer,ndenom), &
                                  ad_derivative(iter,nnumer,ndenom), &
                                  nl_derivative(iter,1,nnumer,ndenom), &
                                  nl_derivative(iter,3,nnumer,ndenom)
                  ENDDO
               ENDDO
            ELSE
               WRITE(90456, FMT='(A,2(A3),2(1x,A17))') '         ','i','f','TL','AD'
               iter = 0
               DO nf = 1, nfinal
                  DO ni = 1, firatio
                     iter = iter+1
                     WRITE(90456, FMT='(I6,2(I5),2(D25.16))')&
                                  iter,ni,nf,&
                                  tl_derivative(iter,nnumer,ndenom), &
                                  ad_derivative(iter,nnumer,ndenom)
                  ENDDO
               ENDDO
            ENDIF        
         ENDDO   
      ENDDO numer_vars_Loop6
   ELSE
      DO ni = 1, ninitial
         iloc_i = locations_i(ni,1)
         jloc_i = locations_i(ni,2)
         kloc_i = locations_i(ni,3)
         iter=(ni-1)*firatio
         DO nf = 1, firatio
            iter = iter+1
            iloc_f = locations_f(iter,1)
            jloc_f = locations_f(iter,2)
            kloc_f = locations_f(iter,3)
            WRITE(90456, FMT='(I6,6(1x,I8))')&
                         iter,&
                         iloc_f,&
                         jloc_f,&
                         kloc_f,&
                         iloc_i,&
                         jloc_i,&
                         kloc_i
         ENDDO
      ENDDO
      numer_vars_Loop7: DO nnumer = 1, vars_count
         IF( ALL(scenario_matrix(nnumer,:) .EQ. 0) ) CYCLE numer_vars_Loop7

         DO ndenom = 1, vars_count
            IF( scenario_matrix(nnumer,ndenom) .EQ. 0) CYCLE
            WRITE(90456, FMT='(5(A))') '        ------check_type = d[',&
                              TRIM(check_name(nnumer)),']/d[',TRIM(check_name(ndenom)),']------'

            IF ( head_grid%check_NL ) THEN
               WRITE(90456, FMT='(A,2(A3),3(1x,A17))') '      ','i','f','TL','NL_BD','NL_FD'
               DO ni = 1, ninitial
                  iter=(ni-1)*firatio
                  DO nf = 1, firatio
                     iter = iter+1
                     WRITE(90456, FMT='(I6,2(I5),3(D25.16))')&
                                  iter,ni,nf,&
                                  tl_derivative(iter,nnumer,ndenom), &
                                  nl_derivative(iter,1,nnumer,ndenom), &
                                  nl_derivative(iter,3,nnumer,ndenom)
                  ENDDO
               ENDDO
            ELSE
               WRITE(90456, FMT='(A,2(A3),2(1x,A17))') '         ','i','f','TL','AD'
               DO ni = 1, ninitial
                  iter=(ni-1)*firatio
                  DO nf = 1, firatio
                     iter = iter+1
                     WRITE(90456, FMT='(I6,2(I5),1(D25.16))')&
                                  iter,ni,nf,&
                                  tl_derivative(iter,nnumer,ndenom)
                  ENDDO
               ENDDO
            ENDIF        
         ENDDO   
      ENDDO numer_vars_Loop7
   ENDIF

   DEALLOCATE(ad_derivative)
   DEALLOCATE(tl_derivative)
   DEALLOCATE(nl_derivative)
  
   DEALLOCATE(check_name)
   DEALLOCATE(scenario_matrix)
   DEALLOCATE(factors)
   DEALLOCATE(locations_f)
   DEALLOCATE(locations_i)

   ! Deallocate all working space
   CALL deallocate_grid()

ENDIF

   ! Release linked list for trajectory
   call xtraj_io_initialize

   ! WRF model clean-up.  This calls MPI_FINALIZE() for DM parallel runs.
   CALL wrf_finalize
   STOP
   END SUBROUTINE wrf_adtl_check_spot

   SUBROUTINE wrf_run_tl( )
!<DESCRIPTION>
!     WRF tangent linear run routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  The start
! and stop times for the domain are set to the start and stop time of the
! model run, and then <a href=integrate.html>integrate</a> is called to
! advance the domain forward through that specified time interval.  On
! return, the simulation is completed.  
! 
!</DESCRIPTION>

   CHARACTER*256 :: timestr
   INTEGER       :: rc
   INTEGER       :: io_auxh17
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_message ( "wrf: calling tangent linear integrate" )
   model_config_rec%dyn_opt = dyn_em_tl
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )
   CALL domain_clock_get( head_grid, start_timestr=timestr )
   CALL domain_clock_set( head_grid, current_timestr=timestr )

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   IF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. SURFDRAGSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
   ELSEIF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. ACMPBLSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
      head_grid%rqcblten = 0.0D0
      head_grid%rqiblten = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. LSCONDSCHEME ) THEN
      head_grid%h_diabatic = 0.0D0
      head_grid%qv_diabatic = 0.0D0
      head_grid%qc_diabatic = 0.0D0
      head_grid%rainnc = 0.0D0
      head_grid%rainncv = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. MKESSLERSCHEME ) THEN
      head_grid%h_diabatic = 0.0
      head_grid%qv_diabatic = 0.0
      head_grid%qc_diabatic = 0.0
      head_grid%rainnc = 0.0
      head_grid%rainncv = 0.0
   ENDIF
   IF ( model_config_rec%cu_physics(head_grid%id) .EQ. DUCUSCHEME ) THEN
      head_grid%rthcuten = 0.0D0
      head_grid%rqvcuten = 0.0D0
      head_grid%rainc = 0.0D0
      head_grid%raincv = 0.0D0
      head_grid%pratec = 0.0
   ENDIF

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_get_io_form_auxhist17( head_grid%id, io_auxh17 )
      CALL nl_set_io_form_auxhist17( head_grid%id, 0 )
   ENDIF

   head_grid%itimestep = 0
   CALL start_domain ( head_grid , .TRUE. )

   IF ( ASSOCIATED(xtraj_tail) ) xtraj_pointer => xtraj_tail
   CALL integrate ( head_grid )

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_set_io_form_auxhist17( head_grid%id, io_auxh17 )
   ENDIF

   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id , 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status , ierr ) 
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags ,&
                           "DATASET=BOUNDARY" )
   ENDIF

   CALL       wrf_message ( "wrf: back from tangent linear integrate" )

   END SUBROUTINE wrf_run_tl


   SUBROUTINE wrf_run_tl_checkpt( )
!<DESCRIPTION>
!     WRF tangent linear run routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  The start
! and stop times for the domain are set to the start and stop time of the
! model run, and then <a href=integrate.html>integrate</a> is called to
! advance the domain forward through that specified time interval.  On
! return, the simulation is completed.  
! 
!</DESCRIPTION>

   CHARACTER*256 :: timestr
   INTEGER       :: rc
   INTEGER       :: io_auxh17
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status

!checkpoint variables
   INTEGER                 :: timestep, ialarm, interval
   integer                           :: icheck
   integer                           :: run_days, run_hours, run_minutes
   character(len=256)                :: STOPtimestr, STARTtimestr, checkstart_timestr
   character(len=256)                :: debugstr
   character(len=1000)               :: msg 
   integer                           :: checkpoint_interval_m, checkpoint_interval_h, checkpoint_interval_d, ncheck
   INTEGER :: start_year,start_month,start_day,start_hour,start_minute
   INTEGER :: end_year,end_month,end_day,end_hour,end_minute
   integer :: start_year0,start_month0,start_day0,start_hour0,start_minute0
   integer :: end_year0,end_month0,end_day0,end_hour0,end_minute0
   integer , dimension(12) :: days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_message ( "wrf: calling tangent linear integrate" )

   IF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. SURFDRAGSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
   ELSEIF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. ACMPBLSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
      head_grid%rqcblten = 0.0D0
      head_grid%rqiblten = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. LSCONDSCHEME ) THEN
      head_grid%h_diabatic = 0.0D0
      head_grid%rainnc = 0.0D0
      head_grid%rainncv = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. MKESSLERSCHEME ) THEN
      head_grid%h_diabatic = 0.0
      head_grid%rainnc = 0.0
      head_grid%rainncv = 0.0
   ENDIF
   IF ( model_config_rec%cu_physics(head_grid%id) .EQ. DUCUSCHEME ) THEN
      head_grid%rthcuten = 0.0D0
      head_grid%rqvcuten = 0.0D0
      head_grid%rainc = 0.0D0
      head_grid%raincv = 0.0D0
      head_grid%pratec = 0.0
   ENDIF

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_get_io_form_auxhist17( head_grid%id, io_auxh17 )
      CALL nl_set_io_form_auxhist17( head_grid%id, 0 )
   ENDIF

   CALL nl_get_time_step ( head_grid%id, timestep)  

   IF(head_grid%checkpoint_interval > 0) then
!BEGIN TLM CHECKPOINTING
      CALL       wrf_message ( "wrf: using multiple checkpoint windows" )
      model_config_rec%rst_inname = "wrf_checkpoint_d<domain>_<date>"
      IF(head_grid%cycle_interval > 0) then
         model_config_rec%input_inname = "wrfinput_d<domain>_<date>"
      ENDIF

      IF (head_grid%checkpoint_interval > 180) then
         CALL wrf_message("WARNING: large checkpoint_interval may &
              &overload memory.  If program crashes, consider reducing &
              &checkpoint_interval and/or number of processors.")
      END IF

      ! Store input boundary form
      CALL nl_get_io_form_boundary( head_grid%id, input_boundary_form )

      IF (.not. head_grid%boundary_io .or. input_boundary_form .eq. 0) xtraj_boundary_form=0
      IF ( model_config_rec%interval_seconds .lt. 60*head_grid%checkpoint_interval .or. &
           mod( model_config_rec%interval_seconds,60*head_grid%checkpoint_interval ) .ne. 0) THEN
         ! Force to turn off boundary input for trajectory calculation
         xtraj_boundary_form = 0
         CALL wrf_message ("WARNING: interval_seconds must be divisible by &
              &60*checkpoint_interval to use boundary updating in a TL or AD &
              &simulation. Turning off boundary updating.")
      ENDIF

      !CALL nl_set_restart_interval_h ( head_grid%id, (ncheck+1)*head_grid%checkpoint_interval/60)
      CALL nl_set_restart_interval_m ( head_grid%id, 0)

      start_year = model_config_rec%start_year(head_grid%id)
      start_month = model_config_rec%start_month(head_grid%id)
      start_day = model_config_rec%start_day(head_grid%id)
      start_hour = model_config_rec%start_hour(head_grid%id)
      start_minute = model_config_rec%start_minute(head_grid%id)

      end_year = model_config_rec%end_year(head_grid%id)
      end_month = model_config_rec%end_month(head_grid%id)
      end_day = model_config_rec%end_day(head_grid%id)
      end_hour = model_config_rec%end_hour(head_grid%id)
      end_minute = model_config_rec%end_minute(head_grid%id)

      ! Store original time keeping variables
      CALL nl_get_run_days  (head_grid%id,run_days)
      CALL nl_get_run_hours (head_grid%id,run_hours)
      CALL nl_get_run_minutes (head_grid%id,run_minutes)
      ncheck = (run_hours + 24 * run_days)

      IF (REAL(ncheck,8) .lt. REAL(head_grid%checkpoint_interval,8)/60.) THEN
         CALL wrf_error_fatal(&
                 'Need to set total run_hours >= checkpoint_interval/60')
      ENDIF
      ncheck = ncheck*60 / head_grid%checkpoint_interval

      end_year0  = end_year
      end_month0 = end_month
      end_day0   = end_day
      end_hour0  = end_hour
      end_minute0  = end_minute

      start_year0 = start_year
      start_month0 = start_month
      start_day0 = start_day
      start_hour0 = start_hour
      start_minute0 = start_minute


      checkpoint_interval_m = mod(head_grid%checkpoint_interval,60)
      checkpoint_interval_h = (head_grid%checkpoint_interval - checkpoint_interval_m) / 60
      IF ( checkpoint_interval_h > 24 ) THEN
         checkpoint_interval_d = (checkpoint_interval_h - mod(checkpoint_interval_h,24)) / 24
         checkpoint_interval_h = mod(checkpoint_interval_h,24)
      ELSE
         checkpoint_interval_d = 0
      ENDIF

      CALL domain_clock_get (head_grid, stop_timestr=STOPtimestr)
      CALL domain_clock_get (head_grid, start_timestr=STARTtimestr)
      checkstart_timestr=STARTtimestr
   
      model_config_rec%run_days = checkpoint_interval_d
      model_config_rec%run_hours = checkpoint_interval_h
      model_config_rec%run_minutes = checkpoint_interval_m

      !!Assuming all forcing has already been calculated read/write settings are in place

! Need to add this subroutine to calculate the forcing for 
! first checkpoint window only.  This will be necessary when the forcing
! isn't calculated outside this subroutine (e.g. in WRF-4DVar or 
! wrf_adtl_check_spot
      !CALL calculate_tl_forcing( head_grid, model_config_rec ) !(Doesn't exist yet.)

      DO icheck = 1,ncheck,1
         head_grid%check_count = icheck
         write(msg,FMT='(A,I5,A,I5)')'WRF TLM: Simulating checkpoint ',head_grid%check_count,' of ',ncheck
         CALL wrf_debug(1,TRIM(msg))

         !Set clock to checkpoint initial condition time stamp for I/O
         CALL domain_clock_set (head_grid, current_timestr=checkstart_timestr )
         CALL domain_clock_set (head_grid, time_step_seconds=60*head_grid%checkpoint_interval)
         CALL domain_clockadvance (head_grid)
         CALL domain_clock_get (head_grid, current_timestr=checkstart_timestr)
         CALL domain_clock_set (head_grid, stop_timestr=checkstart_timestr)
         CALL domain_clock_set (head_grid, time_step_seconds=-60*head_grid%checkpoint_interval)
         CALL domain_clockadvance (head_grid)
         CALL domain_clock_set (head_grid, time_step_seconds=timestep)
         CALL domain_clockprint(150, head_grid, &
                                'get CurrTime from clock,')

         !Read initial conditions
         if ( icheck .eq. 1 .OR. (head_grid%cycle_interval .GT. 0 .AND. &
              mod(head_grid%check_count*head_grid%checkpoint_interval, &
                  head_grid%cycle_interval) .eq. 0 ) ) then
            model_config_rec%restart = .false.
            head_grid%restart = .false.
            call nl_set_restart ( head_grid%id, .false. )
            model_config_rec%io_form_restart = 0
            head_grid%itimestep = 0
         else
            model_config_rec%restart = .true.
            head_grid%restart = .true.
            call nl_set_restart ( head_grid%id, .true. )
            model_config_rec%io_form_restart = 2
            head_grid%itimestep = 1
         endif

         CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )

!#ifdef DM_PARALLEL
!         CALL get_config_as_buffer( configbuf, configbuflen, nbytes )
!         CALL wrf_dm_bcast_bytes( configbuf, nbytes )
!         CALL set_config_as_buffer( configbuf, configbuflen )
!#endif

         CALL med_initialdata_input( head_grid , config_flags )
         model_config_rec%io_form_restart = 0
         head_grid%itimestep = 0
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
         !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
         ! which is necessary to do 4DVAR -JJG 9-23-2014
         IF(.not.model_config_rec%var4d_run) then
            head_grid%scaleant = 1.0D0
            head_grid%scalebb = 1.0D0
         ENDIF
#endif

         ! Set time keeping variables to checkpoint values
         model_config_rec%start_year(head_grid%id) = start_year
         model_config_rec%start_month(head_grid%id) = start_month
         model_config_rec%start_day(head_grid%id) = start_day
         model_config_rec%start_hour(head_grid%id) = start_hour
         model_config_rec%start_minute(head_grid%id) = start_minute


         if (start_minute + checkpoint_interval_m > 59) then
            end_minute = start_minute + checkpoint_interval_m - 60
            end_hour  = start_hour  + 1
         else
            end_minute = start_minute + checkpoint_interval_m
            end_hour  = start_hour
         end if

         if (end_hour + checkpoint_interval_h > 23) then
            end_hour = end_hour + checkpoint_interval_h - 24
            end_day  = start_day  + 1
         else
            end_hour = end_hour + checkpoint_interval_h
            end_day  = start_day
         end if

         end_year  = start_year
         end_month = start_month

         !Account for leap years
         if(start_month .eq. 2 .and. end_day > 24) then  !Assuming the checkpoint interval < 5 days
            if ( mod(end_year,4) .eq. 0 )then
               if ( mod(end_year,100) .eq. 0 ) then
                  if ( mod(end_year,400) .eq. 0 ) then
                     days_in_month(2) = 29
                  else
                     days_in_month(2) = 28
                  end if
               else
                  days_in_month(2) = 29
               end if
            else
               days_in_month(2) = 28
            end if
         end if

         if (end_day + checkpoint_interval_d > days_in_month(start_month)) then
            if(start_month + 1 > 12) then
               end_month = 1
               end_year  = start_year + 1
            else
               end_month = start_month + 1
            end if
            end_day = end_day + checkpoint_interval_d - days_in_month(start_month)
         else
            end_day = end_day + checkpoint_interval_d
         end if

         model_config_rec%end_year(head_grid%id) = end_year
         model_config_rec%end_month(head_grid%id) = end_month
         model_config_rec%end_day(head_grid%id) = end_day
         model_config_rec%end_hour(head_grid%id) = end_hour
         model_config_rec%end_minute(head_grid%id) = end_minute

! Get nonlinear model trajectory
         CALL wrf_message ( "wrf: calculating checkpoint trajectory" )

         ! Set trajectory boundary form
         CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )  

         CALL nl_set_auxhist6_interval_s ( head_grid%id, timestep )
         CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
         CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )

         IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
           IF ( head_grid%domain_clock_created ) THEN
             CALL WRFU_ClockDestroy( head_grid%domain_clock )
             head_grid%domain_clock_created = .FALSE.
           ENDIF
         ENDIF
         IF ( ASSOCIATED( head_grid%alarms ) .AND. &
              ASSOCIATED( head_grid%alarms_created ) ) THEN
           DO alarmid = 1, MAX_WRF_ALARMS
             IF ( head_grid%alarms_created( alarmid ) ) THEN
               CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
               head_grid%alarms_created( alarmid ) = .FALSE.
             ENDIF
           ENDDO
         ENDIF
         CALL Setup_Timekeeping ( head_grid )
         head_grid%start_subtime = domain_get_start_time ( head_grid )
         head_grid%stop_subtime = domain_get_stop_time ( head_grid )

         ! Force to turn off history output in this case
         CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

         model_config_rec%dyn_opt = dyn_em_check

         ! Enable boundary reading for checkpoint intervals shorter than
         ! lbc reading interval
         IF ( xtraj_boundary_form .gt. 0 .and. &
              model_config_rec%interval_seconds .gt. 60*head_grid%checkpoint_interval .and. &
              ( mod( start_hour, model_config_rec%interval_seconds/3600 ) .ne. 0 .or. &
              start_minute .ne. 0 ) ) THEN
            CALL WRFU_AlarmDisable( head_grid%alarms( BOUNDARY_ALARM ), rc=rc )
!            CALL WRFU_AlarmRingerOff( head_grid%alarms( BOUNDARY_ALARM ), rc=rc )
         ENDIF

         !Update time since last boundary read
         head_grid%dtbc = REAL( mod( 60 * (head_grid%check_count-1) * &
                               head_grid%checkpoint_interval, &
                               model_config_rec%interval_seconds ), 8)
!         head_grid%dtbc = REAL( mod( 3600 * start_hour, model_config_rec%interval_seconds ) &
!                          + 60 * start_minute , 8)

         ! Enable auxinput streams for chekpoint_interval < io_intervals
         IF ( icheck .GT. 1 ) THEN
            DO ialarm = first_auxinput, last_auxinput
               CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
               IF ( interval .gt. 60*head_grid%checkpoint_interval ) THEN

                  IF ( mod(start_hour*3600+start_minute*60, interval) .ne. 0 .OR. &
                       mod(interval, 60*head_grid%checkpoint_interval) .ne. 0 ) THEN
                     CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
                  ENDIF

               ELSEIF ( mod(60*head_grid%checkpoint_interval, interval) .ne. 0 ) THEN
                     write(msg,FMT='(A,I2)')&
                       'WARNING: when checkpoint_interval > auxinputXX_interval,&
                       &checkpoint must be divisible by auxinput.  &
                       &Disabling auxinput',&
                       ialarm-first_input
                     CALL wrf_message (TRIM(msg))
                     CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
               ENDIF
            ENDDO
         ENDIF


         CALL domain_clock_get( head_grid, start_timestr=timestr )
         CALL domain_clock_set( head_grid, current_timestr=timestr )

         ! Release linked list for trajectory
         CALL xtraj_io_initialize

         CALL wrf_run

         ! Turn off basic states output
         CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
         CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )

!NOW THE TLM
         CALL wrf_message( "wrf: running TLM" )

         head_grid%itimestep = 0
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
#endif

         !Read initial conditions
!         if ( icheck .eq. 1 ) then
         if ( icheck .eq. 1 .OR. (head_grid%cycle_interval .GT. 0 .AND. &
              mod(head_grid%check_count*head_grid%checkpoint_interval, &
                  head_grid%cycle_interval) .eq. 0 ) ) then
            !Set clock to simulation initial condition time stamp for I/O
            CALL domain_clock_set (head_grid, current_timestr=STARTtimestr )
            model_config_rec%restart = .false.
            head_grid%restart = .false.
            CALL nl_set_restart ( head_grid%id, .false. )
            CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )

!#ifdef DM_PARALLEL
!            CALL get_config_as_buffer( configbuf, configbuflen, nbytes )
!            CALL wrf_dm_bcast_bytes( configbuf, nbytes )
!            CALL set_config_as_buffer( configbuf, configbuflen )
!#endif

            CALL med_initialdata_input( head_grid , config_flags )
#if (WRF_CHEM == 1)
            !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
            ! which is necessary to do 4DVAR -JJG 9-23-2014
            IF(.not.model_config_rec%var4d_run) then
               head_grid%scaleant = 1.0D0
               head_grid%scalebb = 1.0D0
            ENDIF
#endif
!!!         else
!!!            model_config_rec%restart = .true.
         endif

         model_config_rec%dyn_opt = dyn_em_tl

         ! Force to turn off boundary input during derivative model
         CALL nl_set_io_form_boundary( head_grid%id, 0 ) 

         CALL nl_set_auxinput16_interval_s (head_grid%id, timestep )
         CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
         CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

         IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
           IF ( head_grid%domain_clock_created ) THEN
             CALL WRFU_ClockDestroy( head_grid%domain_clock )
             head_grid%domain_clock_created = .FALSE.
           ENDIF
         ENDIF
         IF ( ASSOCIATED( head_grid%alarms ) .AND. &
              ASSOCIATED( head_grid%alarms_created ) ) THEN
           DO alarmid = 1, MAX_WRF_ALARMS
             IF ( head_grid%alarms_created( alarmid ) ) THEN
               CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
               head_grid%alarms_created( alarmid ) = .FALSE.
             ENDIF
           ENDDO
         ENDIF
         CALL Setup_Timekeeping ( head_grid )
         head_grid%start_subtime = domain_get_start_time ( head_grid )
         head_grid%stop_subtime = domain_get_stop_time ( head_grid )

         ! Force to turn off history output in this case
         CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

         IF ( xtraj_boundary_form .gt. 0 .and. &
              model_config_rec%interval_seconds .gt. 60*head_grid%checkpoint_interval ) THEN
            CALL WRFU_AlarmDisable( head_grid%alarms( BOUNDARY_ALARM ), rc=rc )
         ENDIF

         IF ( icheck .GT. 1 ) THEN
            !Disable TL auxinput reading when input interval >= checkpoint_interval
            !or when the checkpoint_interval is not an integer multiple of the
            !input interval
            DO ialarm = first_auxinput, last_auxinput
               CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
!               IF ( interval .gt. 60*head_grid%checkpoint_interval .OR. &
               IF ( interval .ge. 60*head_grid%checkpoint_interval .OR. &
                    ( interval .lt. 60*head_grid%checkpoint_interval .AND. &
                      mod(60*head_grid%checkpoint_interval, interval) .ne. 0 ) )THEN
                  CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
               ENDIF
            ENDDO
         ENDIF

         CALL domain_clock_get( head_grid, start_timestr=timestr )
         CALL domain_clock_set( head_grid, current_timestr=timestr )

         IF ( config_flags%mp_physics .NE. 0 .AND. &
              config_flags%mp_physics .NE. 99 .AND. &
              config_flags%mp_physics .NE. 98 ) &
            CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
         IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. &
              config_flags%bl_pbl_physics .GT. 0 ) &
            CALL nl_set_bl_pbl_physics (head_grid%id, 98)
         IF ( config_flags%cu_physics .GT. 0 ) THEN
            CALL nl_set_cu_physics (head_grid%id, 98)
            head_grid%cudt = 0
         ENDIF
         CALL nl_set_ra_lw_physics (head_grid%id, 0)
         CALL nl_set_ra_sw_physics (head_grid%id, 0)
         IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. &
              config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) &
            CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
         IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. &
              config_flags%sf_surface_physics.NE.PXLSMSCHEME) &
            CALL nl_set_sf_surface_physics (head_grid%id, 0)

         IF (icheck .EQ. 1) THEN
            CALL start_domain ( head_grid , .TRUE. )
         ENDIF
       
         IF ( ASSOCIATED(xtraj_tail) ) xtraj_pointer => xtraj_tail
         CALL integrate ( head_grid )

         CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )

         start_year = end_year
         start_month = end_month
         start_day = end_day
         start_hour = end_hour
         start_minute = end_minute

      ENDDO
      ! Release linked list for trajectory
      CALL xtraj_io_initialize

      model_config_rec%run_days = run_days
      model_config_rec%run_hours = run_hours
      model_config_rec%run_minutes = run_minutes

      model_config_rec%end_year(head_grid%id) = end_year0
      model_config_rec%end_month(head_grid%id) = end_month0
      model_config_rec%end_day(head_grid%id) = end_day0
      model_config_rec%end_hour(head_grid%id) = end_hour0
      model_config_rec%end_minute(head_grid%id) = end_minute0

      model_config_rec%start_year(head_grid%id) = start_year0
      model_config_rec%start_month(head_grid%id) = start_month0
      model_config_rec%start_day(head_grid%id) = start_day0
      model_config_rec%start_hour(head_grid%id) = start_hour0
      model_config_rec%start_minute(head_grid%id) = start_minute0

      model_config_rec%rst_inname = "wrfrst_d<domain>_<date>"
      IF(head_grid%cycle_interval > 0) then
         model_config_rec%input_inname = "wrfinput_d<domain>"
      ENDIF

      !model_config_rec%interval_seconds = 3600*(run_hours + 24*run_days)
      ncheck = 1
      head_grid%check_count = 1

      IF ( xtraj_boundary_form .gt. 0 ) THEN
         ! Close boundary file, as it will be read again from start later
         CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                                                  "DATASET=BOUNDARY" )
      ENDIF
      model_config_rec%restart = .false.
      head_grid%restart = .false.
      call nl_set_restart ( head_grid%id, .false. )
      model_config_rec%io_form_restart = 0

      ! Restore input boundary form
      CALL nl_set_io_form_boundary( head_grid%id, input_boundary_form )  

!END TLM CHECKPOINTING
   ELSE
      ! Setting for TL model
      model_config_rec%dyn_opt = dyn_em_tl

      CALL nl_set_auxinput16_interval_s (head_grid%id, timestep )
      CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
      CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

      IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
        IF ( head_grid%domain_clock_created ) THEN
          CALL WRFU_ClockDestroy( head_grid%domain_clock )
          head_grid%domain_clock_created = .FALSE.
        ENDIF
      ENDIF
      IF ( ASSOCIATED( head_grid%alarms ) .AND. &
           ASSOCIATED( head_grid%alarms_created ) ) THEN
        DO alarmid = 1, MAX_WRF_ALARMS
          IF ( head_grid%alarms_created( alarmid ) ) THEN
            CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
            head_grid%alarms_created( alarmid ) = .FALSE.
          ENDIF
        ENDDO
      ENDIF
      CALL Setup_Timekeeping ( head_grid )
      head_grid%start_subtime = domain_get_start_time ( head_grid )
      head_grid%stop_subtime = domain_get_stop_time ( head_grid )

      ! Force to turn off history output in this case
      CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

      ! Force to turn off boundary input
      CALL nl_set_io_form_boundary( head_grid%id, 0 )

      CALL domain_clock_get( head_grid, start_timestr=timestr )
      CALL domain_clock_set( head_grid, current_timestr=timestr )

      head_grid%itimestep = 0
#if (WRF_CHEM == 1)
      head_grid%ktauc = head_grid%itimestep
#endif

      CALL start_domain ( head_grid , .TRUE. )

      IF ( ASSOCIATED(xtraj_tail) ) xtraj_pointer => xtraj_tail
      CALL integrate ( head_grid )

   ENDIF

   CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_set_io_form_auxhist17( head_grid%id, io_auxh17 )
   ENDIF

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )
   head_grid%start_subtime = domain_get_start_time ( head_grid )
   head_grid%stop_subtime = domain_get_stop_time ( head_grid )


   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id , 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status , ierr ) 
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )
   ENDIF

   CALL       wrf_message ( "wrf: back from tangent linear integrate" )

   END SUBROUTINE wrf_run_tl_checkpt

   SUBROUTINE wrf_run_tl_standalone( )
!<DESCRIPTION>
!     WRF tangent linear code standalone run
!</DESCRIPTION>

!</DESCRIPTION>

   INTEGER :: rc, time_step, id, ierr
   CHARACTER*256 :: timestr

   ! Return immediately if not dyn_em_tl
   IF ( config_flags%dyn_opt .NE. dyn_em_tl ) RETURN

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   IF ( head_grid%trajectory_io ) THEN

   ! Force to read the lateral boundary condition here.
   !CALL med_before_solve_io ( head_grid, config_flags )

   ! Close boundary file, as it will be read again from start later
   !CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )

   CALL init_domain_size ( head_grid, config_flags )

   ! Release linked list for trajectory
   call xtraj_io_initialize

   ! Initialize linked list for adjoint forcing and tl. pertbation
   call adtl_initialize

   CALL       wrf_message ( "wrf: calling nonlinear integrate" )

   ! Set up basic states output
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, time_step )
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL wrf_run

   ! Turn off basic states output
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL       wrf_message ( "wrf: back from nonlinear integrate" )

   ENDIF 

   ! Set the file names and interval for reading basic states.
   model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxinput16_interval_s (head_grid%id, time_step )
   CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

   IF ( config_flags%mp_physics .NE. 0 .and. config_flags%mp_physics .NE. 99 .and. &
        config_flags%mp_physics .NE. 98 ) CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
   IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. config_flags%bl_pbl_physics .GT. 0 ) CALL nl_set_bl_pbl_physics (head_grid%id, 98)
   IF ( config_flags%cu_physics .GT. 0 ) THEN
      CALL nl_set_cu_physics (head_grid%id, 98)
      head_grid%cudt = 0
   ENDIF
   CALL nl_set_ra_lw_physics (head_grid%id, 0)
   CALL nl_set_ra_sw_physics (head_grid%id, 0)
   IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
   IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. config_flags%sf_surface_physics.NE.PXLSMSCHEME) CALL nl_set_sf_surface_physics (head_grid%id, 0)

   head_grid%g_u_1 = 0.0D0
   head_grid%g_v_1 = 0.0D0
   head_grid%g_w_1 = 0.0D0
   head_grid%g_t_1 = 0.0D0
   head_grid%g_ph_1 = 0.0D0
   head_grid%g_mu_1 = 0.0D0
   
   head_grid%g_u_2 = 0.0D0
   head_grid%g_v_2 = 0.0D0
   head_grid%g_w_2 = 0.0D0
   head_grid%g_t_2 = 0.0D0
   head_grid%g_ph_2 = 0.0D0
   head_grid%g_mu_2 = 0.0D0

   head_grid%g_p = 0.0D0
   
   head_grid%g_moist = 0.0D0
   head_grid%g_tracer = 0.0D0
   head_grid%g_scalar = 0.0D0
   head_grid%g_rainnc  = 0.0D0
   head_grid%g_rainncv = 0.0D0
   head_grid%g_rainc  = 0.0D0
   head_grid%g_raincv = 0.0D0
#if (WRF_CHEM == 1)
   head_grid%g_chem = 0.0D0
   head_grid%g_scaleant = 0.0D0
   head_grid%g_scalebb = 0.0D0
   head_grid%g_pm2_5_dry = 0.0D0
   head_grid%g_pm2_5_dry_ec = 0.0D0
   head_grid%g_pm10 = 0.0D0
#endif

   CALL domain_clock_get( head_grid, start_timestr=timestr )
   CALL domain_clock_set( head_grid, current_timestr=timestr )
   config_flags%auxinput9_inname = "init_pert_d01"
   config_flags%io_form_auxinput9 = 2
   CALL nl_set_io_form_auxinput9 ( head_grid%id, 2 )
   config_flags%frames_per_auxinput9 = 1
   CALL med_auxinput_in ( head_grid, auxinput9_alarm, config_flags )
   CALL  wrf_debug ( 0 , 'Read in initial perturbation' )
   config_flags%io_form_auxinput9 = 0

   CALL wrf_run_tl

   ! Release linked list for trajectory
   call xtraj_io_initialize

   END SUBROUTINE wrf_run_tl_standalone

   SUBROUTINE wrf_run_ad( )
!<DESCRIPTION>
!     WRF adjoint run routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  The start
! and stop times for the domain are set to the start and stop time of the
! model run, and then <a href=integrate.html>integrate</a> is called to
! advance the domain forward through that specified time interval.  On
! return, the simulation is completed.  
! 
!</DESCRIPTION>

   CHARACTER*256 :: timestr, timestr1
   INTEGER :: start_year,start_month,start_day,start_hour,start_minute,start_second
   INTEGER :: end_year,end_month,end_day,end_hour,end_minute,end_second
   INTEGER :: rc, runad
   REAL    :: timestep
   TYPE(WRFU_TimeInterval) :: run_interval
   INTEGER       :: io_auxh18
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_message ( "wrf: calling adjoint integrate" )

   ! Seeting for AD model
   model_config_rec%dyn_opt = dyn_em_ad

   model_config_rec%run_days = -1 * model_config_rec%run_days
   model_config_rec%run_hours = -1 * model_config_rec%run_hours
   model_config_rec%run_minutes = -1 * model_config_rec%run_minutes
   model_config_rec%run_seconds = -1 * model_config_rec%run_seconds

   start_year = model_config_rec%start_year(head_grid%id)
   start_month = model_config_rec%start_month(head_grid%id)
   start_day = model_config_rec%start_day(head_grid%id)
   start_hour = model_config_rec%start_hour(head_grid%id)
   start_minute = model_config_rec%start_minute(head_grid%id)
   start_second = model_config_rec%start_second(head_grid%id)

   end_year = model_config_rec%end_year(head_grid%id)
   end_month = model_config_rec%end_month(head_grid%id)
   end_day = model_config_rec%end_day(head_grid%id)
   end_hour = model_config_rec%end_hour(head_grid%id)
   end_minute = model_config_rec%end_minute(head_grid%id)
   end_second = model_config_rec%end_second(head_grid%id)

   model_config_rec%start_year(head_grid%id)  = end_year
   model_config_rec%start_month(head_grid%id)  = end_month
   model_config_rec%start_day(head_grid%id)  = end_day
   model_config_rec%start_hour(head_grid%id)  = end_hour
   model_config_rec%start_minute(head_grid%id)  = end_minute
   model_config_rec%start_second(head_grid%id)  = end_second

   model_config_rec%end_year(head_grid%id)    = start_year
   model_config_rec%end_month(head_grid%id)    = start_month
   model_config_rec%end_day(head_grid%id)    = start_day
   model_config_rec%end_hour(head_grid%id)    = start_hour
   model_config_rec%end_minute(head_grid%id)    = start_minute
   model_config_rec%end_second(head_grid%id)    = start_second

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )
   head_grid%start_subtime = domain_get_start_time ( head_grid )
   head_grid%stop_subtime = domain_get_stop_time ( head_grid )

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   CALL domain_clock_get( head_grid, start_timestr=timestr )
   CALL domain_clock_set( head_grid, current_timestr=timestr )
   CALL domain_clock_set( head_grid, time_step_seconds=-1*model_config_rec%time_step )

   IF ( ASSOCIATED(xtraj_head) ) xtraj_pointer => xtraj_head%next

   ! head_grid%itimestep should be the last
   IF ( head_grid%itimestep .EQ. 0 ) THEN
      timestep=abs(head_grid%dt)
      run_interval = head_grid%stop_subtime - head_grid%start_subtime

      CALL WRFU_TimeIntervalGet( run_interval, S=runad, rc=rc )
      runad = abs(runad)

      head_grid%itimestep = ceiling(real(runad)/timestep)
#if (WRF_CHEM == 1)
      head_grid%ktauc = head_grid%itimestep
#endif
   ENDIF

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_get_io_form_auxhist18( head_grid%id, io_auxh18 )
      CALL nl_set_io_form_auxhist18( head_grid%id, 0 )
   ENDIF

   CALL integrate ( head_grid )

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_set_io_form_auxhist18( head_grid%id, io_auxh18 )
   ENDIF

   CALL start_domain ( head_grid , .TRUE. )

   IF ( .NOT. config_flags%trajectory_io .OR. gradient_out ) THEN
      config_flags%auxhist17_outname = &
          "gradient_wrfplus_d<domain>_<date>"
      config_flags%io_form_auxhist17 = 2
      CALL nl_set_io_form_auxhist17 ( head_grid%id, 2 )
      CALL  med_hist_out ( head_grid , AUXHIST17_ALARM , config_flags )
      CALL nl_set_io_form_auxhist17 ( head_grid%id, 0 )
      CALL  wrf_debug ( 0 , 'Output gradient in WRFPLUS &
                            &(not including LBC gradient)' )
   ENDIF

!  Recover to NL model 
   model_config_rec%dyn_opt = dyn_em

   model_config_rec%run_days = -1 * model_config_rec%run_days
   model_config_rec%run_hours = -1 * model_config_rec%run_hours
   model_config_rec%run_minutes = -1 * model_config_rec%run_minutes
   model_config_rec%run_seconds = -1 * model_config_rec%run_seconds

   model_config_rec%start_year(head_grid%id)  = start_year
   model_config_rec%start_month(head_grid%id)  = start_month
   model_config_rec%start_day(head_grid%id)  = start_day
   model_config_rec%start_hour(head_grid%id)  = start_hour
   model_config_rec%start_minute(head_grid%id)  = start_minute
   model_config_rec%start_second(head_grid%id)  = start_second

   model_config_rec%end_year(head_grid%id)    = end_year
   model_config_rec%end_month(head_grid%id)    = end_month
   model_config_rec%end_day(head_grid%id)    = end_day
   model_config_rec%end_hour(head_grid%id)    = end_hour
   model_config_rec%end_minute(head_grid%id)    = end_minute
   model_config_rec%end_second(head_grid%id)    = end_second

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )
   head_grid%start_subtime = domain_get_start_time ( head_grid )
   head_grid%stop_subtime = domain_get_stop_time ( head_grid )
   CALL domain_clock_set( head_grid, time_step_seconds=model_config_rec%time_step )

   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id , 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status , ierr )
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                           "DATASET=BOUNDARY" )
   ENDIF

   CALL       wrf_message ( "wrf: back from adjoint integrate" )

   END SUBROUTINE wrf_run_ad

   SUBROUTINE wrf_run_ad_checkpt( )
!<DESCRIPTION>
!     WRF adjoint run routine.
!</DESCRIPTION>

!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  The start
! and stop times for the domain are set to the start and stop time of the
! model run, and then <a href=integrate.html>integrate</a> is called to
! advance the domain forward through that specified time interval.  On
! return, the simulation is completed.  
! 
!</DESCRIPTION>

   CHARACTER*256 :: timestr, timestr1
   CHARACTER*128               :: DataSet
   INTEGER                     :: io_form
   INTEGER :: start_year,start_month,start_day,start_hour,start_minute,start_second
   INTEGER :: end_year,end_month,end_day,end_hour,end_minute,end_second
   INTEGER :: rc, runad
   INTEGER    :: timestep, ialarm, interval
   TYPE(WRFU_TimeInterval) :: run_interval
   INTEGER       :: io_auxh18
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status

!checkpoint variables
   integer                           :: icheck
   integer                           :: run_days, run_hours, run_minutes
   character(len=256)                :: STOPtimestr, STARTtimestr, checkend_timestr
   character(len=256)                :: debugstr
   character(len=1000)               :: msg 
   integer                           :: checkpoint_interval_m, checkpoint_interval_h, checkpoint_interval_d, ncheck
   integer :: start_year0,start_month0,start_day0,start_hour0,start_minute0
   integer :: end_year0,end_month0,end_day0,end_hour0,end_minute0
   integer , dimension(12) :: days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_message ( "wrf: calling adjoint integrate" )

   start_year = model_config_rec%start_year(head_grid%id)
   start_month = model_config_rec%start_month(head_grid%id)
   start_day = model_config_rec%start_day(head_grid%id)
   start_hour = model_config_rec%start_hour(head_grid%id)
   start_minute = model_config_rec%start_minute(head_grid%id)
   start_second = model_config_rec%start_second(head_grid%id)

   end_year = model_config_rec%end_year(head_grid%id)
   end_month = model_config_rec%end_month(head_grid%id)
   end_day = model_config_rec%end_day(head_grid%id)
   end_hour = model_config_rec%end_hour(head_grid%id)
   end_minute = model_config_rec%end_minute(head_grid%id)
   end_second = model_config_rec%end_second(head_grid%id)

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_get_io_form_auxhist8( head_grid%id, io_auxh18 )
      CALL nl_set_io_form_auxhist8( head_grid%id, 0 )
   ENDIF

   CALL nl_get_time_step ( head_grid%id, timestep)  

   IF(head_grid%checkpoint_interval > 0) then
!BEGIN ADM CHECKPOINTING
      CALL       wrf_message ( "wrf: using multiple checkpoint windows" )
      model_config_rec%rst_inname = "wrf_checkpoint_d<domain>_<date>"
      IF(head_grid%cycle_interval > 0) then
         model_config_rec%input_inname = "wrfinput_d<domain>_<date>"
      ENDIF

      IF (head_grid%checkpoint_interval > 180) then
         CALL wrf_message("WARNING: large checkpoint_interval may &
              &overload memory.  If program crashes, consider reducing &
              &checkpoint_interval and/or number of processors.")
      END IF

      ! Store input boundary form
      CALL nl_get_io_form_boundary( head_grid%id, input_boundary_form )

      IF (.not. head_grid%boundary_io .or. input_boundary_form .eq. 0) xtraj_boundary_form=0
      IF ( model_config_rec%interval_seconds .lt. 60*head_grid%checkpoint_interval .or. &
           mod( model_config_rec%interval_seconds,60*head_grid%checkpoint_interval ) .ne. 0) THEN
         ! Force to turn off boundary input for trajectory calculation
         xtraj_boundary_form = 0
         CALL wrf_message ("WARNING: interval_seconds must be divisible by &
              &60*checkpoint_interval to use boundary updating in a TL or AD &
              &simulation. Turning off boundary updating.")
      ENDIF

      !CALL nl_set_restart_interval_m ( head_grid%id, (ncheck+1)*head_grid%checkpoint_interval/60)
      CALL nl_set_restart_interval_m ( head_grid%id, 0)

      ! Store original time keeping variables
      CALL nl_get_run_days  (head_grid%id,run_days)
      CALL nl_get_run_hours (head_grid%id,run_hours)
      CALL nl_get_run_minutes (head_grid%id,run_minutes)
      ncheck = (run_hours + 24 * run_days)
      IF (REAL(ncheck,8) .lt. REAL(head_grid%checkpoint_interval,8)/60.) THEN
         CALL wrf_error_fatal(&
                 'Need to set total run_hours >= checkpoint_interval/60')
      ENDIF
      ncheck = ncheck*60 / head_grid%checkpoint_interval

      end_year0  = end_year
      end_month0 = end_month
      end_day0   = end_day
      end_hour0  = end_hour
      end_minute0  = end_minute

      start_year0 = start_year
      start_month0 = start_month
      start_day0 = start_day
      start_hour0 = start_hour
      start_minute0 = start_minute

      checkpoint_interval_m = mod(head_grid%checkpoint_interval,60)
      checkpoint_interval_h = (head_grid%checkpoint_interval - checkpoint_interval_m) / 60
      IF ( checkpoint_interval_h > 24 ) THEN
         checkpoint_interval_d = (checkpoint_interval_h - mod(checkpoint_interval_h,24)) / 24
         checkpoint_interval_h = mod(checkpoint_interval_h,24)
      ELSE
         checkpoint_interval_d = 0
      ENDIF

      CALL domain_clock_get (head_grid, stop_timestr=STOPtimestr)
      CALL domain_clock_get (head_grid, start_timestr=STARTtimestr)
      checkend_timestr=STOPtimestr
 
      model_config_rec%run_days = checkpoint_interval_d
      model_config_rec%run_hours = checkpoint_interval_h
      model_config_rec%run_minutes = checkpoint_interval_m
  
      !!Assuming all forcing has already been calculated read/write settings are in place

      DO icheck = 1,ncheck,1
         head_grid%check_count = ncheck - icheck + 1
         write(msg,FMT='(A,I5,A,I5)')'WRF ADM: Simulating checkpoint ',head_grid%check_count,' of ',ncheck
         CALL wrf_debug(1,TRIM(msg))

         !Set clock to checkpoint initial condition time stamp for I/O
         CALL domain_clock_set (head_grid, current_timestr=checkend_timestr )
         CALL domain_clock_set (head_grid, stop_timestr=checkend_timestr)
         CALL domain_clock_set (head_grid, time_step_seconds=-60*head_grid%checkpoint_interval)
         CALL domain_clockadvance (head_grid)
         CALL domain_clock_get (head_grid, current_timestr=checkend_timestr)
         CALL domain_clock_set (head_grid, time_step_seconds=timestep)
         CALL domain_clockprint(150, head_grid, &
                                'get CurrTime from clock,')

         !Read initial conditions
!         if ( icheck .eq. ncheck ) then
         if ( icheck .eq. ncheck .OR. (head_grid%cycle_interval .GT. 0 .AND. &
              mod(head_grid%check_count*head_grid%checkpoint_interval, &
                  head_grid%cycle_interval) .eq. 0 ) ) then
            model_config_rec%restart = .false.
            head_grid%restart = .false.
            call nl_set_restart ( head_grid%id, .false. )
            model_config_rec%io_form_restart = 0
            head_grid%itimestep = 0
         else
            model_config_rec%restart = .true.
            head_grid%restart = .true.
            call nl_set_restart ( head_grid%id, .true. )
            model_config_rec%io_form_restart = 2
            head_grid%itimestep = 1
         endif

         CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )

!#ifdef DM_PARALLEL
!         CALL get_config_as_buffer( configbuf, configbuflen, nbytes )
!         CALL wrf_dm_bcast_bytes( configbuf, nbytes )
!         CALL set_config_as_buffer( configbuf, configbuflen )
!#endif

         CALL med_initialdata_input( head_grid , config_flags )
         model_config_rec%io_form_restart = 0
         head_grid%itimestep = 0
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
         !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
         ! which is necessary to do 4DVAR -JJG 9-23-2014
         IF(.not.model_config_rec%var4d_run) then
            head_grid%scaleant = 1.0D0
            head_grid%scalebb = 1.0D0
         ENDIF
#endif

         ! Set time keeping variables to checkpoint values
         model_config_rec%end_year(head_grid%id) = end_year
         model_config_rec%end_month(head_grid%id) = end_month
         model_config_rec%end_day(head_grid%id) = end_day
         model_config_rec%end_hour(head_grid%id) = end_hour
         model_config_rec%end_minute(head_grid%id) = end_minute

         if (end_minute - checkpoint_interval_m < 0) then
            start_minute = end_minute - checkpoint_interval_m + 60
            start_hour  = end_hour  - 1
         else
            start_minute = end_minute - checkpoint_interval_m
            start_hour  = end_hour
         end if

         if (start_hour - checkpoint_interval_h < 0) then
            start_hour = start_hour - checkpoint_interval_h + 24
            start_day  = end_day  - 1
         else
            start_hour = start_hour - checkpoint_interval_h
            start_day  = end_day
         end if

         start_year  = end_year
         start_month = end_month

         !Account for leap years
         if(end_month .eq. 3 .and. end_day < 5) then  !Assuming the checkpoint interval < 5 days
            if ( mod(end_year,4) .eq. 0 )then
               if ( mod(end_year,100) .eq. 0 ) then
                  if ( mod(end_year,400) .eq. 0 ) then
                     days_in_month(2) = 29
                  else
                     days_in_month(2) = 28
                  end if
               else
                  days_in_month(2) = 29
               end if
            else
               days_in_month(2) = 28
            end if
         end if

         if (start_day - checkpoint_interval_d < 1) then
            if(end_month - 1 < 1) then
               start_month = 12
               start_year  = end_year - 1
            else
               start_month = end_month - 1
            end if
            start_day = start_day - checkpoint_interval_d + days_in_month(start_month)
         else
            start_day = start_day - checkpoint_interval_d
         end if

         model_config_rec%start_year(head_grid%id) = start_year
         model_config_rec%start_month(head_grid%id) = start_month
         model_config_rec%start_day(head_grid%id) = start_day
         model_config_rec%start_hour(head_grid%id) = start_hour
         model_config_rec%start_minute(head_grid%id) = start_minute

! Get nonlinear model trajectory
         CALL wrf_message ( "wrf: calculating checkpoint trajectory" )

         ! Set trajectory boundary form
         CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )  

         CALL nl_set_auxhist6_interval_s ( head_grid%id, timestep )
         CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
         CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )

         IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
           IF ( head_grid%domain_clock_created ) THEN
             CALL WRFU_ClockDestroy( head_grid%domain_clock )
             head_grid%domain_clock_created = .FALSE.
           ENDIF
         ENDIF
         IF ( ASSOCIATED( head_grid%alarms ) .AND. &
              ASSOCIATED( head_grid%alarms_created ) ) THEN
           DO alarmid = 1, MAX_WRF_ALARMS
             IF ( head_grid%alarms_created( alarmid ) ) THEN
               CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
               head_grid%alarms_created( alarmid ) = .FALSE.
             ENDIF
           ENDDO
         ENDIF
         CALL Setup_Timekeeping ( head_grid )
         head_grid%start_subtime = domain_get_start_time ( head_grid )
         head_grid%stop_subtime = domain_get_stop_time ( head_grid )

         ! Force to turn off history output in this case
         CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

         model_config_rec%dyn_opt = dyn_em_check

         IF ( wrf_dm_on_monitor() ) THEN
            CALL start_timing
         END IF

         ! Enable boundary reading for checkpoint intervals shorter than
         ! lbc reading interval
         IF ( xtraj_boundary_form .gt. 0 .and. &
              model_config_rec%interval_seconds .gt. 60*head_grid%checkpoint_interval ) THEN
            IF ( mod( end_hour, model_config_rec%interval_seconds/3600 ) .eq. 0 .and. &
                 end_minute .eq. 0 ) THEN
               CALL domain_clock_get( head_grid, stop_timestr=timestr )
               CALL domain_clock_set( head_grid, current_timestr=timestr )
               CALL domain_clock_set (head_grid, time_step_seconds=-1*model_config_rec%interval_seconds)
               CALL domain_clockadvance (head_grid)
!               model_config_rec%dyn_opt = dyn_em
               CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
               CALL med_latbound_in ( head_grid, config_flags )
!               CALL med_before_solve_io ( head_grid, config_flags )
!               CALL domain_clock_set (head_grid, time_step_seconds=model_config_rec%interval_seconds)
!               CALL domain_clockadvance (head_grid)
               CALL domain_clock_set (head_grid, time_step_seconds=timestep)
               CALL domain_clockprint(150, head_grid, &
                                      'get CurrTime from clock,')
            ENDIF 
            CALL WRFU_AlarmDisable( head_grid%alarms( BOUNDARY_ALARM ), rc=rc )
         ENDIF

         !Update time since last boundary read
         head_grid%dtbc = REAL( mod( 60 * (head_grid%check_count-1) * &
                               head_grid%checkpoint_interval, &
                               model_config_rec%interval_seconds ), 8)
!         head_grid%dtbc = REAL( mod( 3600 * start_hour, model_config_rec%interval_seconds ) &
!                          + 60 * start_minute , 8)


         ! Enable auxinput streams for chekpoint_interval < io_intervals
         IF ( icheck .LT. ncheck ) THEN
            DO ialarm = first_auxinput, last_auxinput
               CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
               IF ( interval .gt. 60*head_grid%checkpoint_interval ) THEN
                  IF ( mod(end_hour*3600+end_minute*60, interval) .eq. 0 .AND. &
                       mod(interval, 60*head_grid%checkpoint_interval) .eq. 0) THEN
                     CALL domain_clock_get( head_grid, stop_timestr=timestr )
                     CALL domain_clock_set( head_grid, current_timestr=timestr )
                     CALL domain_clock_set (head_grid, time_step_seconds=-1*interval)

                     head_grid%domain_clock%clockint%CurrTime = ESMF_TimeInc (&
                                head_grid%domain_clock%clockint%CurrTime, &
                                head_grid%domain_clock%clockint%TimeStep )

!                     CALL domain_clockadvance (head_grid)
   !               model_config_rec%dyn_opt = dyn_em
                     !CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )
                     CALL domain_clock_set (head_grid, time_step_seconds=timestep)
                     CALL domain_clockprint(1, head_grid, &
                                         'get CurrTime from clock,')


                     IF      ( ialarm-first_input .EQ. 10 ) THEN
                       WRITE(DataSet,'("AUXINPUT10")')
                     ELSE IF ( ialarm-first_input .EQ. 11 ) THEN
                       WRITE(DataSet,'("AUXINPUT11")')
                     ELSE IF ( ialarm-first_input .GE. 10 ) THEN
                       WRITE(DataSet,'("AUXINPUT",I2)')ialarm-first_input
                     ELSE
                       WRITE(DataSet,'("AUXINPUT",I1)')ialarm-first_input
                     ENDIF

                     io_form = io_form_for_dataset( DataSet )

                     write(msg,'("Checkpoint reading for ",A,&
                               &", nframes = ",I2," interval_s = ",I8,&
                               &"io_form = ",I2)') &
                           TRIM(DataSet), head_grid%nframes(ialarm),interval,io_form
                     CALL wrf_debug(1, TRIM(msg))

                     !SST Update
                     IF ( ialarm .EQ. AUXINPUT4_ALARM &
                          .AND. config_flags%io_form_auxinput4 .GT. 0 ) THEN
                        CALL med_auxinput_in ( head_grid, ialarm, config_flags )
#if ( WRF_CHEM == 1 )
                     !Anthropogenic Emissions
                     ELSEIF ( ialarm .EQ. AUXINPUT5_ALARM &
                          .AND. config_flags%io_form_auxinput5 .GT. 0 &
                          .AND. config_flags%emiss_inpt_opt .NE. 0 ) THEN
                        CALL med_read_wrf_chem_emiss ( head_grid, config_flags )

                     !Burning Emissions
                     ELSEIF ( ialarm .EQ. AUXINPUT7_ALARM &
                          .AND. config_flags%io_form_auxinput7 .GT. 0 &
                          .AND. config_flags%biomass_burn_opt .NE. 0 ) THEN
                        CALL med_auxinput_in ( head_grid, ialarm, config_flags )

                     !gocart_bg
                     ELSEIF ( ialarm .EQ. AUXINPUT8_ALARM &
                          .AND. config_flags%io_form_auxinput8 .GT. 0 ) THEN
!                        CALL nl_get_io_form_auxinput8(head_grid%id, io_form)
!                        write(msg,'("IO_FORM_AUXINPUT",I2," = ",I8,I8)') &
!                              ialarm-first_input,config_flags%io_form_auxinput8,io_form
!                        CALL wrf_debug(1, TRIM(msg))
                        CALL med_auxinput_in ( head_grid, ialarm, config_flags )
#endif
                     ENDIF

                     !Include additional streams as necessary
                     !IF ( ialarm .EQ. AUXINPUTXX_ALARM &
                     !     .AND. config_flags%io_form_auxinputXX .GT. 0 ) THEN
                     !   CALL med_auxinput_in ( head_grid, ialarm, config_flags )
                     !ENDIF
                  ENDIF
                  !Disable auxinput reading for all streams with intervals
                  ! greater than checkpoint_interval
                  CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
               ELSEIF ( mod(60*head_grid%checkpoint_interval, interval) .ne. 0 ) THEN
                  write(msg,FMT='(A,I2)')&
                    'WARNING: when checkpoint_interval > auxinputXX_interval,&
                    &checkpoint must be divisible by auxinput.  &
                    &Disabling auxinput',&
                    ialarm-first_input
                  CALL wrf_message (TRIM(msg))
                  CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
               ENDIF
            ENDDO
         ENDIF

         IF ( wrf_dm_on_monitor() ) THEN
            CALL end_timing('AD checkpoint I/O')
         END IF

         CALL domain_clock_get( head_grid, start_timestr=timestr )
         CALL domain_clock_set( head_grid, current_timestr=timestr )

         ! Release linked list for trajectory
         CALL xtraj_io_initialize

         CALL wrf_run

         ! Turn off basic states output
         CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
         CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )


!NOW THE ADM
         CALL wrf_message( "wrf: running ADM" )

! Need to add a subroutine to calculate the forcing for 
! each checkpoint window.  This will be necessary when the forcing
! depends on the trajectory (nonlinear cost function), and the forcing
! is binned for each 4DVAR window (e.g. hourly), and not calculated online (during each time step).  
! That is the case for standard (non chemical) WRFDA-4DVar forcing.
         !CALL calculate_ad_forcing( head_grid, model_config_rec ) !(Doesn't exist yet.)

         ! Reverse timing variables
         model_config_rec%run_days = -1 * model_config_rec%run_days
         model_config_rec%run_hours = -1 * model_config_rec%run_hours
         model_config_rec%run_minutes = -1 * model_config_rec%run_minutes

         model_config_rec%start_year(head_grid%id)  = end_year
         model_config_rec%start_month(head_grid%id)  = end_month
         model_config_rec%start_day(head_grid%id)  = end_day
         model_config_rec%start_hour(head_grid%id)  = end_hour
         model_config_rec%start_minute(head_grid%id)  = end_minute

         model_config_rec%end_year(head_grid%id)    = start_year
         model_config_rec%end_month(head_grid%id)    = start_month
         model_config_rec%end_day(head_grid%id)    = start_day
         model_config_rec%end_hour(head_grid%id)    = start_hour
         model_config_rec%end_minute(head_grid%id)    = start_minute

         model_config_rec%dyn_opt = dyn_em_ad

         ! Force to turn off boundary input during derivative model
         CALL nl_set_io_form_boundary( head_grid%id, 0 ) 

         CALL nl_set_auxinput16_interval_s (head_grid%id, timestep )
         CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
         CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

         !NEED THIS ONCE FORCINGS ARE DEFINED IN 4DVAR BINS (WINDOWS)
         !! Set the file names and interval for reading adjoint forcing.
         !model_config_rec%auxinput17_inname = "af_d<domain>_<date>"
         !call nl_set_auxinput17_interval_s ( head_grid%id, var4d_bin )
         !call nl_set_io_form_auxinput17 ( head_grid%id, 2 )
         !call nl_set_frames_per_auxinput17 ( head_grid%id, 1 )

         IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
           IF ( head_grid%domain_clock_created ) THEN
             CALL WRFU_ClockDestroy( head_grid%domain_clock )
             head_grid%domain_clock_created = .FALSE.
           ENDIF
         ENDIF
         IF ( ASSOCIATED( head_grid%alarms ) .AND. &
              ASSOCIATED( head_grid%alarms_created ) ) THEN
           DO alarmid = 1, MAX_WRF_ALARMS
             IF ( head_grid%alarms_created( alarmid ) ) THEN
               CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
               head_grid%alarms_created( alarmid ) = .FALSE.
             ENDIF
           ENDDO
         ENDIF
         CALL Setup_Timekeeping ( head_grid )
         head_grid%start_subtime = domain_get_start_time ( head_grid )
         head_grid%stop_subtime = domain_get_stop_time ( head_grid )

         ! Force to turn off history output in this case
         CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

         CALL WRFU_AlarmDisable( head_grid%alarms( BOUNDARY_ALARM ), rc=rc )

         !Disable AD auxinput reading when input interval >= checkpoint_interval
         !or when the checkpoint_interval is not an integer multiple of the
         !input interval
         DO ialarm = first_auxinput, last_auxinput
            CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
            IF ( interval .ge. 60*head_grid%checkpoint_interval .OR. &
                 ( interval .lt. 60*head_grid%checkpoint_interval .AND. &
                   mod(60*head_grid%checkpoint_interval, interval) .ne. 0 ) )THEN
               CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
            ENDIF
         ENDDO

         CALL domain_clock_get( head_grid, start_timestr=timestr )
         CALL domain_clock_set( head_grid, current_timestr=timestr )
         CALL domain_clock_set( head_grid, time_step_seconds=-1*timestep )

         IF ( ASSOCIATED(xtraj_head) ) xtraj_pointer => xtraj_head%next

         ! head_grid%itimestep should be the last
         !run_interval = head_grid%stop_subtime - head_grid%start_subtime
         !CALL WRFU_TimeIntervalGet( run_interval, S=runad, rc=rc )
         !runad = abs(runad)
         !head_grid%itimestep = ceiling(real(runad)/real(timestep))
         head_grid%itimestep = real(head_grid%checkpoint_interval*60/timestep)
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
#endif
         IF ( config_flags%mp_physics .NE. 0 .AND. &
              config_flags%mp_physics .NE. 99 .AND. &
              config_flags%mp_physics .NE. 98 ) &
            CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
         IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. &
              config_flags%bl_pbl_physics .GT. 0 ) &
            CALL nl_set_bl_pbl_physics (head_grid%id, 98)
         IF ( config_flags%cu_physics .GT. 0 ) THEN
            CALL nl_set_cu_physics (head_grid%id, 98)
            head_grid%cudt = 0
         ENDIF
         CALL nl_set_ra_lw_physics (head_grid%id, 0)
         CALL nl_set_ra_sw_physics (head_grid%id, 0)
         IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. &
              config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) &
            CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
         IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. &
              config_flags%sf_surface_physics.NE.PXLSMSCHEME) &
            CALL nl_set_sf_surface_physics (head_grid%id, 0)

         CALL integrate ( head_grid )

         CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )

         IF ( icheck .LT. ncheck .AND. ( .NOT. config_flags%trajectory_io .OR. gradient_out ) ) THEN
            config_flags%auxhist17_outname = "adj_sens_d<domain>_<date>"
            config_flags%io_form_auxhist17 = 2
            CALL nl_set_io_form_auxhist17 ( head_grid%id, 2 )
            config_flags%frames_per_auxhist17 = 1
            CALL  med_hist_out ( head_grid , AUXHIST17_ALARM , config_flags )
            config_flags%io_form_auxhist17 = 0
            CALL nl_set_io_form_auxhist17 ( head_grid%id, 0 )
            CALL  wrf_debug ( 0 , 'Output gradient in WRFPLUS' )
         ENDIF

         model_config_rec%run_days = -1 * model_config_rec%run_days
         model_config_rec%run_hours = -1 * model_config_rec%run_hours
         model_config_rec%run_minutes = -1 * model_config_rec%run_minutes

         end_year = start_year
         end_month = start_month
         end_day = start_day
         end_hour = start_hour
         end_minute = start_minute
      ENDDO
      ! Release linked list for trajectory
      CALL xtraj_io_initialize

      model_config_rec%run_days = run_days
      model_config_rec%run_hours = run_hours
      model_config_rec%run_minutes = run_minutes

      model_config_rec%end_year(head_grid%id) = end_year0
      model_config_rec%end_month(head_grid%id) = end_month0
      model_config_rec%end_day(head_grid%id) = end_day0
      model_config_rec%end_hour(head_grid%id) = end_hour0
      model_config_rec%end_minute(head_grid%id) = end_minute0

      model_config_rec%start_year(head_grid%id) = start_year0
      model_config_rec%start_month(head_grid%id) = start_month0
      model_config_rec%start_day(head_grid%id) = start_day0
      model_config_rec%start_hour(head_grid%id) = start_hour0
      model_config_rec%start_minute(head_grid%id) = start_minute0

      model_config_rec%rst_inname = "wrfrst_d<domain>_<date>"
      IF(head_grid%cycle_interval > 0) then
         model_config_rec%input_inname = "wrfinput_d<domain>"
      ENDIF

      !model_config_rec%interval_seconds = 3600*(run_hours + 24*run_days)
      ncheck = 1
      head_grid%check_count = 1
      IF ( xtraj_boundary_form .gt. 0 ) THEN
         ! Close boundary file, as it will be read again from start later
         CALL close_dataset ( head_grid%lbc_fid , config_flags , &
                                                  "DATASET=BOUNDARY" )
      ENDIF
      model_config_rec%restart = .false.
      head_grid%restart = .false.
      call nl_set_restart ( head_grid%id, .false. )
      model_config_rec%io_form_restart = 0

      ! Restore input boundary form
      CALL nl_set_io_form_boundary( head_grid%id, input_boundary_form )  

!END ADM CHECKPOINTING
   ELSE

      ! Setting for AD model
      model_config_rec%dyn_opt = dyn_em_ad

      model_config_rec%run_days = -1 * model_config_rec%run_days
      model_config_rec%run_hours = -1 * model_config_rec%run_hours
      model_config_rec%run_minutes = -1 * model_config_rec%run_minutes
      model_config_rec%run_seconds = -1 * model_config_rec%run_seconds

      model_config_rec%start_year(head_grid%id)  = end_year
      model_config_rec%start_month(head_grid%id)  = end_month
      model_config_rec%start_day(head_grid%id)  = end_day
      model_config_rec%start_hour(head_grid%id)  = end_hour
      model_config_rec%start_minute(head_grid%id)  = end_minute
      model_config_rec%start_second(head_grid%id)  = end_second

      model_config_rec%end_year(head_grid%id)    = start_year
      model_config_rec%end_month(head_grid%id)    = start_month
      model_config_rec%end_day(head_grid%id)    = start_day
      model_config_rec%end_hour(head_grid%id)    = start_hour
      model_config_rec%end_minute(head_grid%id)    = start_minute
      model_config_rec%end_second(head_grid%id)    = start_second

      CALL nl_set_auxinput16_interval_s (head_grid%id, timestep )
      CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
      CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

      IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
        IF ( head_grid%domain_clock_created ) THEN
          CALL WRFU_ClockDestroy( head_grid%domain_clock )
          head_grid%domain_clock_created = .FALSE.
        ENDIF
      ENDIF
      IF ( ASSOCIATED( head_grid%alarms ) .AND. &
           ASSOCIATED( head_grid%alarms_created ) ) THEN
        DO alarmid = 1, MAX_WRF_ALARMS
          IF ( head_grid%alarms_created( alarmid ) ) THEN
            CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
            head_grid%alarms_created( alarmid ) = .FALSE.
          ENDIF
        ENDDO
      ENDIF
      CALL Setup_Timekeeping ( head_grid )
      head_grid%start_subtime = domain_get_start_time ( head_grid )
      head_grid%stop_subtime = domain_get_stop_time ( head_grid )

      ! Force to turn off history output in this case
      CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

      ! Force to turn off boundary input
      CALL nl_set_io_form_boundary( head_grid%id, 0 )

      CALL domain_clock_get( head_grid, start_timestr=timestr )
      CALL domain_clock_set( head_grid, current_timestr=timestr )
      CALL domain_clock_set( head_grid, time_step_seconds=-1*model_config_rec%time_step )

      IF ( ASSOCIATED(xtraj_head) ) xtraj_pointer => xtraj_head%next

      ! head_grid%itimestep should be the last
!      IF ( head_grid%itimestep .EQ. 0 ) THEN
         run_interval = head_grid%stop_subtime - head_grid%start_subtime
         CALL WRFU_TimeIntervalGet( run_interval, S=runad, rc=rc )
         runad = abs(runad)
         head_grid%itimestep = ceiling(real(runad)/real(timestep))
!      ENDIF
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
#endif

      CALL integrate ( head_grid )

      model_config_rec%run_days = -1 * model_config_rec%run_days
      model_config_rec%run_hours = -1 * model_config_rec%run_hours
      model_config_rec%run_minutes = -1 * model_config_rec%run_minutes
      model_config_rec%run_seconds = -1 * model_config_rec%run_seconds

      model_config_rec%start_year(head_grid%id)  = start_year
      model_config_rec%start_month(head_grid%id)  = start_month
      model_config_rec%start_day(head_grid%id)  = start_day
      model_config_rec%start_hour(head_grid%id)  = start_hour
      model_config_rec%start_minute(head_grid%id)  = start_minute
      model_config_rec%start_second(head_grid%id)  = start_second

      model_config_rec%end_year(head_grid%id)    = end_year
      model_config_rec%end_month(head_grid%id)    = end_month
      model_config_rec%end_day(head_grid%id)    = end_day
      model_config_rec%end_hour(head_grid%id)    = end_hour
      model_config_rec%end_minute(head_grid%id)    = end_minute
      model_config_rec%end_second(head_grid%id)    = end_second

   ENDIF

   head_grid%itimestep = 0
   CALL start_domain ( head_grid , .TRUE. )

   IF ( .NOT. config_flags%trajectory_io .OR. gradient_out .OR. head_grid%checkpoint_interval.GT.0) THEN
      config_flags%auxhist17_outname = "adj_sens_d<domain>_<date>"
      !config_flags%auxhist17_outname = "gradient_wrfplus_d<domain>_<date>"
      config_flags%io_form_auxhist17 = 2
      CALL nl_set_io_form_auxhist17 ( head_grid%id, 2 )
      config_flags%frames_per_auxhist17 = 1
      CALL  med_hist_out ( head_grid , AUXHIST17_ALARM , config_flags )
      CALL nl_set_io_form_auxhist17 ( head_grid%id, 0 )
      CALL  wrf_debug ( 0 , 'Output gradient in WRFPLUS (not including LBC gradient)' )
   ENDIF

   IF ( .NOT. config_flags%trajectory_io ) THEN
      CALL nl_set_io_form_auxhist18( head_grid%id, io_auxh18 )
   ENDIF

   CALL nl_set_auxinput16_interval_s (head_grid%id, 0 )

!  Recover to NL model 
   model_config_rec%dyn_opt = dyn_em

   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )
   head_grid%start_subtime = domain_get_start_time ( head_grid )
   head_grid%stop_subtime = domain_get_stop_time ( head_grid )
   CALL domain_clock_set( head_grid, time_step_seconds=model_config_rec%time_step )

   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id , 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status , ierr )
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )
   ENDIF

   CALL       wrf_message ( "wrf: back from adjoint integrate" )

   END SUBROUTINE wrf_run_ad_checkpt

   SUBROUTINE wrf_run_ad_standalone( )
!<DESCRIPTION>
!     WRF adjoint code standalone run
!</DESCRIPTION>

!</DESCRIPTION>

   INTEGER :: rc, time_step, id, ierr
   CHARACTER*256 :: timestr
   CHARACTER (LEN=80) :: bdyname

   ! Return immediately if not dyn_em_ad
   IF ( config_flags%dyn_opt .NE. dyn_em_ad ) RETURN

   ! Force to turn off history output in this case
   CALL WRFU_AlarmRingerOff( head_grid%alarms( HISTORY_ALARM ), rc=rc )

   IF ( head_grid%trajectory_io ) THEN

   ! Force to read the lateral boundary condition here.
   !CALL med_before_solve_io ( head_grid, config_flags )

   ! Close boundary file, as it will be read again from start later
   !CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )

   CALL init_domain_size ( head_grid, config_flags )

   ! Release linked list for trajectory
   call xtraj_io_initialize

   CALL       wrf_message ( "wrf: calling nonlinear integrate" )

   ! Set up basic states output
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, time_step )
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxhist6 ( head_grid%id, 1 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL wrf_run

   ! Turn off basic states output
   CALL nl_set_io_form_auxhist6 ( head_grid%id, 0 )
   CALL nl_set_auxhist6_interval_s ( head_grid%id, 0 )
   IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
     IF ( head_grid%domain_clock_created ) THEN
       CALL WRFU_ClockDestroy( head_grid%domain_clock )
       head_grid%domain_clock_created = .FALSE.
     ENDIF
   ENDIF
   IF ( ASSOCIATED( head_grid%alarms ) .AND. &
        ASSOCIATED( head_grid%alarms_created ) ) THEN
     DO alarmid = 1, MAX_WRF_ALARMS
       IF ( head_grid%alarms_created( alarmid ) ) THEN
         CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
         head_grid%alarms_created( alarmid ) = .FALSE.
       ENDIF
     ENDDO
   ENDIF
   CALL Setup_Timekeeping ( head_grid )

   CALL       wrf_message ( "wrf: back from nonlinear integrate" )

   ENDIF 

   ! Set the file names and interval for reading basic states.
   model_config_rec%auxinput16_inname = "auxhist6_d<domain>_<date>"
   CALL nl_get_time_step ( head_grid%id, time_step )
   CALL nl_set_auxinput16_interval_s (head_grid%id, time_step )
   CALL nl_set_io_form_auxinput16 ( head_grid%id, 2 )
   CALL nl_set_frames_per_auxinput16 ( head_grid%id, 1 )

   IF ( config_flags%mp_physics .NE. 0 .and. config_flags%mp_physics .NE. 99 .and. &
        config_flags%mp_physics .NE. 98 ) CALL nl_set_mp_physics (head_grid%id, config_flags%mp_physics_ad)
   IF ( config_flags%bl_pbl_physics .NE. ACMPBLSCHEME .AND. config_flags%bl_pbl_physics .GT. 0 ) CALL nl_set_bl_pbl_physics (head_grid%id, 98)
   IF ( config_flags%cu_physics .GT. 0 ) THEN
      CALL nl_set_cu_physics (head_grid%id, 98)
      head_grid%cudt = 0
   ENDIF
   CALL nl_set_ra_lw_physics (head_grid%id, 0)
   CALL nl_set_ra_sw_physics (head_grid%id, 0)
   IF ( config_flags%sf_sfclay_physics.NE.SFCLAYSCHEME .AND. config_flags%sf_sfclay_physics.NE.PXSFCSCHEME) CALL nl_set_sf_sfclay_physics (head_grid%id, 0)
   IF ( config_flags%sf_surface_physics.NE.SLABSCHEME .AND. config_flags%sf_surface_physics.NE.PXLSMSCHEME) CALL nl_set_sf_surface_physics (head_grid%id, 0)

   CALL zero_out_ad ( head_grid )

   head_grid%g_u_1 = 0.0D0
   head_grid%g_v_1 = 0.0D0
   head_grid%g_w_1 = 0.0D0
   head_grid%g_t_1 = 0.0D0
   head_grid%g_ph_1 = 0.0D0
   head_grid%g_mu_1 = 0.0D0
   
   head_grid%g_u_2 = 0.0D0
   head_grid%g_v_2 = 0.0D0
   head_grid%g_w_2 = 0.0D0
   head_grid%g_t_2 = 0.0D0
   head_grid%g_ph_2 = 0.0D0
   head_grid%g_mu_2 = 0.0D0

   head_grid%g_p = 0.0D0
   
   head_grid%g_moist = 0.0D0
   head_grid%g_tracer = 0.0D0
   head_grid%g_scalar = 0.0D0
   head_grid%g_rainnc  = 0.0D0
   head_grid%g_rainncv = 0.0D0
   head_grid%g_rainc  = 0.0D0
   head_grid%g_raincv = 0.0D0
#if (WRF_CHEM == 1)
   head_grid%g_chem = 0.0D0
   head_grid%g_scaleant = 0.0D0
   head_grid%g_scalebb = 0.0D0
   head_grid%g_pm2_5_dry = 0.0D0
   head_grid%g_pm2_5_dry_ec = 0.0D0
   head_grid%g_pm10 = 0.0D0

   !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
   ! which is necessary to do 4DVAR -JJG 9-23-2014
   head_grid%scaleant = 1.0D0
   head_grid%scalebb = 1.0D0
#endif

   CALL domain_clock_get( head_grid, stop_timestr=timestr )
   CALL domain_clock_set( head_grid, current_timestr=timestr )
   config_flags%auxinput17_inname = "final_sens_d01"
   config_flags%io_form_auxinput17 = 2
   CALL nl_set_io_form_auxinput17 ( head_grid%id, 2 )
   config_flags%frames_per_auxinput17 = 1
   CALL med_auxinput_in ( head_grid, auxinput17_alarm, config_flags )
   CALL  wrf_debug ( 0 , 'Read in final sensitivuty' )
   config_flags%io_form_auxinput17 = 0
   CALL domain_clock_get( head_grid, start_timestr=timestr )
   CALL domain_clock_set( head_grid, current_timestr=timestr )

   gradient_out = .TRUE.

   CALL wrf_run_ad

!-- Output gradient in boundary
   CALL construct_filename1( bdyname , 'gradient_wrfbdy' , head_grid%id , 2 )
   CALL open_w_dataset ( id, TRIM(bdyname) , head_grid , config_flags , output_boundary , "DATASET=BOUNDARY", ierr )
   IF ( ierr .NE. 0 ) THEN
      CALL wrf_error_fatal( 'real: error opening wrfbdy for writing' )
   END IF
   print *,'Output LBC gradient valid between these times ',current_date, ' ',start_date
   CALL output_boundary ( id, head_grid , config_flags , ierr )

   ! Release linked list for trajectory
   call xtraj_io_initialize

   END SUBROUTINE wrf_run_ad_standalone

   SUBROUTINE wrf_run( )
!<DESCRIPTION>
!     WRF run routine.
!</DESCRIPTION>

#if (WRFPLUS == 1)
   integer :: io_auxh17, io_auxh18
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status
#endif
!<DESCRIPTION>
! Once the top-level domain has been allocated, configured, and
! initialized, the model time integration is ready to proceed.  The start
! and stop times for the domain are set to the start and stop time of the
! model run, and then <a href=integrate.html>integrate</a> is called to
! advance the domain forward through that specified time interval.  On
! return, the simulation is completed.  
! 
!</DESCRIPTION>

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_debug ( 100 , 'wrf: calling integrate' )
#if (WRFPLUS == 1)

#if (WRF_CHEM == 1)
   !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
   ! which is necessary to do 4DVAR -JJG 9-23-2014
   IF(model_config_rec%dyn_opt .ne. dyn_em_check .and. .not.model_config_rec%var4d_run) then
      head_grid%scaleant = 1.0D0
      head_grid%scalebb = 1.0D0
   ENDIF
#endif

   model_config_rec%dyn_opt = dyn_em
   IF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. SURFDRAGSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
   ELSEIF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. ACMPBLSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
      head_grid%rqcblten = 0.0D0
      head_grid%rqiblten = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. LSCONDSCHEME ) THEN
      head_grid%h_diabatic = 0.0D0
      head_grid%qv_diabatic = 0.0D0
      head_grid%qc_diabatic = 0.0D0
      head_grid%rainnc = 0.0D0
      head_grid%rainncv = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. MKESSLERSCHEME ) THEN
      head_grid%h_diabatic = 0.0
      head_grid%qv_diabatic = 0.0
      head_grid%qc_diabatic = 0.0
      head_grid%rainnc = 0.0
      head_grid%rainncv = 0.0
   ENDIF
   IF ( model_config_rec%cu_physics(head_grid%id) .EQ. DUCUSCHEME ) THEN
      head_grid%rthcuten = 0.0D0
      head_grid%rqvcuten = 0.0D0
      head_grid%rainc = 0.0D0
      head_grid%raincv = 0.0D0
   ENDIF

   CALL nl_get_io_form_auxhist17( head_grid%id, io_auxh17 )
   CALL nl_get_io_form_auxhist18( head_grid%id, io_auxh18 )
   CALL nl_set_io_form_auxhist17( head_grid%id, 0 )
   CALL nl_set_io_form_auxhist18( head_grid%id, 0 )

   head_grid%itimestep = 0
#if (WRF_CHEM == 1)
   head_grid%ktauc = head_grid%itimestep
#endif


   CALL start_domain ( head_grid , .TRUE. )
#endif
   CALL integrate ( head_grid )

#if (WRFPLUS == 1)
   CALL nl_set_io_form_auxhist17( head_grid%id, io_auxh17 )
   CALL nl_set_io_form_auxhist18( head_grid%id, io_auxh18 )

   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id , 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status , ierr )
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY" )
   ENDIF
#endif

   CALL       wrf_debug ( 100 , 'wrf: back from integrate' )

   END SUBROUTINE wrf_run


   SUBROUTINE wrf_run_cycle( )
!<DESCRIPTION>
!     WRF with cycling.  All wrfinput files must be generated externally in
!     wrfinput_d01_YYYY-MM0DD_HH:MM:SS format.
!</DESCRIPTION>

   CHARACTER*256 :: timestr
   INTEGER       :: rc
   INTEGER       :: io_auxh17, io_auxh18
   CHARACTER (LEN=80)      :: bdyname
   INTEGER                 :: open_status

!checkpoint variables
   INTEGER                 :: timestep, ialarm, interval
   integer                           :: icycle, i, j, k, ic, ij
   integer                           :: run_days, run_hours, run_minutes
   character(len=256)                :: STOPtimestr, STARTtimestr, cyclestart_timestr
   character(len=256)                :: debugstr
   character(len=1000)               :: msg 
   integer                           :: cycle_interval_m, cycle_interval_h, cycle_interval_d, ncycle
   INTEGER :: start_year,start_month,start_day,start_hour,start_minute
   INTEGER :: end_year,end_month,end_day,end_hour,end_minute
   integer :: start_year0,start_month0,start_day0,start_hour0,start_minute0
   integer :: end_year0,end_month0,end_day0,end_hour0,end_minute0
   integer , dimension(12) :: days_in_month = (/31,28,31,30,31,30,31,31,30,31,30,31/)

   real, allocatable :: chem_temp(:,:,:,:), scaleant_temp(:,:,:,:), scalebb_temp(:,:,:,:)

   !  The forecast integration for the most coarse grid is now started.  The
   !  integration is from the first step (1) to the last step of the simulation.

   CALL       wrf_message ( "wrf: calling nonlinear integrate" )

#if (WRFPLUS == 1)

#if (WRF_CHEM == 1)
   !Need to remove this once scaleant and scalebb are added to wrfinput or another stream,
   ! which is necessary to do 4DVAR -JJG 9-23-2014
   IF(model_config_rec%dyn_opt .ne. dyn_em_check .and. .not.model_config_rec%var4d_run) then
      head_grid%scaleant = 1.0D0
      head_grid%scalebb = 1.0D0
   ENDIF
#endif

   IF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. SURFDRAGSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
   ELSEIF ( model_config_rec%bl_pbl_physics(head_grid%id) .EQ. ACMPBLSCHEME ) THEN
      head_grid%rublten = 0.0D0
      head_grid%rvblten = 0.0D0
      head_grid%rthblten = 0.0D0
      head_grid%rqvblten = 0.0D0
      head_grid%rqcblten = 0.0D0
      head_grid%rqiblten = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. LSCONDSCHEME ) THEN
      head_grid%h_diabatic = 0.0D0
      head_grid%rainnc = 0.0D0
      head_grid%rainncv = 0.0D0
   ENDIF
   IF ( model_config_rec%mp_physics(head_grid%id) .EQ. MKESSLERSCHEME ) THEN
      head_grid%h_diabatic = 0.0
      head_grid%rainnc = 0.0
      head_grid%rainncv = 0.0
   ENDIF
   IF ( model_config_rec%cu_physics(head_grid%id) .EQ. DUCUSCHEME ) THEN
      head_grid%rthcuten = 0.0D0
      head_grid%rqvcuten = 0.0D0
      head_grid%rainc = 0.0D0
      head_grid%raincv = 0.0D0
      head_grid%pratec = 0.0
   ENDIF

   CALL nl_get_io_form_auxhist17( head_grid%id, io_auxh17 )
   CALL nl_set_io_form_auxhist17( head_grid%id, 0 )
   CALL nl_get_io_form_auxhist18( head_grid%id, io_auxh18 )
   CALL nl_set_io_form_auxhist18( head_grid%id, 0 )
#endif

   CALL nl_get_time_step ( head_grid%id, timestep)  

   IF(head_grid%cycle_interval > 0) then
!BEGIN CYCLING
      CALL       wrf_message ( "wrf: using multiple cycle windows" )
      model_config_rec%input_inname = "wrfinput_d<domain>_<date>"

      IF (head_grid%cycle_interval > 1440) then
         CALL wrf_message("WARNING: large cycle_interval may &
              &exhibit unrealiztic physics.")
      END IF

      IF (.not. head_grid%boundary_io ) xtraj_boundary_form=0
      IF ( model_config_rec%interval_seconds/60 .gt. head_grid%cycle_interval .or. &
           mod(head_grid%cycle_interval, model_config_rec%interval_seconds/60 ) .ne. 0) THEN
         ! Force to turn off boundary input
         xtraj_boundary_form = 0
         CALL wrf_message ("WARNING: cycle_interval must be divisible by &
              &interval_seconds/60 to use boundary updating.")
      ENDIF

      start_year = model_config_rec%start_year(head_grid%id)
      start_month = model_config_rec%start_month(head_grid%id)
      start_day = model_config_rec%start_day(head_grid%id)
      start_hour = model_config_rec%start_hour(head_grid%id)
      start_minute = model_config_rec%start_minute(head_grid%id)

      end_year = model_config_rec%end_year(head_grid%id)
      end_month = model_config_rec%end_month(head_grid%id)
      end_day = model_config_rec%end_day(head_grid%id)
      end_hour = model_config_rec%end_hour(head_grid%id)
      end_minute = model_config_rec%end_minute(head_grid%id)

      ! Store original time keeping variables
      CALL nl_get_run_days  (head_grid%id,run_days)
      CALL nl_get_run_hours (head_grid%id,run_hours)
      CALL nl_get_run_minutes (head_grid%id,run_minutes)
      ncycle = (run_hours + 24 * run_days)

      IF (REAL(ncycle,8) .lt. REAL(head_grid%cycle_interval,8)/60.) THEN
         CALL wrf_error_fatal(&
                 'Need to set total run_hours >= cycle_interval/60')
      ENDIF
      ncycle = ncycle*60 / head_grid%cycle_interval

      end_year0  = end_year
      end_month0 = end_month
      end_day0   = end_day
      end_hour0  = end_hour
      end_minute0  = end_minute

      start_year0 = start_year
      start_month0 = start_month
      start_day0 = start_day
      start_hour0 = start_hour
      start_minute0 = start_minute


      cycle_interval_m = mod(head_grid%cycle_interval,60)
      cycle_interval_h = (head_grid%cycle_interval - cycle_interval_m) / 60
      IF ( cycle_interval_h > 24 ) THEN
         cycle_interval_d = (cycle_interval_h - mod(cycle_interval_h,24)) / 24
         cycle_interval_h = mod(cycle_interval_h,24)
      ELSE
         cycle_interval_d = 0
      ENDIF

      CALL domain_clock_get (head_grid, stop_timestr=STOPtimestr)
      CALL domain_clock_get (head_grid, start_timestr=STARTtimestr)
      cyclestart_timestr=STARTtimestr
   
      model_config_rec%run_days = cycle_interval_d
      model_config_rec%run_hours = cycle_interval_h
      model_config_rec%run_minutes = cycle_interval_m

      DO icycle = 1,ncycle,1
         head_grid%cycle_count = icycle
         write(msg,FMT='(A,I5,A,I5)')'WRF: Simulating cycle ',head_grid%cycle_count,' of ',ncycle
         CALL wrf_debug(1,TRIM(msg))

         !Set clock to cycle initial condition time stamp for I/O
         CALL domain_clock_set (head_grid, current_timestr=cyclestart_timestr )
         CALL domain_clock_set (head_grid, time_step_seconds=60*head_grid%cycle_interval)
         CALL domain_clockadvance (head_grid)
         CALL domain_clock_get (head_grid, current_timestr=cyclestart_timestr)
         CALL domain_clock_set (head_grid, stop_timestr=cyclestart_timestr)
         CALL domain_clock_set (head_grid, time_step_seconds=-60*head_grid%cycle_interval)
         CALL domain_clockadvance (head_grid)
         CALL domain_clock_set (head_grid, time_step_seconds=timestep)
         CALL domain_clockprint(150, head_grid, &
                                'get CurrTime from clock,')

#if (WRF_CHEM == 1)
         ! Store scaleant and scalebb in temporary variables to be retrieved
         ! after reading the current wrfinput file
         allocate ( scaleant_temp(head_grid%sm31:head_grid%em31, &
                              1:model_config_rec%num_ant_steps, &
                              head_grid%sm33:head_grid%em33, 1:num_scaleant) )
         allocate ( scalebb_temp(head_grid%sm31:head_grid%em31, &
                              1:model_config_rec%num_bb_steps, &
                              head_grid%sm33:head_grid%em33, 1:num_scalebb) )

         scaleant_temp = 1.D0
         scalebb_temp = 1.D0

         !$OMP PARALLEL DO    &
         !$OMP PRIVATE ( ij, ic, j, k, i )
         DO ij = 1 , head_grid%num_tiles
            DO ic = PARAM_FIRST_SCALAR, num_scaleant
               DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                  DO k = 1, model_config_rec%num_ant_steps
                     DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                        scaleant_temp(i,k,j,ic) = head_grid%scaleant(i,k,j,ic)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO    &
         !$OMP PRIVATE ( ij, ic, j, k, i )
         DO ij = 1 , head_grid%num_tiles
            DO ic = PARAM_FIRST_SCALAR, num_scalebb
               DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                  DO k = 1, model_config_rec%num_bb_steps
                     DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                        scalebb_temp(i,k,j,ic) = head_grid%scalebb(i,k,j,ic)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO

#endif

         !Read initial conditions
         model_config_rec%restart = .false.
         head_grid%restart = .false.
         call nl_set_restart ( head_grid%id, .false. )
!         model_config_rec%io_form_restart = 0
         head_grid%itimestep = 0

         CALL model_to_grid_config_rec ( head_grid%id , model_config_rec , config_flags )

         CALL med_initialdata_input( head_grid , config_flags )
         head_grid%itimestep = 0
#if (WRF_CHEM == 1)
         head_grid%ktauc = head_grid%itimestep
#endif

         ! Set time keeping variables to cycle values
         model_config_rec%start_year(head_grid%id) = start_year
         model_config_rec%start_month(head_grid%id) = start_month
         model_config_rec%start_day(head_grid%id) = start_day
         model_config_rec%start_hour(head_grid%id) = start_hour
         model_config_rec%start_minute(head_grid%id) = start_minute


         if (start_minute + cycle_interval_m > 59) then
            end_minute = start_minute + cycle_interval_m - 60
            end_hour  = start_hour  + 1
         else
            end_minute = start_minute + cycle_interval_m
            end_hour  = start_hour
         end if

         if (end_hour + cycle_interval_h > 23) then
            end_hour = end_hour + cycle_interval_h - 24
            end_day  = start_day  + 1
         else
            end_hour = end_hour + cycle_interval_h
            end_day  = start_day
         end if

         end_year  = start_year
         end_month = start_month

         !Account for leap years
         if(start_month .eq. 2 .and. end_day > 24) then  !Assuming the cycle interval < 5 days
            if ( mod(end_year,4) .eq. 0 )then
               if ( mod(end_year,100) .eq. 0 ) then
                  if ( mod(end_year,400) .eq. 0 ) then
                     days_in_month(2) = 29
                  else
                     days_in_month(2) = 28
                  end if
               else
                  days_in_month(2) = 29
               end if
            else
               days_in_month(2) = 28
            end if
         end if

         if (end_day + cycle_interval_d > days_in_month(start_month)) then
            if(start_month + 1 > 12) then
               end_month = 1
               end_year  = start_year + 1
            else
               end_month = start_month + 1
            end if
            end_day = end_day + cycle_interval_d - days_in_month(start_month)
         else
            end_day = end_day + cycle_interval_d
         end if

         model_config_rec%end_year(head_grid%id) = end_year
         model_config_rec%end_month(head_grid%id) = end_month
         model_config_rec%end_day(head_grid%id) = end_day
         model_config_rec%end_hour(head_grid%id) = end_hour
         model_config_rec%end_minute(head_grid%id) = end_minute

! Get nonlinear model trajectory
         CALL wrf_message ( "wrf: calculating cycle trajectory" )

         ! Set trajectory boundary form
         CALL nl_set_io_form_boundary( head_grid%id, xtraj_boundary_form )  

         IF ( ASSOCIATED( head_grid%domain_clock ) ) THEN
           IF ( head_grid%domain_clock_created ) THEN
             CALL WRFU_ClockDestroy( head_grid%domain_clock )
             head_grid%domain_clock_created = .FALSE.
           ENDIF
         ENDIF
         IF ( ASSOCIATED( head_grid%alarms ) .AND. &
              ASSOCIATED( head_grid%alarms_created ) ) THEN
           DO alarmid = 1, MAX_WRF_ALARMS
             IF ( head_grid%alarms_created( alarmid ) ) THEN
               CALL WRFU_AlarmDestroy( head_grid%alarms( alarmid ) )
               head_grid%alarms_created( alarmid ) = .FALSE.
             ENDIF
           ENDDO
         ENDIF
         CALL Setup_Timekeeping ( head_grid )
         head_grid%start_subtime = domain_get_start_time ( head_grid )
         head_grid%stop_subtime = domain_get_stop_time ( head_grid )


         !Update time since last boundary read
         head_grid%dtbc = REAL( model_config_rec%interval_seconds, 8)

         ! Enable auxinput streams for cycle_interval < io_intervals
         IF ( icycle .GT. 1 ) THEN
            DO ialarm = first_auxinput, last_auxinput
               CALL WRFU_TimeIntervalGet( head_grid%io_intervals( ialarm ), S=interval, rc=rc )
               IF ( interval .gt. 60*head_grid%cycle_interval ) THEN

                  IF ( mod(start_hour*3600+start_minute*60, interval) .ne. 0 .OR. &
                       mod(interval, 60*head_grid%cycle_interval) .ne. 0 ) THEN
                     CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
                  ENDIF

               ELSEIF ( mod(60*head_grid%cycle_interval, interval) .ne. 0 ) THEN
                     write(msg,FMT='(A,I2)')&
                       'WARNING: when cycle_interval > auxinputXX_interval,&
                       &cycle must be divisible by auxinput.  &
                       &Disabling auxinput',&
                       ialarm-first_input
                     CALL wrf_message (TRIM(msg))
                     CALL WRFU_AlarmDisable( head_grid%alarms( ialarm ), rc=rc )
               ENDIF
            ENDDO
         ENDIF


         CALL domain_clock_get( head_grid, start_timestr=timestr )
         CALL domain_clock_set( head_grid, current_timestr=timestr )

#if (WRF_CHEM == 1)
         ! Retrieve chem from end of previous iteration
         IF ( icycle .GT. 1 ) THEN
            head_grid%chem = 0.D0
            !$OMP PARALLEL DO    &
            !$OMP PRIVATE ( ij, ic, j, k, i )
            DO ij = 1 , head_grid%num_tiles
               DO ic = PARAM_FIRST_SCALAR, num_chem
                  DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                     DO k = head_grid%sp32, head_grid%ep32
                        DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                           head_grid%chem(i,k,j,ic) = chem_temp(i,k,j,ic)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
            deallocate ( chem_temp )
         ENDIF
         ! Retrieve scaleant and scalebb from before wrfinput, which overwrites
         ! them.  scaleant and scalebb are needed in wrfinput and restart files
         ! for the AD/TL checkpointing/cycling.
         head_grid%scaleant = 1.D0
         head_grid%scalebb = 1.D0
         !$OMP PARALLEL DO    &
         !$OMP PRIVATE ( ij, ic, j, k, i )
         DO ij = 1 , head_grid%num_tiles
            DO ic = PARAM_FIRST_SCALAR, num_scaleant
               DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                  DO k = 1, model_config_rec%num_ant_steps
                     DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                        head_grid%scaleant(i,k,j,ic) = scaleant_temp(i,k,j,ic)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            DO ic = PARAM_FIRST_SCALAR, num_scalebb
               DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                  DO k = 1, model_config_rec%num_bb_steps
                     DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                        head_grid%scalebb(i,k,j,ic) = scalebb_temp(i,k,j,ic)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         !$OMP END PARALLEL DO
         deallocate ( scaleant_temp )
         deallocate ( scalebb_temp )

         !Need to rewrite wrfinput_d<domain>_<date> with correct chem, scaleant,
         !and scalebb fields for reading in AD/TL models.
         config_flags%input_outname = "wrfinput_d<domain>_<date>"
         call med_filter_out ( head_grid, config_flags )
         config_flags%input_outname = "wrfinput_d<domain>"

#endif

!         CALL wrf_run
         CALL integrate (head_grid)

#if (WRF_CHEM == 1)
         ! Store chem in a temporary variable to be retrieved in next iteration
         ! after reading the next wrfinput file
         IF ( icycle .LT. ncycle ) THEN
            allocate ( chem_temp(head_grid%sm31:head_grid%em31, &
                                 head_grid%sm32:head_grid%em32, &
                                 head_grid%sm33:head_grid%em33, 1:num_chem) )
            chem_temp = 0.D0
            !$OMP PARALLEL DO    &
            !$OMP PRIVATE ( ij, ic, j, k, i )
            DO ij = 1 , head_grid%num_tiles
               DO ic = PARAM_FIRST_SCALAR, num_chem
                  DO j = head_grid%j_start(ij), head_grid%j_end(ij)
                     DO k = head_grid%sp32, head_grid%ep32
                        DO i = head_grid%i_start(ij), head_grid%i_end(ij)
                           chem_temp(i,k,j,ic) = head_grid%chem(i,k,j,ic)
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
            !$OMP END PARALLEL DO
         ENDIF
#endif

         start_year = end_year
         start_month = end_month
         start_day = end_day
         start_hour = end_hour
         start_minute = end_minute

      ENDDO

      model_config_rec%run_days = run_days
      model_config_rec%run_hours = run_hours
      model_config_rec%run_minutes = run_minutes

      model_config_rec%end_year(head_grid%id) = end_year0
      model_config_rec%end_month(head_grid%id) = end_month0
      model_config_rec%end_day(head_grid%id) = end_day0
      model_config_rec%end_hour(head_grid%id) = end_hour0
      model_config_rec%end_minute(head_grid%id) = end_minute0

      model_config_rec%start_year(head_grid%id) = start_year0
      model_config_rec%start_month(head_grid%id) = start_month0
      model_config_rec%start_day(head_grid%id) = start_day0
      model_config_rec%start_hour(head_grid%id) = start_hour0
      model_config_rec%start_minute(head_grid%id) = start_minute0

      ncycle = 1
      head_grid%cycle_count = 1
      model_config_rec%input_inname = "wrfinput_d<domain>"
!END CYCLING
   ELSE
#if (WRFPLUS == 1)
      head_grid%itimestep = 0
#if (WRF_CHEM == 1)
      head_grid%ktauc = head_grid%itimestep
#endif
      CALL start_domain ( head_grid , .TRUE. )
#endif

      CALL integrate ( head_grid )

   ENDIF

#if (WRFPLUS == 1)
   CALL nl_set_io_form_auxhist17( head_grid%id, io_auxh17 )
   CALL nl_set_io_form_auxhist18( head_grid%id, io_auxh18 )

   ! Close boundary file, as it will be read again from start
   CALL construct_filename2a ( bdyname , config_flags%bdy_inname , head_grid%id, 2 , 'dummydate' )
   CALL wrf_inquire_opened(head_grid%lbc_fid , TRIM(bdyname) , open_status ,ierr )
   IF ( open_status .EQ. WRF_FILE_OPENED_FOR_READ ) THEN
      CALL close_dataset ( head_grid%lbc_fid , config_flags , "DATASET=BOUNDARY")
   ENDIF
#endif

   CALL       wrf_message ( "wrf: back from nonlinear integrate" )

   END SUBROUTINE wrf_run_cycle


   SUBROUTINE wrf_finalize( no_shutdown )
!<DESCRIPTION>
!     WRF finalize routine.
!</DESCRIPTION>

!<DESCRIPTION>
! A Mediation Layer-provided
! subroutine, <a href=med_shutdown_io.html>med_shutdown_io</a> is called
! to allow the the model to do any I/O specific cleanup and shutdown, and
! then the WRF Driver Layer routine <a
! href=wrf_shutdown.html>wrf_shutdown</a> (quilt servers would be
! directed to shut down here) is called to properly end the run,
! including shutting down the communications (for example, most comm
! layers would call MPI_FINALIZE at this point if they're using MPI).
! 
!</DESCRIPTION>
     LOGICAL, OPTIONAL, INTENT(IN) :: no_shutdown

   ! shut down I/O
   CALL med_shutdown_io ( head_grid , config_flags )
   CALL       wrf_debug ( 100 , 'wrf: back from med_shutdown_io' )

   CALL       wrf_debug (   0 , 'wrf: SUCCESS COMPLETE WRF' )

   ! Call wrf_shutdown() (which calls MPI_FINALIZE() 
   ! for DM parallel runs).  
   IF ( .NOT. PRESENT( no_shutdown ) ) THEN
     ! Finalize time manager
      IF (coupler_on) THEN 
         CALL cpl_finalize() 
      ELSE
         CALL WRFU_Finalize
         CALL wrf_shutdown
      ENDIF
   ENDIF

   END SUBROUTINE wrf_finalize


   SUBROUTINE wrf_dfi()
!<DESCRIPTION>
! Runs a digital filter initialization procedure.
!</DESCRIPTION>
      IMPLICIT NONE

! #if (EM_CORE == 1)
      ! Initialization procedure
      IF ( config_flags%dfi_opt .NE. DFI_NODFI ) THEN
   
         SELECT CASE ( config_flags%dfi_opt ) 
     
            CASE (DFI_DFL)
               wrf_err_message = 'Initializing with DFL'
               CALL wrf_message(TRIM(wrf_err_message))
   
               wrf_err_message = '   Filtering forward in time'
               CALL wrf_message(TRIM(wrf_err_message))
   
               CALL wrf_dfi_fwd_init()
               CALL wrf_run()
   
               CALL wrf_dfi_array_reset()
   
               CALL wrf_dfi_fst_init()
   
               IF ( config_flags%dfi_write_filtered_input ) THEN
                  CALL wrf_dfi_write_initialized_state()
               END IF
   
            CASE (DFI_DDFI)
               wrf_err_message = 'Initializing with DDFI'
               CALL wrf_message(TRIM(wrf_err_message))
   
               wrf_err_message = '   Integrating backward in time'
               CALL wrf_message(TRIM(wrf_err_message))
   
               CALL wrf_dfi_bck_init()
               CALL wrf_run()
   
               wrf_err_message = '   Filtering forward in time'
               CALL wrf_message(TRIM(wrf_err_message))
   
               CALL wrf_dfi_fwd_init()
               CALL wrf_run()
   
               CALL wrf_dfi_array_reset()
   
               CALL wrf_dfi_fst_init()
   
               IF ( config_flags%dfi_write_filtered_input ) THEN
                  CALL wrf_dfi_write_initialized_state()
               END IF
   
            CASE (DFI_TDFI)
               wrf_err_message = 'Initializing with TDFI'
               CALL wrf_message(TRIM(wrf_err_message))
   
               wrf_err_message = '   Integrating backward in time'
               CALL wrf_message(TRIM(wrf_err_message))
   
               CALL wrf_dfi_bck_init()
               CALL wrf_run()
   
               CALL wrf_dfi_array_reset()
   
               wrf_err_message = '   Filtering forward in time'
               CALL wrf_message(TRIM(wrf_err_message))
   
               CALL wrf_dfi_fwd_init()
               CALL wrf_run()
   
               CALL wrf_dfi_array_reset()
   
               CALL wrf_dfi_fst_init()
   
               IF ( config_flags%dfi_write_filtered_input ) THEN
                  CALL wrf_dfi_write_initialized_state()
               END IF
   
            CASE DEFAULT
               wrf_err_message = 'Unrecognized DFI_OPT in namelist'
               CALL wrf_error_fatal(TRIM(wrf_err_message))
   
         END SELECT
   
      END IF
! #endif

   END SUBROUTINE wrf_dfi

   SUBROUTINE set_derived_rconfigs
!<DESCRIPTION>
! Some derived rconfig entries need to be set based on the value of other,
! non-derived entries before package-dependent memory allocation takes place.
! This might be employed when, for example, we want to allocate arrays in
! a package that depends on the setting of two or more namelist variables.
! In this subroutine, we do just that.
!</DESCRIPTION>

      IMPLICIT NONE

      INTEGER :: i


! #if (EM_CORE == 1)
      IF ( model_config_rec % dfi_opt .EQ. DFI_NODFI ) THEN
        DO i = 1, model_config_rec % max_dom
           model_config_rec % mp_physics_dfi(i) = -1
        ENDDO
      ELSE
        DO i = 1, model_config_rec % max_dom
           model_config_rec % mp_physics_dfi(i) = model_config_rec % mp_physics(i)
        ENDDO
      END IF
! #endif

#if (EM_CORE == 1)
      IF ( model_config_rec % dfi_opt .EQ. DFI_NODFI ) THEN
        DO i = 1, model_config_rec % max_dom
           model_config_rec % bl_pbl_physics_dfi(i) = -1
        ENDDO
      ELSE
        DO i = 1, model_config_rec % max_dom
           model_config_rec % bl_pbl_physics_dfi(i) = model_config_rec % bl_pbl_physics(i)
        ENDDO
      END IF
#endif

#if (DA_CORE == 1)
      IF ( model_config_rec % dyn_opt .EQ. 2 ) THEN
        DO i = 1, model_config_rec % max_dom
           model_config_rec % mp_physics_4dvar(i) = -1
        ENDDO
      ELSE
        DO i = 1, model_config_rec % max_dom
           model_config_rec % mp_physics_4dvar(i) = model_config_rec % mp_physics(i)
        ENDDO
      END IF
#endif

   END SUBROUTINE set_derived_rconfigs

   RECURSIVE SUBROUTINE alloc_doms_for_dfi ( grid )
   
      !  Input variables.

      TYPE (domain) , pointer :: grid

      !  Local variables.

      TYPE (domain) , pointer :: new_nest_loc
      TYPE (grid_config_rec_type) :: parent_config_flags
      INTEGER :: nestid_loc , kid_loc
   
         !  Are there any subdomains from this level.  The output is the nestid (the domain
         !  ID of the nest), and kid (an index to which of the parent's children this new nested
         !  domain represents).
   
         DO WHILE ( nests_to_open( grid , nestid_loc , kid_loc ) )

            !  If we found another child domain, we continue on: allocate, set up time keeping, 
            !  initialize.
   
            CALL alloc_and_configure_domain ( domain_id  = nestid_loc   , &
                                              grid       = new_nest_loc , &
                                              parent     = grid         , &
                                              kid        = kid_loc        )
         
print *,'for parent domain id #',grid%id,', found child domain #',nestid_loc
            !  Since this is a DFI run, set the DFI switches to the same for all domains.

            new_nest_loc%dfi_opt = head_grid%dfi_opt
            new_nest_loc%dfi_stage = DFI_SETUP
         
            !  Set up time keeping for the fine grid space that was just allocated.

            CALL Setup_Timekeeping (new_nest_loc)

            !  With space allocated, and timers set, the fine grid can be initialized with data.

            CALL model_to_grid_config_rec ( grid%id , model_config_rec , parent_config_flags )
            CALL med_nest_initial ( grid , new_nest_loc , config_flags )

            !  Here's the recursive part.  For each of these child domains, we call this same routine.
            !  This will find all of "new_nest_loc" first generation progeny.
   
            CALL alloc_doms_for_dfi ( new_nest_loc )
   
         END DO
   
   END SUBROUTINE alloc_doms_for_dfi

END MODULE module_wrf_top
