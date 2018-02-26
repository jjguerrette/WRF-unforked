  subroutine da_solve ( grid , config_flags)

   !-----------------------------------------------------------------------
   ! Purpose: TBD
   ! 
   ! Edited 09/06/2012: Allow for variable ntmax for each outer loop (Mike Kavulich)
   !-----------------------------------------------------------------------

   implicit none



   type (domain),               intent(inout) :: grid
   type (grid_config_rec_type), intent(inout) :: config_flags

   type (xbx_type)              :: xbx         ! For header & non-grid arrays.
   type (be_type)               :: be          ! Background error structure.
   real*8, allocatable          :: cvt(:)      ! Control variable structure.
   real*8, allocatable          :: xhat(:)     ! Control variable structure.
   real*8, allocatable          :: qhat(:,:)   ! Control variable structure.
   real*8, allocatable          :: eignvec(:,:)
   real*8, allocatable          :: eignval(:)   
!   real, allocatable            :: full_eignvec(:)   
   type (y_type)                :: ob          ! Observation structure.
   type (iv_type)               :: iv          ! Obs. increment structure.
   type (y_type)                :: re          ! Residual (o-a) structure.
   type (y_type)                :: y           ! y = H(x_inc) structure.
   integer                      :: it          ! External loop counter.
   integer                      :: global_it   ! External loop counter.
   integer                      :: neign       
   type (j_type)                :: j           ! Cost function.
   type (j_type)                :: j0_cost     ! Cost function.
   type (y_type)                :: jo_grad_y   ! Grad_y(jo)

   integer                      :: cv_size, i, ichan, k, n
   real                         :: j_grad_norm_target ! Target j norm.

   character (len=3)            :: ci
   character (len=2)            :: outerloop
   character (len=256)          :: timestr
   integer   :: min_yyyy,min_mm,min_dd,min_hh,min_mn,min_ss
   integer   :: max_yyyy,max_mm,max_dd,max_hh,max_mn,max_ss
   character :: s
   real*8    :: time_min, time_max 
   integer   :: jl_start, jl_end
   character(len=256) :: timestr1
   type(x_type) :: shuffle

   real, allocatable :: grid_box_area(:,:), mapfac(:,:)
   character (len=10)    :: variable_name

   integer :: xx, yy, ii, jj,  xy, kk, jjj
   real    :: lat, lon, x_pos, y_pos, dx, dy, dxm, dym
   logical :: outside

   integer   :: cvt_unit, iost
   character(len=8) :: cvtfile
   logical :: ex
#if defined(LAPACK)
   character*10                 :: cproc,cens
   character*10                 :: cit
   character(len=filename_len)  :: cvtfile_rand
   character(len=filename_len)  :: filename
   integer                      :: run_hours0
   integer                      :: iens, nsamp
   logical                      :: rand_exit = .false.
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
!   real*8, allocatable          :: cvt_pert(:)      ! Control variable structure.

   integer                      :: je_start, je_end       ! Start/end indices of Je.
   integer                      :: jp_start, jp_end       ! Start/end indices of Jp.

   integer              :: seed_size
#endif

#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
   integer                           :: mz(13), info, nsstwrk
#else
   integer                           :: mz(7), info, nsstwrk
#endif
#else
#ifdef CLOUD_CV
   integer                           :: mz(15), info, nsstwrk
#else
   integer                           :: mz(9), info, nsstwrk
#endif
#endif

   real :: j0_total
#if (WRF_CHEM == 1)
   real*8, allocatable          :: xhat_prev(:)     ! Control variable structure.
   integer                      :: neign_trunc

   character(len=2) :: itc
   integer          :: cv_s, cv_e, cv_s_ant, cv_e_ant, cv_s_bb, cv_e_bb
#endif
   type(hessian_type) :: hessian_info

   if (trace_use) call da_trace_entry("da_solve")

#ifdef DM_PARALLEL
   call mpi_barrier(comm,ierr)
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif

   if ( config_flags%use_baseparam_fr_nml ) then
      call nl_get_base_pres  ( 1 , base_pres )
      call nl_get_base_temp  ( 1 , base_temp )
      call nl_get_base_lapse ( 1 , base_lapse )
      call nl_get_iso_temp   ( 1 , iso_temp )
      if ( iso_temp .NE. grid%tiso ) THEN
         write(unit=message(1),fmt='(A)') &
           'Namelist iso_temp does not equal iso_temp from fg. Reset namelist value and rerun.'
         write(unit=message(2),fmt='(A,F10.5)')'Namelist iso_temp   = ',iso_temp
         write(unit=message(3),fmt='(A,F10.5)')'Background iso_temp = ',grid%tiso
         call da_error(__FILE__,__LINE__,message(1:3))
      end if
      call nl_get_base_pres_strat  ( 1, base_pres_strat )
      call nl_get_base_lapse_strat ( 1, base_lapse_strat )

      grid%p00   = base_pres
      grid%t00   = base_temp
      grid%tlp   = base_lapse
      grid%tiso  = iso_temp
      grid%p_strat   = base_pres_strat
      grid%tlp_strat = base_lapse_strat
   else
      base_pres  = grid%p00
      base_temp  = grid%t00
      base_lapse = grid%tlp
      iso_temp   = grid%tiso
      base_pres_strat  = grid%p_strat
      base_lapse_strat = grid%tlp_strat
      if ( base_temp < 100.0 .or. base_pres < 10000.0 ) then
         write(unit=message(1),fmt='(A)') &
         'did not find base state parameters in fg. Add use_baseparam_fr_nml = .t. in &dynamics and rerun '
         call da_error(__FILE__,__LINE__,message(1:1))
      end if
   end if

#if (WRF_CHEM == 1)
   !Add all non-chem obs as they become available
   if (use_nonchemobs) then
      use_nonchemobs =  any((/use_synopobs, use_shipsobs, use_metarobs, use_soundobs, &
              use_mtgirsobs, use_tamdarobs, use_bogusobs, use_pilotobs, &
              use_airepobs, use_geoamvobs, use_polaramvobs, use_buoyobs, &
              use_profilerobs, use_satemobs, use_gpspwobs, use_gpsztdobs, &
              use_gpsrefobs, use_ssmiretrievalobs, use_ssmitbobs, &
              use_ssmt1obs, use_ssmt2obs, use_qscatobs, use_hirs2obs, &
              use_hirs3obs, use_hirs4obs, use_mhsobs, use_msuobs, &
              use_amsuaobs, use_amsubobs, use_mwtsobs, use_mwhsobs, &
              use_atmsobs, use_airsobs, use_eos_amsuaobs, use_hsbobs, &
              use_airsretobs, use_rainobs, &
              use_iasiobs, use_radarobs, use_radar_rv, use_radar_rf/))
   end if


   !Add all chem obs as they become available
   use_chemobs = any((/use_chem_surfobs,use_chem_acftobs/))

   !Replace false below with applicable chem obs as they become available (e.g., AOD?)
   use_offlineobs = any((/(use_nonchemobs.and.cv_options.gt.0),&
                          (.false..and.cv_options_chem.gt.0)/)) 

   !Replace false below with applicable non-chem obs as they become available
   use_onlineobs = any((/(.false..and.cv_options.gt.0),&
                          (any((/use_chem_surfobs,use_chem_acftobs/)).and.cv_options_chem.gt.0)/)) 

   init_chem_scale = &
      ( cv_options_chem.eq.1 .and. chem_opt(1).gt.0 .and. .not.osse_chem &
#if defined(LAPACK)
         .and. (rand_outer.eq.1 .or. .not.use_randomblock) &
#endif
         .and. (iouter_restart.eq.1 .or. .not.outer_loop_restart) )

   if ( .not.any((/use_nonchemobs,use_chemobs/)) ) then
      CALL wrf_error_fatal('Must use either meteorological or chemical observations to conduct DA')
   end if
#endif

   ! Calculate the num_fgat_time based on time_window_min, time_window_max
    if ( var4d ) then
       if (time_step == 0) then
          write(unit=message(1),fmt='(A)') &
          'For 4DVAR, in the &domains namelist, "time_step" must be set to a non-zero value'
          call da_error(__FILE__,__LINE__,message(1:1))
       endif
       read(unit=time_window_min,fmt='(i4,5(a1,i2))') min_yyyy,s,min_mm,s,min_dd,s,min_hh,s,min_mn,s,min_ss
       read(unit=time_window_max,fmt='(i4,5(a1,i2))') max_yyyy,s,max_mm,s,max_dd,s,max_hh,s,max_mn,s,max_ss
       call da_get_julian_time(min_yyyy,min_mm,min_dd,min_hh,min_mn,time_min)
       call da_get_julian_time(max_yyyy,max_mm,max_dd,max_hh,max_mn,time_max)
       if ( var4d_bin < time_step ) call nl_set_var4d_bin (1, time_step)
       time_max = (time_max - time_min) * 60   ! unit is : seconds
!       num_fgat_time = NINT(time_max/var4d_bin)
!       if ( NINT(time_max/var4d_bin)*var4d_bin .ne. NINT(time_max) ) then
#if (WRF_CHEM == 1)
       if (.not. use_nonchemobs) then
          num_fgat_time = -1
       else if ( var4d_bin .eq. time_step ) then
          num_fgat_time = 0
       else
#endif
          num_fgat_time = NINT(time_max/var4d_bin)
#if (WRF_CHEM == 1)
       end if
#endif
       if ( NINT(time_max/var4d_bin)*var4d_bin .ne. NINT(time_max) &
#if (WRF_CHEM == 1)
            .and. use_nonchemobs &
#endif
          ) then
          write(unit=message(1),fmt='(A)') &
          '4DVAR assimilation window must be evenly divisible by var4d_bin!'
          write(unit=message(2),fmt='(A,I7)') &
          'var4d_bin       = ',var4d_bin
          write(unit=message(3),fmt='(A,A)') &
          'time_window_max = ',time_window_max
          write(unit=message(4),fmt='(A,A)') &
          'time_window_min = ',time_window_min
          write(unit=message(5),fmt='(A,F10.4)') &
          'time_window_max - time_window_min = ',time_max
          write(unit=message(6),fmt='(A)')'Change var4d_bin, time_window_max, or time_window_min in namelist and rerun'
          call da_error(__FILE__,__LINE__,message(1:6))
       endif
       if ( var4d_bin/time_step*time_step .ne. var4d_bin  &
#if (WRF_CHEM == 1)
            .and. use_nonchemobs &
#endif
          ) then
          write(unit=message(1),fmt='(A)') &
          'var4d_bin must be evenly divisible by time_step!'
          write(unit=message(2),fmt='(A,I7)') &
          'var4d_bin = ',var4d_bin
          write(unit=message(3),fmt='(A,I7)') &
          'time_step = ',time_step
          write(unit=message(4),fmt='(A)')'Change var4d_bin or time_step in namelist and rerun'
          call da_error(__FILE__,__LINE__,message(1:4))
       endif

       num_fgat_time = num_fgat_time + 1
       write(unit=message(1),fmt='(a,i10)') 'num_fgat_time is: ', num_fgat_time
       call da_message(message(1:1))
       if ( use_rainobs ) then
          allocate (fgat_rain_flags(1:num_fgat_time))
          fgat_rain_flags = .false.
          if ( INT(var4d_bin_rain/var4d_bin)*var4d_bin .ne. INT(var4d_bin_rain) ) then
             write(unit=message(1),fmt='(A,A,2I7)') &
             'Please change var4d_bin_rain in namelist and rerun==>', 'var4d_bin_rain, var4d_bin:',var4d_bin_rain,var4d_bin
             call da_error(__FILE__,__LINE__,message(1:1))
          endif
          do n = 1, num_fgat_time, INT(var4d_bin_rain/var4d_bin)
             fgat_rain_flags(n) = .true.
          end do
       end if
    endif

   !---------------------------------------------------------------------------
   ! [1.0] Initial checks
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 1.0 ***")


   if (cv_options_hum /= cv_options_hum_specific_humidity .and. &
       cv_options_hum /= cv_options_hum_relative_humidity) then
      write(unit=message(1),fmt='(A,I3)') &
         'Invalid cv_options_hum = ', cv_options_hum
      call da_error(__FILE__,__LINE__,message(1:1))
   end if

   if (vert_corr == vert_corr_2) then
      if (vertical_ip < vertical_ip_0 .or. vertical_ip > vertical_ip_delta_p) then
         write (unit=message(1),fmt='(A,I3)') &
           'Invalid vertical_ip = ', vertical_ip
         call da_error(__FILE__,__LINE__,message(1:1))
      end if
   end if

   if( use_rf )then
      if (0.5 * real(rf_passes) /= real(rf_passes / 2)) then
         write(unit=stdout,fmt='(A,I4,A)') &
            'rf_passes = ', rf_passes, ' .Should be even.'
         rf_passes = int(real(rf_passes / 2))
         write(unit=stdout,fmt='(A,I4)') 'Resetting rf_passes = ', rf_passes
      end if
   else
      write(stdout,'("da_solve: using wavelet transform")')
   endif

   if ( anal_type_hybrid_dual_res .and. alphacv_method .ne. alphacv_method_xa ) then
      write (unit=message(1),fmt='(A)') &
        'Dual-res hybrid only with alphacv_method = 2'
      call da_error(__FILE__,__LINE__,message(1:1))
   endif

   if (anal_type_randomcv) then
      ntmax = 0
      write(unit=stdout,fmt='(a)') &
         ' Resetting ntmax = 0 for analysis_type = randomcv' 
   end if

#ifdef CLOUD_CV
   if ( cv_options == 3 ) then
      write(unit=message(1),fmt='(A,I3)') &
         'cloud_cv does NOT work with cv_options = ', cv_options
      call da_error(__FILE__,__LINE__,message(1:1))
   end if
#endif

   !---------------------------------------------------------------------------
   ! [2.0] Initialise wrfvar parameters:
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 2.0 ***")

   if ( anal_type_hybrid_dual_res ) then

      !---------------------------------
      ! Get full ensemble grid dimensions
      !---------------------------------
      call da_solve_init(ensemble_grid &
#include "actual_new_args.inc"
)
      ide_ens = ide ! these are unstaggered dimensions of the full ensemble domain
      jde_ens = jde
      kde_ens = kde

      !---------------------------------------
      ! Get "intermediate" grid sizes and tiles
      !---------------------------------------
      call da_solve_init(grid%intermediate_grid &
#include "actual_new_args.inc"
)

      ! these are unstaggered dimensions of the "intermediate" ensemble domain
      !  The intermediate grid is the coarse (ensemble) domain that is co-located with the
      !  hi-resolution (analysis) grid

      ids_int = ids ; jds_int = jds ; kds_int = kds
      ide_int = ide ; jde_int = jde ; kde_int = kde
      
      its_int = its ; ite_int = ite
      jts_int = jts ; jte_int = jte
      kts_int = kts ; kte_int = kte

      ims_int = ims ; ime_int = ime
      jms_int = jms ; jme_int = jme
      kms_int = kms ; kme_int = kme

      ips_int = ips ; ipe_int = ipe
      jps_int = jps ; jpe_int = jpe
      kps_int = kps ; kpe_int = kpe


      grid%imask_xstag = 1   ; grid%imask_ystag = 1
      grid%imask_nostag = 1  ; grid%imask_xystag = 1

      !---------------------------------------------------------------------------
      ! De-allocate parts of grid and replace with grid%intermediate_grid dimensions
      !---------------------------------------------------------------------------
      call reallocate_analysis_grid(grid)

      !----------------------------------------------------------
      ! Allocate and initialize some of grid%intermediate_grid
      !----------------------------------------------------------
      call allocate_intermediate_grid(grid%intermediate_grid)


      !---------------------------------------
      ! Get map projection information for the ensemble
      !---------------------------------------

       call da_setup_firstguess(xbx, ensemble_grid, config_flags, .true. )
       map_info_ens = map_info ! map_info is read in from da_tools.f90, call it something else

   endif

   call da_solve_init(grid &
#include "actual_new_args.inc"
)

   if ( .not. anal_type_hybrid_dual_res ) then
      ide_ens = ide ; jde_ens = jde ; kde_ens = kde

      ids_int = ids ; ide_int = ide 
      jds_int = jds ; jde_int = jde 
      kds_int = kds ; kde_int = kde

      its_int = its ; ite_int = ite 
      jts_int = jts ; jte_int = jte
      kts_int = kts ; kte_int = kte

      ims_int = ims ; ime_int = ime
      jms_int = jms ; jme_int = jme
      kms_int = kms ; kme_int = kme

      ips_int = ips ; ipe_int = ipe
      jps_int = jps ; jpe_int = jpe
      kps_int = kps ; kpe_int = kpe
   endif

   !---------------------------------------------------------------------------
   ! [3.0] Set up first guess field (grid%xb):
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 3.0 ***")

   call da_setup_firstguess(xbx, grid, config_flags, .false.)

   if ( anal_type_hybrid_dual_res ) then

      ! 
      ! Get ensemble grid mapfactor on entire coarse grid
      ! 
      variable_name = 'MAPFAC_M'

      allocate( grid_box_area(1:ide_ens,1:jde_ens), mapfac(1:ide_ens,1:jde_ens) )  
      call da_get_var_2d_real_cdf( input_file_ens, variable_name, mapfac, ide_ens, jde_ens, 1, .false. )
      grid_box_area(:,:) = ( (ensemble_grid%dx)/mapfac(:,:) )**2 
      grid%intermediate_grid%xb%grid_box_area(its_int:ite_int,jts_int:jte_int) = grid_box_area(its_int:ite_int,jts_int:jte_int)
      deallocate(mapfac,grid_box_area)

      !
      ! Make a list of "observations" from the the fine grid lat/lon
      !

      xy = ( ime - ims + 1 ) * ( jme - jms + 1 )

      allocate(ob_locs(1:xy)) ! From da_control

      kk = 0 ; jjj = 0
      do xx = ims, ime !ids, ide 
         do yy = jms, jme !jds, jde 

          outside = .false.

          lat = grid%xb%lat(xx,yy)
          lon = grid%xb%lon(xx,yy)

          x_pos = -1.0 ; y_pos = -1.0
          call da_llxy_wrf(map_info_ens,lat,lon,x_pos,y_pos)
          call da_togrid(x_pos,its_int-2, ite_int+2, ii, dx, dxm )
          call da_togrid(y_pos,jts_int-2, jte_int+2, jj, dy, dym )

          if ((int(x_pos) < ids_int) .or. (int(x_pos) >= ide_int) .or. &
             (int(y_pos) < jds_int) .or. (int(y_pos) >= jde_int)) then
             outside = .true.
          endif

	  if ((ii < ids_int) .or. (ii >= ide_int) .or. &
	     (jj < jds_int) .or. (jj >= jde_int)) then
	     outside     = .true.
	  endif

         if ((ii < its_int-1) .or. (ii > ite_int) .or. &
            (jj < jts_int-1) .or. (jj > jte_int)) then
            outside = .true.
         endif

         if ( .not. outside ) then
              kk = kk + 1
              ob_locs(kk)%x = x_pos
              ob_locs(kk)%y = y_pos
              ob_locs(kk)%i = ii
              ob_locs(kk)%j = jj
              ob_locs(kk)%dx = dx
              ob_locs(kk)%dy = dy
              ob_locs(kk)%dxm = dxm
              ob_locs(kk)%dym = dym
              ob_locs(kk)%xx = xx
              ob_locs(kk)%yy = yy
          else
             jjj = jjj + 1
          endif

        enddo
     enddo

      total_here = kk ! from da_control

   endif 

   !---------------------------------------------------------------------------
   ! [4.0] Set up observations (ob):
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 4.0 ***")

   call da_setup_obs_structures (grid, ob, iv, j)
   if (use_rad) then
      allocate (j % jo % rad(1:iv%num_inst))
      do i=1,iv%num_inst
         allocate (j % jo % rad(i) % jo_ichan(iv%instid(i)%nchan))
         allocate (j % jo % rad(i) % num_ichan(iv%instid(i)%nchan))
      end do
   end if

   !---------------------------------------------------------------------------
   ! [4.1] Observer (ANAL_TYPE="VERIFY")
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 4.1 ***")

   if (anal_type_verify) then
      check_max_iv = .false.
      ntmax=0
      it = 1
      num_qcstat_conv=0

#if defined(RTTOV) || defined(CRTM)
      if (use_rad .and. (use_varbc.or.freeze_varbc)) call da_varbc_init(iv, be)
#endif

      call da_get_innov_vector (it, num_qcstat_conv, ob, iv, grid , config_flags)
      call da_allocate_y (iv, re)
 
      ! write out O-B statistics
      call da_write_diagnostics(it, grid,num_qcstat_conv, ob, iv, re, y, j)

      ! write out Gradient of Jo for adjoint sensitivity
      if (adj_sens) then
         cv_size = 1
         allocate (xhat(cv_size))
         call da_allocate_y (iv, y)
         call da_allocate_y (iv, jo_grad_y)

         call da_calculate_residual(iv, y, re)
         call da_calculate_grady(iv, re, jo_grad_y)
         call da_zero_x(grid%xa)

         call da_transform_xtoy_adj(cv_size, xhat, grid, iv, jo_grad_y, grid%xa)
         call da_transform_xtoxa_adj(grid)
         call da_transfer_wrftltoxa_adj(grid, config_flags, 'fcst', timestr)

         call da_deallocate_y (y)
         call da_deallocate_y (jo_grad_y)
      end if

      call da_deallocate_y(re)
      call da_deallocate_observations (iv)
      if (trace_use) call da_trace_exit ("da_solve")
      return
   end if
   
   !---------------------------------------------------------------------------
   ! [5.0] Set up control variable:
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 5.0 ***")

   be % cv % size_jb = 0
   be % cv % size_je = 0
   be % cv % size_jp = 0
   be % cv % size_js = 0
   be % cv % size_jl = 0
   
   !---------------------------------------------------------------------------
   ! [5.1] Set up background errors (be):
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 5.1 ***")

   if (use_background_errors .and. multi_inc /= 1) then
      call da_setup_background_errors (grid, be)
   else
      be % ne = ensdim_alpha
      be % v1 % mz = 0
      be % v2 % mz = 0
      be % v3 % mz = 0
      be % v4 % mz = 0
      be % v5 % mz = 0
#ifdef CLOUD_CV
      be % v6  % mz = 0
      be % v7  % mz = 0
      be % v8  % mz = 0
      be % v9  % mz = 0
      be % v10 % mz = 0
      be % v11 % mz = 0
#endif
#if (WRF_CHEM == 1)
!      allocate (be%alpha_chem (num_chem) )
!      be % alpha_chem(:) % mz = 0

      allocate (be%alpha_ant_chem (num_scaleant) )
      be % alpha_ant_chem(:) % mz = 0
      allocate (be%alpha_bb_chem (num_scalebb) )
      be % alpha_bb_chem(:) % mz = 0
#endif
   end if

    ! overwrite variables defined in da_setup_cv.inc set in the call to da_setup_background_errors
    if ( anal_type_hybrid_dual_res ) then
       be % cv % size_alphac = (ite_int - its_int + 1) * (jte_int - jts_int + 1) * be % alpha % mz * be % ne
       be % cv % size_je = be % cv % size_alphac
       cv_size_domain_je = (ide_int - ids_int + 1) * (jde_int - jds_int + 1) * be % alpha % mz * be % ne
    endif
   
   !---------------------------------------------------------------------------
   ! [5.2] Set up observation bias correction (VarBC):
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 5.2 ***")
#if defined(RTTOV) || defined(CRTM)
   if (use_rad .and. (use_varbc.or.freeze_varbc)) call da_varbc_init(iv, be)
#endif

   !---------------------------------------------------------------------------
   ! [5.3] Set up satellite control variable:
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 5.3 ***")

#if defined(RTTOV) || defined(CRTM)
   if (ANY(use_satcv)) call da_setup_satcv(iv, be)
#endif
   
   !---------------------------------------------------------------------------
   ! [5.4] Total control variable:
   !---------------------------------------------------------------------------   
   call wrf_message("*** WRF-Var Part 5.4 ***")

   be % cv % size = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_js + be%cv%size_jl
   cv_size = be % cv % size

   !---------------------------------------------------------------------------
   ! [6.0] Set up ensemble perturbation input:
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 6.0 ***")

      grid % ep % ne = be % ne
      if (use_background_errors .and. be % ne > 0) then
!        call da_setup_flow_predictors ( ide, jde, kde, be % ne, grid%ep, &
!                                        its, ite, jts, jte, kts, kte )
         call da_setup_flow_predictors ( ide_ens, jde_ens, kde_ens, be % ne, grid%ep, &
                                         its_int, ite_int, jts_int, jte_int, kts_int, kte_int )
      end if

   !---------------------------------------------------------------------------
   ! [7.0] Setup control variable (cv):
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 7.0 ***")

!  Dynamically allocate the variables which don't rely on ntmax
      allocate (cvt(1:cv_size))
      allocate (xhat(1:cv_size))
!      if (use_lanczos) then
!        allocate (full_eignvec(cv_size))
!      end if

#if (WRF_CHEM != 1)
#ifdef CLOUD_CV
      mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
              be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne /)
#else
      mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne /)
#endif
#else
#ifdef CLOUD_CV
      mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, &
              be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz, be%v10%mz, be%v11%mz, be%alpha%mz, be%ne, &
              be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#else
      mz = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be%ne, &
              be%alpha_ant_chem(PARAM_FIRST_SCALAR)%mz, be%alpha_bb_chem(PARAM_FIRST_SCALAR)%mz /)
#endif
#endif

#if defined(LAPACK)
      if (use_randomblock .and. rand_minimise .and. rand_outer.gt.1) then
         write(cproc,fmt='(i4.4)') myproc
         write(cit,fmt='(i2.2)') rand_outer-1
         cvtfile_rand = 'cvt.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('r', cvtfile_rand, cv_size, cvt(1:cv_size), grid, mz)
         if ( max_ext_its > 1 ) then
           max_ext_its=1
           write(unit=message(1),fmt='(a)') "Re-set max_ext_its = 1 for use_randomblock"
           call da_message(message(1:1))
         end if
      else
#endif
         if ( outer_loop_restart ) then
            global_it = iouter_restart

            call da_get_unit(cvt_unit)
            !cvt_unit=600
            if ( max_ext_its > 1 ) then
              max_ext_its=1
              write(unit=message(1),fmt='(a)') "Re-set max_ext_its = 1 for outer_loop_restart"
              call da_message(message(1:1))
            end if
!            if (use_global_cv_io) then
!               write(cproc,fmt='(i4.4)') myproc
!               cvtfile = 'cvt.p'//trim(adjustl(cproc))
!               if (rootproc) inquire(file=trim(cvtfile), exist=ex)
!!               call da_cv_io ('w', cvtfile, cv_size, omega(1:cv_size,iens), grid, mz)
!            else
!               write(unit=cvtfile,fmt='(a,i4.4)') 'cvt_',myproc
!               inquire(file=trim(cvtfile), exist=ex)
!            end if

            write(unit=cvtfile,fmt='(a,i4.4)') 'cvt_',myproc
            inquire(file=trim(cvtfile), exist=ex)
            if ( ex ) then
              open(unit=cvt_unit,file=trim(cvtfile),iostat=iost,form='UNFORMATTED',status='OLD')
              if (iost /= 0) then
                write(unit=message(1),fmt='(A,I5,A)') &
                   "Error ",iost," opening cvt file "//trim(cvtfile)
                call da_error(__FILE__,__LINE__,message(1:1))
              end if
              write(unit=message(1),fmt='(a)') 'Reading cvt from : '//trim(cvtfile) 
              call da_message(message(1:1))
              read(cvt_unit) cvt
              close(cvt_unit)

!!USED FOR DEBUGGING ONLY
!      if (var4d_inc_out) then
!         write(cproc,fmt='(i4.4)') myproc
!         write(cit,fmt='(i2.2)') global_it-1
!         cvtfile_rand = 'cvt_REF_out.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
!         call da_cv_io ('w', cvtfile_rand, cv_size, cvt(1:cv_size), grid, mz)
!
!         call wrf_message("*** WRF-Var initialized cvt ***")
!         if (trace_use) call da_trace_exit("da_wrfvar_main")
!         if (trace_use) call da_trace_report
!
!         call wrfu_finalize
!         call wrf_shutdown
!         stop
!      endif
!!USED FOR DEBUGGING ONLY

            else
              write(unit=message(1),fmt='(a)') "cvt file '"//trim(cvtfile)//"' does not exists, initiallizing cvt."
              call da_message(message(1:1))
              call da_initialize_cv (cv_size, cvt)

              if (iouter_restart .gt. 1) then
                 write(unit=message(1),fmt='(a)') "Re-set outer_loop_restart=.false. and iouter_restart=1."
                 call da_message(message(1:1))

                 outer_loop_restart = .false.
                 iouter_restart = 1
                 global_it = 1
              end if
            end if
         else
            call da_initialize_cv (cv_size, cvt)
         end if
#if defined(LAPACK)
      end if
#endif

      call da_zero_vp_type (grid%vv)
      call da_zero_vp_type (grid%vp)
  
      if ( var4d ) then
         call da_zero_vp_type (grid%vv6)
         call da_zero_vp_type (grid%vp6)
#if (WRF_CHEM == 1)
         call da_zero_xch_type (grid%vch)
#endif
      end if

#ifdef VAR4D
#if (WRF_CHEM == 1)
   !---------------------------------------------------------------------------
   ! [7.1] Setup CHEM OSSE w/ existing observation locations:
   !---------------------------------------------------------------------------
   if (init_osse_chem) then
      call wrf_message("*** WRF-Var Part 7.1 ***")

      !Generate perturbed first guess and perturbed "true" observations
#if defined(LAPACK)
      if (.not.use_randomblock .or. (rand_stage.eq.0 .and. rand_outer.eq.1) ) then
#endif
      !! Initialize model grid
      disable_traj = .true.
      call da_nl_model(1)

      ! Write OSSE background state vector from model_grid
      call da_init_model_output(filename = "fg_osse")

      call wrf_message("Generating global omega")
      allocate(omega(cv_size,1))

      call da_gen_omega(omega, cv_size, 1, grid, mz)

      !Write omega to file
      use_global_cv_io = .true.
      write(cproc,fmt='(i4.4)') myproc
      filename = 'cvt_true.p'//trim(adjustl(cproc))
      call da_cv_io ('w', filename, cv_size, omega(1:cv_size,iens), grid, mz)

      call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                               omega(1:be%cv%size_jb,1),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                               , grid%vch &
#endif
                               )

      !Transform initial state update to T0 boundary condition update
      call da_transform_xatowrftl_lbc_wrapper( &
              global_it, cv_size, be, grid%ep, &
              omega(:,1), iv, grid%vp6, grid%vv6, xbx, &
              grid, config_flags &
!#if (WRF_CHEM == 1)
!                      , grid%vch6 &
!#endif
              )

      call da_transform_xtoxa (grid)

      call da_transfer_xatoanalysis (it, xbx, grid, config_flags)

      call da_transfer_headtomodel (grid)

      ! Write OSSE truth state vector from model_grid
      call da_init_model_output(filename = "fg_truth")

!Not used in init_osse_chem
!      cvt = cvt + omega(:,1)

      deallocate(omega)
      disable_traj = .false.

      write_checkpoints = 0

      if ( checkpoint_interval .gt. 0) then
         !Turn off all trajectory writing
         write_checkpoints = 3
      end if

      call da_nl_model(0)

      !Restore default trajectory setting
      write_checkpoints = 0

      call da_retrieve_chem_hx( model_grid, iv )
      call da_write_obs_chem( iv, ob, 99 )

!Not necessary with above fg_pert and fg_truth writing
!      ! Re-read model inputs from fg file, then overwrite scaling factors.
!      call da_init_model_input

      ! Write out pseudo obs !!!and perturb head_grid scaling factors
      call da_setup_osse_chem (grid, ob, iv, be, xbx)

!!     !Update xbch
!!      call da_transfer_wrftoxb_chem(grid) !Not necessary if not performing inversion in this execution
!
!      do i = PARAM_FIRST_SCALAR, num_scaleant
!         call kj_swap (grid%scaleant(:,:,:,i), model_grid%scaleant(:,:,:,i), &
!                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, 1, num_ant_steps)
!      enddo
!      do i = PARAM_FIRST_SCALAR, num_scalebb
!         call kj_swap (grid%scalebb(:,:,:,i), model_grid%scalebb(:,:,:,i), &
!                       grid%xp%ims, grid%xp%ime, grid%xp%jms, grid%xp%jme, 1, num_bb_steps)
!      enddo
!
!      ! Write perturbed state vector from model_grid
!      call da_init_model_output(filename = "fg_pert")

      !Clean up and close WRFDA
      call da_deallocate_observations (iv)
      call da_deallocate_y (ob)
      if (use_background_errors) call da_deallocate_background_errors (be)

      call da_wrfvar_finalize
#ifdef VAR4D
      if (var4d) call clean_4dvar
      if (var4d) call da_finalize_model
#endif
      call wrf_message("*** WRF-Var initialized OSSE successfully ***")
      if (trace_use) call da_trace_exit("da_wrfvar_main")
      if (trace_use) call da_trace_report

      call wrfu_finalize
      call wrf_shutdown
      stop
#if defined(LAPACK)
      end if
#endif
   end if !end init_osse_chem
#endif
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: SOLVEINIT')
#endif

   !---------------------------------------------------------------------------
   ! [8] Outerloop and chem Sensitivity
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 8.0 ***")

#if (WRF_CHEM == 1)
   if (.not. chem_sens) then
#endif

   j_grad_norm_target = 1.0
   do it = 1, max_ext_its

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif

!  Dynamically allocate the variables which depend on ntmax
      if (.not.outer_loop_restart) global_it = it
#if defined(LAPACK)
#ifdef VAR4D
      if (use_randomblock .and. var4d) then
         global_it = rand_outer
         use_lanczos = .false.
      end if
#endif
#endif
      nmodes_global=ntmax(global_it)      
      if (use_lanczos) then
         allocate (qhat(1:cv_size, 0:nmodes_global))
         allocate (eignvec(nmodes_global, nmodes_global))
         allocate (eignval(nmodes_global))
#if defined(LAPACK)
         rand_type = 0
         rand_stage = 0
         rand_inner_it = 0
#endif
      end if
#ifdef VAR4D
#if defined(LAPACK)
      if (use_randomblock .and. var4d) then
         nsamp = nmodes_global
         write(randvecdir,fmt='(A,i2.2,A)')'../vectors_', global_it,'/'

         if (rand_type .eq. 3) nmodes_global=nsamp*rand_inner_it
         if ( rand_stage.eq.0 .and. read_omega ) then
            call wrf_message("Generating global omega")
            allocate(omega(cv_size,nsamp))

!            put_rand_seed = .true.
!            seed_array1 = 2343*nsamp
!            seed_array2 = 5677*global_it

            call da_gen_omega(omega, cv_size, nsamp, grid, mz)

            !Write omega to file
            write(cproc,fmt='(i4.4)') myproc
            do iens = 1, nsamp
               write(cens,fmt='(i4.4)') iens

               filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
               call da_cv_io ('w', filename, cv_size, omega(1:cv_size,iens), grid, mz)
            end do
            deallocate(omega)

            if (checkpoint_interval.le.0 .or. rand_type.eq.3) then
               call da_deallocate_observations (iv)
               call da_deallocate_y (ob)
               if (use_background_errors) call da_deallocate_background_errors (be)

#ifdef DM_PARALLEL
               IF ( rootproc ) CALL end_timing('da_end_timing: GEN_OMEGA')
#endif
               call da_wrfvar_finalize

               call wrf_message("*** WRF-Var-RIOT generated OMEGA successfully ***")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: WRFVARRUN')
#endif
               if (trace_use) call da_trace_exit ("da_solve")

               if (trace_use) call da_trace_exit("da_wrfvar_main")
               if (trace_use) call da_trace_report

               call wrfu_finalize
               call wrf_shutdown
               stop
            end if
         end if

         if ( rand_type.eq.3 .and. rand_stage.eq.1 .and. rand_inner_it.eq.1 ) then
            allocate(omega(cv_size,1))

            !Perturb initial state and observations of all but first ensemble member (mean)
            if (ensmember .gt. 1) then
               put_rand_seed = .true.
               seed_array1 = 10*ensmember 
               seed_array2 = 10*nsamp*global_it


               write(unit=message(1),fmt='(a,2i16)') 'Ensemble random number seed_array1, seed_array2 = ', seed_array1, seed_array2
               call da_message(message(1:1))

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: OUTER_it_INIT0')
#endif

               !! Initialize model grid
               disable_traj = .true.
               call da_nl_model(1)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif

               !!-------------------------------------------------------------------
               !! Perturb initial state and T0 boundary conditions with random noise
               !!-------------------------------------------------------------------
               write(cproc,fmt='(i4.4)') myproc
               write(cens,fmt='(i4.4)') ensmember

               if (read_omega) then
                  !Read omega from file
                  filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
                  call da_cv_io ('r', filename, cv_size, omega(:,1), grid, mz)
               else
                  call da_gen_omega(omega, cv_size, 1, grid, mz)
               end if
!               je_start   = be%cv%size_jb + 1
!               je_end     = be%cv%size_jb + be%cv%size_je
!               jp_start   = be%cv%size_jb + be%cv%size_je + 1
!               jp_end     = be%cv%size_jb + be%cv%size_je + be%cv%size_jp
!               jl_start   = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + 1
!               jl_end     = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_jl
!               if ( be%cv%size_jp.gt.0) omega(jp_start:jp_end,1) = 0.D0
!               if ( be%cv%size_je.gt.0) omega(je_start:je_end,1) = 0.D0
!               if ( .not. var4d_lbc .and. be%cv%size_jl.gt.0) omega(jl_start:jl_end,1) = 0.D0

               call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                                        omega(1:be%cv%size_jb,1),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                                        , grid%vch &
#endif
                                        )

               !Transform initial state update to T0 boundary condition update
               call da_transform_xatowrftl_lbc_wrapper( &
                       global_it, cv_size, be, grid%ep, &
                       omega(:,1), iv, grid%vp6, grid%vv6, xbx, &
                       grid, config_flags &
!#if (WRF_CHEM == 1)
!                      , grid%vch6 &
!#endif
                       )

               call da_transform_xtoxa (grid)

               call da_transfer_xatoanalysis (it, xbx, grid, config_flags)

               write(outerloop,'(i2.2)') it
               call da_update_firstguess(grid,'wrfvar_input_'//outerloop//'.e'//cens)

               call da_transfer_headtomodel (grid)

               cvt = cvt + omega(:,1)

               disable_traj = .false.


               !!-------------------------------------------------------------------
               !! Perturb observations with random noise
               !!-------------------------------------------------------------------
               omb_add_noise = .true.

            else
               omega = 0.D0
            end if

            if (.not.read_omega) then
               filename = trim(randvecdir)//'omega.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
               call da_cv_io ('w', filename, cv_size, omega(1:cv_size,1), grid, mz)
            end if

            deallocate(omega)
         end if
         if ((rand_type.eq.3 .and. rand_stage.eq.4) .OR. &
             (rand_type.eq.6 .and. rand_stage.eq.2)) allocate (qhat(1:cv_size, 0:nmodes_global))
         allocate (eignvec(nmodes_global, nmodes_global))
         allocate (eignval(nmodes_global))
      end if
#endif
#endif

! Re-scale the variances and the scale-length for outer-loop > 1:
      if (it > 1 .and. (cv_options == 5 .or. cv_options == 7)) then
         print '(/10X,"===> Re-set BE SCALINGS for outer-loop=",i2)', it
         call  da_scale_background_errors ( be, it )
      else if (it > 1 .and. cv_options == 3) then
         print '(/10X,"===> Re-set CV3 BE SCALINGS for outer-loop=",i2)', it
         call  da_scale_background_errors_cv3 ( grid, be, it )
      endif

      call da_initialize_cv (cv_size, xhat)

      ! [8.1] Calculate nonlinear model trajectory 
      call wrf_message("*** WRF-Var Part 8.1 ***")

!     if (var4d .and. multi_inc /= 2 ) then
      if (var4d) then
#ifdef VAR4D
         if (it > 1) then
            call da_transfer_headtomodel (grid)

            ! Turn off model boundary reading as we already provide a new one.
            if ( checkpoint_interval.eq.0 .or. var4d_lbc .or. .not.boundary_io) &
               call da_model_lbc_off
         endif !end it>1 test

         if (adtl_run_hours .le. 0) adtl_run_hours=run_hours

         call nl_set_var4d_run (head_grid%id, .true.)

         !Restore default trajectory setting
         write_checkpoints = 0

#if defined(LAPACK)
         ! Avoid re-calling nl model in stages where trajectory is not needed
         if (use_randomblock .and. (&
             (rand_type.eq.3 .and. any(rand_stage.eq.(/2,4/)) .OR. &
             (rand_type.eq.6 .and. rand_stage.eq.2)) &

             ) .and. (.not.var4d_lbc .or. checkpoint_interval.gt.0)) disable_traj = .true.

#endif

         if ( checkpoint_interval .gt. 0) then
            write_checkpoints = 1
         end if

#if defined(LAPACK)
         ! Avoid re-calling nl model when obs trajectory values are stored in files for this outer iteration
         ! (currently only chemobs or non-chemobs with checkpoint_interval>0 
         !   have been enabled to read H(x) from files instead of calling NL model)
         if ( use_randomblock .and. &
              ((rand_stage.gt.0 .and. rand_type.ne.3) .or. (rand_stage.gt.1 .and. rand_type.eq.3)) &
              .and. ( checkpoint_interval.gt.0 &
#if (WRF_CHEM == 1)
                .or. (.not.use_offlineobs .and. use_onlineobs .and. use_chemobs) &
#endif
               ) ) then
            write_checkpoints = 2
            disable_traj = .true.
         end if
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: OUTER_it_INIT1')
#endif

         call da_nl_model(it)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif

         ! Use this to exit after writing chem hx values
         ! or after stage 0 of random minimisation
         ii=0

#if (WRF_CHEM == 1)
         if ( use_chemobs ) then
#if defined(LAPACK)
            if (.not.disable_traj) then
#endif
               if (use_onlineobs) then
                  call da_retrieve_chem_hx( model_grid, iv )
                  call da_write_obs_chem( iv, ob, global_it )
               end if
#if defined(LAPACK)
            else
               if (write_checkpoints.eq.2) then
                  call da_read_obs_chem_again( iv, ob, global_it, model_grid)
               else
                  use_onlineobs = .false.
               end if
            end if
#endif
            if (config_flags%calc_hx_only) ii = ii + 1
            call wrf_message("*** WRF-Var generated chem_hx successfully ***")
         end if
#endif
#if defined(LAPACK)
         if (use_randomblock .and. rand_stage.le.0) ii = ii + 1
#endif

#if defined(LAPACK) || (WRF_CHEM == 1)
         if (ii .gt. 0) then
            call da_deallocate_observations (iv)
            call da_deallocate_y (ob)
            if (use_background_errors) call da_deallocate_background_errors (be)

            call da_wrfvar_finalize

#ifdef VAR4D
            if (var4d) call clean_4dvar
            if (var4d) call da_finalize_model
#endif

#ifdef DM_PARALLEL
            IF ( rootproc ) CALL end_timing('da_end_timing: CHEM_HX')
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: WRFVARRUN')
#endif

            if (trace_use) call da_trace_exit("da_wrfvar_main")
            if (trace_use) call da_trace_report

            call wrfu_finalize
            call wrf_shutdown
            stop
         end if
#endif

         !Restore default trajectory setting
         write_checkpoints = 0

!         adtl_run_hours = run_hours0


!      elseif (var4d .and. multi_inc == 2 ) then
#else
         write(unit=message(1),fmt='(A)')'Please re-compile the code with 4dvar option'
         call da_error(__FILE__,__LINE__,message(1:1))
#endif
!end VAR4D test
      end if

      ! [8.2] Calculate innovation vector (O-B):
      call wrf_message("*** WRF-Var Part 8.2 ***")

      num_qcstat_conv=0
#if defined(LAPACK)
      if ((.not.disable_traj &
#if (WRF_CHEM == 1)
           .and. use_offlineobs) .or. (use_onlineobs .and. use_chemobs &
#endif
         )) &
#endif
         call da_get_innov_vector (global_it, num_qcstat_conv, ob, iv, grid , config_flags)
      if ( multi_inc == 1 ) then 
         if (trace_use) call da_trace_exit ("da_solve")
         return
      end if

      if (test_transforms .or. test_gradient) then

         if (test_gradient) then
            call da_allocate_y (iv, re)
            call da_allocate_y (iv, y)
            call da_check_gradient (grid, config_flags, cv_size, xhat, cvt, 1.0e-10, 8, &
                 xbx, be, iv, y, re, j)
            call da_deallocate_y (re)
            call da_deallocate_y (y)
         endif

         if (test_transforms) then
            call da_check (grid, config_flags, cv_size, xbx, be, grid%ep, iv, &
                        grid%vv, grid%vp, y)
         endif

         if (trace_use) call da_trace_exit("da_solve")
         return

      end if


#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: OUTER_it_INIT2')
#endif

      ! [8.4] Minimize cost function:
      call wrf_message("*** WRF-Var Part 8.4 ***")

      call da_allocate_y (iv, re)
      call da_allocate_y (iv, y)
#if (WRF_CHEM == 1)
      call da_allocate_y_chem (iv, re)
      call da_allocate_y_chem (iv, y)
#endif

      !=======================================================
      ! 3 options (use_lanczos, use_randomblock, default - cg)
      !=======================================================
      if (use_lanczos .or. use_randomblock) then
#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif
         if (use_lanczos) then
            if (read_lanczos) then
               call da_lanczos_io('r',cv_size,ntmax(it),neign,eignvec,eignval,qhat)

               call da_kmat_mul(grid,config_flags,it,cv_size,xbx, &
                                be,iv,xhat,qhat,cvt,re,y,j,eignvec,eignval,neign)
             ! Output Cost Function
               call da_calculate_j(it, 1, cv_size, be%cv%size_jb, be%cv%size_je, be%cv%size_jp, &
                                be%cv%size_jl, xbx, be, iv, xhat, cvt, re, y, j, xhat, grid, config_flags )

            else
               if (var4d_inc_out) then
                  hessian_info%nmodes_max = 0
                  allocate(hessian_info%eigs(nmodes_global))
                  allocate(hessian_info%qhat(cv_size,0:nmodes_global))
                  allocate(hessian_info%ghat(cv_size))
               end if

               call da_minimise_lz(grid, config_flags, global_it, cv_size, xbx,& 
                                be, iv, j_grad_norm_target, xhat, qhat, cvt, re, y, j, eignvec, eignval, neign, j0_cost, j0_total, hessian_info=hessian_info)
            end if

            if (write_lanczos) call da_lanczos_io('w',cv_size,ntmax(it),neign,eignvec,eignval,qhat)

            if (adj_sens) call da_sensitivity(grid,config_flags,it,cv_size,xbx, &
                                              be,iv,xhat,qhat,cvt,y,eignvec,eignval,neign )

#ifdef VAR4D
#if defined(LAPACK)
         else if ( use_randomblock .and. var4d) then

            j0_total = 0.D0
            if ( nsamp.gt.0 .and. global_it.gt.0 ) then
               if (rand_stage.eq.1 .and. &
                   ((ensmember.eq.nsamp+1) .or. rand_type.eq.3)) then
                  select case (rand_type)
                     case(3)
                        call da_block_lanczos( &
                              grid, mz, config_flags, cv_size, xbx, &
                              be, iv, cvt, re, y, xhat, nsamp, j_cost = j)

                     case(6)
                        call da_randomise_svd( &
                              grid, mz, config_flags, cv_size, xbx, &
                              be, iv, cvt, re, y, xhat, nsamp, j_cost = j)
                     case default
                        write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: Invalid selection for rand_type', rand_type
                        CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
                  end select
               else
                  select case (rand_type)
                     case(3)
                        if (rand_stage .eq. 4) then
                           if (var4d_inc_out) then
                              allocate(hessian_info%eigs(1))
                              allocate(hessian_info%qhat(cv_size,0:nmodes_global))
                              allocate(hessian_info%ghat(cv_size))
                           end if

                           call da_block_lanczos( &
                                 grid, mz, config_flags, cv_size, xbx, &
                                 be, iv, cvt, re, y, xhat, nsamp, &
                                 qhat=qhat, eignval=eignval, &
                                 eignvec=eignvec, neign=neign , &
                                 hessian_info=hessian_info)
                        else
                           call da_block_lanczos( &
                                 grid, mz, config_flags, cv_size, xbx, &
                                 be, iv, cvt, re, y, xhat, nsamp)
                        end if
                     case(6)
                        if (rand_stage .eq. 2) then
                           if (var4d_inc_out) then
                              allocate(hessian_info%eigs(nmodes_global))
                              allocate(hessian_info%qhat(cv_size,0:nmodes_global))
                              allocate(hessian_info%ghat(cv_size))
                           end if

                           call da_randomise_svd( &
                                 grid, mz, config_flags, cv_size, xbx, &
                                 be, iv, cvt, re, y, xhat, nsamp, &
                                 qhat=qhat, eignval=eignval, &
                                 eignvec=eignvec, neign=neign, &
                                 hessian_info=hessian_info)
                        else
                           call da_randomise_svd( &
                                 grid, mz, config_flags, cv_size, xbx, &
                                 be, iv, cvt, re, y, xhat, nsamp)
                        end if
                     case default
                        write(unit=stdout,fmt='(A,I2)') 'FATAL ERROR: Invalid selection for rand_type', rand_type
                        CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
                  end select
               end if
               if ( &
                   (rand_type.eq.3 .and. (rand_stage.ne.4 .or. rand_inner_it.lt.max_rand_inner)) .OR. &
                   (rand_type.eq.6 .and. rand_stage.ne.2) &
                   ) rand_exit = .true.
            else
               write(unit=stdout,fmt='("nsamp, rand_outer = ",2I4)') &
                   nsamp, global_it

               CALL wrf_error_fatal('In randomised SVD. Check namelist settings.')
            end if

#endif
#endif
         else
            CALL wrf_error_fatal('RIOT can only be used with 4D-Var and requires LAPACK (e.g., from MKL).')
         end if

#ifdef DM_PARALLEL
         IF ( rootproc ) CALL end_timing('da_end_timing: HESS_vector_multiply')
#endif


         if ( var4d .and. var4d_inc_out .and. allocated(hessian_info%ghat) ) then
#ifdef DM_PARALLEL
            IF ( rootproc ) CALL start_timing !Level 1
#endif
            call da_output_increments( grid,config_flags,mz, &
                                       global_it,cv_size,xbx,be, &
                                       cvt,hessian_info)

            if (rand_type .eq. 3) then
               if(allocated(hessian_info%eigs(1)%eignval)) &
                  deallocate(hessian_info%eigs(1)%eignval)
               if(allocated(hessian_info%eigs(1)%eignvec)) &
                  deallocate(hessian_info%eigs(1)%eignvec)
            else
               do ii = 1, hessian_info%nmodes_max
                  if(allocated(hessian_info%eigs(ii)%eignval)) &
                     deallocate(hessian_info%eigs(ii)%eignval)
                  if(allocated(hessian_info%eigs(ii)%eignvec)) &
                     deallocate(hessian_info%eigs(ii)%eignvec)
               end do
            end if

            deallocate(hessian_info%eigs)
            deallocate(hessian_info%qhat)
            deallocate(hessian_info%ghat)
#ifdef DM_PARALLEL
            IF ( rootproc ) CALL end_timing('da_end_timing: increment_diagnostics2')
#endif
         end if

#ifdef VAR4D
#if defined(LAPACK)
         if ( use_randomblock .and. var4d .and. rand_exit) then
            call da_deallocate_observations (iv)
            call da_deallocate_y (ob)
            if (use_background_errors) call da_deallocate_background_errors (be)

            call da_wrfvar_finalize
            if (var4d) call clean_4dvar
            if (var4d) call da_finalize_model

            if (rand_type .ne. 3 .or. rand_stage.lt.4) then
               write(unit=message(1),fmt='(A,I2,A,I2,A)') &
                 "***Completed Stage ",rand_stage," of RIOT in outer iteration",global_it,"***"
               if (rand_type .eq. 3) &
                  write(unit=message(1),fmt='(A,I2,A)') &
                    "***Continuing Block Lanczos inner iteration ",rand_inner_it,"...***"
            else
               write(unit=message(1),fmt='(A,I2,A,I2,A)') &
                 "***Completed Block Lanczos inner iteration ",rand_inner_it," of RIOT in outer iteration",global_it,"***"
            end if
            call wrf_message(message(1))

            if (trace_use) call da_trace_exit("da_wrfvar_main")
            if (trace_use) call da_trace_report

#ifdef DM_PARALLEL
            IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_STAGE')   
#endif
            call wrfu_finalize
            call wrf_shutdown
            stop
         end if
#endif

#if (WRF_CHEM == 1)
#ifdef DM_PARALLEL
         IF ( rootproc ) CALL start_timing
#endif
         if (calc_aminusb) call da_calculate_aminusb(grid,config_flags,global_it,cv_size,xbx, &
                             be,qhat,eignvec,eignval,neign,ntmax(global_it))
#ifdef DM_PARALLEL
         IF ( rootproc ) CALL end_timing('da_end_timing: aminusb')
#endif
#endif

#endif
      else
         call da_minimise_cg( grid, config_flags, it, be % cv % size, & 
              xbx, be, iv, j_grad_norm_target, xhat, cvt, re, y, j, j0_cost, j0_total)
      end if



#ifdef VAR4D
      if (var4d) then
#ifdef DM_PARALLEL
      IF ( rootproc ) CALL start_timing !Level 1
#endif

      if (dgn_search(global_it).eq.1 .and. .not.var4d_lbc) then
         call da_dgn(grid, config_flags, cv_size, xhat, cvt, &
                     xbx, be, iv, ob, &
                     j0_cost, j0_total, global_it, mz)
#if (WRF_CHEM == 1)
      else if ( dgn_search(global_it) .eq. 2 .and. eta_emiss .lt. 1.D0 .and. eta_emiss .gt. 0.D0 ) then
         call da_hdgn(global_it, cv_size, be, xhat)
#endif
      end if

#ifdef DM_PARALLEL
      IF ( rootproc ) CALL end_timing('da_end_timing: LINESEARCH')
#endif
      end if
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif
      ! Update outer-loop control variable
      cvt = cvt + xhat

#if defined(LAPACK)
      if ((use_randomblock .and. rand_minimise)) then
         write(cproc,fmt='(i4.4)') myproc
         write(cit,fmt='(i2.2)') global_it
         cvtfile_rand = 'cvt.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', cvtfile_rand, cv_size, cvt(1:cv_size), grid, mz)
      else
#endif
      if (var4d_inc_out) then
         write(cproc,fmt='(i4.4)') myproc
         write(cit,fmt='(i2.2)') global_it
         cvtfile_rand = 'cvt_REF_out.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
         call da_cv_io ('w', cvtfile_rand, cv_size, cvt(1:cv_size), grid, mz)
      endif

      if ( outer_loop_restart ) then
        open(unit=cvt_unit,status='unknown',file=trim(cvtfile),iostat=iost,form='UNFORMATTED')
        if (iost /= 0) then
           write(unit=message(1),fmt='(A,I5,A)') &
             "Error ",iost," opening cvt file "//trim(cvtfile)
          call da_error(__FILE__,__LINE__,message(1:1))
        end if
        write(unit=message(1),fmt='(a)') 'Writing cvt to : '//trim(cvtfile) 
        call da_message(message(1:1))
        write(cvt_unit) cvt
        close(cvt_unit)
        !call da_free_unit(cvt_unit)
      end if
#if defined(LAPACK)
      end if
#endif
      !------------------------------------------------------------------------
#if (WRF_CHEM == 1)
      if (osse_chem .or. use_iter_metrics) then

#if defined(LAPACK)
         if (use_randomblock .and. rand_minimise) then
            allocate(xhat_prev(1:cv_size))

            if (global_it.gt.1) then 
               write(cproc,fmt='(i4.4)') myproc
               write(cit,fmt='(i2.2)') global_it-1

               cvtfile_rand = 'xhat.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
               call da_cv_io ('r', cvtfile_rand, cv_size, xhat_prev(1:cv_size), grid, mz)
           else
               xhat_prev = 0.D0
            end if
         end if
#endif
         if (it.eq.1.and..not.use_randomblock) then
            allocate(xhat_prev(1:cv_size))
            xhat_prev = 0.D0
         end if

         call da_evaluate_increment(grid, xbx, be, mz, cv_size, &
                                     xhat_prev, xhat)

         deallocate(xhat_prev)

         neign_trunc = neign
#if defined(LAPACK)
         if (use_randomblock .and. rand_minimise) then
            write(cproc,fmt='(i4.4)') myproc
            write(cit,fmt='(i2.2)') global_it
            cvtfile_rand = 'xhat.it'//trim(adjustl(cit))//'.p'//trim(adjustl(cproc))
            call da_cv_io ('w', cvtfile_rand, cv_size, xhat(1:cv_size), grid, mz)
            neign = nsamp
         else if(it .lt. max_ext_its .and. .not.use_randomblock) then
#endif
            allocate(xhat_prev(1:cv_size))
            xhat_prev = xhat
#if defined(LAPACK)
         end if
#endif
         if ((use_lanczos .or. use_randomblock) .and. evaluate_hessian) &
            call da_evaluate_hessian(grid, xbx, be, mz, cv_size, &
                                   neign, &
                                   neign, neign_trunc, &
                                   eignval(1:neign), eignvec(1:neign,1:neign), qhat(:,0:neign))
      end if
#endif
      !------------------------------------------------------------------------

      ! reset cv to random noise
      if (anal_type_randomcv) then
         call da_set_randomcv (cv_size, xhat)
      end if

      ! [8.5] Update latest analysis solution:
      call wrf_message("*** WRF-Var Part 8.5 ***")

      if (.not. var4d) then
         call da_transform_vtox (grid,cv_size,xbx,be,grid%ep,xhat,grid%vv,grid%vp)
      else
         call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,xhat(1:be%cv%size_jb),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
         , grid%vch &
#endif
         )

#ifdef VAR4D
#if defined(LAPACK)
         ! Should always call this subroutine (even for non-random),
         ! and eliminate final call to da_calculate_gradj in both minimise_cg/lz routines

         !Transform initial state update to T0 boundary condition update
         if (use_randomblock) &
            call da_transform_xatowrftl_lbc_wrapper( &
                    global_it, cv_size, be, grid%ep, &
                    xhat, iv, grid%vp6, grid%vv6, xbx, &
                    grid, config_flags &
!#if (WRF_CHEM == 1)
!                   , grid%vch6 &
!#endif
                    )
#endif
#endif
      endif
      call da_transform_xtoxa (grid)

      ! [8.6] Only when use_radarobs = .false. and calc_w_increment =.true.,
      !       the w_increment need to be diagnosed:
      if (calc_w_increment .and. .not. use_radarobs .and. .not. var4d) then
         call wrf_message("*** WRF-Var Part 8.6 ***")

         call da_uvprho_to_w_lin (grid)

#ifdef DM_PARALLEL
#include "HALO_RADAR_XA_W.inc"
#endif
      end if

      ! [8.7] Write out diagnostics
      call wrf_message("*** WRF-Var Part 8.7 ***")

      call da_write_diagnostics (global_it, grid, num_qcstat_conv, ob, iv, re, y, j)

      ! Write "clean" QCed observations if requested:
      if (anal_type_qcobs) then
         ! if (it == 1) then
          if (write_mod_filtered_obs) then
            call da_write_modified_filtered_obs (grid, ob, iv, &
               coarse_ix, coarse_jy, start_x, start_y)
          else
            call da_write_filtered_obs (global_it, grid, ob, iv, &
               coarse_ix, coarse_jy, start_x, start_y)
          end if     
         ! end if     
      end if

      ! [8.7.1] Write Ascii radar OMB and OMA file

      if (use_radarobs) then
         call da_write_oa_radar_ascii (ob,iv,re,global_it)
      end if

      ! [8.3] Interpolate x_g to low resolution grid

      ! [8.8] Write Ascii radiance OMB and OMA file

#if defined(CRTM) || defined(RTTOV)
      if (use_rad .and. write_oa_rad_ascii) then
         call da_write_oa_rad_ascii (global_it,ob,iv,re)
      end if
#endif

      ! [8.9] Update VarBC parameters and write output file
#if defined(CRTM) || defined(RTTOV)      
      if ( use_rad .and. (use_varbc.or.freeze_varbc) ) &
                call da_varbc_update(global_it, cv_size, xhat, iv)
#endif


      !------------------------------------------------------------------------
      ! [8.10] Output WRFVAR analysis and analysis increments:
      !------------------------------------------------------------------------
      call wrf_message("*** WRF-Var Part 8.10 ***")

      call da_transfer_xatoanalysis (global_it, xbx, grid, config_flags)

     if ( it < max_ext_its .and. &
          (print_detail_outerloop &
#ifdef VAR4D
          .or. checkpoint_interval.gt.0 &
          .or. (var4d .and. dgn_search(global_it+1).eq.1 .and. .not.var4d_lbc)&
#endif
           ) ) then
        write(outerloop,'(i2.2)') it
        call da_update_firstguess(grid,'wrfvar_output_'//outerloop)
#ifdef VAR4D
        if (var4d .and. (&
            outer_loop_restart &
            .or. wrfvar_bdy_out &
            .or. (dgn_search(global_it+1).eq.1 .and. .not.var4d_lbc) &
#if defined(LAPACK)
            .or. (use_randomblock .and. rand_minimise) &
#endif
            )) call da_med_initialdata_output_lbc (grid , config_flags, 'wrfvar_bdyout_'//outerloop)
#endif
     end if

     call da_deallocate_y (re)
     call da_deallocate_y (y)


   ! Deallocate arrays which depend on ntmax
     if (use_lanczos &
#if defined(LAPACK)
         .or. use_randomblock &
#endif
         ) then
        deallocate (qhat)
        deallocate (eignvec)
        deallocate (eignval)
     end if

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: OUTER_it_FINISH')
#endif

   end do !END OF OUTER ITERATION

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing !Level 1
#endif

   ! output wrfvar analysis

   if ((config_flags%real_data_init_type == 1) .or. &
       (config_flags%real_data_init_type == 3)) then
      call da_update_firstguess(input_grid)
#ifdef VAR4D
      if (var4d .and. (&
            outer_loop_restart &
            .or. (var4d .and. dgn_search(global_it+1).eq.1 .and. .not.var4d_lbc) &
            .or. wrfvar_bdy_out &
#if defined(LAPACK)
            .or. (use_randomblock .and. rand_minimise) &
#endif
            )) call da_med_initialdata_output_lbc (input_grid , config_flags, 'wrfvar_bdyout')

      !if (var4d) call da_med_initialdata_output_lbc (head_grid , config_flags)
      if ( var4d_lbc ) then
         call domain_clock_get (grid, stop_timestr=timestr1)
         call domain_clock_set( grid, current_timestr=timestr1 )
         call da_med_initialdata_input (grid, config_flags, 'fg02')
         call da_setup_firstguess(xbx, grid, config_flags, .false. )
         shuffle = grid%xa
         jl_start    = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + 1
         jl_end      = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_jl
         grid%xa  = grid%x6a
         call da_transform_vtox(grid, be%cv%size_jl, xbx, be, grid%ep, &
              xhat(jl_start:jl_end), grid%vv6, grid%vp6)
         grid%xa  = shuffle
         call da_transfer_xatoanalysis (it, xbx, grid, config_flags)
         call da_update_firstguess (grid, 'ana02')
         call domain_clock_get (grid, start_timestr=timestr1)
         call domain_clock_set( grid, current_timestr=timestr1 )
      endif
#endif
      call med_shutdown_io (input_grid, config_flags)
   end if

#if (WRF_CHEM == 1)
   else
      if(checkpoint_interval .gt.0 .and. mod(60*run_hours,checkpoint_interval) /= 0) then
         CALL wrf_error_fatal('In checkpoint sensitivity runs, mod(60*run_hours,checkpoint_interval) must be zero.')
      end if

      ! Read chem observations for use in online forcing

      ! Write nonlinear model checkpoint trajectories to disk and evaluate online
      ! observation operators
      if(checkpoint_interval .gt.0 ) write_checkpoints = 1
      call da_nl_model(1)
      if(checkpoint_interval .gt.0 ) write_checkpoints = 0

     if ( use_chemobs ) then
        call da_retrieve_chem_hx( model_grid, iv )

        call da_write_obs_chem( iv, ob, 0 )
     end if

      ! Carry out checkpointed adjoint sensitivity
      if( .not. config_flags%calc_hx_only) then
         !Calculate forcing for a cost function (does not have to be chem-specific)
         ! Chemical 4D-Var cost function
         if ( use_chemobs ) then
            call da_calculate_chem_forcing_ad ( 1, model_grid, iv )
         end if
         call da_ad_model
      end if

   end if
#endif

   !---------------------------------------------------------------------------
   ! [9.0] Tidy up:
   !---------------------------------------------------------------------------
   call wrf_message("*** WRF-Var Part 9.0 ***")


   deallocate (cvt)
   deallocate (xhat)
!   if (use_lanczos) then
!      deallocate (full_eignvec)
!   end if

   ! clean up radiance related arrays
#if defined(RTTOV) || defined(CRTM)
   if (use_rad) then
      call da_deallocate_radiance (ob, iv, j)
      deallocate (time_slots)
#ifdef RTTOV
      if (rtm_option == rtm_option_rttov) then
         deallocate (coefs)
         deallocate (opts)
      end if
#endif
   end if
#endif

   if (var4d .and. use_rainobs) deallocate(fgat_rain_flags)
   call da_deallocate_observations (iv)
   call da_deallocate_y (ob)
   if (use_background_errors) call da_deallocate_background_errors (be)

   if (xbx%pad_num > 0) then
      deallocate (xbx%pad_loc)
      deallocate (xbx%pad_pos)
   end if

   deallocate (xbx % fft_factors_x)
   deallocate (xbx % fft_factors_y)
   deallocate (xbx % fft_coeffs)
   deallocate (xbx % trig_functs_x)
   deallocate (xbx % trig_functs_y)

   if (global) then
      deallocate (xbx%coslat)
      deallocate (xbx%sinlat)
      deallocate (xbx%coslon)
      deallocate (xbx%sinlon)
      deallocate (xbx%int_wgts)
      deallocate (xbx%alp)
      deallocate (xbx%wsave)
      if (jts == jds) then
         deallocate (cos_xls)
         deallocate (sin_xls)
      end if
                                                                                
      if (jte == jde) then
         deallocate (cos_xle)
         deallocate (sin_xle)
      end if
   end if

   if ( anal_type_hybrid_dual_res ) deallocate(ob_locs)



#ifdef DM_PARALLEL
   call mpi_barrier (comm,ierr)
#endif

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: WRFVAR_FINISH')
#endif



   if (trace_use) call da_trace_exit ("da_solve")


contains

#include "da_solve_init.inc"
#include "da_solve_dual_res_init.inc"

end subroutine da_solve

