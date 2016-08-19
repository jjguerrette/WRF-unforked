!WRF:MEDIATION_LAYER:SOLVER

#define BENCH_START(A)
#define BENCH_END(A)

MODULE a_module_first_rk_step_part1

CONTAINS

  SUBROUTINE a_first_rk_step_part1 (   grid , config_flags              &
                             , moist , a_moist , moist_tend , a_moist_tend       &
                             , chem  , a_chem , chem_tend , a_chem_tend          &
#ifdef WRF_CHEM
                             , dvel &
#endif
                             , tracer, a_tracer, tracer_tend, a_tracer_tend      &
                             , scalar , a_scalar , scalar_tend , a_scalar_tend   &
                             , fdda3d, fdda2d                   &
                             , aerod                            &
                             , ru_tendf, a_ru_tendf, rv_tendf, a_rv_tendf        &
                             , rw_tendf, a_rw_tendf, t_tendf, a_t_tendf          &
                             , ph_tendf, a_ph_tendf, mu_tendf, a_mu_tendf        &
                             , tke_tend, a_tke_tend             &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , wspd , a_wspd , gz1oz0 , chklowq &
                             , cu_act_flag , hol , th_phy , a_th_phy      &
                             , pi_phy , a_pi_phy , p_phy , a_p_phy &!, t_phy , a_t_phy &
                             , dz8w , a_dz8w , p8w , a_p8w , t8w , a_t8w   &
                             , ids, ide, jds, jde, kds, kde     &
                             , ims, ime, jms, jme, kms, kme     &
                             , ips, ipe, jps, jpe, kps, kpe     &
                             , imsx,imex,jmsx,jmex,kmsx,kmex    &
                             , ipsx,ipex,jpsx,jpex,kpsx,kpex    &
                             , imsy,imey,jmsy,jmey,kmsy,kmey    &
                             , ipsy,ipey,jpsy,jpey,kpsy,kpey    &
                             , k_start , k_end                  &
                             , f_flux                           &
                            )
    USE module_state_description
    USE module_model_constants
    USE module_domain, ONLY : domain, domain_clock_get, get_ijk_from_subgrid
    USE module_configure, ONLY : grid_config_rec_type, model_config_rec
!   USE module_radiation_driver, ONLY : pre_radiation_driver, radiation_driver
   USE module_surface_driver, ONLY : surface_driver
!   USE module_cumulus_driver, ONLY : cumulus_driver
!   USE module_shallowcu_driver, ONLY : shallowcu_driver
   USE module_pbl_driver, ONLY : pbl_driver
!   USE module_fr_sfire_driver_wrf, ONLY : sfire_driver_em_step
!   USE module_fddagd_driver, ONLY : fddagd_driver
   USE module_em, ONLY : init_zero_tendency
!   USE module_force_scm
!   USE module_convtrans_prep
    USE module_big_step_utilities_em, ONLY : phy_prep
!use module_scalar_tables
#ifdef DM_PARALLEL
    USE module_dm, ONLY : local_communicator, mytask, ntasks, ntasks_x, ntasks_y, local_communicator_periodic, wrf_dm_maxval
    USE module_comm_dm, ONLY : halo_em_phys_a_sub,halo_em_phys_a_ad_sub,halo_pwp_sub,halo_em_phys_a_bl_surf_sub,halo_em_chem_e_3_sub,halo_em_chem_e_5_sub,halo_em_chem_e_3_ad_sub,halo_em_chem_e_5_ad_sub
#endif
#ifdef WRF_CHEM
    USE module_input_chem_data, only: get_last_gas,mozcart_lbc_set
    USE module_dry_dep_vel
    USE a_module_dry_dep_vel
#endif
    USE module_utility
    USE a_module_em, ONLY : a_init_zero_tendency
    USE a_module_big_step_utilities_em, ONLY : a_phy_prep
    USE a_module_pbl_driver, ONLY : a_pbl_driver
    USE a_module_cumulus_driver, ONLY : a_cumulus_driver
    USE a_module_surface_driver, ONLY : a_surface_driver


    IMPLICIT NONE

    TYPE ( domain ), INTENT(INOUT) :: grid
    TYPE ( grid_config_rec_type ), INTENT(IN) :: config_flags
    TYPE(WRFU_Time)                :: currentTime

    INTEGER, INTENT(IN) :: ids, ide, jds, jde, kds, kde,     &
                           ims, ime, jms, jme, kms, kme,     &
                           ips, ipe, jps, jpe, kps, kpe,     &
                           imsx,imex,jmsx,jmex,kmsx,kmex,    &
                           ipsx,ipex,jpsx,jpex,kpsx,kpex,    &
                           imsy,imey,jmsy,jmey,kmsy,kmey,    &
                           ipsy,ipey,jpsy,jpey,kpsy,kpey


    LOGICAL ,INTENT(IN)                        :: adapt_step_flag
    REAL, INTENT(IN)                           :: curr_secs

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: a_moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: a_moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: a_chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: a_chem_tend
    REAL    ,DIMENSION(grid%sm31:grid%em31,1:grid%kdvel,grid%sm33:grid%em33,num_dvel) :: dvel
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: a_tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: a_tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: a_scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: a_scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_fdda3d),INTENT(INOUT)  :: fdda3d
    REAL    ,DIMENSION(ims:ime,1:1,jms:jme,num_fdda2d),INTENT(INOUT)      :: fdda2d
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_aerod),INTENT(INOUT)   :: aerod
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psim
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psih
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: a_wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gz1oz0
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chklowq
    LOGICAL ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cu_act_flag
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: hol

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_p_phy
!    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy
!    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_t_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_t8w

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: tke_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: a_tke_tend

    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: mu_tendf
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: a_mu_tendf

    INTEGER, INTENT(IN)                           ::  k_start, k_end
    LOGICAL, INTENT(IN), OPTIONAL                 ::  f_flux

! Local
    real :: HYDRO_dt
    REAL, DIMENSION( ims:ime, jms:jme ) :: exch_temf  ! 1/7/09 WA

    REAL, DIMENSION( ims:ime, jms:jme ) :: ht_loc, mixht
    INTEGER                             :: ij
    INTEGER  num_roof_layers
    INTEGER  num_wall_layers
    INTEGER  num_road_layers
    INTEGER  iswater
    LOGICAL  :: l_flux
    INTEGER  :: isurban
    INTEGER  rk_step
    INTEGER                         :: yr, month, day, hr, minute, sec, rc
    CHARACTER*80                    :: mesg
    INTEGER                         :: current_month, numgas

   INTEGER                         :: sids , side , sjds , sjde , skds , skde , &
                                      sims , sime , sjms , sjme , skms , skme , &
                                      sips , sipe , sjps , sjpe , skps , skpe

   CHARACTER (LEN=256) :: mminlu
   CHARACTER (LEN=1000) :: message


  CALL get_ijk_from_subgrid (  grid ,                   &
                            sids, side, sjds, sjde, skds, skde,    &
                            sims, sime, sjms, sjme, skms, skme,    &
                            sips, sipe, sjps, sjpe, skps, skpe    )

 ! initialize all tendencies to zero in order to update physics
 ! tendencies first (separate from dry dynamics).

   l_flux=.FALSE.
   if (present(f_flux)) l_flux=f_flux

    rk_step = 1

!BENCH_START(init_zero_tend_tim)
!       !$OMP PARALLEL DO   &
!       !$OMP PRIVATE ( ij )
!
!       DO ij = 1 , grid%num_tiles
!
!         CALL wrf_debug ( 200 , ' call init_zero_tendency' )
!         CALL init_zero_tendency ( ru_tendf, rv_tendf, rw_tendf,     &
!                                   ph_tendf, t_tendf, tke_tend,      &
!                                   mu_tendf,                         &
!                                   moist_tend,chem_tend,scalar_tend, &
!                                   tracer_tend,num_tracer,           &
!                                   num_moist,num_chem,num_scalar,    &
!                                   rk_step,                          &
!                                   ids, ide, jds, jde, kds, kde,     &
!                                   ims, ime, jms, jme, kms, kme,     &
!                                   grid%i_start(ij), grid%i_end(ij), &
!                                   grid%j_start(ij), grid%j_end(ij), &
!                                   k_start, k_end                   )
!
!       END DO
!       !$OMP END PARALLEL DO            
!BENCH_END(init_zero_tend_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_PHYS_A.inc"
#endif

      CALL PUSHREAL8ARRAY(grid%rublten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rvblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rthblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqvblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqcblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqiblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))

!      CALL PUSHREAL8ARRAY(grid%rucuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
!      CALL PUSHREAL8ARRAY(grid%rvcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rthcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqvcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqccuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%rqicuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
!      CALL PUSHREAL8ARRAY(grid%rqscuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call phy_prep' )
        CALL phy_prep ( config_flags,                                    &
                        grid%mut, grid%muu, grid%muv, grid%u_2,          &
                        grid%v_2, grid%p, grid%pb, grid%alt,             &
                        grid%ph_2, grid%phb, grid%t_2, grid%tsk, moist, num_moist,   &
                        grid%rho,th_phy, p_phy, pi_phy, grid%u_phy, grid%v_phy,         &
                        p8w, grid%t_phy, t8w, grid%z, grid%z_at_w, dz8w,      &
                        grid%p_hyd, grid%p_hyd_w, grid%dnw,              &
                        grid%fnm, grid%fnp, grid%znw, grid%p_top,        &
                        grid%rthraten,                                   &
                        grid%rthblten, grid%rublten, grid%rvblten,       &
                        grid%rqvblten, grid%rqcblten, grid%rqiblten,     &
                        grid%rucuten,  grid%rvcuten,  grid%rthcuten,     &
                        grid%rqvcuten, grid%rqccuten, grid%rqrcuten,     &
                        grid%rqicuten, grid%rqscuten,                    &
                        grid%rushten,  grid%rvshten,  grid%rthshten,     &
                        grid%rqvshten, grid%rqcshten, grid%rqrshten,     &
                        grid%rqishten, grid%rqsshten, grid%rqgshten,     &
                        grid%rthften,  grid%rqvften,                     &
                        grid%RUNDGDTEN, grid%RVNDGDTEN, grid%RTHNDGDTEN, &
                        grid%RPHNDGDTEN,grid%RQVNDGDTEN, grid%RMUNDGDTEN,&
!jdf
                        grid%landmask,grid%xland,                 &
!jdf
                        ids, ide, jds, jde, kds, kde,                    &
                        ims, ime, jms, jme, kms, kme,                    &
                        grid%i_start(ij), grid%i_end(ij),                &
                        grid%j_start(ij), grid%j_end(ij),                &
                        k_start, k_end                                   )
      ENDDO
      !$OMP END PARALLEL DO

BENCH_END(phy_prep_tim)
     CALL domain_clock_get( grid, current_time=currentTime, &
                            current_timestr=mesg )
     read(mesg(6:7),FMT='(I2)') current_month

!********* Surface driver
! surface

BENCH_START(surf_driver_tim)

!-----------------------------------------------------------------
! urban related variable are added to arguments of surface_driver
!-----------------------------------------------------------------
      num_roof_layers = grid%num_soil_layers !urban
      num_wall_layers = grid%num_soil_layers !urban
      num_road_layers = grid%num_soil_layers !urban
      CALL nl_get_iswater(grid%id, iswater)
      CALL nl_get_isurban(grid%id, isurban)
      call nl_get_mminlu(grid%id, mminlu)

#ifdef DM_PARALLEL
#     include "HALO_PWP.inc"
#endif

      CALL wrf_debug ( 200 , ' call surface_driver' )

      if( grid%num_nests .lt. 1 )then
          HYDRO_dt = 0
      else
          HYDRO_dt = -1
      endif

      CALL PUSHREAL8ARRAY(grid%ust, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%ustm, (ime-ims+1)*(jme-jms+1))
!      CALL PUSHREAL8ARRAY(grid%rmol, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%mol, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%znt, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%canwat, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%hfx, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%qfx, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%tsk, (ime-ims+1)*(jme-jms+1))
      CALL PUSHREAL8ARRAY(grid%smois, (ime-ims+1)*(jme-jms+1)*grid%num_soil_layers)
      CALL PUSHREAL8ARRAY(grid%tslb, (ime-ims+1)*(jme-jms+1)*grid%num_soil_layers)

      CALL surface_driver(                                                &
     &        HYDRO_dt=HYDRO_dt,  sfcheadrt=grid%sfcheadrt,               &
     &        INFXSRT=grid%INFXSRT,  soldrain=grid%soldrain,              &
     &         ACGRDFLX=grid%acgrdflx  ,ACHFX=grid%achfx        ,ACLHF=grid%aclhf        &
     &        ,ACSNOM=grid%acsnom      ,ACSNOW=grid%acsnow      ,AKHS=grid%akhs          &
     &        ,AKMS=grid%akms          ,ALBBCK=grid%albbck      ,ALBEDO=grid%albedo      &
     &        ,EMBCK=grid%embck                                                          &
     &        ,BR=grid%br              ,CANWAT=grid%canwat      ,CHKLOWQ=chklowq    &
     &        ,CT=grid%ct              ,DT=grid%dt         ,DX=grid%dx         &
     &        ,DZ8W=dz8w          ,DZS=grid%dzs            ,FLHC=grid%flhc          &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &
     &        ,FLQC=grid%flqc          ,GLW=grid%glw            ,GRDFLX=grid%grdflx      &
     &        ,GSW=grid%gsw    ,SWDOWN=grid%swdown        ,GZ1OZ0=gz1oz0      ,HFX=grid%hfx              &
     &        ,HT=grid%ht              ,IFSNOW=config_flags%ifsnow      ,ISFFLX=config_flags%isfflx      &
     &        ,FRACTIONAL_SEAICE=config_flags%fractional_seaice           &
     &        ,SEAICE_ALBEDO_OPT=config_flags%seaice_albedo_opt           &
     &        ,SEAICE_ALBEDO_DEFAULT=config_flags%seaice_albedo_default   &
     &        ,SEAICE_THICKNESS_OPT=config_flags%seaice_thickness_opt     &
     &        ,SEAICE_THICKNESS_DEFAULT=config_flags%seaice_thickness_default     &
     &        ,SEAICE_SNOWDEPTH_OPT=config_flags%seaice_snowdepth_opt     &
     &        ,SEAICE_SNOWDEPTH_MAX=config_flags%seaice_snowdepth_max     &
     &        ,SEAICE_SNOWDEPTH_MIN=config_flags%seaice_snowdepth_min     &
     &        ,TICE2TSK_IF2COLD=config_flags%tice2tsk_if2cold             &
     &        ,IFNDALBSI=grid%ifndalbsi, IFNDICEDEPTH=grid%ifndicedepth   &
     &        ,IFNDSNOWSI=grid%ifndsnowsi                                 &
     &        ,ISLTYP=grid%isltyp      ,ITIMESTEP=grid%itimestep, JULIAN_IN=grid%julian                  &
     &        ,IVGTYP=grid%ivgtyp      ,LH=grid%lh              ,LOWLYR=grid%lowlyr      &
     &        ,MAVAIL=grid%mavail      ,NUM_SOIL_LAYERS=config_flags%num_soil_layers        &
     &        ,P8W=grid%p_hyd_w            ,PBLH=grid%pblh          ,PI_PHY=pi_phy      &
     &        ,PSFC=grid%psfc          ,PSHLTR=grid%pshltr      ,PSIH=psih          &
     &        ,BLDT=grid%bldt     ,CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag  &
     &        ,BLDTACTTIME=grid%bldtacttime                                              & 
     &        ,PSIM=psim          ,P_PHY=grid%p_hyd        ,Q10=grid%q10            &
     &        ,Q2=grid%q2              ,QFX=grid%qfx            ,QSFC=grid%qsfc          &
     &        ,QSHLTR=grid%qshltr      ,QZ0=grid%qz0            ,RAINCV=grid%raincv      &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics         ,RHO=grid%rho            &
     &        ,RMOL=grid%rmol          ,SFCEVP=grid%sfcevp      ,SFCEXC=grid%sfcexc      &
     &        ,SFCRUNOFF=grid%sfcrunoff                                        &
     &        ,opt_thcnd=config_flags%opt_thcnd                                                   &
     &        ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics                        &
     &        ,SF_SURFACE_PHYSICS=config_flags%sf_surface_physics  ,SH2O=grid%sh2o          &
     &        ,SHDMAX=grid%shdmax      ,SHDMIN=grid%shdmin      ,SMOIS=grid%smois        &
     &        ,SMSTAV=grid%smstav      ,SMSTOT=grid%smstot      ,SNOALB=grid%snoalb      &
     &        ,SNOW=grid%snow          ,SNOWC=grid%snowc        ,SNOWH=grid%snowh        &
     &        ,SMCREL=grid%smcrel                                                        &
     &        ,SST=grid%sst            ,SST_INPUT=grid%sst_input,SST_UPDATE=grid%sst_update                  &
     &        ,SSTSK=grid%sstsk        ,DTW=grid%dtw            ,SST_SKIN=grid%sst_skin  &
     &        ,SCM_FORCE_SKINTEMP=grid%scm_force_skintemp                                &
     &        ,SCM_FORCE_FLUX=grid%scm_force_flux                                        &
     &        ,STEPBL=grid%stepbl      ,TH10=grid%th10          ,TH2=grid%th2            &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy      ,TKE_PBL=grid%tke_pbl    &
     &        ,TMN=grid%tmn            ,TSHLTR=grid%tshltr      ,TSK=grid%tsk            &
     &        ,TYR=grid%tyr            ,TYRA=grid%tyra          ,TDLY=grid%tdly          &
     &        ,TLAG=grid%tlag          ,LAGDAY=config_flags%lagday      ,NYEAR=grid%nyear        &
     &        ,NDAY=grid%nday          ,TMN_UPDATE=grid%tmn_update      ,YR=yr            &
     &        ,TSLB=grid%tslb          ,T_PHY=grid%t_phy        ,U10=grid%u10            &
     &        ,URATX=grid%uratx        ,VRATX=grid%vratx   ,TRATX=grid%tratx        &
     &        ,UDRUNOFF=grid%udrunoff  ,UST=grid%ust       ,UZ0=grid%uz0            &
     &        ,U_FRAME=grid%u_frame    ,U_PHY=grid%u_phy   ,V10=grid%v10            &
     &        ,UOCE=grid%uoce          ,VOCE=grid%voce                              &
     &        ,VEGFRA=grid%vegfra      ,VZ0=grid%vz0       ,V_FRAME=grid%v_frame    &
     &        ,V_PHY=grid%v_phy        ,WARM_RAIN=grid%warm_rain                    &
     &        ,WSPD=wspd               ,XICE=grid%xice     ,XLAND=grid%xland        &
     &        ,MAX_EDOM=grid%num_ext_model_couple_dom      ,CPLMASK=grid%cplmask    &
     &        ,Z0=grid%z0              ,Z=grid%z        ,ZNT=grid%znt               &
     &        ,ZS=grid%zs              ,ALBSI=grid%albsi , ICEDEPTH=grid%icedepth   &
     &        ,SNOWSI=grid%snowsi                                                   &
     &        ,XICEM=grid%xicem   ,ISICE=grid%landuse_isice                         &
     &        ,USTM=grid%ustm          ,CK=grid%ck         ,CKA=grid%cka            &
     &                                 ,CD=grid%cd         ,CDA=grid%cda            &
     &        ,ISFTCFLX=config_flags%isftcflx, IZ0TLND=config_flags%iz0tlnd         &
     &        ,SF_OCEAN_PHYSICS=config_flags%sf_ocean_physics                       &
     &        ,OML_HML0=config_flags%oml_hml0 ,OML_GAMMA=config_flags%oml_gamma     &
     &        ,TML=grid%tml, T0ML=grid%t0ml, HML=grid%hml, H0ML=grid%h0ml           &
     &        ,HUML=grid%huml, HVML=grid%hvml, F=grid%f                             &
     &        ,TMOML=grid%TMOML,ISWATER=iswater                                     &
     &        ,OML_RELAXATION_TIME=grid%OML_RELAXATION_TIME                         &
     &        ,lakedepth2d=grid%lakedepth2d,    savedtke12d=grid%savedtke12d        &   
     &        ,snowdp2d=grid%snowdp2d,          h2osno2d=grid%h2osno2d              & !lake
     &        ,snl2d=grid%snl2d,                t_grnd2d=grid%t_grnd2d              &   
     &        ,t_lake3d=grid%t_lake3d,          lake_icefrac3d=grid%lake_icefrac3d  & !lake  
     &        ,z_lake3d=grid%z_lake3d,          dz_lake3d=grid%dz_lake3d            &   
     &        ,t_soisno3d=grid%t_soisno3d,      h2osoi_ice3d=grid%h2osoi_ice3d      & !lake
     &        ,h2osoi_liq3d=grid%h2osoi_liq3d,  h2osoi_vol3d=grid%h2osoi_vol3d      &   
     &        ,z3d=grid%z3d,                    dz3d=grid%dz3d                      & !lake
     &        ,zi3d=grid%zi3d,                  watsat3d=grid%watsat3d              &   
     &        ,csol3d=grid%csol3d,              tkmg3d=grid%tkmg3d                  & !lake     
     &        ,tkdry3d=grid%tkdry3d,            tksatu3d=grid%tksatu3d              &   
     &        ,LakeModel=grid%sf_lake_physics,  lake_min_elev=grid%lake_min_elev    & !lake  
#if ( EM_CORE == 1)
     &        ,LakeMask=grid%LakeMask                                               & ! lake
#endif
! CLM Varaibles
     &        ,NUMC=grid%numc,NUMP=grid%nump,SABV=grid%sabv,SABG=grid%sabg,         &
     &         LWUP=grid%lwup,SNL=grid%snl,                                         &
     &         HISTORY_INTERVAL=config_flags%history_interval ,                     &!ylu add hist inverval for accumulation T max/min
     &         SNOWDP=grid%snowdp, WTC=grid%wtc,WTP=grid%wtp, H2OSNO=grid%h2osno,   &
     &         T_GRND=grid%t_grnd,T_VEG=grid%t_veg,                                 &
     &         H2OCAN=grid%h2ocan, H2OCAN_COL=grid%h2ocan_col,T2M_MAX=grid%t2m_max, &
     &         T2M_MIN=grid%t2m_min,T2CLM=grid%t2clm,                               &
     &         T_REF2M=grid%t_ref2m,H2OSOI_LIQ_S1=grid%h2osoi_liq_s1,               &
     &         H2OSOI_LIQ_S2=grid%h2osoi_liq_s2,                                    &
     &         H2OSOI_LIQ_S3=grid%h2osoi_liq_s3,H2OSOI_LIQ_S4=grid%h2osoi_liq_s4,   &
     &         H2OSOI_LIQ_S5=grid%h2osoi_liq_s5,                                    &
     &         H2OSOI_LIQ1=grid%h2osoi_liq1,H2OSOI_LIQ2=grid%h2osoi_liq2,           &
     &         H2OSOI_LIQ3=grid%h2osoi_liq3,H2OSOI_LIQ4=grid%h2osoi_liq4,           &
     &         H2OSOI_LIQ5=grid%h2osoi_liq5,H2OSOI_LIQ6=grid%h2osoi_liq6,           &
     &         H2OSOI_LIQ7=grid%h2osoi_liq7,H2OSOI_LIQ8=grid%h2osoi_liq8,           &
     &         H2OSOI_LIQ9=grid%h2osoi_liq9, H2OSOI_LIQ10=grid%h2osoi_liq10,        &
     &         H2OSOI_ICE_S1=grid%h2osoi_ice_s1,H2OSOI_ICE_S2=grid%h2osoi_ice_s2,   &
     &         H2OSOI_ICE_S3=grid%h2osoi_ice_s3, H2OSOI_ICE_S4=grid%h2osoi_ice_s4,  &
     &         H2OSOI_ICE_S5=grid%h2osoi_ice_s5,                                    &
     &         H2OSOI_ICE1=grid%h2osoi_ice1, H2OSOI_ICE2=grid%h2osoi_ice2,          &
     &         H2OSOI_ICE3=grid%h2osoi_ice3,H2OSOI_ICE4=grid%h2osoi_ice4,           &
     &         H2OSOI_ICE5=grid%h2osoi_ice5, H2OSOI_ICE6=grid%h2osoi_ice6,          &
     &         H2OSOI_ICE7=grid%h2osoi_ice7,H2OSOI_ICE8=grid%h2osoi_ice8,           &
     &         H2OSOI_ICE9=grid%h2osoi_ice9,H2OSOI_ICE10=grid%h2osoi_ice10,         &
     &         T_SOISNO_S1=grid%t_soisno_s1,T_SOISNO_S2=grid%t_soisno_s2,           &
     &         T_SOISNO_S3=grid%t_soisno_s3,T_SOISNO_S4=grid%t_soisno_s4,           &
     &         T_SOISNO_S5=grid%t_soisno_s5,T_SOISNO1=grid%t_soisno1,               &
     &         T_SOISNO2=grid%t_soisno2,T_SOISNO3=grid%t_soisno3,                   &
     &         T_SOISNO4=grid%t_soisno4,T_SOISNO5=grid%t_soisno5,                   &
     &         T_SOISNO6=grid%t_soisno6,T_SOISNO7=grid%t_soisno7,                   &
     &         T_SOISNO8=grid%t_soisno8,T_SOISNO9=grid%t_soisno9,                   &
     &         T_SOISNO10=grid%t_soisno10,DZSNOW1=grid%dzsnow1,DZSNOW2=grid%dzsnow2,& 
     &         DZSNOW3=grid%dzsnow3,DZSNOW4=grid%dzsnow4, DZSNOW5=grid%dzsnow5,     &
     &         SNOWRDS1=grid%snowrds1,SNOWRDS2=grid%snowrds2,                       &
     &         SNOWRDS3=grid%snowrds3 ,SNOWRDS4=grid%snowrds4,                      &
     &         SNOWRDS5=grid%snowrds5,                                              &
     &         T_LAKE1=grid%t_lake1,T_LAKE2=grid%t_lake2,T_LAKE3=grid%t_lake3,      &
     &         T_LAKE4=grid%t_lake4,                                                &
     &         T_LAKE5=grid%t_lake5,T_LAKE6=grid%t_lake6, T_LAKE7=grid%t_lake7,     &
     &         T_LAKE8=grid%t_lake8, T_LAKE9=grid%t_lake9,T_LAKE10=grid%t_lake10,   &
     &         H2OSOI_VOL1=grid%h2osoi_vol1,H2OSOI_VOL2=grid%h2osoi_vol2,           &
     &         H2OSOI_VOL3=grid%h2osoi_vol3,H2OSOI_VOL4=grid%h2osoi_vol4,           &
     &         H2OSOI_VOL5=grid%h2osoi_vol5,                                        &
     &         H2OSOI_VOL6=grid%h2osoi_vol6,H2OSOI_VOL7=grid%h2osoi_vol7,           &
     &         H2OSOI_VOL8=grid%h2osoi_vol8,                                        &
     &         H2OSOI_VOL9=grid%h2osoi_vol9,H2OSOI_VOL10=grid%h2osoi_vol10,         &
     &         MAXPATCH=config_flags%maxpatch,                                      &
     &         INEST=grid%id,ALBEDOsubgrid=grid%ALBEDOsubgrid,                      &
     &         LHsubgrid=grid%LHsubgrid,                                            &
     &         HFXsubgrid=grid%HFXsubgrid,LWUPsubgrid=grid%LWUPsubgrid,             &
     &         Q2subgrid=grid%Q2subgrid,SABVsubgrid=grid%SABVsubgrid,               &
     &         SABGsubgrid=grid%SABGsubgrid,NRAsubgrid=grid%NRAsubgrid,             &
     &         SWUPsubgrid=grid%SWUPsubgrid,LHsoi=grid%LHsoi,                       &
     &         LHveg=grid%LHveg, LHtran=grid%LHtran                                 &
! end of CLM variables
     &        ,SLOPE_RAD=config_flags%slope_rad,TOPO_SHADING=config_flags%topo_shading & ! solar
     &        ,SHADOWMASK=grid%shadowmask,DIFFUSE_FRAC=grid%diffuse_frac               & ! solar
     &        ,SLOPE=grid%slope, SLP_AZI=grid%slp_azi, SWNORM=grid%swnorm              & ! solar
     &        ,DECLIN=grid%declin ,SOLCON=grid%solcon ,COSZEN=grid%coszen ,HRANG=grid%hrang    &
     &        ,xlat_urb2d=grid%XLAT                                       & !I urban
     &        ,NUM_ROOF_LAYERS=num_roof_layers                            & !I urban
     &        ,NUM_WALL_LAYERS=num_wall_layers                            & !I urban
     &        ,NUM_ROAD_LAYERS=num_road_layers                            &
     &        ,DZR=grid%dzr ,DZB=grid%dzb ,DZG=grid%dzg                   & !I urban
     &        ,TR_URB2D=grid%tr_urb2d ,TB_URB2D=grid%tb_urb2d             &
     &        ,TG_URB2D=grid%tg_urb2d                                     & !H urban
     &        ,TC_URB2D=grid%tc_urb2d ,QC_URB2D=grid%qc_urb2d             & !H urban
     &        ,UC_URB2D=grid%uc_urb2d                                     & !H urban
     &        ,XXXR_URB2D=grid%xxxr_urb2d                                 &
     &        ,XXXB_URB2D=grid%xxxb_urb2d                                 & !H urban
     &        ,XXXG_URB2D=grid%xxxg_urb2d                                 &
     &        ,XXXC_URB2D=grid%xxxc_urb2d                                 & !H urban
     &        ,CMCR_URB2D=grid%cmcr_urb2d,TGR_URB2D=grid%tgr_urb2d        & !H urban
     &        ,TGRL_URB3D=grid%tgrl_urb3d,SMR_URB3D=grid%smr_urb3d        & !H urban
     &        ,JULIAN=grid%julday, JULYR=grid%julyr                       & !I urban
     &        ,DRELR_URB2D=grid%drelr_urb2d,DRELB_URB2D=grid%drelb_urb2d  & !H urban
     &        ,DRELG_URB2D=grid%drelg_urb2d                               & !H urban
     &        ,FLXHUMR_URB2D=grid%flxhumr_urb2d                           & !H urban
     &        ,FLXHUMB_URB2D=grid%flxhumb_urb2d                           & !H urban
     &        ,FLXHUMG_URB2D=grid%flxhumg_urb2d                           & !H urban
     &        ,TRL_URB3D=grid%trl_urb3d   ,TBL_URB3D=grid%tbl_urb3d       & !H urban
     &        ,TGL_URB3D=grid%tgl_urb3d                                   & !H urban
     &        ,SH_URB2D=grid%sh_urb2d     ,LH_URB2D=grid%lh_urb2d         &
     &        ,G_URB2D=grid%g_urb2d                                       & !H urban
     &        ,RN_URB2D=grid%rn_urb2d     , TS_URB2D=grid%ts_urb2d        & !H urban
     &        ,FRC_URB2D=grid%frc_urb2d                                   & !H urban
     &        ,UTYPE_URB2D=grid%utype_urb2d                               & !H urban
         ! Optional urban for BEP scheme
     &        ,SF_URBAN_PHYSICS=config_flags%sf_urban_physics             &
     &        ,NUM_URBAN_LAYERS=config_flags%num_urban_layers             & !multi-layer urban
     &        ,NUM_URBAN_HI=config_flags%num_urban_hi                     & !multi-layer urban
     &        ,TRB_URB4D=grid%trb_urb4d,TW1_URB4D=grid%tw1_urb4d          & !multi-layer urban
     &        ,TW2_URB4D=grid%tw2_urb4d,TGB_URB4D=grid%tgb_urb4d          & !multi-layer urban
     &        ,TLEV_URB3D=grid%tlev_urb3d                               & !multi-layer urban
     &        ,QLEV_URB3D=grid%qlev_urb3d                               & !multi-layer urban
     &        ,TW1LEV_URB3D=grid%tw1lev_urb3d                           & !multi-layer urban
     &        ,TW2LEV_URB3D=grid%tw2lev_urb3d                           & !multi-layer urban
     &        ,TGLEV_URB3D=grid%tglev_urb3d                             & !multi-layer urban
     &        ,TFLEV_URB3D=grid%tflev_urb3d                             & !multi-layer urban
     &        ,SF_AC_URB3D=grid%sf_ac_urb3d                             & !multi-layer urban
     &        ,LF_AC_URB3D=grid%lf_ac_urb3d                             & !multi-layer urban
     &        ,CM_AC_URB3D=grid%cm_ac_urb3d                             & !multi-layer urban
     &        ,SFVENT_URB3D=grid%sfvent_urb3d                           & !multi-layer urban
     &        ,LFVENT_URB3D=grid%lfvent_urb3d                           & !multi-layer urban
     &        ,SFWIN1_URB3D=grid%sfwin1_urb3d                           & !multi-layer urban
     &        ,SFWIN2_URB3D=grid%sfwin2_urb3d                           & !multi-layer urban
     &        ,SFW1_URB3D=grid%sfw1_urb3d,SFW2_URB3D=grid%sfw2_urb3d      & !multi-layer urban
     &        ,SFR_URB3D=grid%sfr_urb3d,SFG_URB3D=grid%sfg_urb3d          & !multi-layer urban
     &        ,LP_URB2D=grid%lp_urb2d,HI_URB2D=grid%hi_urb2d              & !multi-layer urban
     &        ,LB_URB2D=grid%lb_urb2d,HGT_URB2D=grid%hgt_urb2d            & !multi-layer urban
     &        ,MH_URB2D=grid%mh_urb2d,STDH_URB2D=grid%stdh_urb2d          & !SLUCM
     &        ,LF_URB2D=grid%lf_urb2d                                           & 
     &        ,GMT=grid%gmt,XLAT=grid%xlat,XLONG=grid%xlong,JULDAY=grid%julday  &
     &        ,A_U_BEP=grid%a_u_bep,A_V_BEP=grid%a_v_bep,A_T_BEP=grid%a_t_bep   &
     &        ,A_Q_BEP=grid%a_q_bep                                             &
     &        ,B_U_BEP=grid%b_u_bep,B_V_BEP=grid%b_v_bep,B_T_BEP=grid%b_t_bep   &
     &        ,B_Q_BEP=grid%b_q_bep                                             &
     &        ,SF_BEP=grid%sf_bep,VL_BEP=grid%vl_bep                            &
     &        ,A_E_BEP=grid%a_e_bep,B_E_BEP=grid%b_e_bep,DLG_BEP=grid%dlg_bep   &
     &        ,DL_U_BEP=grid%dl_u_bep                                           &
     &        ,CMR_SFCDIF=grid%cmr_sfcdif, CHR_SFCDIF=grid%chr_sfcdif     & !I/O urban
     &        ,CMC_SFCDIF=grid%cmc_sfcdif, CHC_SFCDIF=grid%chc_sfcdif     & !I/O urban
     &        ,CMGR_SFCDIF=grid%cmgr_sfcdif, CHGR_SFCDIF=grid%chgr_sfcdif & !I/O urban
           ! P-X LSM Variables
     &        ,LANDUSEF=grid%landusef, SOILCTOP=grid%soilctop             &   ! P-X LSM
     &        ,SOILCBOT=grid%soilcbot                                     &   ! P-X LSM
     &        ,RA=grid%ra, RS=grid%rs, LAI=grid%lai, IMPERV=grid%imperv   &   ! P-X LSM
     &        ,CANFRA=grid%canfra, NLCAT=grid%num_land_cat                &   ! P-X LSM
     &        ,NSCAT=grid%num_soil_cat                                    &   ! P-X LSM
     &        ,VEGF_PX=grid%vegf_px, SNOWNCV=grid%snowncv                 &   ! P-X LSM
     &        ,ANAL_INTERVAL=config_flags%auxinput9_interval_s+config_flags%auxinput9_interval_m*60  &   ! P-X LSM
     &        ,PXLSM_SMOIS_INIT=config_flags%pxlsm_smois_init             &   ! P-X LSM
     &        ,PXLSM_SOIL_NUDGE=config_flags%pxlsm_soil_nudge             &   ! P-X LSM
           ! SSiB LSM variables (fds 06/2010)
     &        ,alswvisdir=grid%alswvisdir, alswvisdif=grid%alswvisdif    & !ssib
     &        ,alswnirdir=grid%alswnirdir, alswnirdif=grid%alswnirdif    & !ssib
     &        ,swvisdir=grid%swvisdir, swvisdif=grid%swvisdif            & !ssib
     &        ,swnirdir=grid%swnirdir, swnirdif=grid%swnirdif            & !ssib 
     &        ,ssib_br=grid%ssib_br,  ssib_fm=grid%ssib_fm               & !ssib
     &        ,ssib_fh=grid%ssib_fh,  ssib_cm=grid%ssib_cm               & !ssib
     &        ,ssibxdd=grid%ssibxdd,  ssib_lhf=grid%ssib_lhf             & !ssib
     &        ,ssib_shf=grid%ssib_shf, ssib_ghf=grid%ssib_ghf            & !ssib
     &        ,ssib_egs=grid%ssib_egs, ssib_eci=grid%ssib_eci            & !ssib
     &        ,ssib_ect=grid%ssib_ect, ssib_egi=grid%ssib_egi            & !ssib
     &        ,ssib_egt=grid%ssib_egt, ssib_sdn=grid%ssib_sdn            & !ssib
     &        ,ssib_sup=grid%ssib_sup, ssib_ldn=grid%ssib_ldn            & !ssib
     &        ,ssib_lup=grid%ssib_lup, ssib_wat=grid%ssib_wat            & !ssib
     &        ,ssib_shc=grid%ssib_shc, ssib_shg=grid%ssib_shg            & !ssib
     &        ,ssib_lai=grid%ssib_lai, ssib_vcf=grid%ssib_vcf            & !ssib
     &        ,ssib_z00=grid%ssib_z00, ssib_veg=grid%ssib_veg            & !ssib
     &        ,cldfra=grid%cldfra                                        & !ssib
     &        ,ISNOW=grid%isnow, SWE=grid%swe, SNOWDEN=grid%snowden      & !ssib snow
     &        ,SNOWDEPTH=grid%snowdepth, TKAIR=grid%tkair                & !ssib snow
     &        ,DZO1=grid%dzo1, WO1=grid%wo1, TSSN1=grid%tssn1, TSSNO1=grid%tssno1     & !ssib snow
     &        ,BWO1=grid%bwo1, BTO1=grid%bto1, CTO1=grid%cto1, FIO1=grid%fio1         & !ssib snow
     &        ,FLO1=grid%flo1, BIO1=grid%bio1, BLO1=grid%blo1, HO1=grid%ho1           & !ssib snow
     &        ,DZO2=grid%dzo2, WO2=grid%wo2, TSSN2=grid%tssn2, TSSNO2=grid%tssno2     & !ssib snow
     &        ,BWO2=grid%bwo2, BTO2=grid%bto2, CTO2=grid%cto2, FIO2=grid%fio2         & !ssib snow
     &        ,FLO2=grid%flo2, BIO2=grid%bio2, BLO2=grid%blo2, HO2=grid%ho2           & !ssib snow
     &        ,DZO3=grid%dzo3, WO3=grid%wo3, TSSN3=grid%tssn3, TSSNO3=grid%tssno3     & !ssib snow
     &        ,BWO3=grid%bwo3, BTO3=grid%bto3, CTO3=grid%cto3, FIO3=grid%fio3         & !ssib snow
     &        ,FLO3=grid%flo3, BIO3=grid%bio3, BLO3=grid%blo3, HO3=grid%ho3           & !ssib snow
     &        ,DZO4=grid%dzo4, WO4=grid%wo4, TSSN4=grid%tssn4, TSSNO4=grid%tssno4     & !ssib snow
     &        ,BWO4=grid%bwo4, BTO4=grid%bto4, CTO4=grid%cto4, FIO4=grid%fio4         & !ssib snow
     &        ,FLO4=grid%flo4, BIO4=grid%bio4, BLO4=grid%blo4, HO4=grid%ho4           & !ssib snow
     &        ,RA_SW_PHYSICS=config_flags%ra_sw_physics                  & !ssib
!------------------------------------------------------------------------------
           ! Optional PX LSM nudging
     &        ,t2_ndg_old=grid%t2_ndg_old                  &
     &        ,q2_ndg_old=grid%q2_ndg_old                  &
     &        ,t2_ndg_new=grid%t2_ndg_new                  &
     &        ,q2_ndg_new=grid%q2_ndg_new                  &
     &        ,sn_ndg_old=grid%sn_ndg_old                  &
     &        ,sn_ndg_new=grid%sn_ndg_new                  &
! for Noah-MP LSM
     &        ,idveg=config_flags%dveg,         iopt_crs=config_flags%opt_crs &
     &        ,iopt_btr=config_flags%opt_btr,   iopt_run=config_flags%opt_run &
     &        ,iopt_sfc=config_flags%opt_sfc,   iopt_frz=config_flags%opt_frz &
     &        ,iopt_inf=config_flags%opt_inf,   iopt_rad=config_flags%opt_rad &
     &        ,iopt_alb=config_flags%opt_alb,   iopt_snf=config_flags%opt_snf &
     &        ,iopt_tbot=config_flags%opt_tbot, iopt_stc=config_flags%opt_stc & 
     &        ,iopt_gla=config_flags%opt_gla,   iopt_rsf=config_flags%opt_rsf & 
     &        , isnowxy=grid%isnowxy  ,    tvxy=grid%tvxy    ,    tgxy=grid%tgxy     &
     &        ,canicexy=grid%canicexy ,canliqxy=grid%canliqxy,   eahxy=grid%eahxy    &
     &        ,   tahxy=grid%tahxy    ,    cmxy=grid%cmxy    ,    chxy=grid%chxy     &
     &        ,  fwetxy=grid%fwetxy   ,sneqvoxy=grid%sneqvoxy,alboldxy=grid%alboldxy &
     &        , qsnowxy=grid%qsnowxy  ,wslakexy=grid%wslakexy,   zwtxy=grid%zwtxy    &
     &        ,    waxy=grid%waxy     ,    wtxy=grid%wtxy    ,  tsnoxy=grid%tsnoxy   &
     &        , zsnsoxy=grid%zsnsoxy  , snicexy=grid%snicexy , snliqxy=grid%snliqxy  &
     &        ,lfmassxy=grid%lfmassxy ,rtmassxy=grid%rtmassxy,stmassxy=grid%stmassxy &
     &        ,  woodxy=grid%woodxy   ,stblcpxy=grid%stblcpxy,fastcpxy=grid%fastcpxy &
     &        , grainxy=grid%grainxy  ,   gddxy=grid%gddxy                           &
     &        ,  xsaixy=grid%xsaixy   , taussxy=grid%taussxy                         &
     &        ,  t2mvxy=grid%t2mvxy   ,  t2mbxy=grid%t2mbxy                          &
     &        ,  q2mvxy=grid%q2mvxy   ,  q2mbxy=grid%q2mbxy                          & 
     &        ,  tradxy=grid%tradxy   ,   neexy=grid%neexy   ,   gppxy=grid%gppxy    &
     &        ,   nppxy=grid%nppxy    ,  fvegxy=grid%fvegxy  , runsfxy=grid%runsfxy  &
     &        , runsbxy=grid%runsbxy  ,  ecanxy=grid%ecanxy  ,  edirxy=grid%edirxy   &
     &        , etranxy=grid%etranxy  ,   fsaxy=grid%fsaxy   ,  firaxy=grid%firaxy   &
     &        ,  aparxy=grid%aparxy   ,   psnxy=grid%psnxy   ,   savxy=grid%savxy    &
     &        ,   sagxy=grid%sagxy    , rssunxy=grid%rssunxy , rsshaxy=grid%rsshaxy  &
     &        ,  bgapxy=grid%bgapxy   ,  wgapxy=grid%wgapxy  ,   tgvxy=grid%tgvxy    &
     &        ,   tgbxy=grid%tgbxy    ,   chvxy=grid%chvxy   ,   chbxy=grid%chbxy    &
     &        ,   shgxy=grid%shgxy    ,   shcxy=grid%shcxy   ,   shbxy=grid%shbxy    &
     &        ,   evgxy=grid%evgxy    ,   evbxy=grid%evbxy   ,   ghvxy=grid%ghvxy    &
     &        ,   ghbxy=grid%ghbxy    ,   irgxy=grid%irgxy   ,   ircxy=grid%ircxy    &
     &        ,   irbxy=grid%irbxy    ,    trxy=grid%trxy    ,   evcxy=grid%evcxy    &
     &        ,chleafxy=grid%chleafxy ,  chucxy=grid%chucxy                          &                       
     &        ,  chv2xy=grid%chv2xy   ,  chb2xy=grid%chb2xy , chstarxy=grid%chstarxy &                       
           !Optional hydro variables in NOAHMP
     &        ,smcwtdxy=grid%smcwtdxy ,rechxy=grid%rechxy ,deeprechxy=grid%deeprechxy &
     &        ,fdepthxy=grid%fdepthxy, areaxy=grid%areaxy, rivercondxy=grid%rivercondxy  &
     &        ,riverbedxy=grid%riverbedxy, eqzwt=grid%eqzwt, pexpxy=grid%pexpxy &
     &        ,qrfxy=grid%qrfxy, qspringxy=grid%qspringxy, qslatxy=grid%qslatxy, qrfsxy=grid%qrfsxy &
     &        ,qspringsxy=grid%qspringsxy, smoiseq=grid%smoiseq, wtddt=config_flags%wtddt, stepwtd=grid%stepwtd &
           ! Noah UA changes
     &        ,ua_phys=config_flags%ua_phys,flx4=grid%flx4,fvb=grid%fvb       &
     &        ,fbur=grid%fbur,fgsn=grid%fgsn                                  &
           ! Indexes
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe          &
     &        , I_START=grid%i_start,I_END=min(grid%i_end, ide-1)         &
     &        , J_START=grid%j_start,J_END=min(grid%j_end, jde-1)         &
     &        , KTS=k_start, KTE=min(k_end,kde-1)                         &
     &        , NUM_TILES=grid%num_tiles                                  &
          ! Variables required by TEMF PBL - WA 1/7/09
              ,te_temf=grid%te_temf,hd_temf=grid%hd_temf                  &
              ,fCor=grid%f,exch_temf=exch_temf,wm_temf=grid%wm_temf       &
          ! Variables required by IDEAL SCM sfc scheme - WA 1/6/10
              ,hfx_force=grid%hfx_force,lh_force=grid%lh_force            &
              ,tsk_force=grid%tsk_force                                   &
              ,hfx_force_tend=grid%hfx_force_tend                         &
              ,lh_force_tend=grid%lh_force_tend                           &
              ,tsk_force_tend=grid%tsk_force_tend                         &
           ! Optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV                 &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC                 &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR                 &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI                 &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,CAPG=grid%capg, EMISS=grid%emiss, HOL=hol,MOL=grid%mol     &
     &        ,T2OBS=grid%t2obs, Q2OBS=grid%q2obs                         &
     &        ,RAINBL=grid%rainbl,SR=grid%sr,RAINSHV=grid%rainshv         &
     &        ,GRAUPELNCV=grid%graupelncv, HAILNCV=grid%hailncv           &
     &        ,RAINNCV=grid%rainncv,REGIME=grid%regime,T2=grid%t2,THC=grid%thc            &
     &        ,QSG=grid%qsg,QVG=grid%qvg,QCG=grid%qcg,SOILT1=grid%soilt1,TSNAV=grid%tsnav & ! ruc lsm
     &        ,SMFR3D=grid%smfr3d,KEEPFR3DFLAG=grid%keepfr3dflag,DEW=grid%dew             & ! ruc lsm
     &        ,POTEVP=grid%POTEVP, SNOPCX=grid%SNOPCX, SOILTB=grid%SOILTB                 & ! ruc lsm
     &        ,rhosnf=grid%rhosnf ,precipfr=grid%precipfr                 &   ! RUC LSM
     &        ,snowfallac=grid%snowfallac                                 &   ! RUC LSM
     &        ,MOSAIC_LU=config_flags%mosaic_lu                           & ! RUC LSM
     &        ,MOSAIC_SOIL=config_flags%mosaic_soil                       &   ! RUC LSM
     &        ,ISURBAN=isurban, MMINLU=TRIM(mminlu)                       &
     &        ,SNOTIME = grid%SNOTIME                                     &
     &        ,RDLAI2D=config_flags%rdlai2d                               &
     &        ,usemonalb=config_flags%usemonalb                           &
     &        ,NOAHRES=grid%noahres                                       &
     &        ,TSK_SAVE=grid%tsk_save                                     &
     &        ,ch=grid%ch,tsq=grid%tsq,qsq=grid%qsq,cov=grid%cov          & !MYNN - MP
     &        ,Sh3d=grid%sh3d,EL_PBL=grid%el_pbl                          & !JOE- MYNN cloudpdf
     &        ,bl_mynn_cloudpdf=config_flags%bl_mynn_cloudpdf             & !JOE- MYNN cloudpdf
     &        ,icloud_bl=config_flags%icloud_bl                           & !JOE- subgrid cloud
     &        ,qc_bl=grid%qc_bl,cldfra_bl=grid%cldfra_bl                  & !JOE- subgrid cloud
     &        ,fgdp=grid%fgdp,dfgdp=grid%dfgdp,vdfg=grid%vdfg             & !Katata - fogdes
     &        ,grav_settling=config_flags%grav_settling                   & !Katata - fogdes
     &        ,OM_TMP=grid%om_tmp, OM_S=grid%om_s, OM_U=grid%om_u, OM_V=grid%om_v   & !cyl:3DPWP
     &        ,OM_DEPTH=grid%om_depth, OM_ML=grid%OM_ML, OM_LON=grid%om_lon         & !cyl:3DPWP
     &        ,OM_LAT=grid%om_lat                                                  & !cy:3DPWP
     &        , okms = 1, okme=config_flags%ocean_levels                  & ! cyl:3DPWP
     &        ,rdx=grid%rdx, rdy=grid%rdy,msfu=grid%msfu,msfv=grid%msfv,msft=grid%msft &!cyl: 3DPWP
     &        ,XTIME=grid%xtime,OM_TINI=grid%om_tini,OM_SINI=grid%om_sini,id=grid%id,omdt=config_flags%omdt &!cyl: 3DPWP 
     &        ,sf_surface_mosaic=config_flags%sf_surface_mosaic,  mosaic_cat=config_flags%mosaic_cat    &
     &        ,mosaic_cat_index=grid%mosaic_cat_index                                                   &  !danli mosaic
     &        ,landusef2=grid%landusef2,TSK_mosaic=grid%TSK_mosaic,QSFC_mosaic=grid%QSFC_mosaic         &
     &        ,TSLB_mosaic=grid%TSLB_mosaic,SMOIS_mosaic=grid%SMOIS_mosaic,SH2O_mosaic=grid%SH2O_mosaic &  !danli mosaic
     &        ,CANWAT_mosaic=grid%CANWAT_mosaic,SNOW_mosaic=grid%SNOW_mosaic                            &
     &        ,SNOWH_mosaic=grid%SNOWH_mosaic,SNOWC_mosaic=grid%SNOWC_mosaic                            &  !danli mosaic
     &        ,ALBEDO_mosaic=grid%ALBEDO_mosaic,ALBBCK_mosaic=grid%ALBBCK_mosaic                        &
     &        ,EMISS_mosaic=grid%EMISS_mosaic, EMBCK_mosaic=grid%EMBCK_mosaic                           &
     &        ,ZNT_mosaic=grid%ZNT_mosaic, Z0_mosaic=grid%Z0_mosaic                                     &  !danli mosaic
     &        ,HFX_mosaic=grid%HFX_mosaic,QFX_mosaic=grid%QFX_mosaic, LH_mosaic=grid%LH_mosaic          &
     &        ,GRDFLX_mosaic=grid%GRDFLX_mosaic,SNOTIME_mosaic=grid%SNOTIME_mosaic                      &  !danli mosaic
     &        ,TR_URB2D_mosaic=grid%TR_URB2D_mosaic,TB_URB2D_mosaic=grid%TB_URB2D_mosaic                &  !danli mosaic 
     &        ,TG_URB2D_mosaic=grid%TG_URB2D_mosaic,TC_URB2D_mosaic=grid%TC_URB2D_mosaic                &  !danli mosaic 
     &        ,QC_URB2D_mosaic=grid%QC_URB2D_mosaic,UC_URB2D_mosaic=grid%UC_URB2D_mosaic                &  !danli mosaic
     &        ,TRL_URB3D_mosaic=grid%TRL_URB3D_mosaic,TBL_URB3D_mosaic=grid%TBL_URB3D_mosaic            &  !danli mosaic 
     &        ,TGL_URB3D_mosaic=grid%TGL_URB3D_mosaic                                                   &  !danli mosaic 
     &        ,SH_URB2D_mosaic=grid%SH_URB2D_mosaic,LH_URB2D_mosaic=grid%LH_URB2D_mosaic                &  !danli mosaic 
     &        ,G_URB2D_mosaic=grid%G_URB2D_mosaic,RN_URB2D_mosaic=grid%RN_URB2D_mosaic                  &  !danli mosaic 
     &        ,TS_URB2D_mosaic=grid%TS_URB2D_mosaic                                                     &  !danli mosaic 
     &        ,TS_RUL2D_mosaic=grid%TS_RUL2D_mosaic                                                     &  !danli mosaic
     &        ,ZOL=grid%ZOL                                                                             &
     &        ,SDA_HFX=grid%SDA_HFX, SDA_QFX=grid%SDA_QFX,HFX_BOTH=grid%HFX_BOTH                        &  !fasdas
     &        ,QFX_BOTH=grid%QFX_BOTH,QNORM=grid%QNORM,fasdas=config_flags%fasdas                       &  !fasdas
     &                                                           )

#ifdef WRF_HYDRO
        if(HYDRO_dt .gt. 1 ) call wrf_drv_HYDRO(HYDRO_dt, grid,  &
     &     grid%i_start(1),min(grid%i_end(1), ide-1),      &
     &     grid%j_start(1),min(grid%j_end(1), jde-1) )
#endif

BENCH_END(surf_driver_tim)

!*********
! pbl
      CALL wrf_debug ( 200 , ' call pbl_driver' )
BENCH_START(pbl_driver_tim)
      CALL PUSHREAL8ARRAY(grid%pblh, (ime-ims+1)*(jme-jms+1))

#if ( WRF_CHEM == 1)
      numgas = get_last_gas(config_flags%chem_opt)
      if (config_flags%vertmix_onoff>0 &
          .and. grid%num_vert_mix.gt.0 &
          .and. config_flags%chem_opt.gt.0 &
         ) then

         if (grid%itimestep.gt.2) then

            if(config_flags%cu_diag .eq. 0 ) grid%raincv_b(:,:) = grid%raincv(:,:)

            call wrf_debug(15,'calling dry_dep_vel')
            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )

            DO ij = 1 , grid%num_tiles

               call dry_dep_vel(grid%id,curr_secs,grid%itimestep,grid%dt,config_flags,          &
                    grid%gmt,ifix(grid%julian),grid%t_phy,moist,p8w,t8w,                &
                    grid%alt,p_phy,chem,grid%rho,dz8w,grid%hfx,grid%dx,      & 
                    grid%raincv_b,            &
                    grid%ivgtyp,grid%tsk,grid%gsw,grid%vegfra,grid%pblh,   &
                    grid%rmol,grid%ust,grid%znt,grid%xlat,grid%xlong,                    &
                    grid%z,grid%z_at_w,grid%xland,                                &
                    grid%h2oaj,grid%h2oai,grid%nu3,grid%ac3,grid%cor3,grid%asulf,        &
                    grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,grid%cvalk1,grid%cvole1,&
                    grid%cvapi1,grid%cvapi2,grid%cvlim1,grid%cvlim2,     &
                    numgas,current_month,dvel,grid%snowh,  &
                    grid%dustdrydep_1,grid%dustdrydep_2,grid%dustdrydep_3,               &
                    grid%dustdrydep_4,grid%dustdrydep_5,                                 &
                    grid%depvelocity,                                                    &
                    grid%is_CAMMGMP_used,                  & !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
                    grid%dep_vel,grid%num_vert_mix,                      &
                    ids,ide, jds,jde, kds,kde,                                           &
                    ims,ime, jms,jme, kms,kme,                                           &
                    grid%i_start(ij), grid%i_end(ij), &
                    grid%j_start(ij), grid%j_end(ij), &
                    k_start, k_end                   )
                 
               if( config_flags%chem_opt == MOZART_KPP .or. &
                 config_flags%chem_opt == MOZCART_KPP .or. &
                 config_flags%chem_opt == MOZART_MOSAIC_4BIN_KPP .or. &
                 config_flags%chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP) then
                  call mozcart_lbc_set( chem, num_chem, grid%id, &
                                        ims, ime, jms, jme, kms, kme,    & 
                                        grid%i_start(ij), grid%i_end(ij), &
                                        grid%j_start(ij), grid%j_end(ij), &
                                        k_start )
               end if
            END DO
            !$OMP END PARALLEL DO
         else
            grid%dep_vel(:,:,:,:) = 0.D0
         end if
      else
         grid%dep_vel(:,:,:,:) = 0.D0
      endif
#endif

#ifdef DM_PARALLEL
      IF ( config_flags%bl_pbl_physics .EQ. 98 ) THEN
#     include "HALO_EM_PHYS_A_BL_SURF.inc"
      ENDIF
#endif

      CALL PUSHREAL8ARRAY(chem, (ime-ims+1)*(kme-kms+1)*(jme-jms+1)*num_chem)
      CALL PUSHREAL8ARRAY(grid%pblh, (ime-ims+1)*(jme-jms+1))

      CALL pbl_driver(                                                              &
     &         AKHS=grid%akhs          ,AKMS=grid%akms                              &
     &        ,BL_PBL_PHYSICS=config_flags%bl_pbl_physics                           &
     &        ,WINDFARM_OPT=config_flags%windfarm_opt,power=grid%power              &
     &        ,BLDT=grid%bldt, CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag &
     &        ,BLDTACTTIME=grid%bldtacttime                               &  
     &        ,BR=grid%br         ,CHKLOWQ=chklowq    ,CT=grid%ct                   &
     &        ,DT=grid%dt         ,DX=grid%dx         ,DY=grid%dy                   &
     &        ,DZ8W=dz8w          &
     &        ,EXCH_H=grid%exch_h     ,EXCH_M=grid%exch_m                           &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &
     &        ,F=grid%f  ,GRDFLX=grid%grdflx      &
     &        ,GZ1OZ0=gz1oz0      ,HFX=grid%hfx            ,HT=grid%ht              &
     &        ,ID=grid%id         ,ITIMESTEP=grid%itimestep     ,KPBL=grid%kpbl          &
     &        ,LH=grid%lh              ,LOWLYR=grid%lowlyr      ,P8W=grid%p_hyd_w   &
     &        ,PBLH=grid%pblh          ,PI_PHY=pi_phy      ,PSIH=psih          &
     &        ,PSIM=psim          ,P_PHY=grid%p_hyd        ,QFX=grid%qfx            &
     &        ,QSFC=grid%qsfc          ,QZ0=grid%qz0  ,MIXHT=mixht                  &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics                   &
     &        ,RHO=grid%rho            ,RQCBLTEN=grid%rqcblten  ,RQIBLTEN=grid%rqiblten  &
     &        ,RQVBLTEN=grid%rqvblten  ,RTHBLTEN=grid%rthblten  ,RUBLTEN=grid%rublten    &
     &        ,RVBLTEN=grid%rvblten    ,SNOW=grid%snow          ,STEPBL=grid%stepbl      &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy                               &
     &        ,TSK=grid%tsk            ,T_PHY=grid%t_phy        ,UST=grid%ust       &
     &        ,U10=grid%u10 ,UZ0=grid%uz0      ,U_FRAME=grid%u_frame    ,U_PHY=grid%u_phy   &
     &        ,V10=grid%v10 ,VZ0=grid%vz0      ,V_FRAME=grid%v_frame    ,V_PHY=grid%v_phy   &
     &        ,W=grid%w_2   ,UOCE=grid%uoce    ,VOCE=grid%voce                              &
              ,T2=grid%t2 &
     &        ,WARM_RAIN=grid%warm_rain                    ,WSPD=wspd          &
     &        ,XICE=grid%xice          ,XLAND=grid%xland        ,Z=grid%z                &
     &        ,ZNT=grid%znt                                                    &
     &        ,ysu_topdown_pblmix=config_flags%ysu_topdown_pblmix              &
     &        ,shinhong_tke_diag=config_flags%shinhong_tke_diag               &
! paj: topo_wind
     &        ,CTOPO=grid%ctopo,CTOPO2=grid%ctopo2                                   &
! variables added for BEP
     &          ,FRC_URB2D=grid%frc_urb2d                                                  &
     &          ,A_U_BEP=grid%a_u_bep,A_V_BEP=grid%a_v_bep,A_T_BEP=grid%a_t_bep          &
     &          ,A_Q_BEP=grid%a_q_bep                                                    &
     &          ,B_U_BEP=grid%b_u_bep,B_V_BEP=grid%b_v_bep,B_T_BEP=grid%b_t_bep          &
     &          ,B_Q_BEP=grid%b_q_bep                                                    &
     &          ,SF_BEP=grid%sf_bep,VL_BEP=grid%vl_bep                                   &
     &          ,A_E_BEP=grid%a_e_bep,B_E_BEP=grid%b_e_bep,DLG_BEP=grid%dlg_bep          &
     &          ,DL_U_BEP=grid%dl_u_bep                                                  &
     &          ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics                        &
     &          ,SF_URBAN_PHYSICS=config_flags%sf_urban_physics                          &
! Bep changes end
! add tke_pbl, and turbulent fluxes
     &          ,TKE_PBL=grid%tke_pbl,EL_PBL=grid%el_pbl,WU_TUR=grid%wu_tur              &
     &          ,WV_tur=grid%wv_tur,WT_tur=grid%wt_tur,WQ_tur=grid%wq_tur                &
! end add tke_pbl, and turbulent fluxes
! GBMPBL change: add exch_tke, rthraten
     &        ,EXCH_TKE=grid%exch_tke, RTHRATEN=grid%rthraten             &
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,I_START=grid%i_start,I_END=min(grid%i_end, ide-1)          &
     &        ,J_START=grid%j_start,J_END=min(grid%j_end, jde-1)          &
     &        ,KTS=k_start, KTE=min(k_end,kde-1)                          &
     &        ,NUM_TILES=grid%num_tiles                                   &
          ! Variables Required by ACM PBL  - jp
     &        ,ZNU=grid%znu,ZNW=grid%znw,MUT=grid%mut,P_TOP=grid%p_top    &
          ! Variables required by TEMF PBL - WA 9/9/08
              ,te_temf=grid%te_temf                                 &
              ,kh_temf=grid%kh_temf,km_temf=grid%km_temf            &
              ,shf_temf=grid%shf_temf,qf_temf=grid%qf_temf          &
              ,uw_temf=grid%uw_temf,vw_temf=grid%vw_temf            &
              ,hd_temf=grid%hd_temf,lcl_temf=grid%lcl_temf          &
              ,wupd_temf=grid%wupd_temf,mf_temf=grid%mf_temf          &
              ,thup_temf=grid%thup_temf,qtup_temf=grid%qtup_temf          &
              ,qlup_temf=grid%qlup_temf          &
              ,cf3d_temf=grid%cf3d_temf,cfm_temf=grid%cfm_temf          &
              ,hct_temf=grid%hct_temf                               &
              ,flhc=grid%flhc,flqc=grid%flqc                        &
              ,exch_temf=exch_temf                      &
          ! optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV                 &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC                 &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR                 &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI                 &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,HOL=HOL, MOL=grid%mol, REGIME=grid%REGIME                  &
!mynn mp@
     &        ,QKE=grid%qke, Sh3d=grid%sh3d                               &
     &        ,QKE_ADV=scalar(ims,kms,jms,P_qke_adv)                      & !ACF-QKE advection
     &        ,bl_mynn_tkeadvect=config_flags%bl_mynn_tkeadvect           & !ACF-QKE advection
     &        ,tsq=grid%tsq, qsq=grid%qsq, cov=grid%cov                   &
     &        ,DQKE=grid%dqke,QWT=grid%qWT                                &
     &        ,QSHEAR=grid%qSHEAR,QBUOY=grid%qBUOY,QDISS=grid%qDISS       &
     &        ,bl_mynn_tkebudget=config_flags%bl_mynn_tkebudget           &
     &        ,bl_mynn_cloudpdf=config_flags%bl_mynn_cloudpdf             &
     &        ,bl_mynn_mixlength=config_flags%bl_mynn_mixlength           &
     &        , icloud_bl=config_flags%icloud_bl                        & !JOE: subgrid cloud
     &        , qc_bl=grid%qc_bl,cldfra_bl=grid%cldfra_bl                 & !JOE: subgrid cloud
     &        ,bl_mynn_edmf=config_flags%bl_mynn_edmf                     & !JOE- MYNN edmf
     &        ,bl_mynn_edmf_mom=config_flags%bl_mynn_edmf_mom             & !JOE- MYNN edmf
     &        ,bl_mynn_edmf_tke=config_flags%bl_mynn_edmf_tke             & !JOE- MYNN edmf
     &        ,bl_mynn_edmf_part=config_flags%bl_mynn_edmf_part           & !JOE- MYNN edmf
     &        ,bl_mynn_cloudmix=config_flags%bl_mynn_cloudmix             & !JOE- MYNN cloud mixing
     &        ,bl_mynn_mixqt=config_flags%bl_mynn_mixqt                   & !JOE- MYNN tendency method
     &        ,edmf_a=grid%edmf_a,edmf_w=grid%edmf_w                      &
     &        ,edmf_thl=grid%edmf_thl,edmf_qt=grid%edmf_qt                &
     &        ,edmf_ent=grid%edmf_ent,edmf_qc=grid%edmf_qc                & !JOE- MYNN edmf
     &        ,rmol=grid%rmol, ch=grid%ch                                 &
     &        ,qcg=grid%qcg, grav_settling=config_flags%grav_settling     &
!    &        ,K_m=grid%K_m, K_h=grid%K_h, K_q=grid%K_q                   &
     &        ,vdfg=grid%vdfg                                             & ! Katata - fogdep
!GWD for ARW
     &        ,GWD_OPT=config_flags%gwd_opt &
     &        ,DTAUX3D=grid%dtaux3d,DTAUY3D=grid%dtauy3d &
     &        ,DUSFCG=grid%dusfcg,DVSFCG=grid%dvsfcg &
     &        ,VAR2D=grid%var2d,OC12D=grid%oc12d     &
     &        ,OA1=grid%oa1,OA2=grid%oa2,OA3=grid%oa3,OA4=grid%oa4        &
     &        ,OL1=grid%ol1,OL2=grid%ol2,OL3=grid%ol3,OL4=grid%ol4        &
     &        ,MFSHCONV=config_flags%mfshconv                             &
     &        ,MASSFLUX_EDKF=grid%massflux_EDKF                           & 
     &        ,ENTR_EDKF=grid%entr_EDKF, DETR_EDKF=grid%detr_EDKF         &
     &        ,THL_UP=grid%thl_up                                         &
     &        ,THV_UP=grid%thv_up, RT_UP=grid%rt_up ,RV_UP=grid%rv_up     &
     &        ,RC_UP=grid%rc_up, U_UP=grid% u_up, V_UP=grid%v_up          &
     &        ,FRAC_UP=grid%frac_up, RC_MF=grid%RC_MF                     &
! For Wind Turbine Drag Parameterizations
     &        ,phb=grid%phb                                               &
     &        ,XLAT_U=grid%xlat_u,XLONG_U=grid%xlong_u                    &
!Variables required for camuwpbl scheme
     &        ,Z_AT_W=grid%z_at_w,CLDFRA_OLD_MP=grid%cldfra_old_mp        &
     &        ,CLDFRA=grid%cldfra                                        &
     &        ,RTHRATENLW=grid%rthratenlw,TAURESX2D=grid%tauresx2d        &
     &        ,TAURESY2D=grid%tauresy2d                                   &
     &        ,TPERT2D=grid%tpert2d,QPERT2D=grid%qpert2d                  &
     &        ,WPERT2D=grid%wpert2d,WSEDL3D=grid%wsedl3d                  &
     &        ,TURBTYPE3D=grid%turbtype3d,SMAW3D=grid%smaw3d              &
     &        ,QNC_CURR=scalar(ims,kms,jms,P_QNC), F_QNC=f_qnc            &
     &        ,QNI_CURR=scalar(ims,kms,jms,P_QNI), F_QNI=f_qni            &
     &        ,RQNIBLTEN=grid%rqniblten                                   &
     &        ,XLAT_V=grid%xlat_v,XLONG_V=grid%xlong_v,FNM=grid%fnm       &
     &        ,FNP=grid%fnp, IS_CAMMGMP_USED = grid%is_CAMMGMP_used       &
! for grims shallow convection with ysupbl
     &        ,WSTAR=grid%wstar_ysu,DELTA=grid%delta_ysu                  &
! for pbl mixing of scalars and tracers
     &        ,SCALAR=scalar,SCALAR_TEND=scalar_tend,NUM_SCALAR=num_scalar&
     &        ,TRACER=tracer,TRACER_TEND=tracer_tend,NUM_TRACER=num_tracer&
     &        ,SCALAR_PBLMIX=config_flags%scalar_pblmix                  &
     &        ,TRACER_PBLMIX=config_flags%tracer_pblmix                  &
#if (WRF_CHEM == 1)
     &        ,CHEM=chem,VD=grid%dep_vel                                 &
     &        ,NCHEM=num_chem,kdvel=config_flags%kdepvel                 &
     &        ,ndvel=config_flags%ndepvel                                &
     &        ,NUM_VERT_MIX=grid%num_vert_mix                            &
#endif
     &        ,QNORM=grid%QNORM, fasdas=config_flags%fasdas              & !fasdas
     &        )

BENCH_END(pbl_driver_tim)

#ifdef WRF_CHEM
# ifdef DM_PARALLEL
   IF ( num_chem >= PARAM_FIRST_SCALAR .AND. config_flags%bl_pbl_physics == ACMPBLSCHEME) then
     CALL wrf_debug ( 200 , ' call HALO CHEM AFTER PBL' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_CHEM_E_3.inc"
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_CHEM_E_5.inc"

     ELSE
       WRITE(message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(message))
     ENDIF
   ENDIF
# endif
# endif



!  Adjoint part begin
!#if ( WRF_CHEM == 1 )
!        grid%raincv_b(:,:) = grid%raincv(:,:) 
!#endif
! cumulus para.

      CALL wrf_debug ( 200 , ' call a_cumulus_driver' )

BENCH_START(a_cu_driver_tim)
      CALL a_cumulus_driver(grid                                             &
                 ! Prognostic variables
!     &             ,U=grid%u_phy, V=grid%v_phy, TH=th_phy, a_TH=a_th_phy, T=grid%t_phy, a_T=grid%a_t_phy  &
!     &             ,W=grid%w_2, P=grid%p_hyd, a_P=grid%a_p_hyd, PI=pi_phy, RHO=grid%rho, a_RHO=grid%a_rho &
     &             ,U=grid%u_phy, V=grid%v_phy, TH=th_phy, THB=a_th_phy, T=grid%t_phy, TB=grid%a_t_phy  &
     &             ,W=grid%w_2, P=grid%p_hyd, PB=grid%a_p_hyd, PI=pi_phy, RHO=grid%rho, RHOB=grid%a_rho &
                 ! Other arguments
     &             ,ITIMESTEP=grid%itimestep ,DT=grid%dt      ,DX=grid%dx              &
     &             ,CUDT=grid%cudt,CURR_SECS=curr_secs,ADAPT_STEP_FLAG=adapt_step_flag &
     &             ,CUDTACTTIME=grid%cudtacttime                                       &
     &             ,RAINC=grid%rainc   ,RAINCV=grid%raincv   ,PRATEC=grid%pratec       &
     &             ,RAINCVB=grid%a_raincv, PRATECB=grid%a_pratec                       &
     &             ,NCA=grid%nca                                                       &
     &             ,HTOP=grid%cutop     ,HBOT=grid%cubot       ,KPBL=grid%kpbl         &
     &             ,Z=grid%z ,Z_AT_W=grid%z_at_w ,MAVAIL=grid%mavail ,PBLH=grid%pblh   &
     &             ,DZ8W=dz8w     ,P8W=grid%p_hyd_w, PSFC=grid%psfc, TSK=grid%tsk      &
     &             ,ZB=grid%a_z, DZ8WB=a_dz8w                                        &
     &             ,TKE_PBL=grid%tke_pbl, UST=grid%ust                                 &
     &             ,W0AVG=grid%w0avg   ,STEPCU=grid%stepcu                             &
     &             ,CLDEFI=grid%cldefi ,LOWLYR=grid%lowlyr ,XLAND=grid%xland           &
     &             ,APR_GR=grid%apr_gr ,APR_W=grid%apr_w   ,APR_MC=grid%apr_mc         &
     &             ,APR_ST=grid%apr_st ,APR_AS=grid%apr_as ,APR_CAPMA=grid%apr_capma   &
     &             ,APR_CAPME=grid%apr_capme          ,APR_CAPMI=grid%apr_capmi     &
     &             ,MASS_FLUX=grid%mass_flux          ,XF_ENS=grid%xf_ens           &
     &             ,PR_ENS=grid%pr_ens ,HT=grid%ht,EDT_OUT=grid%edt_out             &
     &             ,imomentum=grid%imomentum,clos_choice=grid%clos_choice           &
     &             ,ishallow=config_flags%ishallow                                  &
     &             ,cugd_tten=grid%cugd_tten,cugd_qvten=grid%cugd_qvten,cugd_qcten=grid%cugd_qcten   &
     &             ,cugd_ttens=grid%cugd_ttens,cugd_qvtens=grid%cugd_qvtens   &
     &             ,ENSDIM=config_flags%ensdim ,MAXIENS=config_flags%maxiens ,MAXENS=config_flags%maxens         &
     &             ,MAXENS2=config_flags%maxens2 ,MAXENS3=config_flags%maxens3       &
     &             ,CU_ACT_FLAG=cu_act_flag   ,WARM_RAIN=grid%warm_rain   &
     &             ,HFX=grid%hfx, QFX=grid%qfx                            &
     &             ,CLDFRA=grid%cldfra, TPERT2D=grid%tpert2d              &
     &             ,GSW=grid%gsw,cugd_avedx=config_flags%cugd_avedx       &
     &             ,k22_shallow=grid%k22_shallow,kbcon_shallow=grid%kbcon_shallow &
     &             ,ktop_shallow=grid%ktop_shallow,xmb_shallow=grid%xmb_shallow   &
     &             ,PERIODIC_X=(config_flags%polar .OR. config_flags%periodic_x)  &
     &             ,PERIODIC_Y=config_flags%periodic_y                    &
                 ! Zhang-McFarlane outputs
     &             ,CAPE=grid%cape ,ZMMU=grid%zmmu ,ZMMD=grid%zmmd        &
     &             ,ZMDT=grid%zmdt ,ZMDQ=grid%zmdq                        &
     &             ,DLF=grid%dlf, RLIQ=grid%rliq                          &
     &             ,PCONVB=grid%pconvb, PCONVT=grid%pconvt                &
     &             ,EVAPTZM=grid%evaptzm, FZSNTZM=grid%fzsntzm, EVSNTZM=grid%evsntzm     &
     &             ,EVAPQZM=grid%evapqzm, ZMFLXPRC=grid%zmflxprc          &
     &             ,ZMFLXSNW=grid%zmflxsnw, ZMNTPRPD=grid%zmntprpd        &
     &             ,ZMNTSNPD=grid%zmntsnpd, ZMEIHEAT=grid%zmeiheat        &
     &             ,CMFMC=grid%cmfmc, CMFMCDZM=grid%cmfmcdzm              &
     &             ,PRECCDZM=grid%preccdzm, PRECZ=grid%precz              &
     &             ,ZMMTU=grid%zmmtu, ZMMTV=grid%zmmtv  &
     &             ,ZMUPGU=grid%zmupgu, ZMUPGD=grid%zmupgd                &
     &             ,ZMVPGU=grid%zmvpgu, ZMVPGD=grid%zmvpgd                &
     &             ,ZMICUU=grid%zmicuu, ZMICUD=grid%zmicud                &
     &             ,ZMICVU=grid%zmicvu, ZMICVD=grid%zmicvd                &
     &             ,ZMDICE=grid%zmdice, ZMDLIQ=grid%zmdliq                &
                 ! Selection flag
     &             ,pgcon=config_flags%sas_pgcon                          &
     &             ,CU_PHYSICS=config_flags%cu_physics                    &
     &             ,BL_PBL_PHYSICS=config_flags%bl_pbl_physics            &
     &             ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics      &
     &             ,KFETA_TRIGGER=config_flags%kfeta_trigger              &
                 ! Dimension arguments
     &             ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde     &
     &             ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme     &
     &             ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe     &
     &             ,I_START=grid%i_start,I_END=min(grid%i_end, ide-1)     &
     &             ,J_START=grid%j_start,J_END=min(grid%j_end, jde-1)     &
     &             ,KTS=k_start, KTE=min(k_end,kde-1)                     &
     &             ,NUM_TILES=grid%num_tiles                              &
                 ! Moisture tendency arguments
     &             ,RQVCUTEN=grid%rqvcuten , RQCCUTEN=grid%rqccuten       &
     &             ,RQVCUTENB=grid%a_rqvcuten                            &
     &             ,RQSCUTEN=grid%rqscuten , RQICUTEN=grid%rqicuten       &
     &             ,RQRCUTEN=grid%rqrcuten                                &
     &             ,RQVBLTEN=grid%rqvblten , RQVFTEN=grid%rqvften         &
                 ! Other tendency arguments
     &             ,RTHRATEN=grid%rthraten , RTHBLTEN=grid%rthblten       &
     &             ,RUCUTEN=grid%rucuten   , RVCUTEN=grid%rvcuten         &
     &             ,RTHCUTEN=grid%rthcuten , RTHFTEN=grid%rthften         &
     &             ,RTHCUTENB=grid%a_rthcuten                            &
                 ! Moisture tracer arguments
     &             ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV            &
     &             ,QV_CURRB=a_moist(ims,kms,jms,P_QV)                   &
     &             ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC            &
     &             ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR            &
     &             ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI            &
     &             ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS            &
     &             ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG            &
! Variables for Tiedtke and NSAS schemes
     &             ,ZNU=grid%znu                                          &
     &             ,MP_PHYSICS=config_flags%mp_physics                    &
     &             ,GD_CLOUD=grid%GD_CLOUD,GD_CLOUD2=grid%GD_CLOUD2       &
     &             ,cfu1=grid%cfu1,cfd1=grid%cfd1,dfu1=grid%dfu1,efu1=grid%efu1,dfd1=grid%dfd1,efd1=grid%efd1,f_flux=l_flux)
BENCH_END(a_cu_driver_tim)
!


#if (WRF_CHEM == 1)
# ifdef DM_PARALLEL
   IF ( num_chem >= PARAM_FIRST_SCALAR .AND. config_flags%bl_pbl_physics == ACMPBLSCHEME) then
     CALL wrf_debug ( 200 , ' call HALO CHEM AFTER PBL' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_CHEM_E_3_AD.inc"
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_CHEM_E_5_AD.inc"

     ELSE
       WRITE(message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(message))
     ENDIF
   ENDIF
# endif
# endif



!*********
! pbl

      CALL wrf_debug ( 200 , ' call a_pbl_driver' )
BENCH_START(a_pbl_driver_tim)
      CALL POPREAL8ARRAY(grid%pblh, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(chem, (ime-ims+1)*(kme-kms+1)*(jme-jms+1)*num_chem)
      CALL a_pbl_driver(                                                            &
     &         AKHS=grid%akhs          ,AKMS=grid%akms                              &
     &        ,BL_PBL_PHYSICS=config_flags%bl_pbl_physics                           &
     &        ,BLDT=grid%bldt, CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag &
     &        ,BLDTACTTIME=grid%bldtacttime                                         &
     &        ,BR=grid%br              ,CHKLOWQ=chklowq    ,CT=grid%ct                   &
     &        ,DT=grid%dt              ,DX=grid%dx    ,DZ8W=dz8w, A_DZ8W=a_dz8w     &
     &        ,EXCH_H=grid%exch_h     ,EXCH_M=grid%exch_m                           &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &
     &        ,F=grid%f  ,GRDFLX=grid%grdflx      &
     &        ,GZ1OZ0=gz1oz0      ,HFX=grid%hfx ,a_HFX=grid%a_hfx  ,HT=grid%ht      &
     &        ,ID=grid%id         ,ITIMESTEP=grid%itimestep     ,KPBL=grid%kpbl     &
     &        ,LH=grid%lh ,LOWLYR=grid%lowlyr ,P8W=grid%p_hyd_w, a_P8W=grid%a_p_hyd_w   &
     &        ,PBLH=grid%pblh, A_PBLH=grid%a_pblh ,PI_PHY=pi_phy,a_PI_PHY=a_pi_phy, PSIH=psih  &
     &        ,PSIM=psim ,P_PHY=grid%p_hyd, QFX=grid%qfx  ,a_QFX=grid%a_qfx         &
     &        ,QSFC=grid%qsfc          ,QZ0=grid%qz0  ,MIXHT=mixht                  &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics                   &
     &        ,RHO=grid%rho,     RQCBLTEN=grid%rqcblten  ,RQIBLTEN=grid%rqiblten  &
     &        ,A_RHO=grid%a_rho, A_RQCBLTEN=grid%a_rqcblten  ,A_RQIBLTEN=grid%a_rqiblten &
     &        ,RQVBLTEN=grid%rqvblten  ,a_RQVBLTEN=grid%a_rqvblten                  &
     &        ,RTHBLTEN=grid%rthblten  ,a_RTHBLTEN=grid%a_rthblten                  &
     &        ,RUBLTEN=grid%rublten, a_RUBLTEN=grid%a_rublten    &
     &        ,RVBLTEN=grid%rvblten ,a_RVBLTEN=grid%a_rvblten   ,SNOW=grid%snow   ,STEPBL=grid%stepbl  &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy     ,A_TH_PHY=a_th_phy            &
     &        ,TSK=grid%tsk            ,T_PHY=grid%t_phy  ,A_T_PHY=grid%a_t_phy                   &
     &        ,UST=grid%ust ,a_UST=grid%a_ust &
     &        ,U10=grid%u10 ,UZ0=grid%uz0 ,U_FRAME=grid%u_frame ,U_PHY=grid%u_phy ,a_U_PHY=grid%a_u_phy &
     &        ,V10=grid%v10 ,VZ0=grid%vz0 ,V_FRAME=grid%v_frame ,V_PHY=grid%v_phy ,a_V_PHY=grid%a_v_phy &
     &        ,T2=grid%t2 &
     &        ,WARM_RAIN=grid%warm_rain  ,WSPD=wspd    ,a_WSPD=a_wspd          &
     &        ,XICE=grid%xice   ,XLAND=grid%xland  ,Z=grid%z ,a_Z=grid%a_z      &
     &        ,ZNT=grid%znt                                                    &
! paj: topo_wind
     &        ,CTOPO=grid%ctopo,CTOPO2=grid%ctopo2                                   &
! variables added for BEP
     &          ,FRC_URB2D=grid%frc_urb2d                                                  &
     &          ,A_U_BEP=grid%a_u_bep,A_V_BEP=grid%a_v_bep,A_T_BEP=grid%a_t_bep          &
     &          ,A_Q_BEP=grid%a_q_bep                                                    &
     &          ,B_U_BEP=grid%b_u_bep,B_V_BEP=grid%b_v_bep,B_T_BEP=grid%b_t_bep          &
     &          ,B_Q_BEP=grid%b_q_bep                                                    &
     &          ,SF_BEP=grid%sf_bep,VL_BEP=grid%vl_bep                                   &
     &          ,A_E_BEP=grid%a_e_bep,B_E_BEP=grid%b_e_bep,DLG_BEP=grid%dlg_bep          &
     &          ,DL_U_BEP=grid%dl_u_bep                                                  &
     &          ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics                        &
     &          ,SF_URBAN_PHYSICS=config_flags%sf_urban_physics                          &
! Bep changes end
! add tke_pbl, and turbulent fluxes
     &          ,TKE_PBL=grid%tke_pbl,EL_PBL=grid%el_pbl,WU_TUR=grid%wu_tur,WV_tur=grid%wv_tur,WT_tur=grid%wt_tur,WQ_tur=grid%wq_tur   &
! end add tke_pbl, and turbulent fluxes
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        ,I_START=grid%i_start,I_END=min(grid%i_end, ide-1)          &
     &        ,J_START=grid%j_start,J_END=min(grid%j_end, jde-1)          &
     &        ,KTS=k_start, KTE=min(k_end,kde-1)                          &
     &        ,NUM_TILES=grid%num_tiles                                   &
          ! Variables Required by ACM PBL  - jp
     &        ,ZNU=grid%znu,ZNW=grid%znw,MUT=grid%mut,a_MUT=grid%a_mut,P_TOP=grid%p_top    &
          ! Variables required by TEMF PBL - WA 9/9/08
              ,te_temf=grid%te_temf                                 &
              ,kh_temf=grid%kh_temf,km_temf=grid%km_temf            &
              ,shf_temf=grid%shf_temf,qf_temf=grid%qf_temf          &
              ,uw_temf=grid%uw_temf,vw_temf=grid%vw_temf            &
              ,hd_temf=grid%hd_temf,lcl_temf=grid%lcl_temf          &
              ,wupd_temf=grid%wupd_temf,mf_temf=grid%mf_temf          &
              ,thup_temf=grid%thup_temf,qtup_temf=grid%qtup_temf          &
              ,qlup_temf=grid%qlup_temf          &
              ,cf3d_temf=grid%cf3d_temf,cfm_temf=grid%cfm_temf          &
              ,hct_temf=grid%hct_temf                               &
              ,flhc=grid%flhc,flqc=grid%flqc                        &
              ,exch_temf=exch_temf                      &
          ! optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV                 &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC), F_QC=F_QC                 &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR), F_QR=F_QR                 &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI), F_QI=F_QI                 &
     &        ,a_QV_CURR=a_moist(ims,kms,jms,P_QV)                        &
     &        ,a_QC_CURR=a_moist(ims,kms,jms,P_QC)                        &
     &        ,a_QI_CURR=a_moist(ims,kms,jms,P_QI)                        &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,HOL=HOL, MOL=grid%mol, REGIME=grid%REGIME                  &
!mynn mp
     &        ,QKE=scalar(ims,kms,jms,P_qke_adv), tsq=grid%tsq, qsq=grid%qsq  &
     &        ,cov=grid%cov                                               &
     &        ,DQKE=grid%dqke,QWT=grid%qWT                                &
     &        ,QSHEAR=grid%qSHEAR,QBUOY=grid%qBUOY,QDISS=grid%qDISS &
     &        ,bl_mynn_tkebudget=config_flags%bl_mynn_tkebudget           &
     &        ,rmol=grid%rmol, ch=grid%ch                                 &
     &        ,qcg=grid%qcg, grav_settling=config_flags%grav_settling             &
!    &        ,K_m=grid%K_m, K_h=grid%K_h, K_q=grid%K_q                   &
!GWD for ARW
     &        ,GWD_OPT=config_flags%gwd_opt &
     &        ,DTAUX3D=grid%dtaux3d,DTAUY3D=grid%dtauy3d &
     &        ,A_DTAUX3D=grid%a_dtaux3d,A_DTAUY3D=grid%a_dtauy3d &
     &        ,DUSFCG=grid%dusfcg,DVSFCG=grid%dvsfcg &
     &        ,A_DUSFCG=grid%a_dusfcg,A_DVSFCG=grid%a_dvsfcg &
     &        ,VAR2D=grid%var2d,OC12D=grid%oc12d     &
     &        ,OA1=grid%oa1,OA2=grid%oa2,OA3=grid%oa3,OA4=grid%oa4        &
     &        ,OL1=grid%ol1,OL2=grid%ol2,OL3=grid%ol3,OL4=grid%ol4        &
     &        ,MFSHCONV=grid%mfshconv                                     & 
     &        ,MASSFLUX_EDKF=grid%massflux_EDKF                           & 
     &        ,ENTR_EDKF=grid%entr_EDKF, DETR_EDKF=grid%detr_EDKF         &
     &        ,THL_UP=grid%thl_up                                         &
     &        ,THV_UP=grid%thv_up, RT_UP=grid%rt_up ,RV_UP=grid%rv_up     &
     &        ,RC_UP=grid%rc_up, U_UP=grid% u_up, V_UP=grid%v_up          &
     &        ,FRAC_UP=grid%frac_up, RC_MF=grid%RC_MF                     &
! For Wind Turbine Drag Parameterizations
     &        ,phb=grid%phb                                               &
     &        ,XLAT_U=grid%xlat_u,XLONG_U=grid%xlong_u                    &
!Variables required for camuwpbl scheme
     &        ,Z_AT_W=grid%z_at_w,CLDFRA=grid%cldfra                      &
     &        ,RTHRATENLW=grid%rthratenlw,TAURESX2D=grid%tauresx2d        &
     &        ,TAURESY2D=grid%tauresy2d                                   &
     &        ,TPERT2D=grid%tpert2d,QPERT2D=grid%qpert2d                  &
     &        ,WPERT2D=grid%wpert2d                                       &
     &        ,XLAT_V=grid%xlat_v,XLONG_V=grid%xlong_v                    &
#if (WRF_CHEM == 1)
     &        ,CHEM=chem  ,A_CHEM=a_chem  ,VD=grid%dep_vel ,A_VD=grid%a_dep_vel   &
     &        ,NCHEM=num_chem,kdvel=config_flags%kdepvel                 &
     &        ,ndvel=config_flags%ndepvel                                &
     &        ,NUM_VERT_MIX=grid%num_vert_mix                     &
#endif
     &        )


      CALL POPREAL8ARRAY(grid%pblh, (ime-ims+1)*(jme-jms+1))
#if ( WRF_CHEM == 1)
      numgas = get_last_gas(config_flags%chem_opt)
      if (config_flags%vertmix_onoff>0 &
          .and. grid%num_vert_mix.gt.0 &
          .and. config_flags%chem_opt.gt.0 &
         ) then

         if (grid%itimestep.gt.2) then

!            if(config_flags%cu_diag .eq. 0 ) grid%raincv_b(:,:) = grid%raincv(:,:)

            call wrf_debug(15,'calling dry_dep_vel')
            !$OMP PARALLEL DO   &
            !$OMP PRIVATE ( ij )

            DO ij = 1 , grid%num_tiles

               call a_dry_dep_vel(grid%id,curr_secs,grid%itimestep,grid%dt,config_flags,          &
                    grid%gmt,ifix(grid%julian),grid%t_phy,grid%a_t_phy,moist,a_moist,p8w,t8w,                &
                    grid%alt,p_phy,a_p_phy,chem,grid%rho,dz8w,grid%hfx,grid%dx,      & 
                    grid%raincv_b,            &
                    grid%ivgtyp,grid%tsk,grid%gsw,grid%vegfra,grid%pblh,grid%a_pblh,   &
                    grid%rmol,grid%a_rmol,grid%ust,grid%a_ust,grid%znt,grid%a_znt,     &
                    grid%xlat,grid%xlong,                    &
                    grid%z,grid%z_at_w,grid%a_z_at_w,grid%xland,                                &
                    grid%h2oaj,grid%h2oai,grid%nu3,grid%ac3,grid%cor3,grid%asulf,        &
                    grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,grid%cvalk1,grid%cvole1,&
                    grid%cvapi1,grid%cvapi2,grid%cvlim1,grid%cvlim2,     &
                    numgas,current_month,dvel,grid%snowh,  &
                    grid%dustdrydep_1,grid%dustdrydep_2,grid%dustdrydep_3,               &
                    grid%dustdrydep_4,grid%dustdrydep_5,                                 &
                    grid%depvelocity,                                                    &
                    grid%is_CAMMGMP_used,                  & !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
                    grid%dep_vel,grid%a_dep_vel,grid%num_vert_mix,                      &
                    ids,ide, jds,jde, kds,kde,                                           &
                    ims,ime, jms,jme, kms,kme,                                           &
                    grid%i_start(ij), grid%i_end(ij), &
                    grid%j_start(ij), grid%j_end(ij), &
                    k_start, k_end                   )

!               if( config_flags%chem_opt == MOZART_KPP .or. &
!                 config_flags%chem_opt == MOZCART_KPP .or. &
!                 config_flags%chem_opt == MOZART_MOSAIC_4BIN_KPP .or. &
!                 config_flags%chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP) then
!                  call mozcart_lbc_set( chem, num_chem, grid%id, &
!                                        ims, ime, jms, jme, kms, kme,    & 
!                                        grid%i_start(ij), grid%i_end(ij), &
!                                        grid%j_start(ij), grid%j_end(ij), &
!                                        k_start )
!               end if
            END DO
            !$OMP END PARALLEL DO
         else
            grid%a_dep_vel(:,:,:,:) = 0.D0
         end if
      else
         grid%a_dep_vel(:,:,:,:) = 0.D0
      endif
#endif
BENCH_END(a_pbl_driver_tim)

!!#ifdef DM_PARALLEL
!!      IF ( config_flags%bl_pbl_physics .EQ. 98 ) THEN
!!#     include "HALO_EM_PHYS_A_BL_SURF.inc"
!!      ENDIF
!!#endif

BENCH_START(a_surf_driver_tim)

      CALL POPREAL8ARRAY(grid%tslb, (ime-ims+1)*(jme-jms+1)*grid%num_soil_layers)
      CALL POPREAL8ARRAY(grid%smois, (ime-ims+1)*(jme-jms+1)*grid%num_soil_layers)
      CALL POPREAL8ARRAY(grid%tsk, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%qfx, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%hfx, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%canwat, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%znt, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%mol, (ime-ims+1)*(jme-jms+1))
!      CALL POPREAL8ARRAY(grid%rmol, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%ustm, (ime-ims+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%ust, (ime-ims+1)*(jme-jms+1))
CALL a_surface_driver(                                                &
     &        HYDRO_dt=HYDRO_dt,  sfcheadrt=grid%sfcheadrt,               &
     &        INFXSRT=grid%INFXSRT,  soldrain=grid%soldrain,              &
     &         ACGRDFLX=grid%acgrdflx  ,ACHFX=grid%achfx        ,ACLHF=grid%aclhf        &       
     &        ,ACSNOM=grid%acsnom      ,ACSNOW=grid%acsnow      ,AKHS=grid%akhs          &
     &        ,AKMS=grid%akms          ,ALBBCK=grid%albbck      ,ALBEDO=grid%albedo      &
     &        ,EMBCK=grid%embck                                                          &
     &        ,BR=grid%br              ,CANWAT=grid%canwat,a_CANWAT=grid%a_canwat      ,CHKLOWQ=chklowq    &
     &        ,CT=grid%ct              ,DT=grid%dt         ,DX=grid%dx         &
     &        ,DZ8W=dz8w     ,a_DZ8W=a_dz8w    ,DZS=grid%dzs            ,FLHC=grid%flhc          &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &    
     &        ,FLQC=grid%flqc          ,GLW=grid%glw            ,GRDFLX=grid%grdflx      &
     &        ,GSW=grid%gsw    ,SWDOWN=grid%swdown        ,GZ1OZ0=gz1oz0     &
     &        ,HFX=grid%hfx            ,a_HFX=grid%a_hfx              &
     &        ,HT=grid%ht              ,IFSNOW=config_flags%ifsnow      ,ISFFLX=config_flags%isfflx      &      
     &        ,FRACTIONAL_SEAICE=config_flags%fractional_seaice           &
     &        ,SEAICE_ALBEDO_OPT=config_flags%seaice_albedo_opt           &         
     &        ,SEAICE_ALBEDO_DEFAULT=config_flags%seaice_albedo_default   &         
     &        ,SEAICE_THICKNESS_OPT=config_flags%seaice_thickness_opt     &         
     &        ,SEAICE_THICKNESS_DEFAULT=config_flags%seaice_thickness_default     & 
     &        ,SEAICE_SNOWDEPTH_OPT=config_flags%seaice_snowdepth_opt     &         
     &        ,SEAICE_SNOWDEPTH_MAX=config_flags%seaice_snowdepth_max     &
     &        ,SEAICE_SNOWDEPTH_MIN=config_flags%seaice_snowdepth_min     &
     &        ,TICE2TSK_IF2COLD=config_flags%tice2tsk_if2cold             &
     &        ,IFNDALBSI=grid%ifndalbsi, IFNDICEDEPTH=grid%ifndicedepth   &            
     &        ,IFNDSNOWSI=grid%ifndsnowsi                                 &            
     &        ,ISLTYP=grid%isltyp      ,ITIMESTEP=grid%itimestep, JULIAN_IN=grid%julian                  &   
     &        ,IVGTYP=grid%ivgtyp      ,LH=grid%lh              ,LOWLYR=grid%lowlyr      &
     &        ,MAVAIL=grid%mavail      ,NUM_SOIL_LAYERS=config_flags%num_soil_layers        &
     &        ,P8W=grid%p_hyd_w        ,a_P8W=grid%a_p_hyd_w   &          
     &        ,PBLH=grid%pblh          ,a_PBLH=grid%a_pblh    ,PI_PHY=pi_phy  &
     &        ,PSFC=grid%psfc          ,PSHLTR=grid%pshltr      ,PSIH=psih          &
     &        ,BLDT=grid%bldt     ,CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag  &
     &        ,BLDTACTTIME=grid%bldtacttime                                              &
     &        ,PSIM=psim               ,P_PHY=grid%p_hyd        ,a_P_PHY=grid%a_p_hyd    &
     &        ,Q10=grid%q10            &
     &        ,Q2=grid%q2              ,QFX=grid%qfx    ,a_QFX=grid%a_qfx  ,QSFC=grid%qsfc &
     &        ,QSHLTR=grid%qshltr      ,QZ0=grid%qz0            ,RAINCV=grid%raincv      &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics     ,RHO=grid%rho   ,a_RHO=grid%a_rho      &
     &        ,RMOL=grid%rmol ,a_RMOL=grid%a_rmol ,SFCEVP=grid%sfcevp      ,SFCEXC=grid%sfcexc  &
     &        ,SFCRUNOFF=grid%sfcrunoff                                        &
     &        ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics                        &
     &        ,SF_SURFACE_PHYSICS=config_flags%sf_surface_physics  ,SH2O=grid%sh2o          &
     &        ,SHDMAX=grid%shdmax      ,SHDMIN=grid%shdmin      &
     &        ,SMOIS=grid%smois        ,a_SMOIS=grid%a_smois        &
     &        ,SMSTAV=grid%smstav      ,SMSTOT=grid%smstot      ,SNOALB=grid%snoalb      &
     &        ,SNOW=grid%snow          ,SNOWC=grid%snowc        ,SNOWH=grid%snowh        &
     &        ,SMCREL=grid%smcrel                                                        &
     &        ,SST=grid%sst            ,SST_UPDATE=grid%sst_update                  &
     &        ,SSTSK=grid%sstsk        ,DTW=grid%dtw            ,SST_SKIN=grid%sst_skin  &
     &        ,SCM_FORCE_SKINTEMP=grid%scm_force_skintemp                                &
     &        ,SCM_FORCE_FLUX=grid%scm_force_flux                                        &
     &        ,STEPBL=grid%stepbl      ,TH10=grid%th10          ,TH2=grid%th2            &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy  ,a_TH_PHY=a_th_phy                &
     &        ,TKE_PBL=grid%tke_pbl    &
     &        ,TMN=grid%tmn            ,TSHLTR=grid%tshltr      ,TSK=grid%tsk            &
     &        ,TYR=grid%tyr            ,TYRA=grid%tyra          ,TDLY=grid%tdly          &
     &        ,TLAG=grid%tlag          ,LAGDAY=config_flags%lagday      ,NYEAR=grid%nyear        &
     &        ,NDAY=grid%nday          ,TMN_UPDATE=grid%tmn_update      ,YR=yr            &
     &        ,TSLB=grid%tslb  ,a_TSLB=grid%a_tslb ,T_PHY=grid%t_phy  ,a_T_PHY=grid%a_t_phy   &
     &        ,U10=grid%u10            &
     &        ,URATX=grid%uratx        ,VRATX=grid%vratx   ,TRATX=grid%tratx        &
     &        ,UDRUNOFF=grid%udrunoff  ,UST=grid%ust ,a_UST=grid%a_ust  ,UZ0=grid%uz0    &
     &        ,U_FRAME=grid%u_frame    ,U_PHY=grid%u_phy  ,a_U_PHY=grid%a_u_phy   ,V10=grid%v10    &
     &        ,VEGFRA=grid%vegfra      ,VZ0=grid%vz0       ,V_FRAME=grid%v_frame    &
     &        ,V_PHY=grid%v_phy             ,a_V_PHY=grid%a_v_phy    ,WARM_RAIN=grid%warm_rain     &
     &        ,WSPD=wspd       ,a_WSPD=a_wspd        ,XICE=grid%xice     ,XLAND=grid%xland        &
     &        ,Z0=grid%z0              ,Z=grid%z        ,ZNT=grid%znt, a_ZNT=grid%a_znt &
     &        ,ZS=grid%zs              ,ALBSI=grid%albsi , ICEDEPTH=grid%icedepth   &
     &        ,SNOWSI=grid%snowsi                                                   &
     &        ,XICEM=grid%xicem   ,ISICE=grid%landuse_isice                         &
     &        ,USTM=grid%ustm     ,a_USTM=grid%a_ustm      ,CK=grid%ck    ,CKA=grid%cka   &
     &                                 ,CD=grid%cd         ,CDA=grid%cda            &
     &        ,ISFTCFLX=config_flags%isftcflx, IZ0TLND=config_flags%iz0tlnd         &
     &        ,SF_OCEAN_PHYSICS=config_flags%sf_ocean_physics                       &
     &        ,OML_HML0=config_flags%oml_hml0 ,OML_GAMMA=config_flags%oml_gamma     &
     &        ,TML=grid%tml, T0ML=grid%t0ml                                         &
     &        ,HML=grid%hml, H0ML=grid%h0ml, HUML=grid%huml, HVML=grid%hvml, F=grid%f &
     &        ,TMOML=grid%TMOML,ISWATER=iswater &
     &        ,SLOPE_RAD=config_flags%slope_rad,TOPO_SHADING=config_flags%topo_shading & ! solar
     &        ,SHADOWMASK=grid%shadowmask                                              & ! solar
     &        ,SLOPE=grid%slope, SLP_AZI=grid%slp_azi, SWNORM=grid%swnorm              & ! solar
     &        ,DECLIN=grid%declin ,SOLCON=grid%solcon ,COSZEN=grid%coszen ,HRANG=grid%hrang    &
           ! P-X LSM Variables
     &        ,LANDUSEF=grid%landusef, SOILCTOP=grid%soilctop             &   ! P-X LSM
     &        ,SOILCBOT=grid%soilcbot                                     &   ! P-X LSM
     &        ,RA=grid%ra, RS=grid%rs, LAI=grid%lai                       &   ! P-X LSM
     &        ,NLCAT=grid%num_land_cat,  NSCAT=grid%num_soil_cat          &   ! P-X LSM
     &        ,VEGF_PX=grid%vegf_px, SNOWNCV=grid%snowncv                 &   ! P-X LSM
     &        ,ANAL_INTERVAL=config_flags%auxinput9_interval_s+config_flags%auxinput9_interval_m*60  &   ! P-X LSM
     &        ,PXLSM_SMOIS_INIT=config_flags%pxlsm_smois_init             &   ! P-X LSM
     &        ,PXLSM_SOIL_NUDGE=config_flags%pxlsm_soil_nudge             &   ! P-X LSM
!------------------------------------------------------------------------------
           ! Optional PX LSM nudging
     &        ,t2_ndg_old=grid%t2_ndg_old                  &
     &        ,q2_ndg_old=grid%q2_ndg_old                  &
     &        ,t2_ndg_new=grid%t2_ndg_new                  &
     &        ,q2_ndg_new=grid%q2_ndg_new                  &
     &        ,sn_ndg_old=grid%sn_ndg_old                  &
     &        ,sn_ndg_new=grid%sn_ndg_new                  &
! for Noah-MP LSM
     &        ,idveg=config_flags%dveg,         iopt_crs=config_flags%opt_crs &
     &        ,iopt_btr=config_flags%opt_btr,   iopt_run=config_flags%opt_run &
     &        ,iopt_sfc=config_flags%opt_sfc,   iopt_frz=config_flags%opt_frz &
     &        ,iopt_inf=config_flags%opt_inf,   iopt_rad=config_flags%opt_rad &
     &        ,iopt_alb=config_flags%opt_alb,   iopt_snf=config_flags%opt_snf &
     &        ,iopt_tbot=config_flags%opt_tbot, iopt_stc=config_flags%opt_stc &
           ! Noah UA changes
     &        ,ua_phys=config_flags%ua_phys,flx4=grid%flx4,fvb=grid%fvb       &
     &        ,fbur=grid%fbur,fgsn=grid%fgsn                                  &
           ! Indexes
     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
     &        , I_START=grid%i_start,I_END=min(grid%i_end, ide-1)         &
     &        , J_START=grid%j_start,J_END=min(grid%j_end, jde-1)         &
     &        , KTS=k_start, KTE=min(k_end,kde-1)                         &
     &        , NUM_TILES=grid%num_tiles                                  &
           ! Optional
     &        ,QV_CURR=moist(ims,kms,jms,P_QV)                            &
     &        ,a_QV_CURR=a_moist(ims,kms,jms,P_QV)                            &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC)                            &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR)                            &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI)                            &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS)                            &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG)                            &
     &        ,CAPG=grid%capg, EMISS=grid%emiss, HOL=hol,MOL=grid%mol,a_MOL=grid%a_mol      &
     &        ,T2OBS=grid%t2obs, Q2OBS=grid%q2obs                         &
     &        ,RAINBL=grid%rainbl,SR=grid%sr,RAINSHV=grid%rainshv         &
     &        ,RAINNCV=grid%rainncv,REGIME=grid%regime,T2=grid%t2,THC=grid%thc            &
     &        ,QSG=grid%qsg,QVG=grid%qvg,QCG=grid%qcg,SOILT1=grid%soilt1,TSNAV=grid%tsnav & ! ruc lsm
     &        ,SMFR3D=grid%smfr3d,KEEPFR3DFLAG=grid%keepfr3dflag,DEW=grid%dew             & ! ruc lsm
     &        ,POTEVP=grid%POTEVP, SNOPCX=grid%SNOPCX, SOILTB=grid%SOILTB                 & ! ruc lsm
     &        ,MOSAIC_LU=config_flags%mosaic_lu                           &   ! RUC LSM
     &        ,MOSAIC_SOIL=config_flags%mosaic_soil                       &   ! RUC LSM
     &        ,ISURBAN=isurban, MMINLU=TRIM(mminlu)                       &
     &        ,SNOTIME = grid%SNOTIME                                     &
     &        ,RDLAI2D=config_flags%rdlai2d                               &
     &        ,usemonalb=config_flags%usemonalb                           &
     &        ,NOAHRES=grid%noahres                                       &
     &        ,TSK_SAVE=grid%tsk_save                                     &
!mynn mp
     &         ,ch=grid%ch,tsq=grid%tsq,qsq=grid%qsq,cov=grid%cov &
     &                                                              )


BENCH_END(a_surf_driver_tim)
!      CALL POPREAL8ARRAY(grid%rqscuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rqicuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rqccuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rqvcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rthcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
!      CALL POPREAL8ARRAY(grid%rvcuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
!      CALL POPREAL8ARRAY(grid%rucuten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))

      CALL POPREAL8ARRAY(grid%rqiblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rqcblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rqvblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rthblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rvblten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))
      CALL POPREAL8ARRAY(grid%rublten, (ime-ims+1)*(kme-kms+1)*(jme-jms+1))


      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call a_phy_prep' )
        CALL a_phy_prep ( config_flags,                           &
                        grid%mut,grid%a_mut, grid%muu, grid%a_muu, grid%muv, grid%a_muv, &
                        grid%u_2,grid%a_u_2, grid%v_2,grid%a_v_2,             &
                        grid%p,grid%a_p, grid%pb, grid%alt,grid%a_alt,        &
                        grid%ph_2,grid%a_ph_2, grid%phb, grid%t_2,grid%a_t_2, &
                        grid%tsk, moist,a_moist, num_moist,        &
                        grid%rho,grid%a_rho, th_phy,a_th_phy, p_phy,a_p_phy, pi_phy,a_pi_phy, &
                        grid%u_phy,grid%a_u_phy, grid%v_phy,grid%a_v_phy,                         &
                        p8w,a_p8w, grid%t_phy,grid%a_t_phy, t8w,a_t8w, grid%z,grid%a_z, &
                        grid%z_at_w,grid%a_z_at_w, dz8w,a_dz8w,               &
                        grid%p_hyd,grid%a_p_hyd, grid%p_hyd_w,grid%a_p_hyd_w, &
                        grid%dnw,                                             &
                        grid%fnm, grid%fnp, grid%znw, grid%p_top,             &
                        grid%rthraten, grid%a_rthraten,                       &
                        grid%rthblten, grid%a_rthblten,                       &
                        grid%rublten, grid%a_rublten, grid%rvblten, grid%a_rvblten,            &
                        grid%rqvblten, grid%a_rqvblten, grid%rqcblten, grid%a_rqcblten, grid%rqiblten, grid%a_rqiblten,         &
                        grid%rucuten,  grid%a_rucuten , grid%rvcuten, grid%a_rvcuten,  grid%rthcuten, grid%a_rthcuten,    &
                        grid%rqvcuten, grid%a_rqvcuten, grid%rqccuten, grid%a_rqccuten, grid%rqrcuten, grid%a_rqrcuten,    &
                        grid%rqicuten, grid%a_rqicuten, grid%rqscuten, grid%a_rqscuten,                    &
                        grid%rushten,  grid%a_rushten, grid%rvshten,  grid%a_rvshten, grid%rthshten, grid%a_rthshten,    &
                        grid%rqvshten, grid%a_rqvshten, grid%rqcshten, grid%a_rqcshten, grid%rqrshten, grid%a_rqrshten,    &
                        grid%rqishten, grid%a_rqishten, grid%rqsshten, grid%a_rqsshten, grid%rqgshten, grid%a_rqgshten,    &
                        grid%rthften,  grid%a_rthften, grid%rqvften, grid%a_rqvften,                    &
                        grid%RUNDGDTEN, grid%a_RUNDGDTEN, grid%RVNDGDTEN, grid%a_RVNDGDTEN, grid%RTHNDGDTEN, grid%a_RTHNDGDTEN, &
                        grid%RPHNDGDTEN,grid%a_RPHNDGDTEN,grid%RQVNDGDTEN, grid%a_RQVNDGDTEN,grid%RMUNDGDTEN,&
!jdf
                        grid%landmask,grid%xland,                 &
!jdf
                        ids, ide, jds, jde, kds, kde,           &
                        ims, ime, jms, jme, kms, kme,           &
                        grid%i_start(ij), grid%i_end(ij),       &
                        grid%j_start(ij), grid%j_end(ij),       &
                        k_start, k_end                         )
      ENDDO
      !$OMP END PARALLEL DO

BENCH_END(a_phy_prep_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_PHYS_A_AD.inc"
#endif

BENCH_START(a_init_zero_tend_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call init_zero_tendency' )
         CALL a_init_zero_tendency(a_ru_tendf, &
                                   a_rv_tendf, &
                                   a_rw_tendf, &
                                   a_ph_tendf, &
                                   a_t_tendf,  &
                                   a_tke_tend, &
                                   a_mu_tendf,  &
                                   a_moist_tend,  &
                                   a_chem_tend,   &
                                   a_scalar_tend, &
                                   a_tracer_tend, &
                                   num_tracer,                   &
                                   num_moist,num_chem,num_scalar, &
                                   rk_step,                          &
                                   ids, ide, jds, jde, kds, kde,     &
                                   ims, ime, jms, jme, kms, kme,     &
                                   grid%i_start(ij), grid%i_end(ij), &
                                   grid%j_start(ij), grid%j_end(ij), &
                                   k_start, k_end                   )

       END DO
       !$OMP END PARALLEL DO
BENCH_END(a_init_zero_tend_tim)

  END SUBROUTINE a_first_rk_step_part1

END MODULE a_module_first_rk_step_part1
