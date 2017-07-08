!WRF:MEDIATION_LAYER:SOLVER

#define BENCH_START(A)
#define BENCH_END(A)

MODULE g_module_first_rk_step_part1

CONTAINS

  SUBROUTINE g_first_rk_step_part1 (   grid , config_flags              &
                             , moist , g_moist, moist_tend, g_moist_tend             &
                             , chem  , g_chem, chem_tend, g_chem_tend                &
#if (WRF_CHEM == 1)
                             , dvel &
#endif
                             , tracer, g_tracer, tracer_tend, g_tracer_tend          &
                             , scalar , g_scalar, scalar_tend, g_scalar_tend         &
                             , fdda3d, fdda2d                   &
                             , aerod                            &
                             , ru_tendf, g_ru_tendf, rv_tendf, g_rv_tendf               &
                             , rw_tendf, g_rw_tendf, t_tendf , g_t_tendf               &
                             , ph_tendf, g_ph_tendf, mu_tendf, g_mu_tendf               &
                             , tke_tend, g_tke_tend                         &
                             , adapt_step_flag , curr_secs      &
                             , psim , psih , wspd , g_wspd , gz1oz0 , chklowq &
                             , cu_act_flag , hol , th_phy , g_th_phy        &
                             , pi_phy, g_pi_phy, p_phy ,g_p_phy & !, t_phy ,g_t_phy  &
                             , dz8w , g_dz8w, p8w ,g_p8w, t8w ,g_t8w  &
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
   USE module_radiation_driver, ONLY : pre_radiation_driver, radiation_driver
!   USE module_surface_driver, ONLY : surface_driver
!   USE module_cumulus_driver, ONLY : cumulus_driver
!   USE module_shallowcu_driver, ONLY : shallowcu_driver
!   USE module_pbl_driver, ONLY : pbl_driver
!   USE module_fr_sfire_driver_wrf, ONLY : sfire_driver_em_step
!   USE module_fddagd_driver, ONLY : fddagd_driver
!   USE module_em, ONLY : init_zero_tendency
!   USE module_force_scm
   USE module_convtrans_prep
!   USE module_big_step_utilities_em, ONLY : phy_prep
!use module_scalar_tables
#ifdef DM_PARALLEL
    USE module_dm, ONLY : local_communicator, mytask, ntasks, ntasks_x, ntasks_y, local_communicator_periodic, wrf_dm_maxval
    USE module_comm_dm, ONLY : halo_em_phys_a_tl_sub,halo_em_phys_a_bl_surf_sub,halo_pwp_sub,halo_em_chem_e_3_tl_sub,halo_em_chem_e_5_tl_sub
#endif
#if (WRF_CHEM == 1)
    USE module_input_chem_data, only: get_last_gas,mozcart_lbc_set
    USE g_module_dry_dep_vel
#endif

    USE module_utility

    USE g_module_em, ONLY : g_init_zero_tendency
    USE g_module_big_step_utilities_em, ONLY : g_phy_prep
    USE g_module_pbl_driver, ONLY : g_pbl_driver
    USE g_module_cumulus_driver, ONLY : g_cumulus_driver
    USE g_module_surface_driver, ONLY : g_surface_driver

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
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: g_moist
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_moist),INTENT(INOUT)   :: g_moist_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: g_chem
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: chem_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_chem),INTENT(INOUT)   :: g_chem_tend
#if (WRF_CHEM == 1)
    REAL    ,DIMENSION(grid%sm31:grid%em31,1:grid%kdvel,grid%sm33:grid%em33,num_dvel) :: dvel
#endif
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: g_tracer
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_tracer),INTENT(INOUT)   :: g_tracer_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: g_scalar
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_scalar),INTENT(INOUT)   :: g_scalar_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_fdda3d),INTENT(INOUT)  :: fdda3d
    REAL    ,DIMENSION(ims:ime,1:1,jms:jme,num_fdda2d),INTENT(INOUT)      :: fdda2d
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme,num_aerod),INTENT(INOUT)   :: aerod
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psim
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: psih
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: g_wspd
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: gz1oz0
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: chklowq
    LOGICAL ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: cu_act_flag
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT)         :: hol

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_th_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_pi_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_p_phy
!    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_phy
!    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_t_phy
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_dz8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_p8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t8w
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_t8w

    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_ru_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_rv_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_rw_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_ph_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_t_tendf
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: tke_tend
    REAL    ,DIMENSION(ims:ime,kms:kme,jms:jme), INTENT(INOUT) :: g_tke_tend

    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: mu_tendf
    REAL    ,DIMENSION(ims:ime,jms:jme), INTENT(INOUT) :: g_mu_tendf

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

#if ( WRF_DFI_RADAR == 1 )
    INTEGER  do_capsupress   ! =1 do CAP supress, other = don't
#endif

  CALL get_ijk_from_subgrid (  grid ,                   &
                            sids, side, sjds, sjde, skds, skde,    &
                            sims, sime, sjms, sjme, skms, skme,    &
                            sips, sipe, sjps, sjpe, skps, skpe    )

 ! initialize all tendencies to zero in order to update physics
 ! tendencies first (separate from dry dynamics).

   l_flux=.FALSE.
   if (present(f_flux)) l_flux=f_flux

    rk_step = 1

BENCH_START(g_init_zero_tend_tim)
       !$OMP PARALLEL DO   &
       !$OMP PRIVATE ( ij )

       DO ij = 1 , grid%num_tiles

         CALL wrf_debug ( 200 , ' call init_zero_tendency' )
         CALL g_init_zero_tendency(ru_tendf,  g_ru_tendf, &
                                   rv_tendf,  g_rv_tendf, &
                                   rw_tendf,  g_rw_tendf, &
                                   ph_tendf,  g_ph_tendf, &
                                   t_tendf,   g_t_tendf,  &
                                   tke_tend, g_tke_tend, &
                                   mu_tendf,  g_mu_tendf,  &
                                   moist_tend, g_moist_tend,  &
                                   chem_tend, g_chem_tend,  &
                                   scalar_tend, g_scalar_tend, &
                                   tracer_tend, g_tracer_tend, &
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
BENCH_END(g_init_zero_tend_tim)

#ifdef DM_PARALLEL
#     include "HALO_EM_PHYS_A_TL.inc"
#endif

      !$OMP PARALLEL DO   &
      !$OMP PRIVATE ( ij )
      DO ij = 1 , grid%num_tiles

        CALL wrf_debug ( 200 , ' call g_phy_prep' )
        CALL g_phy_prep ( config_flags,                           &
                        grid%mut,grid%g_mut, grid%muu, grid%g_muu, grid%muv, grid%g_muv, &
                        grid%u_2,grid%g_u_2, grid%v_2,grid%g_v_2,             &
                        grid%p,grid%g_p, grid%pb, grid%alt,grid%g_alt,        &
                        grid%ph_2,grid%g_ph_2, grid%phb, grid%t_2,grid%g_t_2, &
                        grid%tsk, moist,g_moist, num_moist,        &
                        grid%rho,grid%g_rho, th_phy,g_th_phy, p_phy,g_p_phy, pi_phy,g_pi_phy, &
                        grid%u_phy,grid%g_u_phy, grid%v_phy,grid%g_v_phy,                         &
                        p8w,g_p8w, grid%t_phy,grid%g_t_phy, t8w,g_t8w, grid%z,grid%g_z, &
                        grid%z_at_w,grid%g_z_at_w, dz8w,g_dz8w,               &
                        grid%p_hyd,grid%g_p_hyd, grid%p_hyd_w,grid%g_p_hyd_w, &
                        grid%dnw,                                             &
                        grid%fnm, grid%fnp, grid%znw, grid%p_top,             &
                        grid%rthraten, grid%g_rthraten,                       &
                        grid%rthblten, grid%g_rthblten,                       &
                        grid%rublten, grid%g_rublten, grid%rvblten, grid%g_rvblten,            &
                        grid%rqvblten, grid%g_rqvblten, grid%rqcblten, grid%g_rqcblten, grid%rqiblten, grid%g_rqiblten,         &
                        grid%rucuten,  grid%g_rucuten , grid%rvcuten, grid%g_rvcuten,  grid%rthcuten, grid%g_rthcuten,    &
                        grid%rqvcuten, grid%g_rqvcuten, grid%rqccuten, grid%g_rqccuten, grid%rqrcuten, grid%g_rqrcuten,    &
                        grid%rqicuten, grid%g_rqicuten, grid%rqscuten, grid%g_rqscuten,                    &
                        grid%rushten,  grid%g_rushten, grid%rvshten,  grid%g_rvshten, grid%rthshten, grid%g_rthshten,    &
                        grid%rqvshten, grid%g_rqvshten, grid%rqcshten, grid%g_rqcshten, grid%rqrshten, grid%g_rqrshten,    &
                        grid%rqishten, grid%g_rqishten, grid%rqsshten, grid%g_rqsshten, grid%rqgshten, grid%g_rqgshten,    &
                        grid%rthften,  grid%g_rthften, grid%rqvften, grid%g_rqvften,                    &
                        grid%RUNDGDTEN, grid%g_RUNDGDTEN, grid%RVNDGDTEN, grid%g_RVNDGDTEN, grid%RTHNDGDTEN, grid%g_RTHNDGDTEN, &
                        grid%RPHNDGDTEN,grid%g_RPHNDGDTEN,grid%RQVNDGDTEN, grid%g_RQVNDGDTEN,grid%RMUNDGDTEN,&
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

BENCH_END(g_phy_prep_tim)

! radiation
     CALL domain_clock_get( grid, current_time=currentTime, &
                            current_timestr=mesg )
     read(mesg(6:7),FMT='(I2)') current_month

     CALL WRFU_TimeGet( currentTime, YY=yr, dayOfYear=day, H=hr, M=minute, S=sec, rc=rc)
         IF( rc/= WRFU_SUCCESS)THEN
         CALL wrf_error_fatal('WRFU_TimeGet failed')
         ENDIF

! GSW and GLW (in future RTHRATEN) are needed for tangent calls to pbl and surface drivers (and module_first_rk_step_part2_tl).
! Two options:
! (1) call radiation_driver radt (seconds) or stepra (# timesteps) inline here, since timing is correct in forward integration
! (2) store/retrieve these variables from trajectory and do not overwrite here
! We are using (2) above since it corresponds to adjoint version of this code. The downside is that the radiation_driver outputs will be offset by one time step from the appropriate forward model variables, but the AD/TL will be symmetric.
!-JJ Guerrette 

!!! this driver is only needed to handle non-local shadowing effects
!!      CALL pre_radiation_driver ( grid, config_flags                        &
!!     &        ,itimestep=grid%itimestep, ra_call_offset=grid%ra_call_offset    &
!!     &        ,XLAT=grid%xlat, XLONG=grid%xlong, GMT=grid%gmt                  &
!!     &        ,julian=grid%julian, xtime=grid%xtime, RADT=grid%radt            &
!!     &        ,STEPRA=grid%stepra                                              &
!!     &        ,ht=grid%ht,dx=grid%dx,dy=grid%dy,sina=grid%sina,cosa=grid%cosa  &
!!     &        ,shadowmask=grid%shadowmask,slope_rad=config_flags%slope_rad     &
!!     &        ,topo_shading=config_flags%topo_shading                          &
!!     &        ,shadlen=config_flags%shadlen,ht_shad=grid%ht_shad,ht_loc=ht_loc &
!!     &        ,ht_shad_bxs=grid%ht_shad_bxs, ht_shad_bxe=grid%ht_shad_bxe      &
!!     &        ,ht_shad_bys=grid%ht_shad_bys, ht_shad_bye=grid%ht_shad_bye      &
!!     &        ,nested=config_flags%nested, min_ptchsz=grid%min_ptchsz          &
!!     &        ,spec_bdy_width=config_flags%spec_bdy_width                      &
!!            ! indexes
!!     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
!!     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
!!     &        ,IPS=ips,IPE=ipe, JPS=jps,JPE=jpe, KPS=kps,KPE=kpe          &
!!     &        ,i_start=grid%i_start,i_end=min(grid%i_end, ide-1)          &
!!     &        ,j_start=grid%j_start,j_end=min(grid%j_end, jde-1)          &
!!     &        ,kts=k_start, kte=min(k_end,kde-1)                          &
!!     &        ,num_tiles=grid%num_tiles                                   )
!!
!!BENCH_START(rad_driver_tim)
!!      !Placeholder until g_radiation_driver present
!!      grid%g_rthraten = 0.D0
!!
!!      CALL radiation_driver(                                                  &
!!     &         p_top=grid%p_top & !DJW 140312 added p_top for vertical nesting
!!     &        ,ACFRCV=grid%acfrcv      ,ACFRST=grid%acfrst      ,ALBEDO=grid%albedo  &
!!     &        ,CFRACH=grid%cfrach      ,CFRACL=grid%cfracl      ,CFRACM=grid%cfracm  &
!!     &        ,CUPPT=grid%cuppt        ,CZMEAN=grid%czmean      ,DT=grid%dt          &
!!     &        ,DZ8W=dz8w               ,EMISS=grid%emiss        ,GLW=grid%glw        &
!!     &        ,GMT=grid%gmt            ,GSW=grid%gsw            ,HBOT=grid%hbot      &
!!     &        ,HTOP=grid%htop          ,HBOTR=grid%hbotr        ,HTOPR=grid%htopr    &
!!     &        ,ICLOUD=config_flags%icloud                                            &
!!     &        ,ITIMESTEP=grid%itimestep,JULDAY=grid%julday      , JULIAN=grid%julian &
!!     &        ,JULYR=grid%julyr        ,LW_PHYSICS=config_flags%ra_lw_physics        &
!!     &        ,NCFRCV=grid%ncfrcv      ,NCFRST=grid%ncfrst      ,NPHS=1              &
!!     &        ,o3input=config_flags%o3input     ,O3rad=grid%o3rad                    &
!!     &        ,aer_opt=config_flags%aer_opt ,aerod=aerod(:,:,:,P_ocarbon:P_upperaer) &
!!     &        ,swint_opt=config_flags%swint_opt                                      &
!!     &        ,P8W=grid%p_hyd_w        ,P=grid%p_hyd            ,PI=pi_phy           &
!!     &        ,RADT=grid%radt          ,RA_CALL_OFFSET=grid%ra_call_offset           &
!!     &        ,RHO=grid%rho            ,RLWTOA=grid%rlwtoa                           &
!!     &        ,RSWTOA=grid%rswtoa      ,RTHRATEN=grid%rthraten                       &
!!     &        ,RTHRATENLW=grid%rthratenlw       ,RTHRATENSW=grid%rthratensw          &
!!     &        ,SNOW=grid%snow          ,STEPRA=grid%stepra      ,SWDOWN=grid%swdown  &
!!     &        ,SWDOWNC=grid%swdownc    ,SW_PHYSICS=config_flags%ra_sw_physics        &
!!     &        ,T8W=t8w                 ,T=grid%t_phy           ,TAUCLDC=grid%taucldc &
!!     &        ,TAUCLDI=grid%taucldi    ,TSK=grid%tsk            ,VEGFRA=grid%vegfra  &
!!     &        ,WARM_RAIN=grid%warm_rain ,XICE=grid%xice         ,XLAND=grid%xland    &
!!     &        ,XLAT=grid%xlat          ,XLONG=grid%xlong        ,YR=yr               &
!!           ! SSiB LSM radiation components (fds 06/2010)
!!     &        ,ALSWVISDIR=grid%alswvisdir ,ALSWVISDIF=grid%alswvisdif      &  !ssib 
!!     &        ,ALSWNIRDIR=grid%alswnirdir ,ALSWNIRDIF=grid%alswnirdif      &  !ssib
!!     &        ,SWVISDIR=grid%swvisdir ,SWVISDIF=grid%swvisdif              &  !ssib
!!     &        ,SWNIRDIR=grid%swnirdir ,SWNIRDIF=grid%swnirdif              &  !ssib
!!     &        ,SF_SURFACE_PHYSICS=config_flags%sf_surface_physics          &  !ssib
!!! WRF-solar and aerosol variables from jararias 2013/8 and 2013/11
!!     &       ,SWDDIR=grid%swddir,SWDDNI=grid%swddni,SWDDIF=grid%swddif                              & 
!!     &       ,Gx=grid%Gx,Bx=grid%Bx,gg=grid%gg,bb=grid%bb                                           &
!!     &       ,swdown_ref=grid%swdown_ref,swddir_ref=grid%swddir_ref                                 &
!!     &       ,coszen_ref=grid%coszen_ref                                                            &
!!     &       ,aer_type=config_flags%aer_type                                                        &
!!     &       ,aer_aod550_opt=config_flags%aer_aod550_opt,aer_aod550_val=config_flags%aer_aod550_val &
!!     &       ,aer_angexp_opt=config_flags%aer_angexp_opt,aer_angexp_val=config_flags%aer_angexp_val &
!!     &       ,aer_ssa_opt=config_flags%aer_ssa_opt,aer_ssa_val=config_flags%aer_ssa_val             &
!!     &       ,aer_asy_opt=config_flags%aer_asy_opt,aer_asy_val=config_flags%aer_asy_val             &
!!     &       ,aod5502d=grid%aod5502d,angexp2d=grid%angexp2d,aerssa2d=grid%aerssa2d                  &
!!     &       ,aerasy2d=grid%aerasy2d,aod5503d=grid%aod5503d                                         &
!!     &       ,taod5502d=grid%taod5502d,taod5503d=grid%taod5503d                                     & ! Trude
!!!Optional solar variables
!!     &        ,DECLINX=grid%declin ,SOLCONX=grid%solcon ,COSZEN=grid%coszen ,HRANG=grid%hrang    &
!!     &        , CEN_LAT=grid%cen_lat                                      &
!!     &        ,Z=grid%z                                                   &
!!     &        ,ALEVSIZ=grid%alevsiz, no_src_types=grid%no_src_types       &
!!     &        ,LEVSIZ=grid%levsiz, N_OZMIXM=num_ozmixm                    &
!!     &        ,N_AEROSOLC=num_aerosolc                                    &
!!     &        ,PAERLEV=grid%paerlev   ,ID=grid%id                         &
!!     &        ,CAM_ABS_DIM1=grid%cam_abs_dim1, CAM_ABS_DIM2=grid%cam_abs_dim2 &
!!     &        ,CAM_ABS_FREQ_S=grid%cam_abs_freq_s                         &
!!     &        ,XTIME=grid%xtime                                                &
!!              ,CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag       &
!!!BSINGH - For WRFCuP scheme
!!     &        ,CU_PHYSICS=config_flags%cu_physics                         & !CuP, wig 5-Oct-2006
!!     &        ,SHALLOWCU_FORCED_RA=config_flags%shallowcu_forced_ra       & !CuP, wig
!!     &        ,CUBOT=grid%cubot, CUTOP=grid%cutop                         & !CuP, wig 9-Oct-2006
!!     &        ,CLDFRA_CUP=grid%cldfra_cup                                 & !CuP, wig 1-Oct-2006
!!     &        ,SHALL=grid%shall                                           & !CuP, wig 4-Feb-2008
!!!BSINGH - ENDS
!!            ! indexes
!!     &        ,IDS=ids,IDE=ide, JDS=jds,JDE=jde, KDS=kds,KDE=kde          &
!!     &        ,IMS=ims,IME=ime, JMS=jms,JME=jme, KMS=kms,KME=kme          &
!!     &        ,i_start=grid%i_start,i_end=min(grid%i_end, ide-1)          &
!!     &        ,j_start=grid%j_start,j_end=min(grid%j_end, jde-1)          &
!!     &        ,kts=k_start, kte=min(k_end,kde-1)                          &
!!     &        ,num_tiles=grid%num_tiles                                   &
!!            ! Optional
!!!JJS 20101020 vvvvv
!!     &        , TLWDN=grid%tlwdn, TLWUP=grid%tlwup                        & ! goddard schemes
!!     &        , SLWDN=grid%slwdn, SLWUP=grid%slwup                        & ! goddard schemes
!!     &        , TSWDN=grid%tswdn, TSWUP=grid%tswup                        & ! goddard schemes
!!     &        , SSWDN=grid%sswdn, SSWUP=grid%sswup                        & ! goddard schemes
!!!JJS 20101020 ^^^^^
!!!ZCX+ cloud fraction for CLWRF
!!     &        , CLDT=grid%cldt, ZNU=grid%znu                              &
!!!ZCX-
!!     &        , CLDFRA=grid%cldfra, CLDFRA_MP_ALL=grid%cldfra_mp_all      &
!!     &        , LRADIUS=grid%LRADIUS,IRADIUS=grid%IRADIUS                 & !BSINGH(01/22/2014)
!!     &        , CLDFRA_DP=grid%cldfra_dp                                  & ! ckay for subgrid cloud
!!     &        , CLDFRA_SH=grid%cldfra_sh                                  &
!!     &        , icloud_bl=config_flags%icloud_bl                        & !JOE: subgrid BL clouds
!!     &        , qc_bl=grid%qc_bl,cldfra_bl=grid%cldfra_bl                 & !JOE: subgrid bl clouds
!!     &        , re_cloud=grid%re_cloud, re_ice=grid%re_ice, re_snow=grid%re_snow & ! G. Thompson
!!     &        , has_reqc=grid%has_reqc, has_reqi=grid%has_reqi, has_reqs=grid%has_reqs & ! G. Thompson
!!     &        , PB=grid%pb                                                &
!!     &        , F_ICE_PHY=grid%f_ice_phy,F_RAIN_PHY=grid%f_rain_phy       &
!!     &        , QV=moist(ims,kms,jms,P_QV), F_QV=F_QV                     &
!!     &        , QC=moist(ims,kms,jms,P_QC), F_QC=F_QC                     &
!!     &        , QR=moist(ims,kms,jms,P_QR), F_QR=F_QR                     &
!!     &        , QI=moist(ims,kms,jms,P_QI), F_QI=F_QI                     &
!!     &        , QS=moist(ims,kms,jms,P_QS), F_QS=F_QS                     &
!!     &        , QG=moist(ims,kms,jms,P_QG), F_QG=F_QG                     &
!!     &        , QNDROP=scalar(ims,kms,jms,P_QNDROP), F_QNDROP=F_QNDROP    &
!!     &        ,QNIFA=scalar(ims,kms,jms,P_QNIFA),F_QNIFA=F_QNIFA          & !Trude
!!     &        ,QNWFA=scalar(ims,kms,jms,P_QNWFA),F_QNWFA=F_QNWFA          & !Trude
!!     &        ,ACSWUPT=grid%acswupt    ,ACSWUPTC=grid%acswuptc            &
!!     &        ,ACSWDNT=grid%acswdnt    ,ACSWDNTC=grid%acswdntc            &
!!     &        ,ACSWUPB=grid%acswupb    ,ACSWUPBC=grid%acswupbc            &
!!     &        ,ACSWDNB=grid%acswdnb    ,ACSWDNBC=grid%acswdnbc            &
!!     &        ,ACLWUPT=grid%aclwupt    ,ACLWUPTC=grid%aclwuptc            &
!!     &        ,ACLWDNT=grid%aclwdnt    ,ACLWDNTC=grid%aclwdntc            &
!!     &        ,ACLWUPB=grid%aclwupb    ,ACLWUPBC=grid%aclwupbc            &
!!     &        ,ACLWDNB=grid%aclwdnb    ,ACLWDNBC=grid%aclwdnbc            &
!!     &        ,SWUPT=grid%swupt    ,SWUPTC=grid%swuptc                    &
!!     &        ,SWDNT=grid%swdnt    ,SWDNTC=grid%swdntc                    &
!!     &        ,SWUPB=grid%swupb    ,SWUPBC=grid%swupbc                    &
!!     &        ,SWDNB=grid%swdnb    ,SWDNBC=grid%swdnbc                    &
!!     &        ,LWUPT=grid%lwupt    ,LWUPTC=grid%lwuptc                    &
!!     &        ,LWDNT=grid%lwdnt    ,LWDNTC=grid%lwdntc                    &
!!     &        ,LWUPB=grid%lwupb    ,LWUPBC=grid%lwupbc                    &
!!     &        ,LWDNB=grid%lwdnb    ,LWDNBC=grid%lwdnbc                    &
!!     &        ,LWCF=grid%lwcf                                                  &
!!     &        ,SWCF=grid%swcf                                                  &
!!     &        ,OLR=grid%olr                                                    &
!!     &        ,AERODM=grid%aerodm, PINA=grid%pina, AODTOT=grid%aodtot          &
!!     &        ,OZMIXM=grid%ozmixm, PIN=grid%pin                                &
!!     &        ,M_PS_1=grid%m_ps_1, M_PS_2=grid%m_ps_2, AEROSOLC_1=grid%aerosolc_1        &
!!     &        ,AEROSOLC_2=grid%aerosolc_2, M_HYBI0=grid%m_hybi                      &
!!     &        ,ABSTOT=grid%abstot, ABSNXT=grid%absnxt, EMSTOT=grid%emstot                &
!!     &        ,RADTACTTIME=grid%radtacttime                                         &  
!!     &        ,ICLOUD_CU=config_flags%ICLOUD_CU                            &
!!     &        ,QC_CU=grid%QC_CU , QI_CU=grid%QI_CU                         &
!!#if (WRF_CHEM == 1)
!!     &        ,AER_RA_FEEDBACK=config_flags%aer_ra_feedback                &
!!     &        ,PM2_5_DRY=grid%pm2_5_dry, PM2_5_WATER=grid%pm2_5_water               &
!!     &        ,PM2_5_DRY_EC=grid%pm2_5_dry_ec                                  &
!!     &        ,TAUAER300=grid%tauaer1, TAUAER400=grid%tauaer2 & ! jcb
!!     &        ,TAUAER600=grid%tauaer3, TAUAER999=grid%tauaer4 & ! jcb
!!     &        ,GAER300=grid%gaer1, GAER400=grid%gaer2, GAER600=grid%gaer3, GAER999=grid%gaer4 & ! jcb
!!     &        ,WAER300=grid%waer1, WAER400=grid%waer2, WAER600=grid%waer3, WAER999=grid%waer4 & ! jcb
!!     &        ,TAUAERlw1  =grid%tauaerlw1,  TAUAERlw2=grid%tauaerlw2    & 
!!     &        ,TAUAERlw3  =grid%tauaerlw3,  TAUAERlw4=grid%tauaerlw4    & 
!!     &        ,TAUAERlw5  =grid%tauaerlw5,  TAUAERlw6=grid%tauaerlw6    & 
!!     &        ,TAUAERlw7  =grid%tauaerlw7,  TAUAERlw8=grid%tauaerlw8    & 
!!     &        ,TAUAERlw9  =grid%tauaerlw9,  TAUAERlw10=grid%tauaerlw10    & 
!!     &        ,TAUAERlw11 =grid%tauaerlw11, TAUAERlw12=grid%tauaerlw12    & 
!!     &        ,TAUAERlw13 =grid%tauaerlw13, TAUAERlw14=grid%tauaerlw14    & 
!!     &        ,TAUAERlw15 =grid%tauaerlw15, TAUAERlw16=grid%tauaerlw16   & 
!!     &        ,progn=config_flags%progn                                            &
!!#endif
!!     &         ,slope_rad=config_flags%slope_rad,topo_shading=config_flags%topo_shading     &
!!     &         ,shadowmask=grid%shadowmask,ht=grid%ht,dx=grid%dx,dy=grid%dy &           
!!     &         ,diffuse_frac=grid%diffuse_frac &           
!!     &         ,IS_CAMMGMP_USED = grid%is_CAMMGMP_used    &
!!     &         ,MP_PHYSICS=CONFIG_FLAGS%MP_PHYSICS )
!!
!!BENCH_END(rad_driver_tim)

BENCH_START(surf_driver_tim)
    
! !gmm halo of wtd and riverflow for leafhydro
!#ifdef DM_PARALLEL
!  IF ( config_flags%sf_surface_physics.eq.NOAHMPSCHEME ) THEN
!       IF ( config_flags%opt_run.eq.5.and.mod(grid%itimestep,grid%STEPWTD).eq.0 )  THEN
!#     include "HALO_EM_HYDRO_NOAHMP.inc"
!       ENDIF
!  ENDIF
!#endif

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

CALL g_surface_driver(                                                &
     &        HYDRO_dt=HYDRO_dt,  sfcheadrt=grid%sfcheadrt,               &
     &        INFXSRT=grid%INFXSRT,  soldrain=grid%soldrain,              &
     &         ACGRDFLX=grid%acgrdflx  ,ACHFX=grid%achfx        ,ACLHF=grid%aclhf        &
     &        ,ACSNOM=grid%acsnom      ,ACSNOW=grid%acsnow      ,AKHS=grid%akhs          &
     &        ,AKMS=grid%akms          ,ALBBCK=grid%albbck      ,ALBEDO=grid%albedo      &
     &        ,EMBCK=grid%embck                                                          &
     &        ,BR=grid%br              ,CANWAT=grid%canwat,g_CANWAT=grid%g_canwat      ,CHKLOWQ=chklowq    &
     &        ,CT=grid%ct              ,DT=grid%dt         ,DX=grid%dx         &
     &        ,DZ8W=dz8w     ,g_DZ8W=g_dz8w    ,DZS=grid%dzs            ,FLHC=grid%flhc          &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &
     &        ,FLQC=grid%flqc          ,GLW=grid%glw            ,GRDFLX=grid%grdflx      &
     &        ,GSW=grid%gsw    ,SWDOWN=grid%swdown        ,GZ1OZ0=gz1oz0     &
     &        ,HFX=grid%hfx            ,g_HFX=grid%g_hfx              &
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
     &        ,P8W=grid%p_hyd_w        ,g_P8W=grid%g_p_hyd_w   &
     &        ,PBLH=grid%pblh          ,g_PBLH=grid%g_pblh    ,PI_PHY=pi_phy  &
     &        ,PSFC=grid%psfc          ,PSHLTR=grid%pshltr      ,PSIH=psih          &
     &        ,BLDT=grid%bldt     ,CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag  &
     &        ,BLDTACTTIME=grid%bldtacttime                                              &
     &        ,PSIM=psim               ,P_PHY=grid%p_hyd        ,g_P_PHY=grid%g_p_hyd    &
     &        ,Q10=grid%q10            &
     &        ,Q2=grid%q2              ,QFX=grid%qfx    ,g_QFX=grid%g_qfx  ,QSFC=grid%qsfc &
     &        ,QSHLTR=grid%qshltr      ,QZ0=grid%qz0            ,RAINCV=grid%raincv      &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics     ,RHO=grid%rho   ,g_RHO=grid%g_rho      &
     &        ,RMOL=grid%rmol ,g_RMOL=grid%g_rmol ,SFCEVP=grid%sfcevp   ,SFCEXC=grid%sfcexc  &
     &        ,SFCRUNOFF=grid%sfcrunoff                                        &
     &        ,SF_SFCLAY_PHYSICS=config_flags%sf_sfclay_physics                        &
     &        ,SF_SURFACE_PHYSICS=config_flags%sf_surface_physics  ,SH2O=grid%sh2o          &
     &        ,SHDMAX=grid%shdmax      ,SHDMIN=grid%shdmin      &
     &        ,SMOIS=grid%smois        ,g_SMOIS=grid%g_smois        &
     &        ,SMSTAV=grid%smstav      ,SMSTOT=grid%smstot      ,SNOALB=grid%snoalb      &
     &        ,SNOW=grid%snow          ,SNOWC=grid%snowc        ,SNOWH=grid%snowh        &
     &        ,SMCREL=grid%smcrel                                                        &
     &        ,SST=grid%sst            ,SST_UPDATE=grid%sst_update                  &
     &        ,SSTSK=grid%sstsk        ,DTW=grid%dtw            ,SST_SKIN=grid%sst_skin  &
     &        ,SCM_FORCE_SKINTEMP=grid%scm_force_skintemp                                &
     &        ,SCM_FORCE_FLUX=grid%scm_force_flux                                        &
     &        ,STEPBL=grid%stepbl      ,TH10=grid%th10          ,TH2=grid%th2            &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy  ,g_TH_PHY=g_th_phy                &
     &        ,TKE_PBL=grid%tke_pbl    &
     &        ,TMN=grid%tmn            ,TSHLTR=grid%tshltr      ,TSK=grid%tsk            &
     &        ,TYR=grid%tyr            ,TYRA=grid%tyra          ,TDLY=grid%tdly          &
     &        ,TLAG=grid%tlag          ,LAGDAY=config_flags%lagday      ,NYEAR=grid%nyear        &
     &        ,NDAY=grid%nday          ,TMN_UPDATE=grid%tmn_update      ,YR=yr            &
     &        ,TSLB=grid%tslb  ,g_TSLB=grid%g_tslb ,T_PHY=grid%t_phy  ,g_T_PHY=grid%g_t_phy   &
     &        ,U10=grid%u10            &
     &        ,URATX=grid%uratx        ,VRATX=grid%vratx   ,TRATX=grid%tratx        &
     &        ,UDRUNOFF=grid%udrunoff  ,UST=grid%ust ,g_UST=grid%g_ust  ,UZ0=grid%uz0    &
     &        ,U_FRAME=grid%u_frame    ,U_PHY=grid%u_phy  ,g_U_PHY=grid%g_u_phy   ,V10=grid%v10    &
     &        ,VEGFRA=grid%vegfra      ,VZ0=grid%vz0       ,V_FRAME=grid%v_frame    &
     &        ,V_PHY=grid%v_phy        ,g_V_PHY=grid%g_v_phy    ,WARM_RAIN=grid%warm_rain     &
     &        ,WSPD=wspd      ,g_WSPD=g_wspd         ,XICE=grid%xice     ,XLAND=grid%xland        &
     &        ,Z0=grid%z0              ,Z=grid%z        ,ZNT=grid%znt, g_ZNT=grid%g_znt &
     &        ,ZS=grid%zs              ,ALBSI=grid%albsi , ICEDEPTH=grid%icedepth   &
     &        ,SNOWSI=grid%snowsi                                                   &
     &        ,XICEM=grid%xicem   ,ISICE=grid%landuse_isice                         &
     &        ,USTM=grid%ustm     ,g_USTM=grid%g_ustm      ,CK=grid%ck    ,CKA=grid%cka   &
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
     &        ,g_QV_CURR=g_moist(ims,kms,jms,P_QV)                            &
     &        ,QC_CURR=moist(ims,kms,jms,P_QC)                            &
     &        ,QR_CURR=moist(ims,kms,jms,P_QR)                            &
     &        ,QI_CURR=moist(ims,kms,jms,P_QI)                            &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS)                            &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG)                            &
     &        ,CAPG=grid%capg, EMISS=grid%emiss, HOL=hol,MOL=grid%mol,g_MOL=grid%g_mol      &
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

BENCH_END(g_surf_driver_tim)

!*********
! pbl


      CALL wrf_debug ( 200 , ' call g_pbl_driver' )
BENCH_START(g_pbl_driver_tim)
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
               call g_dry_dep_vel(grid%id,curr_secs,grid%itimestep,grid%dt,config_flags,          &
                    grid%gmt,ifix(grid%julian),grid%t_phy,grid%g_t_phy,moist,g_moist,p8w,t8w,                &
                    grid%alt,p_phy,g_p_phy,chem,grid%rho,dz8w,grid%hfx,grid%dx,      & 
                    grid%raincv_b,            &
                    grid%ivgtyp,grid%tsk,grid%gsw,grid%vegfra,grid%pblh,grid%g_pblh,   &
                    grid%rmol,grid%g_rmol,grid%ust,grid%g_ust,grid%znt,grid%g_znt,     &
                    grid%xlat,grid%xlong,                    &
                    grid%z,grid%z_at_w,grid%g_z_at_w,grid%xland,                                &
                    grid%h2oaj,grid%h2oai,grid%nu3,grid%ac3,grid%cor3,grid%asulf,        &
                    grid%ahno3,grid%anh3,grid%cvaro1,grid%cvaro2,grid%cvalk1,grid%cvole1,&
                    grid%cvapi1,grid%cvapi2,grid%cvlim1,grid%cvlim2,     &
                    numgas,current_month,dvel,grid%snowh,  &
                    grid%dustdrydep_1,grid%dustdrydep_2,grid%dustdrydep_3,               &
                    grid%dustdrydep_4,grid%dustdrydep_5,                                 &
                    grid%depvelocity,                                                    &
                    grid%is_CAMMGMP_used,                  & !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
                    grid%dep_vel,grid%g_dep_vel,grid%num_vert_mix,                      &
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
!                                        k_start, k_end                   )
!               end if
            END DO
            !$OMP END PARALLEL DO
         else
            grid%dep_vel(:,:,:,:) = 0.D0
            grid%g_dep_vel(:,:,:,:) = 0.D0
         end if
      else
         grid%dep_vel(:,:,:,:) = 0.D0
         grid%g_dep_vel(:,:,:,:) = 0.D0
      endif
#endif

#ifdef DM_PARALLEL
      IF ( config_flags%bl_pbl_physics .EQ. 98 ) THEN
#     include "HALO_EM_PHYS_A_BL_SURF.inc"
      ENDIF
#endif

      CALL g_pbl_driver(                                                              &
     &         AKHS=grid%akhs          ,AKMS=grid%akms                              &
     &        ,BL_PBL_PHYSICS=config_flags%bl_pbl_physics                           &
     &        ,BLDT=grid%bldt, CURR_SECS=curr_secs, ADAPT_STEP_FLAG=adapt_step_flag &
     &        ,BLDTACTTIME=grid%bldtacttime                                         & 
     &        ,BR=grid%br              ,CHKLOWQ=chklowq    ,CT=grid%ct                   &
     &        ,DT=grid%dt              ,DX=grid%dx    ,DZ8W=dz8w, G_DZ8W=g_dz8w     &
     &        ,EXCH_H=grid%exch_h     ,EXCH_M=grid%exch_m                           &
     &        ,FM=grid%fm         ,FHH=grid%fh                                      &
     &        ,F=grid%f  ,GRDFLX=grid%grdflx      &
     &        ,GZ1OZ0=gz1oz0      ,HFX=grid%hfx ,g_HFX=grid%g_hfx  ,HT=grid%ht      &
     &        ,ID=grid%id         ,ITIMESTEP=grid%itimestep     ,KPBL=grid%kpbl          &
     &        ,LH=grid%lh ,LOWLYR=grid%lowlyr ,P8W=grid%p_hyd_w, g_P8W=grid%g_p_hyd_w   &
     &        ,PBLH=grid%pblh     ,G_PBLH=grid%g_pblh ,PI_PHY=pi_phy ,g_PI_PHY=g_pi_phy ,PSIH=psih &
     &        ,PSIM=psim   ,P_PHY=grid%p_hyd, QFX=grid%qfx ,g_QFX=grid%g_qfx        &
     &        ,QSFC=grid%qsfc          ,QZ0=grid%qz0  ,MIXHT=mixht                  &
     &        ,RA_LW_PHYSICS=config_flags%ra_lw_physics                   &
     &        ,RHO=grid%rho            ,RQCBLTEN=grid%rqcblten  ,RQIBLTEN=grid%rqiblten  &
     &        ,G_RHO=grid%g_rho        ,G_RQCBLTEN=grid%g_rqcblten  ,G_RQIBLTEN=grid%g_rqiblten  &
     &        ,RQVBLTEN=grid%rqvblten  ,g_RQVBLTEN=grid%g_rqvblten                  &
     &        ,RTHBLTEN=grid%rthblten ,g_RTHBLTEN=grid%g_rthblten                   &
     &        ,RUBLTEN=grid%rublten, g_RUBLTEN=grid%g_rublten                 &
     &        ,RVBLTEN=grid%rvblten ,g_RVBLTEN=grid%g_rvblten   ,SNOW=grid%snow ,STEPBL=grid%stepbl &
     &        ,THZ0=grid%thz0          ,TH_PHY=th_phy     ,G_TH_PHY=g_th_phy  &
     &        ,TSK=grid%tsk            ,T_PHY=grid%t_phy  ,G_T_PHY=grid%g_t_phy         &
     &        ,UST=grid%ust  ,g_UST=grid%g_ust                                &
     &        ,U10=grid%u10 ,UZ0=grid%uz0 ,U_FRAME=grid%u_frame ,U_PHY=grid%u_phy ,g_U_PHY=grid%g_u_phy &
     &        ,V10=grid%v10 ,VZ0=grid%vz0 ,V_FRAME=grid%v_frame ,V_PHY=grid%v_phy ,g_V_PHY=grid%g_v_phy &
     &        ,T2=grid%t2 &
     &        ,WARM_RAIN=grid%warm_rain ,WSPD=wspd        ,g_WSPD=g_wspd        &
     &        ,XICE=grid%xice   ,XLAND=grid%xland  ,Z=grid%z ,g_Z=grid%g_z      &
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
     &        ,ZNU=grid%znu,ZNW=grid%znw,MUT=grid%mut,g_MUT=grid%g_mut,P_TOP=grid%p_top    &
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
     &        ,g_QV_CURR=g_moist(ims,kms,jms,P_QV)                        &
     &        ,g_QC_CURR=g_moist(ims,kms,jms,P_QC)                        &
     &        ,g_QI_CURR=g_moist(ims,kms,jms,P_QI)                        &
     &        ,QS_CURR=moist(ims,kms,jms,P_QS), F_QS=F_QS                 &
     &        ,QG_CURR=moist(ims,kms,jms,P_QG), F_QG=F_QG                 &
     &        ,HOL=HOL, MOL=grid%mol, REGIME=grid%REGIME                  &
!mynn mp
     &        ,QKE=scalar(ims,kms,jms,P_qke_adv), tsq=grid%tsq, qsq=grid%qsq  &
     &        ,cov=grid%cov                                               &
     &        ,DQKE=grid%dqke,QWT=grid%qWT                                &
     &        ,QSHEAR=grid%qSHEAR,QBUOY=grid%qBUOY,QDISS=grid%qDISS       &
     &        ,bl_mynn_tkebudget=config_flags%bl_mynn_tkebudget           &
     &        ,rmol=grid%rmol, ch=grid%ch                                 &
     &        ,qcg=grid%qcg, grav_settling=config_flags%grav_settling             &
!    &        ,K_m=grid%K_m, K_h=grid%K_h, K_q=grid%K_q                   &
!GWD for ARW
     &        ,GWD_OPT=config_flags%gwd_opt &
     &        ,DTAUX3D=grid%dtaux3d,DTAUY3D=grid%dtauy3d &
     &        ,G_DTAUX3D=grid%g_dtaux3d,G_DTAUY3D=grid%g_dtauy3d &
     &        ,DUSFCG=grid%dusfcg,DVSFCG=grid%dvsfcg &
     &        ,G_DUSFCG=grid%g_dusfcg,G_DVSFCG=grid%g_dvsfcg &
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
     &        ,CHEM=chem  ,G_CHEM=g_chem  ,VD=grid%dep_vel ,G_VD=grid%g_dep_vel   &
     &        ,NCHEM=num_chem,kdvel=config_flags%kdepvel                 &
     &        ,ndvel=config_flags%ndepvel                                &
     &        ,NUM_VERT_MIX=grid%num_vert_mix                     &
#endif
     &        )
                                                         

BENCH_END(g_pbl_driver_tim)

#if (WRF_CHEM == 1)
# ifdef DM_PARALLEL
   IF ( num_chem >= PARAM_FIRST_SCALAR .AND. config_flags%bl_pbl_physics == ACMPBLSCHEME) then
     CALL wrf_debug ( 200 , ' call HALO CHEM AFTER PBL' )
     IF      ( config_flags%h_mom_adv_order <= 4 ) THEN
#      include "HALO_EM_CHEM_E_3_TL.inc"
     ELSE IF ( config_flags%h_mom_adv_order <= 6 ) THEN
#      include "HALO_EM_CHEM_E_5_TL.inc"

     ELSE
       WRITE(message,*)'solve_em: invalid h_mom_adv_order = ',config_flags%h_mom_adv_order
       CALL wrf_error_fatal(TRIM(message))
     ENDIF
   ENDIF
# endif
# endif


! cumulus para.

      CALL wrf_debug ( 200 , ' call g_cumulus_driver' )

BENCH_START(g_cu_driver_tim)
      CALL g_cumulus_driver(grid                                             &
                 ! Prognostic variables
!     &             ,U=grid%u_phy, V=grid%v_phy, TH=th_phy, g_TH=g_th_phy, T=grid%t_phy, g_T=grid%g_t_phy  &
!     &             ,W=grid%w_2, P=grid%p_hyd, g_P=grid%g_p_hyd, PI=pi_phy, RHO=grid%rho, g_RHO=grid%g_rho &
     &             ,U=grid%u_phy, V=grid%v_phy, TH=th_phy, THD=g_th_phy, T=grid%t_phy, TD=grid%g_t_phy  &
     &             ,W=grid%w_2, P=grid%p_hyd, PD=grid%g_p_hyd, PI=pi_phy, RHO=grid%rho, RHOD=grid%g_rho &
                 ! Other arguments
     &             ,ITIMESTEP=grid%itimestep ,DT=grid%dt      ,DX=grid%dx              &
     &             ,CUDT=grid%cudt,CURR_SECS=curr_secs,ADAPT_STEP_FLAG=adapt_step_flag &
     &             ,CUDTACTTIME=grid%cudtacttime                                       &
     &             ,RAINC=grid%rainc   ,RAINCV=grid%raincv   ,PRATEC=grid%pratec       &
     &             ,RAINCVD=grid%g_raincv, PRATECD=grid%g_pratec                       &
     &             ,NCA=grid%nca                                                       &
     &             ,HTOP=grid%cutop     ,HBOT=grid%cubot       ,KPBL=grid%kpbl         &
     &             ,Z=grid%z ,Z_AT_W=grid%z_at_w ,MAVAIL=grid%mavail ,PBLH=grid%pblh   &
     &             ,DZ8W=dz8w     ,P8W=grid%p_hyd_w, PSFC=grid%psfc, TSK=grid%tsk      &
     &             ,ZD=grid%g_z, DZ8WD=g_dz8w                                        &
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
     &             ,RQVCUTEND=grid%g_rqvcuten                            &
     &             ,RQSCUTEN=grid%rqscuten , RQICUTEN=grid%rqicuten       &
     &             ,RQRCUTEN=grid%rqrcuten                                &
     &             ,RQVBLTEN=grid%rqvblten , RQVFTEN=grid%rqvften         &
                 ! Other tendency arguments
     &             ,RTHRATEN=grid%rthraten , RTHBLTEN=grid%rthblten       &
     &             ,RUCUTEN=grid%rucuten   , RVCUTEN=grid%rvcuten         &
     &             ,RTHCUTEN=grid%rthcuten , RTHFTEN=grid%rthften         &
     &             ,RTHCUTEND=grid%g_rthcuten                            &
                 ! Moisture tracer arguments
     &             ,QV_CURR=moist(ims,kms,jms,P_QV), F_QV=F_QV            &
     &             ,QV_CURRD=g_moist(ims,kms,jms,P_QV)                   &
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
BENCH_END(g_cu_driver_tim)
!

  END SUBROUTINE g_first_rk_step_part1

END MODULE g_module_first_rk_step_part1
