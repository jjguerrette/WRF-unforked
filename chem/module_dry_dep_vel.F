!
!   WRF-chem V3.2 : Original version of dry_dep_driver written by Georg Grell (ESRL/GSD)
!                   VERTMX was originally developed by Mariusz Pagowski and modified by
!                   Richard C. Easter (PNNL)
! 10/11/2011 - Ravan Ahmadov (NOAA) updated to include the RACM_SOA_VBS option
!
! 10/14/2013 - Jonathan Guerrette separated dep velocity calc from dry_dep_driver
!               only for ACM2 PBL scheme
!
!WRF:MODEL_LAYER:CHEMICS
!
MODULE module_dry_dep_vel
  IMPLICIT NONE

CONTAINS

    subroutine dry_dep_vel(id,curr_secs,ktau,dtstep,config_flags,      &
               gmt,julday,t_phy,moist,p8w,t8w,alt,               &
               p_phy,chem,rho_phy,dz8w,hfx,dx,           &  
               raincv,                   &
               ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
               xland,h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,       &
               anh3,cvaro1,cvaro2,   &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,      &
               numgas,current_month,dvel,snowh,                &
               dustdrydep_1,dustdrydep_2,dustdrydep_3,                    &
               dustdrydep_4,dustdrydep_5,                                 &
               depvelocity,                                               &
               is_CAMMGMP_used,                                     & !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
               dep_vel,num_vert_mix,                                      &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!----------------------------------------------------------------------
  USE module_model_constants
  USE module_configure
  USE module_state_description
  USE module_domain_type, only : domain
  USE module_dep_simple
!  USE module_vertmx_wrf
  USE module_data_sorgam
  USE module_aerosols_sorgam
!  USE module_gocart_settling
!  USE module_vash_settling
  USE module_gocart_drydep
  USE module_mosaic_drydep, only:  mosaic_drydep_driver
!  USE module_mixactivate_wrappers, only: mosaic_mixactivate, sorgam_mixactivate
  USE module_aer_drydep
  USE module_aerosols_soa_vbs, only: soa_vbs_depdriver  
!!  use module_cam_support, only: pcnst => pcnst_runtime
!!  USE module_data_cam_mam_asect, only: lptr_chem_to_q, lptr_chem_to_qqcw !Balwinder.Singh@pnnl.gov: Added to avoid mixing of CHEM array constituents multiple times
  USE modal_aero_data,         only: numptr_amode, lmassptr_amode, ntot_amode, nspec_amode !Added by Balwinder.Singh@pnnl.gov to avoid mixing of CHEM array constituents multiple times
  USE module_cam_mam_drydep, only:  cam_mam_drydep_driver
  use module_scalar_tables,     only: chem_dname_table !Balwinder.Singh@pnnl.gov:Added for MAM aerosols dry deposition
  USE module_aerosols_sorgam_vbs, only: sorgam_vbs_depdriver

  IMPLICIT NONE

   TYPE(grid_config_rec_type),  INTENT(IN   )    :: config_flags
   LOGICAL, INTENT(IN)                           :: is_CAMMGMP_used !BSINGH:01/31/2013: Added is_CAMMGMP_used for MAM drydep
   INTEGER,      INTENT(IN   ) :: id,julday,                    &
                                  numgas,                       &
                                  current_month,                &
                                  ids,ide, jds,jde, kds,kde,    &
                                  ims,ime, jms,jme, kms,kme,    &
                                  its,ite, jts,jte, kts,kte
   INTEGER,      INTENT(IN   ) :: ktau
   integer l
   REAL(KIND=8), INTENT(IN   ) :: curr_secs
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_moist ),        &
         INTENT(IN ) ::                                   moist
   REAL, DIMENSION( ims:ime, kms:kme, jms:jme, num_chem ),         &
         INTENT(INOUT ) ::                                 chem

   REAL, DIMENSION( ims:ime, config_flags%kdepvel, jms:jme, config_flags%ndepvel ), &
         INTENT(INOUT ) ::                                   dep_vel

   REAL, DIMENSION( ims:ime, config_flags%kdvel, jms:jme, num_dvel ), &
         INTENT(INOUT ) ::                                   dvel


   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(IN   ) ::                                      &
                                                      t_phy,    &
                                                        alt,    &
                                                      p_phy,    &
                                                      dz8w,     &
                                              t8w,p8w,z_at_w ,  &
                                              rho_phy,z
   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )         ,    &
          INTENT(INOUT) ::                                      &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2,    &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2
   INTEGER,DIMENSION( ims:ime , jms:jme )                  ,    &
          INTENT(IN   ) ::                                      &
                                                     ivgtyp
   REAL,  DIMENSION( ims:ime , jms:jme )                   ,    &
          INTENT(INOUT) ::                                      &
                                                     tsk,       &
                                                     gsw,       &
                                                  vegfra,       &
                                                     pbl,       &
                                                     rmol,       &
                                                     ust,       &
                                                     hfx,       &
                                                     xlat,      &
                                                     xlong,     &
                                                     snowh,     &
                                          xland,znt,raincv

   REAL,  DIMENSION( ims:ime , jms:jme , num_chem )         ::  &
                    ddflx    !dry deposition flux

      REAL,      INTENT(IN   ) ::                               &
                             dtstep,gmt,dx

      INTEGER,  INTENT(INOUT) :: num_vert_mix

   REAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) ::        &
                dustdrydep_1, dustdrydep_2, dustdrydep_3,       &
                dustdrydep_4, dustdrydep_5, depvelocity
!--- deposition and emissions stuff
! .. Parameters ..
! ..
! .. Local Scalars ..
      REAL ::  clwchem,  dvfog, dvpart,  &
        rad, rhchem, ta, ustar, vegfrac, z1,zntt

      INTEGER :: iland, iprt, iseason, jce, jcs,  &
                 n, nr, ipr, jpr, nvr,   &
                 idrydep_onoff, aer_mech_id
      INTEGER :: l2,m,lnum,lmass

      LOGICAL :: highnh3, rainflag, vegflag, wetflag
! ..
! .. Local Arrays ..
      REAL :: p(kts:kte)
   REAL, DIMENSION( its:ite, jts:jte, num_chem ) ::   ddvel

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme ) :: dryrho_phy
   REAL,  DIMENSION( kms:kme ) :: dryrho_1d

! turbulent transport
   real :: pblst(kts:kte),ekmfull(kts:kte+1),zzfull(kts:kte+1),zz(kts:kte)
   integer :: ii,jj,kk,i,j,k,nv
   integer :: ll   
!
! necessary for aerosols (module dependent)
!
   REAL, DIMENSION( its:ite, jts:jte ) ::   aer_res, aer_res_def, aer_res_zcen

!!   INTEGER, DIMENSION( pcnst )    ::   lptr_q_to_chem !Balwinder.Singh@pnnl.gov: pointer from Q array to CHEM array
   LOGICAL, DIMENSION( num_chem ) ::   vertMixAero    !Balwinder.Singh@pnnl.gov: Flag (true or false) which decides whether to  vertical mix a constituent or not

   real, parameter :: m2cm = 100.

   REAL RSI                   ! gas constant in SI units (J/mol-K)
   PARAMETER (RSI = 8.314510) ! RSI is the same as RGASUNIV in CONST.EXT


   integer :: k_a, k_c, kmax, m_mam
   real, dimension( its:ite, jts:jte ) :: frac_removed

   REAL,  DIMENSION( ims:ime , kms:kme , jms:jme )  :: rh


! ..
! .. Intrinsic Functions ..
      INTRINSIC max, min

      !Balwinder.Singh@pnnl.gov: Following code sets vertMixAero(:) flag to true or false based on whether
      !a particular CHEM array constituent should be mixed (vertically) by dry_dep_driver or not
      !The constituents which are  mixed in CAMMGMP microphysics are NOT mixed in dry_dep_driver. 
!!      if(is_CAMMGMP_used) then ! If CAMMGMP microphysics is selected
!!         vertMixAero(:)    = .FALSE. 
!!         lptr_q_to_chem(:) = -999888777 ! pointer which maps Q array (CAM's STATE%Q array) to CHEM array
!!         
!!         do nv = 2, num_chem
!!            l2 = lptr_chem_to_q(nv)
!!            if (l2 >= 0) then 
!!               vertMixAero(nv) = .TRUE. !Set default true
!!               lptr_q_to_chem(l2) = nv
!!            end if
!!         enddo
!!         
!!         !Set constitutents, which are mixed in CAMMGMP, False
!!         do m = 1, ntot_amode
!!            lnum = numptr_amode(m)
!!            if( lnum > 0 ) then
!!               vertMixAero(lptr_q_to_chem(lnum)) = .FALSE.
!!            endif
!!            do l = 1, nspec_amode(m)
!!               lmass = lmassptr_amode(l,m)
!!               vertMixAero(lptr_q_to_chem(lmass)) = .FALSE.
!!            enddo
!!         enddo
!!      endif

!
! compute dry deposition velocities = ddvel
!
! 28-jun-2005 rce - initialize ddvel=0; call aerosol drydep routine
!           only when gas_drydep_opt == WESELY
!       the wesely_driver routine computes aer_res, and currently
!	    you cannot compute aerosol drydep without it !!
! 08-jul-2005 rce - pass idrydep_onoff to mixactivate routines
!
   ddvel(:,:,:) = 0.0
   idrydep_onoff = 0

   drydep_select: SELECT CASE(config_flags%gas_drydep_opt)

     CASE ( WESELY )
!
! gas_drydep_opt == WESELY means 
!     wesely for gases 
!     other (appropriate) routine for aerosols
!
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES WITH WESELY METHOD')

       IF( config_flags%chem_opt /= CHEM_TRACER    .and.                  &
           config_flags%chem_opt /= CHEM_TRACE2    .and.                  &
           config_flags%chem_opt /= CO2_TRACER     .and.                  & 
           config_flags%chem_opt /= GHG_TRACER     .and.                  &
           config_flags%chem_opt /= CHEM_VASH      .and.                  &
           config_flags%chem_opt /= CHEM_VOLC_4BIN .and.                  &
           config_flags%chem_opt /= DUST           .and.                  &
           config_flags%chem_opt /= GOCART_SIMPLE  .and.                  &
#if ( WRFPLUS == 1 )
           config_flags%chem_opt /= GOCART_BCCO .and.                  &
#endif
           config_flags%chem_opt /= GOCARTRACM_KPP )THEN
          call wesely_driver(id,ktau,dtstep,                              &
               config_flags,current_month,                                &
               gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
               p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
               ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
               snowh,numgas,                                              &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!-----------------------------------------------------------------
! this does aerosol species (dust,seas, bc,oc) for gocart only
!-----------------------------------------------------------------
          IF ( config_flags%chem_opt == MOZCART_KPP ) then
            call gocart_drydep_driver( dtstep,                            &
                  config_flags, numgas,                                   &
                  t_phy, moist, p8w, t8w, rmol,aer_res_def,               &
                  p_phy, chem, rho_phy, dz8w, ddvel, xland, hfx,          &
                  ivgtyp, tsk, vegfra, pbl, ust, znt, xlat, xlong,        &
                  dustdrydep_1,dustdrydep_2,dustdrydep_3,                 &
                  dustdrydep_4,dustdrydep_5,                              &
                  depvelocity,                                            &
                  ids,ide, jds,jde, kds,kde,                              &
                  ims,ime, jms,jme, kms,kme,                              &
                  its,ite, jts,jte, kts,kte )
          ENDIF
          if( config_flags%diagnostic_chem == DEPVEL1 .and. &
              (config_flags%chem_opt == MOZCART_KPP .or. &
              config_flags%chem_opt == MOZART_KPP   .or. &
              config_flags%chem_opt == MOZART_MOSAIC_4BIN_KPP .or. &
              config_flags%chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP ) ) then
               do j = jts,jte
                  dvel(its:ite,1,j,p_dvel_o3) = m2cm*ddvel(its:ite,j,p_o3)
                  dvel(its:ite,1,j,p_dvel_no) = m2cm*ddvel(its:ite,j,p_no)
                  dvel(its:ite,1,j,p_dvel_no2) = m2cm*ddvel(its:ite,j,p_no2)
                  dvel(its:ite,1,j,p_dvel_nh3) = m2cm*ddvel(its:ite,j,p_nh3)
                  dvel(its:ite,1,j,p_dvel_hno3) = m2cm*ddvel(its:ite,j,p_hno3)
                  dvel(its:ite,1,j,p_dvel_hno4) = m2cm*ddvel(its:ite,j,p_hno4)
                  dvel(its:ite,1,j,p_dvel_h2o2) = m2cm*ddvel(its:ite,j,p_h2o2)
                  dvel(its:ite,1,j,p_dvel_co) = m2cm*ddvel(its:ite,j,p_co)
                  dvel(its:ite,1,j,p_dvel_ch3ooh) = m2cm*ddvel(its:ite,j,p_ch3ooh)
                  dvel(its:ite,1,j,p_dvel_hcho) = m2cm*ddvel(its:ite,j,p_hcho)
                  dvel(its:ite,1,j,p_dvel_ch3oh) = m2cm*ddvel(its:ite,j,p_ch3oh)
                  dvel(its:ite,1,j,p_dvel_eo2) = m2cm*ddvel(its:ite,j,p_eo2)
                  dvel(its:ite,1,j,p_dvel_ald) = m2cm*ddvel(its:ite,j,p_ald)
                  dvel(its:ite,1,j,p_dvel_ch3cooh) = m2cm*ddvel(its:ite,j,p_ch3cooh)
                  dvel(its:ite,1,j,p_dvel_acet) = m2cm*ddvel(its:ite,j,p_acet)
                  dvel(its:ite,1,j,p_dvel_mgly) = m2cm*ddvel(its:ite,j,p_mgly)
                  dvel(its:ite,1,j,p_dvel_paa) = m2cm*ddvel(its:ite,j,p_paa)
                  dvel(its:ite,1,j,p_dvel_pooh) = m2cm*ddvel(its:ite,j,p_c3h6ooh)
                  dvel(its:ite,1,j,p_dvel_mpan) = m2cm*ddvel(its:ite,j,p_mpan)
                  dvel(its:ite,1,j,p_dvel_mco3) = m2cm*ddvel(its:ite,j,p_mco3)
                  dvel(its:ite,1,j,p_dvel_mvkooh) = m2cm*ddvel(its:ite,j,p_mvkooh)
                  dvel(its:ite,1,j,p_dvel_c2h5oh) = m2cm*ddvel(its:ite,j,p_c2h5oh)
                  dvel(its:ite,1,j,p_dvel_etooh) = m2cm*ddvel(its:ite,j,p_etooh)
                  dvel(its:ite,1,j,p_dvel_prooh) = m2cm*ddvel(its:ite,j,p_prooh)
                  dvel(its:ite,1,j,p_dvel_acetp) = m2cm*ddvel(its:ite,j,p_acetp)
                  dvel(its:ite,1,j,p_dvel_onit) = m2cm*ddvel(its:ite,j,p_onit)
                  dvel(its:ite,1,j,p_dvel_onitr) = m2cm*ddvel(its:ite,j,p_onitr)
                  dvel(its:ite,1,j,p_dvel_isooh) = m2cm*ddvel(its:ite,j,p_isooh)
                  dvel(its:ite,1,j,p_dvel_acetol) = m2cm*ddvel(its:ite,j,p_acetol)
                  dvel(its:ite,1,j,p_dvel_glyald) = m2cm*ddvel(its:ite,j,p_glyald)
                  dvel(its:ite,1,j,p_dvel_hydrald) = m2cm*ddvel(its:ite,j,p_hydrald)
                  dvel(its:ite,1,j,p_dvel_alkooh) = m2cm*ddvel(its:ite,j,p_alkooh)
                  dvel(its:ite,1,j,p_dvel_mekooh) = m2cm*ddvel(its:ite,j,p_mekooh)
                  dvel(its:ite,1,j,p_dvel_tolooh) = m2cm*ddvel(its:ite,j,p_tolooh)
                  dvel(its:ite,1,j,p_dvel_xooh) = m2cm*ddvel(its:ite,j,p_xooh)
                  dvel(its:ite,1,j,p_dvel_so2) = m2cm*ddvel(its:ite,j,p_so2)
                  dvel(its:ite,1,j,p_dvel_so4) = m2cm*ddvel(its:ite,j,p_sulf)
                  dvel(its:ite,1,j,p_dvel_pan) = m2cm*ddvel(its:ite,j,p_pan)
                  dvel(its:ite,1,j,p_dvel_terpooh) = m2cm*ddvel(its:ite,j,p_terpooh)
               enddo
          endif

          if ( config_flags%chem_opt == MOZART_MOSAIC_4BIN_AQ_KPP ) then
               do j = jts,jte
                   dvel(its:ite,1,j,p_dvel_cvasoaX) = 0.0
                   dvel(its:ite,1,j,p_dvel_cvasoa1) = m2cm*ddvel(its:ite,j,p_cvasoa1)
                   dvel(its:ite,1,j,p_dvel_cvasoa2) = m2cm*ddvel(its:ite,j,p_cvasoa2)
                   dvel(its:ite,1,j,p_dvel_cvasoa3) = m2cm*ddvel(its:ite,j,p_cvasoa3)
                   dvel(its:ite,1,j,p_dvel_cvasoa4) = m2cm*ddvel(its:ite,j,p_cvasoa4)
                   dvel(its:ite,1,j,p_dvel_cvbsoaX) = 0.0
                   dvel(its:ite,1,j,p_dvel_cvbsoa1) = m2cm*ddvel(its:ite,j,p_cvbsoa1)
                   dvel(its:ite,1,j,p_dvel_cvbsoa2) = m2cm*ddvel(its:ite,j,p_cvbsoa2)
                   dvel(its:ite,1,j,p_dvel_cvbsoa3) = m2cm*ddvel(its:ite,j,p_cvbsoa3)
                   dvel(its:ite,1,j,p_dvel_cvbsoa4) = m2cm*ddvel(its:ite,j,p_cvbsoa4)
                enddo
           endif
       ELSEIF ( config_flags%chem_opt == GOCART_SIMPLE &
#if ( WRFPLUS == 1 )
                .OR. config_flags%chem_opt == GOCART_BCCO &
#endif
              ) then
          call wesely_driver(id,ktau,dtstep,                              &
               config_flags,current_month,                                &
               gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
               p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    & 
               ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
               snowh, numgas,                                             &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
         call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                       &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
               dustdrydep_1,dustdrydep_2,dustdrydep_3,                 &
               dustdrydep_4,dustdrydep_5,                              &
               depvelocity,                                            &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
       ELSEIF ( config_flags%chem_opt == DUST) then
!
! this does aerosol species (dust) for gocart only
!
         call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                       &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
               dustdrydep_1,dustdrydep_2,dustdrydep_3,                 &
               dustdrydep_4,dustdrydep_5,                              &
               depvelocity,                                            &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )


       ELSEIF ( config_flags%chem_opt == GOCARTRACM_KPP) then
!
! this does aerosol species (dust,seas, bc,oc) for gocart only
!
          call wesely_driver(id,ktau,dtstep,                              &
               config_flags,current_month,                                &
               gmt,julday,t_phy,moist,p8w,t8w,raincv,                     &
               p_phy,chem,rho_phy,dz8w,ddvel,aer_res_def,aer_res_zcen,    &
               ivgtyp,tsk,gsw,vegfra,pbl,rmol,ust,znt,xlat,xlong,z,z_at_w,&
               snowh, numgas,                                             &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )
         call gocart_drydep_driver(dtstep,                                &
               config_flags,numgas,                                       &
               t_phy,moist,p8w,t8w,rmol,aer_res_def,                      &
               p_phy,chem,rho_phy,dz8w,ddvel,xland,hfx,                   &
               ivgtyp,tsk,vegfra,pbl,ust,znt,xlat,xlong,                  &
               dustdrydep_1,dustdrydep_2,dustdrydep_3,                 &
               dustdrydep_4,dustdrydep_5,                              &
               depvelocity,                                            &
               ids,ide, jds,jde, kds,kde,                                 &
               ims,ime, jms,jme, kms,kme,                                 &
               its,ite, jts,jte, kts,kte                                  )

       ELSE
          !Set dry deposition velocity to zero when using the
          !chemistry tracer mode.
          ddvel(:,:,:) = 0.
       END IF

! diagnose output dry deposition flux
       if( config_flags%diagnostic_dep == 1) then
       ddflx(:,:,:) = 0.D0
       do i = its, ite
       do j = jts, jte
       ddflx(i, j,1:numgas)=ddflx(i,j,1:numgas)+chem(i,kts,j,1:numgas)*p_phy(i,kts,j)/(RSI*t_phy(i,kts,j))*ddvel(i,j,1:numgas)*dtstep*1.E-6
       enddo
       enddo
       end if


       if (config_flags%aer_aerodynres_opt == 2) then
          ! use aerodynamic resistance from center of layer kts to surface
          aer_res(:,:) = aer_res_zcen(:,:)
       else
          ! this is the default -- use aerodynamic resistance from 
          ! "default reference height" (currently 2 m) to surface
          aer_res(:,:) = aer_res_def(:,:)
       end if

       idrydep_onoff = 1
       aer_mech_id_select: SELECT CASE(config_flags%chem_opt)
!!! TUCCELLA
          CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_AQCHEM,RADM2SORG_KPP, &
                RACM_ESRLSORG_KPP,RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP, &
                CBMZSORG,CBMZSORG_AQ, &
                CB05_SORG_AQ_KPP,CB05_SORG_VBS_AQ_KPP)
             aer_mech_id = 1
          CASE (RACMSORG_AQ,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP,RACMSORG_KPP)
             aer_mech_id = 2
          CASE ( CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_KPP, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, &
                 CBMZ_MOSAIC_8BIN_AQ,SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
                 CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
                 MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
                 CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP, &!BSINGH(12/03/2013): Added SAPRC 8 bin
                 SAPRC99_MOSAIC_8BIN_VBS2_KPP)!BSINGH(04/04/2014): Added SAPRC 8 bin non-aq
             aer_mech_id = 3
          CASE ( CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ )
             aer_mech_id = 4
          CASE DEFAULT
             aer_mech_id = 0
       END SELECT aer_mech_id_select


!  config_flags%aer_drydep_opt <= 0   -- aerosol depositon velocities are set to zero
!  config_flags%aer_drydep_opt >= 100 -- aerosol depvels are calculated in subr aer_depvel_driver
!                                        (see module_aer_drydep.F for details)
!  config_flags%aer_drydep_opt == 1   -- SORGAM chem packages use subr vdvg_2 of module_aerosols_sorgam.F
!                                        MOSAIC chem packages subr aerosol_depvel_2 of module_mosaic_drydep.F
!  config_flags%aer_drydep_opt == 11  -- SORGAM chem packages use subr vdvg of module_aerosols_sorgam.F

       if ((config_flags%aer_drydep_opt <= 0) .or. (aer_mech_id <= 0)) then
          CALL wrf_debug(15,'AEROSOL DRY DEP VELOCITIES  = 0.0')

       else if (config_flags%aer_drydep_opt <= 99) then

   adrydep_select: SELECT CASE(config_flags%chem_opt)
     CASE (RADM2SORG,RADM2SORG_AQ,RADM2SORG_AQCHEM,RADM2SORG_KPP,RACM_ESRLSORG_KPP,CBMZSORG,CBMZSORG_AQ)
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR AEROSOLS/RADM')
       call sorgam_depdriver (id,config_flags,ktau,dtstep,              &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,               &
               aer_res,ddvel(:,:,numgas+1:num_chem),                    &
               numgas,ddflx(:,:,numgas+1:num_chem), &
               num_chem-numgas,                                         &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
     CASE (RACMSORG_AQ,RACMSORG_AQCHEM_KPP,RACM_ESRLSORG_AQCHEM_KPP,RACMSORG_KPP)
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR AEROSOLS/RACM')
       call sorgam_depdriver (id,config_flags,ktau,dtstep,              &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,               &
               aer_res,ddvel(:,:,numgas+1:num_chem),                    &
               numgas,ddflx(:,:,numgas+1:num_chem), &
               num_chem-numgas,                                         &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
     CASE (CB05_SORG_AQ_KPP)
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR AEROSOLS/CB05')
       call sorgam_depdriver (id,config_flags,ktau,dtstep,              &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,z,z_at_w,                    &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,cvaro1,cvaro2, &
               cvalk1,cvole1,cvapi1,cvapi2,cvlim1,cvlim2,               &
               aer_res,ddvel(:,:,numgas+1:num_chem),                    &
               numgas,ddflx(:,:,numgas+1:num_chem), &
               num_chem-numgas,                                         &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
     CASE (CB05_SORG_VBS_AQ_KPP)
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR AEROSOLS/CB05')
       call sorgam_vbs_depdriver (id,config_flags,ktau,dtstep,          &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,rh,z,z_at_w,                 &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,               &
               aer_res,ddvel(:,:,numgas+1:num_chem),                    &
               numgas,ddflx(:,:,numgas+1:num_chem), &
               num_chem-numgas,                                         &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )

     CASE ( CBMZ_MOSAIC_4BIN, CBMZ_MOSAIC_KPP, CBMZ_MOSAIC_8BIN, CBMZ_MOSAIC_4BIN_AQ, CBMZ_MOSAIC_8BIN_AQ, &
          CBMZ_MOSAIC_DMS_4BIN, CBMZ_MOSAIC_DMS_8BIN, CBMZ_MOSAIC_DMS_4BIN_AQ, CBMZ_MOSAIC_DMS_8BIN_AQ, &
          SAPRC99_MOSAIC_4BIN_VBS2_KPP, &
          MOZART_MOSAIC_4BIN_KPP, MOZART_MOSAIC_4BIN_AQ_KPP, &
          CRI_MOSAIC_8BIN_AQ_KPP, CRI_MOSAIC_4BIN_AQ_KPP, &
          SAPRC99_MOSAIC_8BIN_VBS2_AQ_KPP,SAPRC99_MOSAIC_8BIN_VBS2_KPP)!BSINGH(12/03/2013): Added SAPRC 8 bin and non-aq on 04/04/2014
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR MOSAIC AEROSOLS')
       call mosaic_drydep_driver(                                       &
               id, curr_secs, ktau, dtstep, config_flags,               &
               gmt, julday,                                             &
               t_phy, rho_phy, p_phy,                                   &
               ust, aer_res,                                            &
               moist, chem, ddvel,                                      &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
!!! TUCCELLA
     CASE ( RACM_SOA_VBS_KPP,RACM_SOA_VBS_AQCHEM_KPP )
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR SOA_VBS AEROSOLS')
       call soa_vbs_depdriver (id,config_flags,ktau,dtstep,             &
               ust,t_phy,moist,p8w,t8w,rmol,znt,pbl,                    &
               alt,p_phy,chem,rho_phy,dz8w,rh,z,z_at_w,                 &
               h2oaj,h2oai,nu3,ac3,cor3,asulf,ahno3,anh3,               &
               aer_res,ddvel(:,:,numgas+1:num_chem),                    &
               num_chem-numgas,                                         &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
    CASE (CBMZ_CAM_MAM3_NOAQ, CBMZ_CAM_MAM3_AQ, CBMZ_CAM_MAM7_NOAQ, CBMZ_CAM_MAM7_AQ)
       CALL wrf_debug(15,'DOING DRY DEP VELOCITIES FOR CAM_MAM AEROSOLS')
       call cam_mam_drydep_driver(                                      &
               id, curr_secs, ktau, dtstep, config_flags,               &
               gmt, julday,                                             &
               t_phy, rho_phy, p_phy,                                   &
               ust, aer_res,                                            &
               moist, chem, ddvel,                                      &
               ids,ide, jds,jde, kds,kde,                               &
               ims,ime, jms,jme, kms,kme,                               &
               its,ite, jts,jte, kts,kte                                )
!
     CASE DEFAULT

     END SELECT adrydep_select
   else ! (config_flags%aer_drydep_opt > 99)
              CALL wrf_debug(15,'DOING DRY DEP VELOCITIES THRU AER_DRYDEP_DRIVER')
              call aer_drydep_driver(                                          &
                      id, ktau, dtstep, config_flags, aer_mech_id,           &
                      gmt, julday,                                             &
                      t_phy, rho_phy, p_phy,                                   &
                      alt, p8w, t8w, dz8w, z, z_at_w,                          &
                      ust, aer_res, ivgtyp, vegfra, pbl, rmol, znt,            &
                      moist, chem, ddvel,                                      &
                      h2oai, h2oaj, numgas,                                    &
                      ids,ide, jds,jde, kds,kde,                               &
                      ims,ime, jms,jme, kms,kme,                               &
                      its,ite, jts,jte, kts,kte                                )
   end if
       if (config_flags%aer_drydep_opt > 0) then
          if ((aer_mech_id > 0) .and. (aer_mech_id <= 4)) then
             ! limit aerosol ddvels to <= 0.5 m/s
             ! drydep routines occasionally produce unrealistically-large particle
             !    diameter leading to unrealistically-large sedimentation velocity
             ddvel(:,:,numgas+1:num_chem) = min( 0.50, ddvel(:,:,numgas+1:num_chem) )
          end if
       end if

       !For MAM aersols dry deposition
       if ((aer_mech_id == 4)) then
          do m_mam = 1, num_chem
             !Look for species with _a or _c and adjust these species for dry deposition
             k_a   = index(trim(adjustl(chem_dname_table(1,m_mam))),'_a') !To include *_a* species
             k_c   = index(trim(adjustl(chem_dname_table(1,m_mam))),'_c') !To include *_c* species
             kmax = max(k_a, k_c)
             if(kmax > 0 ) then !kmax>0 means a species has _a or _c
                frac_removed(its:ite,jts:jte) = max( 0.0, min( 1.0, ddvel(its:ite,jts:jte,m_mam)*dtstep/dz8w(its:ite,1,jts:jte) ) )
                chem(its:ite,1,jts:jte,m_mam) = chem(its:ite,1,jts:jte,m_mam)*(1.0 - frac_removed(its:ite,jts:jte)) !update lowest layer
             endif
          enddo
       endif
       
    CASE DEFAULT 
                                                     
   END SELECT drydep_select                              

! Add dep_vel here for all aerosol and gas species in the dvel array 
      ll = max( 1, min( config_flags%ndepvel, num_vert_mix ) )
      dep_vel(:,:,:,:) = 0.
      do l=1,ll
      do j=jts,jte
      do k=1,config_flags%kdepvel
      do i=its,ite
        dep_vel(i,k,j,l) = ddvel(i,j,l)
      enddo
      enddo
      enddo
      enddo

!      if( config_flags%diagnostic_dep == 1) then
!       ddlen(its:ite,jts:jte,:)=ddlen(its:ite,jts:jte,:)+ddvel(its:ite,jts:jte,:)*m2cm
!      end if

       CALL wrf_debug(15,'end of dry_dep_vel')

END SUBROUTINE dry_dep_vel

END MODULE module_dry_dep_vel
