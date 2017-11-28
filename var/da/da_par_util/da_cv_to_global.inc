subroutine da_cv_to_global(cv_size, cv_size_global, x, grid, mzs, xg )


   !-----------------------------------------------------------------------
   ! Purpose: Gathers local cv-array x into domain cv-array xg(:).  
   ! Global cv-array xg will only be valid on the "monitor" task.  
   !
   ! Must be called by all MPI tasks.  
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in)    :: cv_size        ! Size of local cv-array
   integer,          intent(in)    :: cv_size_global ! Size of domain cv-array
   real,             intent(in)    :: x(1:cv_size)   ! local cv-array
   type(domain),     intent(in)    :: grid
   integer,          intent(in)    :: mzs(:)  ! mz for each variable
                                         ! (to identify empty or 2D arrays)
!   real,             intent(inout) :: xg(1:cv_size_global) ! global cv-array
   real,             intent(inout) :: xg(:) ! global cv-array (needs to be dimension cv_size_global on rootproc, else dimension 1)

#ifdef DM_PARALLEL

   ! Local declarations
   type (vp_type) :: vv_x    ! Grdipt/EOF cv_array (local)
   type (vp_type) :: vv_xg   ! Grdipt/EOF cv_array (global)
   type (xpose_type) :: xpg  ! global dimensions
#if (WRF_CHEM == 1)
   type (xch_type) :: vvc_x  ! Grdipt/EOF cv_array (local)
   type (xch_type) :: vvc_xg ! Grdipt/EOF cv_array (global)
   integer :: ic, mz_ant, mz_bb
#endif
!   integer :: ids, ide, jds, jde, kds, &
!              ims, ime, jms, jme, kms, &
!              ips, ipe, jps, jpe, kps 

   integer :: n

   if (trace_use) call da_trace_entry("da_cv_to_global")      

   !
   ! Gather to mimic serial summation order.  
   !

   ! k?e varies with variable v1 - v5
!   ids = ids; ide = ide; jds = jds; jde = jde; kds = kds
!   ims = ims; ime = ime; jms = jms; jme = jme; kms = grid%xp%kms
   kms = grid%xp%kms
   ips = grid%xp%ips; ipe = grid%xp%ipe; jps = grid%xp%jps; jpe = grid%xp%jpe; kps = grid%xp%kps

   ! TOdo:  encapsulate this crap!  
   ! allocate vv_x
#if (WRF_CHEM == 1)
   if (cv_options .gt. 0) then
#endif
   allocate(vv_x%v1(ims:ime,jms:jme,mzs(1)))
   allocate(vv_x%v2(ims:ime,jms:jme,mzs(2)))
   allocate(vv_x%v3(ims:ime,jms:jme,mzs(3)))
   allocate(vv_x%v4(ims:ime,jms:jme,mzs(4)))
   allocate(vv_x%v5(ims:ime,jms:jme,mzs(5)))
#ifdef CLOUD_CV
   allocate(vv_x%v6(ims:ime,jms:jme,mzs(6)))
   allocate(vv_x%v7(ims:ime,jms:jme,mzs(7)))
   allocate(vv_x%v8(ims:ime,jms:jme,mzs(8)))
   allocate(vv_x%v9(ims:ime,jms:jme,mzs(9)))
   allocate(vv_x%v10(ims:ime,jms:jme,mzs(10)))
   allocate(vv_x%v11(ims:ime,jms:jme,mzs(11)))   
   allocate(vv_x%alpha(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int,mzs(13)))
#if (WRF_CHEM == 1)
   end if

   mz_ant=mzs(14)
   mz_bb =mzs(15)
#endif
#else
   allocate(vv_x%alpha(ims_int:ime_int,jms_int:jme_int,kms_int:kme_int,mzs(7)))
#if (WRF_CHEM == 1)
   end if

   mz_ant=mzs(8)
   mz_bb =mzs(9)
#endif
#endif

#if (WRF_CHEM == 1)
   if (cv_options_chem .gt. 0) then
      allocate(vvc_x%alpha_ant(ims:ime,jms:jme,mz_ant,num_scaleant))
      allocate(vvc_x%alpha_bb(ims:ime,jms:jme,mz_bb,num_scalebb))
   end if
#endif

   call da_cv_to_vv (cv_size, x, mzs, vv_x &
#if (WRF_CHEM == 1)
      , vvc_x &
#endif
      )

   if (rootproc) then
#if (WRF_CHEM == 1)
      if (cv_options .gt. 0) then
#endif
      ! allocate vv_xg
      allocate(vv_xg%v1(ids:ide,jds:jde,mzs(1)))
      allocate(vv_xg%v2(ids:ide,jds:jde,mzs(2)))
      allocate(vv_xg%v3(ids:ide,jds:jde,mzs(3)))
      allocate(vv_xg%v4(ids:ide,jds:jde,mzs(4)))
      allocate(vv_xg%v5(ids:ide,jds:jde,mzs(5)))
#ifdef CLOUD_CV
      allocate(vv_xg%v6(ids:ide,jds:jde,mzs(6)))
      allocate(vv_xg%v7(ids:ide,jds:jde,mzs(7)))
      allocate(vv_xg%v8(ids:ide,jds:jde,mzs(8)))
      allocate(vv_xg%v9(ids:ide,jds:jde,mzs(9)))
      allocate(vv_xg%v10(ids:ide,jds:jde,mzs(10)))
      allocate(vv_xg%v11(ids:ide,jds:jde,mzs(11)))
!     allocate(vv_xg%alpha(ids:ide,jds:jde,kds:kde,mzs(13)))
      allocate(vv_xg%alpha(ids_int:ide_int,jds_int:jde_int,kds_int:kde_int,mzs(13)))
#else
!     allocate(vv_xg%alpha(ids:ide,jds:jde,kds:kde,mzs(7)))
      allocate(vv_xg%alpha(ids_int:ide_int,jds_int:jde_int,kds_int:kde_int,mzs(7)))
#endif
#if (WRF_CHEM == 1)
      end if

      if (cv_options_chem .gt. 0) then
         allocate(vvc_xg%alpha_ant(ids:ide,jds:jde,mz_ant,num_scaleant))
         allocate(vvc_xg%alpha_bb(ids:ide,jds:jde,mz_bb,num_scalebb))
      end if
#endif
   else
      ! Allocate dummy array for non-monitor process to keep Fortran
      ! CICO happy...
#if (WRF_CHEM == 1)
      if (cv_options .gt. 0) then
#endif
      allocate(vv_xg%v1(1,1,1))
      allocate(vv_xg%v2(1,1,1))
      allocate(vv_xg%v3(1,1,1))
      allocate(vv_xg%v4(1,1,1))
      allocate(vv_xg%v5(1,1,1))
#ifdef CLOUD_CV
      allocate(vv_xg%v6(1,1,1))
      allocate(vv_xg%v7(1,1,1))
      allocate(vv_xg%v8(1,1,1))
      allocate(vv_xg%v9(1,1,1))
      allocate(vv_xg%v10(1,1,1))
      allocate(vv_xg%v11(1,1,1))
#endif
      allocate(vv_xg%alpha(1,1,1,1))
#if (WRF_CHEM == 1)
      end if

      if (cv_options_chem .gt. 0) then
         allocate(vvc_xg%alpha_ant(1,1,1,1))
         allocate(vvc_xg%alpha_bb(1,1,1,1))
      end if
#endif
   end if

   ! TOdo:  encapsulate this crap!  
   ! gather to global data structures
   ! it is possible to gather straight into a globally-sized cv-array
   ! TOdo:  add this optimization later
#if (WRF_CHEM == 1)
   if (cv_options .gt. 0) then
#endif
   call da_patch_to_global(grid, vv_x%v1,    vv_xg%v1,    mzs(1))
   call da_patch_to_global(grid, vv_x%v2,    vv_xg%v2,    mzs(2))
   call da_patch_to_global(grid, vv_x%v3,    vv_xg%v3,    mzs(3))
   call da_patch_to_global(grid, vv_x%v4,    vv_xg%v4,    mzs(4))
   call da_patch_to_global(grid, vv_x%v5,    vv_xg%v5,    mzs(5))
#ifdef CLOUD_CV
   call da_patch_to_global(grid, vv_x%v6,    vv_xg%v6,    mzs(6))
   call da_patch_to_global(grid, vv_x%v7,    vv_xg%v7,    mzs(7))
   call da_patch_to_global(grid, vv_x%v8,    vv_xg%v8,    mzs(8))
   call da_patch_to_global(grid, vv_x%v9,    vv_xg%v9,    mzs(9))
   call da_patch_to_global(grid, vv_x%v10,    vv_xg%v10,    mzs(10))
   call da_patch_to_global(grid, vv_x%v11,    vv_xg%v11,    mzs(11))
#endif
#ifdef CLOUD_CV
   if ( mzs(13) > 0 ) then
      do n = 1, mzs(13) ! Ensemble size
         if ( anal_type_hybrid_dual_res ) then
            call da_patch_to_global_dual_res(grid%intermediate_grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(12))
         else
            call da_patch_to_global(grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(12))
         endif
      end do
   end if
#else
   if ( mzs(7) > 0 ) then
      do n = 1, mzs(7) ! Ensemble size
         if ( anal_type_hybrid_dual_res ) then
            call da_patch_to_global_dual_res(grid%intermediate_grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(10))
         else
            call da_patch_to_global(grid, vv_x%alpha(:,:,:,n), vv_xg%alpha(:,:,:,n), mzs(6))
         endif
      end do
   end if
#endif
#ifdef CLOUD_CV
   ! deallocate vv_x
   deallocate (vv_x%v1, vv_x%v2, vv_x%v3, vv_x%v4, vv_x%v5, vv_x%v6, vv_x%v7, vv_x%v8, vv_x%v9, vv_x%v10, vv_x%v11, vv_x%alpha)
#else
   ! deallocate vv_x
   deallocate (vv_x%v1, vv_x%v2, vv_x%v3, vv_x%v4, vv_x%v5, vv_x%alpha)
#endif
#if (WRF_CHEM == 1)
   end if

   if (cv_options_chem .gt. 0) then
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         call da_patch_to_global(grid, vvc_x%alpha_ant(:,:,:,ic),    vvc_xg%alpha_ant(:,:,:,ic),    mz_ant)
      end do
      do ic = PARAM_FIRST_SCALAR, num_scalebb
         call da_patch_to_global(grid, vvc_x%alpha_bb(:,:,:,ic),    vvc_xg%alpha_bb(:,:,:,ic),    mz_bb)
      end do

      ! deallocate vvc_x
      deallocate (vvc_x%alpha_ant, vvc_x%alpha_bb)
   end if
#endif
   if (rootproc) then
      ! finally, collapse data back into a globally-sized cv-array
      xpg%ids = ids; xpg%ide = ide
      xpg%ims = ids; xpg%ime = ide
      xpg%its = ids; xpg%ite = ide
      xpg%jds = jds; xpg%jde = jde
      xpg%jms = jds; xpg%jme = jde
      xpg%jts = jds; xpg%jte = jde
      xpg%kds = kds; xpg%kde = kde
      xpg%kms = kds; xpg%kme = kde
      xpg%kts = kds; xpg%kte = kde
      call da_vv_to_cv(vv_xg, xpg, mzs, cv_size_global, xg &
#if (WRF_CHEM == 1)
         , vvc_xg &
#endif
      )
   end if
#if (WRF_CHEM == 1)
   if (cv_options .gt. 0) then
#endif
#ifdef CLOUD_CV
   ! deallocate vv_xg
   deallocate(vv_xg%v1, vv_xg%v2, vv_xg%v3, vv_xg%v4, vv_xg%v5, vv_xg%v6, vv_xg%v7,vv_xg%v8,vv_xg%v9, vv_xg%v10, vv_xg%v11, vv_xg%alpha)
#else
   ! deallocate vv_xg
   deallocate(vv_xg%v1, vv_xg%v2, vv_xg%v3, vv_xg%v4, vv_xg%v5, vv_xg%alpha)
#endif
#if (WRF_CHEM == 1)
   end if

   ! deallocate vvc_xg
   if (cv_options_chem .gt. 0) then
      deallocate (vvc_xg%alpha_ant, vvc_xg%alpha_bb)
   end if
#endif

   if (trace_use) call da_trace_exit("da_cv_to_global")    

#endif

end subroutine da_cv_to_global


