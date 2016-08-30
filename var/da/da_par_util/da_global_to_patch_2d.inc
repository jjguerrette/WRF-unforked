subroutine da_global_to_patch_2d (grid, vlocal, vglobal)

   !---------------------------------------------------------------------
   ! Purpose: Reverse of Gathers local 2D array vlocal into global array vglobal. 
   !
   ! Must be called by all MPI tasks.
   !---------------------------------------------------------------------  

   implicit none

   type(domain), intent(in)  :: grid
   real,         intent(out) :: vlocal(:,:)
   real,         intent(in)  :: vglobal(:,:)

   real, allocatable :: vlocal3d(:,:,:), vglobal3d(:,:,:)

   if (trace_use_frequent) call da_trace_entry("da_global_to_patch_2d")

   allocate(vlocal3d (ims:ime, jms:jme, 1:1))
   allocate(vglobal3d(ids:ide, jds:jde, 1:1))

   if (rootproc) then
       vglobal3d(:,:,1) = vglobal(:,:)
   end if

   call da_global_to_patch_3d(grid, vlocal3d, vglobal3d, 1)

   vlocal(:,:) = vlocal3d(:,:,1)

   deallocate(vlocal3d)
   deallocate(vglobal3d)

   if (trace_use_frequent) call da_trace_exit("da_global_to_patch_2d")

end subroutine da_global_to_patch_2d


