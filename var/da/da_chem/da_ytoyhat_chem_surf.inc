subroutine da_ytoyhat_chem_surf(iv, y, yhat)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on y-vector
   !----------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector.
   type (y_type),     intent(inout) :: y           ! Residual vector.
   real, allocatable, intent(inout) :: yhat(:)

   integer :: n, ichem, iobs
   
   if (trace_use) call da_trace_entry("da_ytoyhat_chem_surf")       

   if (allocated(yhat)) deallocate(yhat)

   if ((num_chem_surf-1) .gt. 0 ) then

      allocate(yhat( iv%info(chem_surf)%nlocal * (num_chem_surf-1) ))

      yhat(:) = 0.D0

      iobs = 1
      do n=1, iv%info(chem_surf)%nlocal
         do ichem = 1, num_chem_surf-1
            if (iv%chem_surf(n)%chem(ichem)%qc < obs_qc_pointer) &
               y%chem_surf(n)%chem(ichem) = 0.0

            yhat(iobs) = -y%chem_surf(n)%chem(ichem) / iv%chem_surf(n)%chem(ichem)%error

            iobs = iobs + 1
         end do
      end do
   end if

   if (trace_use) call da_trace_exit("da_ytoyhat_chem_surf")  
     
end subroutine da_ytoyhat_chem_surf


