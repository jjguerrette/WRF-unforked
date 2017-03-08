subroutine da_ytoyhat_chem_acft(iv, y, yhat)

   !----------------------------------------------------------------------
   ! Purpose: Applies (obs inverse)^1/2 on yhat
   !----------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector.
   type (y_type),     intent(inout) :: y           ! Residual vector.
   real, allocatable, intent(inout) :: yhat(:)

   integer :: n, ichem, iplt, iobs
   
   if (trace_use) call da_trace_entry("da_ytoyhat_chem_acft")       

   if (allocated(yhat)) deallocate(yhat)

   if ((num_chem_acft-1) .gt. 0 .and. &
       num_platform .gt. 0) then

      allocate(yhat( iv%info(chem_acft)%nlocal * (num_chem_acft-1) * num_platform ))

      yhat(:) = 0.D0

      iobs = 1
      do iplt = 1, num_platform
         do n=1, iv%info(chem_acft)%nlocal
            do ichem = 1, num_chem_acft-1

!!!JJG - !!write(unit=stdout,fmt='("DEBUG_01:",3I5,2D25.14)') n,ichem,iplt,y%chem_acft(n)%chem(ichem,iplt),iv%chem_acft(n)%chem(ichem,iplt)%error

               if (iv%chem_acft(n)%chem(ichem,iplt)%qc < obs_qc_pointer) &
                  y%chem_acft(n)%chem(ichem,iplt) = 0.0

               yhat(iobs) = -y%chem_acft(n)%chem(ichem,iplt) / iv%chem_acft(n)%chem(ichem,iplt)%error

               iobs = iobs + 1
            end do
         end do
      end do
   end if
   
   if (trace_use) call da_trace_exit("da_ytoyhat_chem_acft")  
     
end subroutine da_ytoyhat_chem_acft


