subroutine da_ytoyhat_chem_acft_adj(iv, qhat, jo_grad_y)

   !----------------------------------------------------------------------
   ! Purpose: Applies (obs inverse)^1/2 on qhat
   !----------------------------------------------------------------------

   implicit none

   type (iv_type),    intent(in)    :: iv          ! Innovation vector.
   real, allocatable, intent(inout) :: qhat(:)
   type (y_type),     intent(inout) :: jo_grad_y   ! Grad_y(Jo)

   integer :: n, ichem, iplt, iobs
   
   if (trace_use) call da_trace_entry("da_ytoyhat_chem_acft_adj")       

   if (allocated(qhat)) then
      iobs = 1
      do n=1, iv%info(chem_acft)%nlocal
         do ichem = 1, num_chem_acft-1
            do iplt = 1, num_platform
               if (iv%chem_acft(n)%chem(ichem,iplt)%qc < obs_qc_pointer) &
                  qhat(iobs) = 0.0

               jo_grad_y%chem_acft(n)%chem(ichem,iplt) = qhat(iobs) / iv%chem_acft(n)%chem(ichem,iplt)%error

               iobs = iobs + 1
            end do
         end do
      end do
   end if
   
   if (trace_use) call da_trace_exit("da_ytoyhat_chem_acft_adj")  
     
end subroutine da_ytoyhat_chem_acft_adj


