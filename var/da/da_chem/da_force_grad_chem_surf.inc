subroutine da_force_grad_chem_surf(iv, y, jo_grad_y, ensmember, obscount)

   !----------------------------------------------------------------------
   ! Purpose: Applies obs inverse on y-vector for Hessian calculation
   !----------------------------------------------------------------------

   implicit none

   type (iv_type), intent(in)    :: iv          ! Innovation vector.
   type (y_type),  intent(inout) :: y          ! Residual vector.
   type (y_type),  intent(inout) :: jo_grad_y   ! Grad_y(Jo)
   integer,        intent(in)    :: ensmember
   integer,        intent(out)   :: obscount

   integer :: n, ichem
   
   if (trace_use) call da_trace_entry("da_force_grad_chem_surf")       

   obscount = 0

   do n=1, iv%info(chem_surf)%nlocal
      do ichem = 1, num_chem_surf-1
          if (iv%chem_surf(n)%chem(ichem)%qc < obs_qc_pointer) y%chem_surf(n)%chem(ichem) = 0.0
          if (abs(iv%chem_surf(n)%chem(ichem)%inv) .gt. 1.D-12) obscount = obscount + 1
          if (obscount .eq. ensmember) jo_grad_y%chem_surf(n)%chem(ichem) = &
             1.D0 / (iv%chem_surf(n)%chem(ichem)%error)
      end do
   end do
   
   if (trace_use) call da_trace_exit("da_force_grad_chem_surf")  
     
end subroutine da_force_grad_chem_surf


