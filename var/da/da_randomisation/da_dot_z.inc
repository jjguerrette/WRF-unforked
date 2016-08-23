complex function da_dot_z(n,x,y)

   !-----------------------------------------------------------------------
   ! Purpose: forms the dot product of two vectors.
   ! uses unrolled loops for increments equal to one.
   !-----------------------------------------------------------------------

   implicit none

   integer, intent(in) :: n
   complex*16,    intent(in) :: x(n)
   complex*16,    intent(in) :: y(n)

   complex*16    :: dtemp1
   integer :: i,m,mp1

   da_dot_z = (0.D0, 0.D0)
   if (n <= 0) return

   if (trace_use) call da_trace_entry("da_dot_z")    

   dtemp1 = 0.0

   ! code for both increments equal to 1

   if (n > 0) then
      m = mod(n,5)
      if (m /= 0) then
         do i = 1,m
            dtemp1 = dtemp1 + x(i)*y(i)
         end do
      end if
      if (n >= 5) then
         mp1 = m + 1
         do i = mp1,n,5
            dtemp1 = dtemp1 + x(i   )*y(i   ) + x(i + 1)*y(i + 1) + &
                              x(i + 2)*y(i + 2) + x(i + 3)*y(i + 3) + &
                              x(i + 4)*y(i + 4)
         end do
      end if
   end if

   da_dot_z = dtemp1

   if (trace_use) call da_trace_exit("da_dot_z")    


end function da_dot_z


