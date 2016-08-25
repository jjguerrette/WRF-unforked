real function da_dot_obs(nobs, x, y)

   !-----------------------------------------------------------------------
   ! Purpose: Forms the dot product of two vectors of length nobs
   !
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in) :: nobs           ! Size of array (tile).
   real,             intent(in) :: x(:)        ! 1st vector.
   real,             intent(in) :: y(:)        ! 1st vector.

   real, pointer :: xx(:), yy(:)                     ! Temporary vectors.
   real, pointer :: xg(:), yg(:)                     ! Serial data arrays.
   real          :: dtemp1(1), dtmp                  ! Temporary.

   if (trace_use) call da_trace_entry("da_dot_obs")


   if (nobs .gt. 0) then
      allocate(xx(1:nobs))
      allocate(yy(1:nobs))
      xx = x(1:nobs)
      yy = y(1:nobs)

      dtemp1(1) = da_dot(nobs, xx, yy)
   else
      dtemp1(1) = 0
   end if

   dtmp = dtemp1(1)
   ! summation across processors:
   dtemp1(1) = wrf_dm_sum_real(dtmp)
 
   if (allocated(xx)) deallocate(xx,yy)

   da_dot_obs = dtemp1(1)

   if (trace_use) call da_trace_exit("da_dot_obs")

end function da_dot_obs
