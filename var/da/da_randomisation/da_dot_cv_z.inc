complex function da_dot_cv_z(cv_size, x, y, grid, mzs, jp_start, jp_end)

   !-----------------------------------------------------------------------
   ! Purpose: Forms the dot product of two vectors that are organized in the 
   ! format of a "cv_type".  
   !
   ! Capable of producing bitwise-exact results for distributed-memory 
   ! parallel runs for testing.  This feature is very slow and consumes 
   ! lots of memory. 
   !-----------------------------------------------------------------------

   implicit none

   integer,          intent(in) :: cv_size           ! Size of array (tile).
   complex*16,             intent(in) :: x(cv_size)        ! 1st vector.
   complex*16,             intent(in) :: y(cv_size)        ! 1st vector.
   type(domain),     intent(in) :: grid              ! decomposed dimensions
   integer,          intent(in) :: mzs(:)            ! mz for each variable
                                                     ! (to identify 2D arrays)

   integer, optional,intent(in) :: jp_start, jp_end

   logical       :: lvarbc
   complex*16, pointer :: xx(:), yy(:)                     ! Temporary vectors.
   complex*16, pointer :: xg(:), yg(:)                     ! Serial data arrays.
   complex*16          :: dtemp1(1)                        ! Temporary.
   real                :: dtmpr, dtmpi                     ! Temporary.

   if (trace_use) call da_trace_entry("da_dot_cv_z")

   allocate(xx(1:cv_size))
   allocate(yy(1:cv_size))
   xx = x
   yy = y

 ! VarBC parameters are global (no summation across processors)
 !-------------------------------------------------------------
   lvarbc = present(jp_start) .and. present(jp_end)
   if (lvarbc) lvarbc = lvarbc .and. (jp_end >= jp_start)   
   if (lvarbc .and. .not. rootproc) then
      xx(jp_start:jp_end) = 0.
      yy(jp_start:jp_end) = 0.
   end if
      
      dtemp1(1) = da_dot_z(cv_size, xx, yy)
      
      !if (.not. global) then
         dtmpr = REAL( dtemp1(1) ,8)
         dtmpi = REAL(IMAG( dtemp1(1)) ,8)

         ! summation across processors:
         dtemp1(1) = dcmplx(wrf_dm_sum_real(dtmpr) , wrf_dm_sum_real(dtmpi))

      !end if
 
   deallocate(xx,yy)

   da_dot_cv_z = dtemp1(1)

   if (trace_use) call da_trace_exit("da_dot_cv_z")

end function da_dot_cv_z
