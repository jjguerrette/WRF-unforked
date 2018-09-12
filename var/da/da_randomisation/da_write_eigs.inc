#if defined(LAPACK)
subroutine da_write_eigs(    eignval, neign, nsamp, suffix)

   !-------------------------------------------------------------------------
   ! Purpose:   Write eigs for Pa and increment approximations
   !
   ! Called from da_compare_decomp_methods
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   real*8, intent(in)        :: eignval(neign)
   integer, intent(in)       :: neign, nsamp
   character(*), intent(in)  :: suffix

   real*8                       :: vals(nsamp)
   character(len=filename_len)  :: filename
   character(len=5)             :: cnsamp
   integer                      :: ounit, iostat

   if (trace_use) call da_trace_entry("da_write_eigs")
   vals = -9999.9
   vals(1:neign) = eignval
   filename = "EIGNVALUES_"//trim(suffix)

   write(cnsamp,fmt='(I0)') nsamp

   if(rootproc) then
      ! =========
      ! open file
      ! =========
      call da_get_unit(ounit)
      if (neign .eq. 1) then
         open(unit   = ounit,     &
            FILE   = trim(filename), &
            FORM   = 'FORMATTED',  &
            ACCESS = 'SEQUENTIAL', &
            iostat =  iostat,     &
            STATUS = 'REPLACE')
      else
         open(unit   = ounit,     &
            FILE   = trim(filename), &
            FORM   = 'FORMATTED',  &
            POSITION = 'APPEND', &
            iostat =  iostat,     &
            STATUS = 'UNKNOWN')
      end if

      if (iostat /= 0) then
         write(unit=message(1),fmt='(A,I5, I5, A)') "Error message: ",iostat, trim(filename)
         call da_warning(__FILE__,__LINE__,message(1:1))
         call da_wrf_error_fatal('Error opening NORM file')
      end if
      write(ounit,fmt='(I6,'//trim(adjustl(cnsamp))//'E27.16E3)') neign, vals
      close(ounit)
      call da_free_unit(ounit)
   end if
#ifdef DM_PARALLEL
   call mpi_barrier(comm, ierr)
#endif

   if (trace_use) call da_trace_exit("da_write_eigs")

end subroutine da_write_eigs
#endif
