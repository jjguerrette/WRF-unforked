subroutine da_hessian_io (io_config, cv_size, it, nens, neign, eignvec, eignval, qhat, ghat)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write Hessian eigenpairs
   !
   ! Called from da_solve
   !
   ! History: 04/19/2017  Modified from da_lanczos_io (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none
   
   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   integer,          intent(in)    :: cv_size
   integer,          intent(in)    :: it
   integer,          intent(in)    :: nens
   integer,          intent(inout) :: neign
   real*8, optional, intent(inout) :: eignvec(nens, nens)
   real*8, optional, intent(inout) :: eignval(nens)
   real*8, optional, intent(inout) :: qhat(1:cv_size, 0:nens)
   real*8, optional, intent(inout) :: ghat(1:cv_size)

   character(len=filename_len)     :: filename               ! I/O filename
   character*10                    :: cproc
   character*10                    :: cit
   integer                         :: ep_unit
   integer                         :: i, dummy1, iostat

   if (trace_use) call da_trace_entry("da_hessian_io")
   
   write(cproc,fmt='(i4.4)') myproc
   write(cit,fmt='(i2.2)') it
   filename = '../hessian_eigenpairs.it'//trim(adjustl(cit))//'.'//trim(adjustl(cproc))

   call da_get_unit (ep_unit)

   if (io_config == 'r') then
      write(*,*) 'Reading Hessian eigenpairs: ',trim(filename)
      open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'old')  
      if (iostat .ne. 0) then
         write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
         call da_warning(__FILE__,__LINE__,message(1:1))
         call da_wrf_error_fatal('Error opening hessian file')
      end if

      read(unit=ep_unit) neign, dummy1
      if (dummy1.ne.cv_size) then
         write(unit=message(1),fmt='(A,2I7)') "Error message: ",dummy1,cv_size
         call da_warning(__FILE__,__LINE__,message(1:1))
         call da_wrf_error_fatal('Incorrect cv_size in hessian_eigenpairs file')
      end if
      if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
         do i = 1, nens
            read(unit=ep_unit) eignval(i)
            read(unit=ep_unit) eignvec(1:nens,i)
            read(unit=ep_unit) qhat(1:cv_size,i)
         end do
         if (present(ghat)) read(unit=ep_unit) ghat(1:cv_size)
      end if
      close(unit=ep_unit)
   else if (io_config == 'w') then
      write(*,*) 'Writing Hessian eigenpairs: ',trim(filename)
      open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'replace')  
      if (iostat .ne. 0) then
         write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
         call da_warning(__FILE__,__LINE__,message(1:1))
         call da_wrf_error_fatal('Error opening hessian file')
      end if

      write(unit=ep_unit) neign, cv_size
      if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
         do i = 1, nens
            write(unit=ep_unit) eignval(i)
            write(unit=ep_unit) eignvec(1:nens,i)
            write(unit=ep_unit) qhat(1:cv_size,i)
         end do
         if (present(ghat)) write(unit=ep_unit) ghat(1:cv_size)
      end if
      close(unit=ep_unit)
   else
      write(*,*) 'Unknow configuration for Hessian I/O routine'
   end if

   call da_free_unit (ep_unit)

   if (trace_use) call da_trace_exit ("da_hessian_io")

end subroutine da_hessian_io
