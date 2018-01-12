subroutine da_cv_io_global (io_config, grid, filename, cv_size, rcv, mz, ascii_bool)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write yhat, omega, qhat, and cv vectors
   !
   ! Called from da_solve and da_randomise_svd
   !
   ! History: 09/01/2016  Creation (JJ Guerrette), modifed from da_lanczos_io.inc
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)     :: grid
   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   character(*),     intent(in)    :: filename               ! I/O filename
   integer,          intent(in)    :: cv_size
   real*8,           intent(inout) :: rcv(1:cv_size)
   integer,          intent(in)    :: mz(:)
   logical,      optional, intent(in)    :: ascii_bool

   integer                         :: cv_unit, dummy

   real*8, allocatable       :: cv_global(:)

   if (trace_use) call da_trace_entry("da_cv_io_global")
#ifdef DM_PARALLEL
   cv_size_domain = wrf_dm_sum_integer(cv_size)

   if (.not.rootproc) cv_size_domain = 1
#else
   cv_size_domain = cv_size
#endif
   if (rootproc) then
      allocate(cv_global(1:cv_size_domain))
   else
      allocate(cv_global(1))
   end if

   call da_get_unit (cv_unit)

   if (io_config == 'r') then
      if (rootproc) then
         write(*,*) 'Reading cv vector from '//trim(filename)
         if (present(ascii_bool)) then
            open (unit=cv_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
            read(unit=cv_unit,fmt='(I10)') dummy
         else
            open (unit=cv_unit, file = trim(filename), form = 'unformatted', status = 'old')
            read(unit=cv_unit) dummy
         end if
         if ( dummy .ne. cv_size_domain ) then
            write(unit=message(1),fmt='(A,2I7)') "Error message: ",dummy,cv_size_domain
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Incorrect cv_size_domain in global cv I/O')
         end if
         if (present(ascii_bool)) then
            read(unit=cv_unit,fmt='(E27.16E3)') cv_global(1:cv_size_domain)
         else
            read(unit=cv_unit) cv_global(1:cv_size_domain)
         end if
         close(unit=cv_unit)
      else
         cv_global = 0.D0
      end if

#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
      call da_global_to_cv(cv_size, cv_size_domain, rcv, grid, mz, cv_global)
#else
      rcv = cv_global
#endif

   else if (io_config == 'w') then
#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
      call da_cv_to_global(cv_size, cv_size_domain, rcv, grid, mz, cv_global)
#else
      cv_global = rcv
#endif

      if (rootproc) then
         write(*,*) 'Writing cv vector to '//trim(filename)
         if (present(ascii_bool)) then
            open (unit=cv_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
            write(unit=cv_unit,fmt='(I10)') cv_size_domain
            write(unit=cv_unit,fmt='(E27.16E3)') cv_global(1:cv_size_domain)
         else
            open (unit=cv_unit, file = trim(filename), form = 'unformatted', status = 'replace')
            write(unit=cv_unit) cv_size_domain
            write(unit=cv_unit) cv_global(1:cv_size_domain)
         end if
         close(unit=cv_unit)
      end if

#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
#endif

   else
      write(*,*) 'Unknown configuration for CV I/O routine'
   end if

   deallocate(cv_global)

   call da_free_unit (cv_unit)

   if (trace_use) call da_trace_exit ("da_cv_io_global")

end subroutine da_cv_io_global
