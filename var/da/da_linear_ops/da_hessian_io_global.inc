subroutine da_hessian_io_global (io_config, grid, cv_size, nens, neign, filename, eignvec, eignval, qhat, mz, ghat)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write Hessian eigenpairs
   !
   ! Called from da_hessian_io
   !
   ! History: 04/19/2017  Modified from da_lanczos_io (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none
   type(domain), intent(inout)     :: grid
   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   integer,          intent(in)    :: cv_size
   integer,          intent(in)    :: nens
   integer,          intent(inout) :: neign
   character(len=filename_len), intent(in)  :: filename               ! I/O filename

   real*8, optional, intent(inout) :: eignvec(nens, nens)
   real*8, optional, intent(inout) :: eignval(nens)
   real*8, optional, intent(inout) :: qhat(1:cv_size, 0:nens)
   integer, optional,intent(in)    :: mz(:)
   real*8, optional, intent(inout) :: ghat(1:cv_size)

   character*10                    :: cproc
   character*10                    :: cit
   integer                         :: ep_unit
   integer                         :: i, dummy, iostat

   real*8, allocatable             :: dumcv(:)
   real*8, allocatable             :: dumcv_g(:)
   real*8, dimension(1:nens)       :: temp
   real*8, dimension(1:nens,1:nens):: tempa


   if (trace_use) call da_trace_entry("da_hessian_io_global")
   
#ifdef DM_PARALLEL
   cv_size_domain = wrf_dm_sum_integer(cv_size)

   if (.not.rootproc) cv_size_domain = 1
#else
   cv_size_domain = cv_size
#endif
   if (rootproc) then
      allocate(dumcv_g(1:cv_size_domain))
   else
      allocate(dumcv_g(1))
   end if
   allocate(dumcv(cv_size))

   call da_get_unit (ep_unit)

   if (io_config == 'r') then

      neign = 0

      if (rootproc) then
         write(*,*) 'Reading Hessian eigenpairs: ',trim(filename)
         open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'old')  
         if (iostat .ne. 0) then
            write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening hessian file')
         end if

         read(unit=ep_unit) neign, dummy

         if (dummy.ne.cv_size_domain) then
            write(unit=message(1),fmt='(A,2I7)') "Error message: ",dummy,cv_size_domain
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Incorrect cv_size_domain in hessian_eigenpairs file')
         end if
      end if

      call mpi_barrier(comm, ierr)
      call wrf_dm_bcast_integer(neign,1)

      if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
         if (rootproc) read(unit=ep_unit) dummy

         do i = 1, nens
            if (rootproc) then
               read(unit=ep_unit) temp(i)
               read(unit=ep_unit) tempa(:,i)
               read(unit=ep_unit) dumcv_g
            else
               temp = 0.D0
               tempa = 0.D0
               dumcv_g = 0.D0
            end if
#ifdef DM_PARALLEL
            call mpi_barrier(comm, ierr)
            call da_global_to_cv(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)
            qhat(1:cv_size,i) = dumcv
#else
            qhat(1:cv_size,i) = dumcv_g
#endif
            eignvec(:,i) = tempa(:,i)
         end do
         eignval = temp

         if (present(ghat)) then
            if (rootproc) then
               read(unit=ep_unit) dumcv_g
            else
               dumcv_g = 0.D0
            end if
#ifdef DM_PARALLEL
            call mpi_barrier(comm, ierr)
            call da_global_to_cv(cv_size, cv_size_domain, ghat, grid, mz, dumcv_g)
#else
            ghat = dumcv_g
#endif
         end if
!      else
!         if (rootproc) read(unit=ep_unit) nens
!#ifdef DM_PARALLEL
!         call mpi_barrier(comm, ierr)
!         call wrf_dm_bcast_integer(nens,1)
!#endif
      end if
      if (rootproc) close(unit=ep_unit)
   else if (io_config == 'w') then
      if (rootproc) write(*,*) 'Writing Hessian eigenpairs: ',trim(filename)
      if (rootproc) then
         open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'replace')  
         if (iostat .ne. 0) then
            write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening hessian file')
         end if

         write(unit=ep_unit) neign, cv_size_domain
         write(unit=ep_unit) nens
      end if

      if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
         do i = 1, nens
#ifdef DM_PARALLEL
            dumcv = qhat(1:cv_size,i)
            call mpi_barrier(comm, ierr)
            call da_cv_to_global(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)
#else
            dumcv_g = qhat(1:cv_size,i)
#endif
            if (rootproc) then
               write(unit=ep_unit) eignval(i)
               write(unit=ep_unit) eignvec(1:nens,i)
               write(unit=ep_unit) dumcv_g
            end if
         end do
         if (present(ghat)) then
#ifdef DM_PARALLEL
            call mpi_barrier(comm, ierr)
            call da_cv_to_global(cv_size, cv_size_domain, ghat, grid, mz, dumcv_g)
#else
            dumcv_g = ghat
#endif
            if (rootproc) write(unit=ep_unit) dumcv_g
         end if
      end if
      if (rootproc) close(unit=ep_unit)
   else
      write(*,*) 'Unknow configuration for Hessian I/O routine'
   end if

   deallocate(dumcv_g)
   deallocate(dumcv)

   call da_free_unit (ep_unit)

   if (trace_use) call da_trace_exit ("da_hessian_io_global")

end subroutine da_hessian_io_global
