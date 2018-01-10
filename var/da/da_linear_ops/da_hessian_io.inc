subroutine da_hessian_io (io_config, grid, cv_size, it, nens, neign, eignvec, eignval, qhat, mz, ghat, hess_file)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write Hessian eigenpairs
   !
   ! Called from da_minimise_lz, da_randomise_svd, da_spectral_precon, and others
   !
   ! History: 04/19/2017  Modified from da_lanczos_io (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none
   type(domain), intent(inout)     :: grid
   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   integer,          intent(in)    :: cv_size
   integer,          intent(in)    :: it
   integer,          intent(inout) :: nens
   integer,          intent(inout) :: neign

   real*8, optional, intent(inout) :: eignvec(nens, nens)
   real*8, optional, intent(inout) :: eignval(nens)
   real*8, optional, intent(inout) :: qhat(1:cv_size, 0:nens)
   integer, optional,intent(in)    :: mz(:)
   real*8, optional, intent(inout) :: ghat(1:cv_size)
   character(len=filename_len), optional, intent(in)  :: hess_file  ! I/O filename

   character(len=filename_len)     :: filename               ! I/O filename
   character*10                    :: cproc
   character*10                    :: cit
   integer                         :: ep_unit
   integer                         :: i, dummy, iostat

   if (trace_use) call da_trace_entry("da_hessian_io")

   if ( present(hess_file) ) then
      filename = hess_file
   else
      write(cproc,fmt='(i4.4)') myproc
      write(cit,fmt='(i2.2)') it
      filename = '../hessian_eigenpairs.it'//trim(adjustl(cit))//'.'//trim(adjustl(cproc))
   end if

   if ( io_config == 'r' .and. .not.any((/present(eignvec),present(eignval),present(qhat)/)) ) then
      if (rootproc) then
         write(*,*) 'Reading nens and neign from Hessian eigenpair file: ',trim(filename)
         call da_get_unit (ep_unit)
         open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'old')
         if (iostat .ne. 0) then
            write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening hessian file')
         end if

         read(unit=ep_unit) neign, dummy
         read(unit=ep_unit) nens
      end if
#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
      call wrf_dm_bcast_integer(neign,1)
      call wrf_dm_bcast_integer(nens,1)
#endif
      if (rootproc) then
         close(unit=ep_unit)
         call da_free_unit (ep_unit)
      end if
   else
      if (use_global_cv_io) then
         if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
            if (present(ghat)) then
               call da_hessian_io_global(io_config, grid, cv_size, nens, neign, filename, &
                                         eignvec,eignval,qhat,mz,ghat)
            else
               call da_hessian_io_global(io_config, grid, cv_size, nens, neign, filename, &
                                         eignvec,eignval,qhat,mz)
            end if
!         else
!            call da_hessian_io_global(io_config, grid, cv_size, nens, neign, filename)
         end if
      else

         call da_get_unit (ep_unit)

         if (io_config == 'r') then
            write(*,*) 'Reading Hessian eigenpairs: ',trim(filename)
            open (unit=ep_unit, file = trim(filename), form = 'unformatted', iostat =  iostat, status = 'old')  
            if (iostat .ne. 0) then
               write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
               call da_warning(__FILE__,__LINE__,message(1:1))
               call da_wrf_error_fatal('Error opening hessian file')
            end if

            read(unit=ep_unit) neign, dummy
            if (dummy.ne.cv_size) then
               write(unit=message(1),fmt='(A,2I7)') "Error message: ",dummy,cv_size
               call da_warning(__FILE__,__LINE__,message(1:1))
               call da_wrf_error_fatal('Incorrect cv_size in hessian_eigenpairs file')
            end if
            if (present(eignvec) .and. present(eignval) .and. present(qhat)) then
               read(unit=ep_unit) dummy
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
            write(unit=ep_unit) nens

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

      end if
   end if
   if (trace_use) call da_trace_exit ("da_hessian_io")

end subroutine da_hessian_io
