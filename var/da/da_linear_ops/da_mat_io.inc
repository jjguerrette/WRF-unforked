subroutine da_mat_io (io_config, filename, n, m, MAT, ascii_bool)

   !-------------------------------------------------------------------------------
   ! Purpose:        Read / Write arbitrary matrix, stored identically on all cores
   !
   ! Called from da_block_lanczos
   !
   ! History: 12/10/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------------

   implicit none

   character, intent(in)         :: io_config    ! 'r' = Read; 'w' = Write
   character(*), intent(in)      :: filename     ! I/O filename
   integer, intent(in)           :: n,m
   real*8, intent(inout)         :: MAT(n,m)
   logical, optional, intent(in) :: ascii_bool

   integer   :: i
   integer   :: mat_unit, dummy
   real*8    :: temp(m)
   character(5) :: cm

   if (trace_use) call da_trace_entry("da_mat_io")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

   call da_get_unit (mat_unit)

   write(cm,fmt='(I0)') m

   if (io_config == 'r') then
      write(*,*) 'Reading matrix from '//trim(filename)
      if (rootproc) then
         if (present(ascii_bool)) then
            open (unit=mat_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
         else
            open (unit=mat_unit, file = trim(filename), form = 'unformatted', status = 'old')
         end if
      end if
      do i = 1, n
         if (rootproc) then
            if (present(ascii_bool)) then
               read(unit=mat_unit,fmt='('//trim(cm)//'E29.20E3)') temp
            else
               read(unit=mat_unit) temp
            end if
         end if
         call mpi_barrier(comm, ierr)
#ifdef DM_PARALLEL
         call wrf_dm_bcast_real(temp,m)
#endif
         MAT(i,:) = temp
      end do
   else if (io_config == 'w') then
      write(*,*) 'Writing matrix to '//trim(filename)
      if (rootproc) then
         if (present(ascii_bool)) then
            open (unit=mat_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
            do i = 1, n
               write(unit=mat_unit,fmt='('//trim(cm)//'E29.20E3)') MAT(i,:)
            end do
         else
            open (unit=mat_unit, file = trim(filename), form = 'unformatted', status = 'replace')
            do i = 1, n
               write(unit=mat_unit) MAT(i,:)
            end do
         end if
      end if
      call mpi_barrier(comm, ierr)
   else
      write(*,*) 'Unknow configuration for Matrix I/O routine'
   end if

   call da_free_unit (mat_unit)

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: RIOT_IO_mat')
#endif

   if (trace_use) call da_trace_exit ("da_mat_io")

end subroutine da_mat_io
