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

   if (trace_use) call da_trace_entry("da_mat_io")

   call da_get_unit (mat_unit)

   if (io_config == 'r') then
      write(*,*) 'Reading matrix from '//trim(filename)
      if (rootproc) then
         if (present(ascii_bool)) then
            open (unit=mat_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
         else
            open (unit=mat_unit, file = trim(filename), form = 'unformatted', status = 'old')
         end if
      end if
      do i = 1, m
         if (rootproc) then
            if (present(ascii_bool)) then
               read(unit=mat_unit,fmt='(E27.18E3)') temp_mat(:,i)
            else
               read(unit=mat_unit) temp_mat(:,i)
            end if
         end if
         call mpi_barrier(comm, ierr)
#ifdef DM_PARALLEL
         call wrf_dm_bcast_real(temp_mat(:,i),n)
#endif
         MAT(:,i) = temp_mat(:,i)
      end do
   else if (io_config == 'w') then
      write(*,*) 'Writing matrix to '//trim(filename)
      if (rootproc) then
         if (present(ascii_bool)) then
            open (unit=mat_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
            do i = 1, m
                  write(unit=mat_unit,fmt='(E27.18E3)') MAT(:,i)
            end do
         else
            open (unit=mat_unit, file = trim(filename), form = 'unformatted', status = 'replace')
            do i = 1, m
                  write(unit=mat_unit) MAT(:,i)
            end do
         end if
      end if
      call mpi_barrier(comm, ierr)
   else
      write(*,*) 'Unknow configuration for Matrix I/O routine'
   end if

   call da_free_unit (mat_unit)

   if (trace_use) call da_trace_exit ("da_mat_io")

end subroutine da_mat_io
