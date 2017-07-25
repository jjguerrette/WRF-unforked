subroutine da_cv_io (io_config, filename, cv_size, rcv, grid, mz)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write yhat, omega, qhat, and cv vectors
   !
   ! Called from da_solve and da_randomise_svd
   !
   ! History: 09/01/2016  Creation (JJ Guerrette), modifed from da_lanczos_io.inc
   !
   !-------------------------------------------------------------------------

   implicit none

   character,        intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   character(*),     intent(in)    :: filename               ! I/O filename
   integer,          intent(in)    :: cv_size
   real*8,           intent(inout) :: rcv(1:cv_size)
   type(domain), optional, intent(inout) :: grid
   integer,      optional, intent(in)    :: mz(:)

   integer                         :: svd_unit, dummy

   if (trace_use) call da_trace_entry("da_cv_io")

   if (use_global_cv_io .and. present(grid) .and. present(mz)) then
      call da_cv_io_global (io_config, grid, filename, cv_size, rcv, mz)
   else 
      call da_get_unit (svd_unit)

      ! Unformatted
      if (io_config == 'r') then
         write(*,*) 'Reading randomized SVD cv vector from '//trim(filename)
         open (unit=svd_unit, file = trim(filename), form = 'unformatted', status = 'old')
         read(unit=svd_unit) dummy
         read(unit=svd_unit) rcv(1:cv_size)
         close(unit=svd_unit)
      else if (io_config == 'w') then
         write(*,*) 'Writing randomized SVD cv vector to '//trim(filename)
         open (unit=svd_unit, file = trim(filename), form = 'unformatted', status = 'replace')
         write(unit=svd_unit) cv_size
         write(unit=svd_unit) rcv(1:cv_size)
         close(unit=svd_unit)
      else
         write(*,*) 'Unknow configuration for randomized SVD I/O routine'
      end if

!! Use 'formatted' when cv file contents need to be ascii-type for debugging
!      if (io_config == 'r') then
!         write(*,*) 'Reading randomized SVD cv vector from '//trim(filename)
!         open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
!         read(unit=svd_unit,fmt='(I10)') dummy
!         read(unit=svd_unit,fmt='(E27.16E3)') rcv(1:cv_size)
!         close(unit=svd_unit)
!      else if (io_config == 'w') then
!         write(*,*) 'Writing randomized SVD cv vector to '//trim(filename)
!         open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
!         write(unit=svd_unit,fmt='(I10)') cv_size
!         write(unit=svd_unit,fmt='(E27.16E3)') rcv(1:cv_size)
!         close(unit=svd_unit)
!      else
!         write(*,*) 'Unknow configuration for randomized SVD I/O routine'
!      end if

      call da_free_unit (svd_unit)

   end if

   if (trace_use) call da_trace_exit ("da_cv_io")

end subroutine da_cv_io
