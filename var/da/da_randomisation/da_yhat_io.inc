subroutine da_yhat_io (io_config, filename, nall, nobs, obs)

   !-------------------------------------------------------------------------
   ! Purpose:        Read / Write yhat_obs and qhat_obs
   !
   ! Called from da_randomise_svd_51
   !
   ! History: 09/01/2016  Creation (JJ Guerrette), modifed from da_lanczos_io.inc
   !
   !-------------------------------------------------------------------------

   implicit none
   
   character,         intent(in)    :: io_config              ! 'r' = Read; 'w' = Write
   character(*),      intent(in)    :: filename               ! I/O filename
   integer,           intent(inout) :: nall(num_ob_indexes)
   integer,           intent(inout) :: nobs
   real, allocatable, intent(inout) :: obs(:)

   integer                         :: svd_unit, dummy

   if (trace_use) call da_trace_entry("da_yhat_io")
   
   call da_get_unit (svd_unit)

!   ! Unformatted
!   if (io_config == 'r') then
!      write(*,*) 'Reading randomized SVD obs vector from '//trim(filename)
!      open (unit=svd_unit, file = trim(filename), form = 'unformatted', status = 'old')
!      read(unit=svd_unit) nobs
!      if (nobs .gt. 0) then
!         if (allocated(obs)) deallocate(obs)
!         allocate(obs(1:nobs))
!         read(unit=svd_unit) obs(1:nobs)
!         read(unit=svd_unit) nall
!      else
!         nall(:) = 0
!      end if
!      close(unit=svd_unit)
!   else if (io_config == 'w') then
!      write(*,*) 'Writing randomized SVD obs vector to '//trim(filename)
!      open (unit=svd_unit, file = trim(filename), form = 'unformatted', status = 'replace')
!      write(unit=svd_unit) nobs
!      if (nobs .gt. 0 .and. allocated(obs) ) then
!         if (size(obs) .eq. nobs) then
!            write(unit=svd_unit) obs(1:nobs)
!            write(unit=svd_unit) nall
!         else
!            write(*,*) 'ERROR in da_yhat_io: obs vector is wrong size', size(obs), nobs
!         end if
!      end if

   if (io_config == 'r') then
      write(*,*) 'Reading randomized SVD obs vector from '//trim(filename)
      open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'old',access='sequential')
      read(unit=svd_unit,fmt='(I10)') nobs
      if (nobs .gt. 0) then
         if (allocated(obs)) deallocate(obs)
         allocate(obs(1:nobs))
         read(unit=svd_unit,fmt='(E27.16E3)') obs(1:nobs)
         read(unit=svd_unit,fmt='(I10)') nall
      else
         nall(:) = 0
      end if
      close(unit=svd_unit)
   else if (io_config == 'w') then
      write(*,*) 'Writing randomized SVD obs vector to '//trim(filename)
      open (unit=svd_unit, file = trim(filename), form = 'formatted', status = 'replace',access='sequential')
      write(unit=svd_unit,fmt='(I10)') nobs
      if (nobs .gt. 0 .and. allocated(obs) ) then
         if (size(obs) .eq. nobs) then
            write(unit=svd_unit,fmt='(E27.16E3)') obs(1:nobs)
            write(unit=svd_unit,fmt='(I10)') nall
         else
            write(*,*) 'da_yhat_io: obs vector is wrong size', size(obs), nobs
         end if
      end if
      close(unit=svd_unit)
   else
      write(*,*) 'Unknow configuration for randomized SVD I/O routine'
   end if

   call da_free_unit (svd_unit)

   if (trace_use) call da_trace_exit ("da_yhat_io")

end subroutine da_yhat_io
