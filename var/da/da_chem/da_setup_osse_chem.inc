subroutine da_setup_osse_chem( grid, ob, iv, be, xbx)

   !---------------------------------------------------------------------------
   ! Purposes: 
   !  (1) Calculate sigma (relative + absolute variances) from truth (iv)
   !  (2) Calculate perturbed observation using sigma and random vector epsilon_y (local variable iv_pert = iv + pert)
   !  (3) Assign y+dy to ob (call da_fill_obs_structures_chem(iv_pert, ob))
   !  (4) Perturb emissions using random vector epsilon_x and U from B=UU^T
   !---------------------------------------------------------------------------

   implicit none

   type (domain),                  intent(inout) :: grid
   type (y_type),                  intent(inout) :: ob     ! Observation structure.
   type (iv_type),                 intent(inout) :: iv     ! O-B structure.
   type(be_type),                  intent(in)    :: be     ! background error structure.
   type(xbx_type),                 intent(inout) :: xbx    ! For header & non-grid arrays.

   real                        :: omega(1:be%cv%size_jb)   ! random CV (local).

!!   type (iv_type)              :: iv_pert          ! O-B structure.

   integer                     :: i, j, k, ic, itime, ichem, iplt, iobs, &
                                  icount, iostat, iunit
   character*10                :: cproc
   character(len=filename_len) :: filename
   real, allocatable           :: e_error(:), o_error(:), dummy_out(:)
   real                        :: sigmaval, chemval
   character*256               :: o_file = "O_ERROR.txt"
   character*256               :: e_file = "E_ERROR.txt"

   if (trace_use) call da_trace_entry("da_setup_osse_chem")


!!   iv_pert%info(chem_surf)%ntotal = 0
!!   iv_pert%info(chem_surf)%nlocal = 0
!!
!!   iv_pert%info(chem_acft)%ntotal = num_ts
!!   iv_pert%info(chem_acft)%nlocal = num_ts
!!
!!   call da_allocate_observations_chem (iv_pert)
!!
!!   allocate( iv_pert%chem_acft (num_ts) )
!================================================================================


   ! [1] Initialize random number generator
   call da_random_seed

   ! [2] Initialize obs gaussian error
   allocate( o_error (num_ts * num_platform * (num_chem_acft-1)) )
   o_error(:) = 0.D0

!!   if ( use_obs_error ) then
      ! [3] Calculate obs gaussian error on root node
      if(rootproc) then
         icount = 0
         do iplt = 1,num_platform
            do itime = 1, num_ts
               do ichem = 1, num_chem_acft-1
                  if ( grid%acft_obs_count(itime,iplt) .gt. 0 ) then
                     icount = icount + 1
                     call da_gauss_noise(o_error(icount))
                  end if
               end do
            end do
         end do

         ! [4] Write obs gaussian error to file
         ! =========
         ! open file
         ! =========
         call da_get_unit(iunit)
         open(unit   = iunit,     &
            FILE   = trim(o_file), &
            FORM   = 'FORMATTED',  &
            ACCESS = 'SEQUENTIAL', &
            iostat =  iostat,     &
            STATUS = 'REPLACE')

         if (iostat /= 0) then
            write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening emission error file')
         end if

         do iobs=1, num_ts * num_platform * (num_chem_acft-1)
            write (iunit, fmt = '(F20.15)')   o_error(iobs)
            if (iostat /= 0) then
               write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
               call da_warning(__FILE__,__LINE__,message(1:1))
               call da_wrf_error_fatal('Error reading emission error file')
            end if
         end do
         close(iunit)
         call da_free_unit(iunit)
      end if

      ! [5] Write obs gaussian errors to other nodes
      allocate( dummy_out (num_ts * num_platform * (num_chem_acft-1)) )
#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
      call wrf_dm_sum_reals(o_error,  dummy_out)
      o_error = dummy_out
#endif
      deallocate( dummy_out )

!!   end if

   ! [6] Calculate uncertainty and perturbed observation   
   icount = 0
   do iplt = 1,num_platform
      do itime = 1, num_ts
         do ichem = 1, num_chem_acft-1
            if ( grid%acft_obs_count(itime,iplt) .gt. 0 ) then
               icount = icount + 1
               chemval = grid%acft_hx(itime,iplt,ichem+1)
               sigmaval = osse_chem_rel * chemval + osse_chem_abs
               iv%chem_acft(itime)%chem(ichem,iplt)%inv = max(chemval + sigmaval * o_error(icount), osse_chem_abs)
               iv%chem_acft(itime)%chem(ichem,iplt)%error = sigmaval
               iv%chem_acft(itime)%chem(ichem,iplt)%qc = 0
            else
               iv%chem_acft(itime)%chem(ichem,iplt)%inv = missing_r
               iv%chem_acft(itime)%chem(ichem,iplt)%error = missing_r
!               iv%chem_acft(itime)%chem(ichem,iplt)%error = 0
               iv%chem_acft(itime)%chem(ichem,iplt)%qc = -77
               iv%chem_acft(itime)%weight(ichem,iplt) = missing_r
            end if
         end do
      end do
   end do

   deallocate( o_error )

!   ! [7] Transfer perturbed iv to ob structure
!   call da_deallocate_y(ob)
!   call da_allocate_y(iv, ob)
!   call da_allocate_y_chem(iv, ob)
!
!   call da_fill_obs_structures_chem(iv, ob)

   ! Write perturbed observations and assoc. variances
   call da_read_obs_chem_multiplat(grid, iv, 'w')


!================================================================================
!================================================================================
!NEXT FOR EMISSIONS
   write(cproc,fmt='(i4.4)') myproc
   if (read_omega) then
      ! Read in omega vector portion for each processor
      write(unit=stdout,fmt='("Reading omega ensemble member : ")') 
      filename = 'e_error.p'//trim(adjustl(cproc))
      call da_cv_io ('r', filename, be%cv%size_jb, omega(1:be%cv%size_jb))
   else
      call da_set_randomcv(be%cv%size_jb, omega)

      !Write omega to file portion for each processor
      filename = 'e_error.p'//trim(adjustl(cproc))
      call da_cv_io ('w', filename, be%cv%size_jb, omega(1:be%cv%size_jb))
   end if


!NEED ...

   call da_transform_vtox(grid, be%cv%size_jb, xbx, be, &
                          grid%ep, omega(1:be%cv%size_jb), &
                          grid%vv, grid%vp, grid%vch)

!!   if (random_emiss_error) then

      do ic = PARAM_FIRST_SCALAR, num_scaleant
         do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
                  grid%scaleant(i,j,k,ic) = exp( grid%xach%alpha_ant(i,j,k,ic) )
               end do
            end do
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         do k = 1, num_bb_steps
            do j=jts,jte
               do i=its,ite
                  grid%scalebb(i,j,k,ic) = exp( grid%xach%alpha_bb(i,j,k,ic) )
               end do
            end do
         end do
      end do

!!   else
!!      do i = PARAM_FIRST_SCALAR, num_scaleant
!!         grid%scaleant(:,:,:,i) = REAL ( init_scale * scale_scale, 8 )
!!      end do
!!
!!      do i = PARAM_FIRST_SCALAR, num_scalebb
!!         grid%scalebb(:,:,:,i) = REAL ( init_scale * scale_scale, 8 )
!!      end do
!!   end if


   if (trace_use) call da_trace_exit("da_setup_osse_chem")

end subroutine da_setup_osse_chem
