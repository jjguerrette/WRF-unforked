subroutine da_read_obs_chem_again ( iv, ob, it, model_grid )
   !-------------------------------------------------------------------------
   ! Purpose: Read modeled observations from ascii files that were previously
   !          calculated (specifically modeled k_at_p)
   !
   ! Called from da_solve
   !
   ! History: 09/2016  Creation (Jonathan Guerrette)
   !
   !-------------------------------------------------------------------------

   type(domain),      intent(inout) :: model_grid
   type (iv_type),    intent(inout) :: iv
   type (y_type),     intent(inout) :: ob
   integer, intent(in)              :: it
   integer                                    :: iobs, ichem, itime, iplt, iloc, jloc, kloc, &
                                                 iunit, iostat, ij, ilocal, dumi
   real, dimension(:,:), allocatable :: dummy_surf
   real, dimension(:,:), allocatable :: dummy_acft

   real :: k_at_p, dumr
   real, allocatable :: conc(:)

   character*256         :: surf_file! = "SURFACE_Hx_y.txt"
   character*256         :: acft_file! = "AIRCRAFT_Hx_y.txt"  

   character*3         :: dummychar1="", dummychar2=""

   if (num_surf_obs > 0) then
      call da_message((/"Reading in chem SURFACE observations and model concentrations to avoid NL model integration"/))

      write(surf_file,'("SURFACE_Hx_y_i",I2.2,".txt")') it
      if(rootproc) then
         ! =========
         ! open file
         ! =========
         call da_get_unit(iunit)  
         open(unit   = iunit,     &
            FILE   = trim(surf_file), &
            FORM   = 'FORMATTED',  &
            ACCESS = 'SEQUENTIAL', &
            iostat =  iostat,     &
            STATUS = 'old', &
            POSITION = 'rewind')
   
         if (iostat /= 0) then
            write(unit=message(1),fmt='(A,I5,2A)') "Error message: ",iostat,"; ",trim(surf_file)
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening surface Hx_y file')
         end if
      end if

      allocate ( dummy_surf (2,num_chem_surf-1) )
      allocate ( conc (num_surf_hx-1) )

      ilocal = 0
      do iobs = 1, iv%info(chem_surf)%ntotal !num_surf_obs
         if ( .NOT. (head_grid%ti_surf(iobs) .GE. 1 .AND. head_grid%tf_surf(iobs) .LE. num_ts) ) cycle
         ilocal = ilocal + 1

         if (rootproc) then
            write(dummychar1,fmt='(I0)') num_chem_surf-1
            write(dummychar2,fmt='(I0)') num_surf_hx-1
            read(iunit,fmt='(3I6,D14.5,' // trim(dummychar1) // 'D14.5,' // trim(dummychar1) // 'D23.14,' // trim(dummychar2) // 'D23.14)') &
               iloc, jloc, kloc, &
               dumr, &
               dummy_surf(1,:), &
               dummy_surf(2,:), &
               conc !!head_grid%surf_hx(iobs,2:num_surf_hx)
            if (iostat /= 0 ) then
               write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
               call da_warning(__FILE__,__LINE__,message(1:1))
               call da_wrf_error_fatal('Error reading surface chem Hx_y')
            end if
         end if

#ifdef DM_PARALLEL
         call mpi_barrier(comm, ierr)

         call wrf_dm_bcast_real(conc, num_surf_hx-1)

         model_grid%surf_hx(iobs,2:num_surf_hx) = conc
         head_grid%surf_hx(iobs,2:num_surf_hx) = conc
#endif

      end do

      deallocate ( dummy_surf )
      deallocate ( conc )

      if (rootproc) then
         close(iunit)
         call da_free_unit(iunit)
      end if
   end if
#ifdef DM_PARALLEL
   call mpi_barrier(comm, ierr)
#endif

   if (num_acft_obs(1) > 0) then  

      call da_message((/"Reading in chem AIRCRAFT observations and model concentrations to avoid NL model integration"/))

      write(acft_file,'("AIRCRAFT_Hx_y_i",I2.2,".txt")') it

      if(rootproc) then
         ! =========
         ! open file
         ! =========
         call da_get_unit(iunit)  
         open(unit   = iunit,     &
            FILE   = trim(acft_file), &
            FORM   = 'FORMATTED',  &
            ACCESS = 'SEQUENTIAL', &
            iostat =  iostat,     &
            STATUS = 'old', &
            POSITION = 'rewind')
      
         if (iostat /= 0) then
            write(unit=message(1),fmt='(A,I5,2A)') "Error message: ",iostat,"; ",trim(acft_file)
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening aircraft Hx_y file')
         end if   
      end if

      allocate ( dummy_acft (2,num_chem_acft-1) )
      allocate ( conc (num_acft_hx-1) )

      do iplt = 1, num_platform
         ilocal = 0
         do itime = 1, iv%info(chem_acft)%ntotal
            if( .NOT. head_grid%acft_obs_count(itime,iplt) > 0 ) cycle
            ilocal = ilocal + 1

            if(rootproc) then
               write(dummychar1,fmt='(I0)') num_chem_acft-1
               write(dummychar2,fmt='(I0)') num_acft_hx-1
               read(iunit,fmt='(3I6,F12.6,3I6,' // trim(dummychar1) // 'D14.5,' // trim(dummychar1) // 'D23.14,' // trim(dummychar2) // 'D23.14)') &

                  iloc, jloc, kloc, k_at_p, dumi, dumi, &
                  dumi, &
                  dummy_acft(1,:), &
                  dummy_acft(2,:), &
                  conc !!head_grid%acft_hx(itime,iplt,2:num_acft_hx)
   
               if (iostat /= 0 ) then
                  write(unit=message(1),fmt='(A,I5)') "Error message: ",iostat
                  call da_warning(__FILE__,__LINE__,message(1:1))
                  call da_wrf_error_fatal('Error reading aircraft chem Hx_y')
               end if
            end if
#ifdef DM_PARALLEL
            call mpi_barrier(comm, ierr)

            call wrf_dm_bcast_real(k_at_p, 1)
            head_grid%k_acft(itime,iplt) = k_at_p
            model_grid%k_acft(itime,iplt) = k_at_p

            call wrf_dm_bcast_real(conc, num_acft_hx-1)
            model_grid%acft_hx(itime,iplt,2:num_acft_hx) = conc
            head_grid%acft_hx(itime,iplt,2:num_acft_hx) = conc
#endif

         end do
      end do

      deallocate ( dummy_acft )
      deallocate ( conc )

      if(rootproc) then
         close(iunit)
         call da_free_unit(iunit)
      end if !rootproc
   end if

#ifdef DM_PARALLEL
   call mpi_barrier(comm, ierr)
#endif

end subroutine da_read_obs_chem_again

