subroutine da_dgn(grid, config_flags, cv_size, xhat, cv, &
                  xbx, be, iv, ob, &
                  j_cost, ja_total, it, mz)

   !---------------------------------------------------------------------------
   ! Purpose: Uses line search to scale the cv increments; improves nonlinear convergence
   ! Reference: Guerrette and Henze (2017)
   ! Creation: 07/25/2016           JJ Guerrette 
   ! Modified: 12/10/2017           JJ Guerrette
   !  - Added independent treatment of met and chem variables
   !---------------------------------------------------------------------------

   implicit none

   type (domain),       intent(inout) :: grid
   integer, intent(in)                :: cv_size    ! Total cv size.
   type (xbx_type),intent(inout)      :: xbx    ! For header & non-grid arrays.
   type (be_type), intent(in)         :: be     ! background error structure.
   type (iv_type), intent(inout)      :: iv     ! innovation vector (o-b).
   real*8, intent(inout)              :: xhat(1:cv_size) ! control variables.
   real*8, intent(in)                 :: cv(1:cv_size)   ! control variables.
   type (y_type) , intent(inout)      :: ob          ! Observation structure.
   type(grid_config_rec_type), intent(inout) :: config_flags
   type (j_type), intent(in)          :: j_cost   ! cost function
   real, intent(in)                   :: ja_total ! Cost function at eta = 0
   integer, intent(in)                :: it 
   integer, intent(in)                :: mz(:)

   ! Need to allow ncv_type to be +1 when var4d_lbc
#if (WRF_CHEM == 1)
   integer, parameter    :: ncv_type = 2
#else
   integer, parameter    :: ncv_type = 1
#endif
   real*8  :: grat, jo_partial
   integer :: isearch, i
   real*8, dimension(1:ncv_type) :: &
             eta, eta_a, eta_b, eta_c, eta_d, &
             jout, ja, jb, jc, jd
   integer, dimension(1:ncv_type) :: &
             cv_s, cv_e, ncv, branch
   logical, dimension(1:ncv_type) :: &
             exit_array
   character(1) :: ncv_char

   integer :: jl_start, jl_end
   character*256 :: timestr, timestr1
   type(x_type) :: shuffle

   if (trace_use) call da_trace_entry("da_dgn")

   write(unit=stdout,fmt='("Starting Damped Gauss Newton (DGN) line search")')

   ! Determine the MET CV indices
   cv_s(1) = 1
   cv_e(1) = cv_s(1) &
          + be % cv % size1 &
          + be % cv % size2 &
          + be % cv % size3 &
          + be % cv % size4 &
          + be % cv % size5 &
#ifdef CLOUD_CV
          + be % cv % size6 &
          + be % cv % size7 &
          + be % cv % size8 &
          + be % cv % size9 &
#endif
          + be%cv%size_je &
          - 1

#if (WRF_CHEM == 1)
   ! Determine the CHEM CV indices
   cv_s(2) = cv_e(1) + 1
   cv_e(2) = cv_s(2) &
             + sum(be % cv % sizeant (:)) &
             + sum(be % cv % sizebb (:)) &
             - 1
#endif

!   ! Still need to figure this out, probably group MET with MET lbc and CHEM with CHEM lbc
!   if (var4d_lbc) then
!      jl_start    = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + 1
!      jl_end      = be%cv%size_jb + be%cv%size_je + be%cv%size_jp + be%cv%size_jl
!      cv_s(ncv_type+1) = jl_start
!      cv_e(ncv_type+1) = jl_end
!!#if (WRF_CHEM == 1)
!!      cv_e(ncv_type+1) = jl_end - sum(be % cv % sizechem (:))
!!      cv_s(ncv_type+2) = cv_e(ncv_type+1) + 1
!!      cv_e(ncv_type+2) = jl_end
!!#endif
!   end if

   write(unit=stdout,fmt='("cv_s, cv_e, cv_s_chem, cv_e_chem, cv_size, size_jb: ",6I8)') cv_s(1), cv_e(1), cv_s(2), cv_e(2), cv_size, be % cv % size_jb

   !---------------------------------------------------
   ! Perform line search using golden section algorithm
   !---------------------------------------------------
   grat = (sqrt(5.D0) - 1.D0) / 2.D0 ! Golden Ratio

   do i = 1, ncv_type
      ncv(i) = cv_e(i) - cv_s(i) + 1
      eta_a(i) = 0.D0
      eta_b(i) = 1.D0
      eta_c(i) = eta_b(i) - grat * (eta_b(i) - eta_a(i))
      eta_d(i) = eta_a(i) + grat * (eta_b(i) - eta_a(i))

      eta(i) = eta_d(i)

      exit_array(i) = .false.
      if ( ncv(i) .le. 0 ) then
         exit_array(i) = .true.
         eta(i) = 1.D0
      end if
   end do

   ! The value of ja is nonessential, but informative if available
   ja = 0.D0
   if (ja_total .gt. 0.D0) then
      ja(1) = ja_total
#if (WRF_CHEM == 1)
      if ( ncv(2) .gt. 0 ) then
         jo_partial = &
                    j_cost%jo%chem_acft &
                  + j_cost%jo%chem_surf 
         ja(2) = wrf_dm_sum_real(jo_partial)

         ja(2) = ja(2) + jb_factor * 0.5 * &
              da_dot_cv(ncv(2),  cv(cv_s(2):cv_e(2)), cv(cv_s(2):cv_e(2)), grid, mz)

         ja(1) = ja(1) - ja(2)
      end if
#endif
   end if
   jb = 0.D0
   jc = 0.D0
   jd = 0.D0

   ! Don't need to calculate initial jb (@ eta=1.) except to demonstrate need for DGN
   ! Could use jb/ja > 1 as criteria for needing line search
   call da_evaluate_j(grid, config_flags, cv_size, be%cv%size_jb, &
                         cv_s, cv_e, eta_b, jb, &
                         xbx, be, iv, xhat, cv, ob, &
                         it, mz, 1)

   write(unit=stdout,fmt='("DGN:   i  eta_b          jb")')
   do i = 1, ncv_type
      if ( ncv(i) .le. 0 ) cycle
      write(unit=stdout,fmt='("DGN:",I4,1F7.3,1F11.3)') i,eta_b(i),jb(i)
   end do

   call da_evaluate_j(grid, config_flags, cv_size, be%cv%size_jb, &
                         cv_s, cv_e, eta_c, jc, &
                         xbx, be, iv, xhat, cv, ob, &
                         it, mz, 2)

   write(unit=stdout,fmt='("DGN:   i  eta_c          jc")')
   do i = 1, ncv_type
      if ( ncv(i) .le. 0 ) cycle
      write(unit=stdout,fmt='("DGN:",I4,1F7.3,1F11.3)') i,eta_c(i),jc(i)
   end do

   branch = 2
   isearch = 2
   do while ( .not.all(exit_array) )
      isearch = isearch + 1
      call da_evaluate_j(grid, config_flags, cv_size, be%cv%size_jb, &
                            cv_s, cv_e, eta, jout, &
                            xbx, be, iv, xhat, cv, ob, &
                            it, mz, 2+isearch)

      write( unit=stdout,fmt=&
       '("DGN: icv   i  eta_a  eta_c  eta_d  eta_b       ja         jc         jd         jb")' )

      do i = 1, ncv_type
         if ( ncv(i) .le. 0 ) cycle

         if (branch(i) .eq. 1) then
            jc(i) = jout(i)
         else
            jd(i) = jout(i)
         end if
         write(unit=stdout,fmt='("DGN:",2I4,4F7.3,4F11.3)') &
            i, isearch, eta_a(i), eta_c(i), eta_d(i), eta_b(i), ja(i), jc(i), jd(i), jb(i)

         if (jc(i) .lt. jd(i)) then
            eta_b(i) = eta_d(i)
            jb(i) = jd(i)
            eta_d(i) = eta_c(i)
            jd(i) = jc(i)
            eta_c(i) = eta_b(i) - grat * (eta_b(i) - eta_a(i))
            eta(i) = eta_c(i)
            branch(i) = 1
         else
            eta_a(i) = eta_c(i)
            ja(i) = jc(i)
            eta_c(i) = eta_d(i)
            jc(i) = jd(i)
            eta_d(i) = eta_a(i) + grat * (eta_b(i) - eta_a(i))
            eta(i) = eta_d(i)
            branch(i) = 2
         end if

         if (.not.exit_array(i)) then
            if (ja(i).gt.1e-6 &
                   .and. abs(jb(i)-ja(i))/min(ja(i),jb(i),jd(i)).lt.dgn_tol(it)/10. &
                   .and. abs(jb(i)-jd(i))/min(ja(i),jb(i),jd(i)).lt.dgn_tol(it)/10.) then
               if (all(jb(i).lt.(/ja(i),jc(i),jd(i)/))) eta_a(i) = eta_b(i)
               exit_array(i) = .true.
            end if

            if (eta(i).gt.(1.d0 - dgn_tol(it))) then
               eta(i) = 1.D0
               exit_array(i) = .true.
            end if

            if ((eta_b(i) - eta_a(i)) .le. dgn_tol(it)) exit_array(i) = .true.
         end if
      end do
   end do

   do i = 1, ncv_type
      if ( ncv(i) .le. 0 ) cycle
      if (eta(i) .lt. 1.d0 ) eta(i) = (eta_a(i) + eta_b(i)) / 2.D0
      xhat(cv_s(i):cv_e(i)) = eta(i) * xhat(cv_s(i):cv_e(i))
   end do
   write(ncv_char,fmt='(I1)') ncv_type
   write(unit=stdout,fmt='(A,'//ncv_char//'F8.5)') 'Final eta from DGN: ',eta


   call da_transform_vtox (grid,be%cv%size_jb,xbx,be,grid%ep,&
                            xhat(1:be%cv%size_jb),grid%vv,grid%vp &
#if (WRF_CHEM == 1)
                            , grid%vch &
#endif
                            )

   call da_transform_xatowrftl_lbc_wrapper( &
           it, cv_size, be, grid%ep, &
           xhat, iv, grid%vp6, grid%vv6, xbx, &
           grid, config_flags &
!#if (WRF_CHEM == 1)
!          , grid%vch6 &
!#endif
           )

!Outputting fg_02 here would require restoring xb at initial time afterward
!   if ((config_flags%real_data_init_type == 1) .or. &
!       (config_flags%real_data_init_type == 3)) then
!#ifdef VAR4D
!      if ( var4d_lbc ) then
!         call domain_clock_get (grid, stop_timestr=timestr1)
!         call domain_clock_set( grid, current_timestr=timestr1 )
!         if (it .eq. 1) then
!            filename = 'fg02'
!         else
!            write(filename,fmt='("fg02_it",I2.2)') it-1
!         end if
!         call da_med_initialdata_input (grid, config_flags, trim(filename))
!         call da_setup_firstguess(xbx, grid, config_flags, .false. )
!         shuffle = grid%xa
!         grid%xa  = grid%x6a
!         call da_transform_vtox(grid, be%cv%size_jl, xbx, be, grid%ep, &
!              xhat(jl_start:jl_end), grid%vv6, grid%vp6)
!         grid%xa  = shuffle
!         call da_transfer_xatoanalysis (it, xbx, grid, config_flags)
!         write(filename,fmt='("fg02_it",I2.2)') it
!         call da_update_firstguess (grid, trim(filename))
!         call domain_clock_get (grid, start_timestr=timestr1)
!         call domain_clock_set( grid, current_timestr=timestr1 )
!      endif
!#endif
!   end if


   if (trace_use) call da_trace_exit("da_dgn")

end subroutine da_dgn
