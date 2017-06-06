subroutine da_calculate_aminusb(grid, config_flags,   &
                                it, cv_size, xbx, be, &
                                qhat, eignvec, eignval, neign,nmodes, lcin)

   !--------------------------------------------------------------------------------
   ! Purpose:  Calculate and write out difference between prior and posterior
   !              variance using Tridiagonal eigenpairs and Lanczos orthonormal
   !              vectors (qhat)
   !              Note: Process is general for all variable types, but
   !                    netcdf write is only written for exponential 
   !                    emission scaling factors --> need non-chem still
   !
   ! Called from da_solve
   !
   ! History: 12/2015  Creation (JJ Guerrette)
   !
   !--------------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   type(grid_config_rec_type), intent(inout) :: config_flags
   integer, intent(in)               :: it                           ! external iteration.
   integer, intent(in)               :: cv_size                      ! Total cv size
   integer, intent(in)               :: neign
   integer, intent(in)               :: nmodes

   type (xbx_type),intent(inout)     :: xbx                          ! Header & non-gridded vars.
   type (be_type), intent(in)        :: be                           ! background error structure.
!   type (iv_type), intent(inout)     :: iv                           ! ob. increment vector.
!   real, intent(out)                 :: xhat(1:cv_size)              ! Output vector: xhat=K.d
   real, intent(in)                  :: qhat(1:cv_size, 0:nmodes) ! Ritz vectors
!   real, intent(in)                  :: cv(1:cv_size)                ! control variable (local).
!   type (y_type), intent(inout)      :: re                           ! residual (o-a) structure.
!   type (y_type), intent(inout)      :: y                            ! y = H(x_inc) structure.
!   type (j_type), intent(out)        :: j_cost                       ! cost function
   real*8, intent(in)                :: eignvec(nmodes, nmodes)
   real*8, intent(in)                :: eignval(nmodes)
   character(len=3), intent(in), optional :: lcin


   real*8                            :: nu(1:be%cv%size_jb)  ! Premultiplied eigenvectors, Q_l * e_lj
#if (WRF_CHEM == 1)
   real*8, allocatable               :: Pa_ant(:,:,:,:), Pa_bb(:,:,:,:)  ! Posterior Variance
#endif

   real*8                            :: lam                          ! (A-B) eigenvalues
   integer :: jj,ii,i,j,k,ic, ij

   character(len=2) :: itc
   character(len=3) :: kc
   character(len=3) :: lc

   if (trace_use) call da_trace_entry("da_calculate_aminusb")

   write(*,*) 'Computing Analysis Variance Reduction: v = A-B = (B^-1+H^T.R^-1.H)^-1' !!! ~ SUM(lam^-1 - 1)(U*nu)(U*nu)^T'

   WRITE(unit=itc,fmt='(I2.2)') it
   if (present(lcin)) then
      lc = lcin
   else
      WRITE(unit=lc,fmt='(I3.3)') neign
   end if

#if (WRF_CHEM == 1)
   allocate(Pa_ant(its:ite,jts:jte,1:num_ant_steps,1:num_scaleant))
   allocate(Pa_bb(its:ite,jts:jte,1:num_bb_steps,1:num_scalebb))

   ! Initialize Variance Reduction to zero.
   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, ic )
   do ij = 1 , grid%num_tiles
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         grid%g_scaleant(its:ite,jts:jte,:,ic) = 0.D0
         Pa_ant(its:ite,jts:jte,:,ic) = 0.D0
      end do
   end do
   !$OMP END PARALLEL DO

   !$OMP PARALLEL DO &
   !$OMP PRIVATE ( ij, ic )
   do ij = 1 , grid%num_tiles
      do ic = PARAM_FIRST_SCALAR, num_scalebb
         grid%g_scalebb(its:ite,jts:jte,:,ic) = 0.D0
         Pa_bb(its:ite,jts:jte,:,ic) = 0.D0
      end do
   end do
   !$OMP END PARALLEL DO
#endif

   ! Calculate Variance Reduction
   do jj = 1, neign
      nu(:) = 0.D0
      ! Calculate nu_j = Q_l * e_lj
      ! l = neign = inner loop count
      ! Q_l = qhat(:,1:neign)
      ! e_lj = eignvec(1:neign,jj)
      ! lambda_j = eignval(jj)
      do ii = 1, neign
         nu = nu + qhat(1:be%cv%size_jb, ii) * eignvec(ii,jj)
      end do

      ! Calculate U * nu_j and store in grid%xa,grid%xach
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                             nu(1:be%cv%size_jb), grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                             , grid%vch &
#endif
                             )
      lam = ( 1.D0 / eignval(jj)  - 1.D0 )

#if (WRF_CHEM == 1)
      ! Store Uv in corresponding I/O capable variables
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
                     grid%g_scaleant(i,j,k,ic) = grid%xach%alpha_ant(i,j,k,ic)
               end do
            end do
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         do k = 1, num_bb_steps
            do j=jts,jte
               do i=its,ite
                     grid%g_scalebb(i,j,k,ic) = grid%xach%alpha_bb(i,j,k,ic)
               end do
            end do
         end do
      end do

      WRITE(unit=kc,fmt='(I3.3)') jj

      config_flags%auxhist18_outname = "Uv_k"//kc//"_l"//lc//"_it"//itc
      config_flags%io_form_auxhist18 = 2
      CALL nl_set_io_form_auxhist18 ( grid%id, 2 )
      config_flags%frames_per_auxhist18 = 1
      CALL  med_hist_out ( grid , AUXHIST18_ALARM , config_flags )

      ! Store variances in Pa arrays
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
                  Pa_ant(i,j,k,ic) = Pa_ant(i,j,k,ic) + lam * grid%xach%alpha_ant(i,j,k,ic) ** 2
               end do
            end do
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         do k = 1, num_bb_steps
            do j=jts,jte
               do i=its,ite
                  Pa_bb(i,j,k,ic) = Pa_bb(i,j,k,ic) + lam * grid%xach%alpha_bb(i,j,k,ic) ** 2
               end do
            end do
         end do
      end do
#endif
   end do

#if (WRF_CHEM == 1)
   ! Store result in corresponding I/O capable variables
   do ic = PARAM_FIRST_SCALAR, num_scaleant
      do k = 1, num_ant_steps
         do j=jts,jte
            do i=its,ite
                  grid%g_scaleant(i,j,k,ic) = Pa_ant(i,j,k,ic)
            end do
         end do
      end do
   end do

   do ic = PARAM_FIRST_SCALAR, num_scalebb
      do k = 1, num_bb_steps
         do j=jts,jte
            do i=its,ite
                  grid%g_scalebb(i,j,k,ic) = Pa_bb(i,j,k,ic)
            end do
         end do
      end do
   end do

   deallocate(Pa_ant)
   deallocate(Pa_bb)

   config_flags%auxhist18_outname = "AminusB_diag_it"//itc
   config_flags%io_form_auxhist18 = 2
   CALL nl_set_io_form_auxhist18 ( grid%id, 2 )
   config_flags%frames_per_auxhist18 = 1
   CALL  med_hist_out ( grid , AUXHIST18_ALARM , config_flags )
#endif

   CALL nl_set_io_form_auxhist18 ( grid%id, 0 )
   config_flags%auxhist18_outname = "auxhist18_d<domain>_<date>"


   if (trace_use) call da_trace_exit ("da_calculate_aminusb")

end subroutine da_calculate_aminusb
