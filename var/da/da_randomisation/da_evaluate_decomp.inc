#if defined(LAPACK)
subroutine da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                              rank_rsvd_Pa, rankA, neign, ntrunc, &
                              xhat_hess, ghat0, &
                              omega_Pav, &
                              yhat_Pav_LRU, yhat_Pav_LRA, &
                              yhat_Pa_LRU, yhat_Pa_LRA, &
                              La_reduced, &
                              eignval, eignvec, qhat, &
                              eignval_hess, eignvec_hess, qhat_hess, &
                              dvnorm_hess, dxnorm_hess, &
                              EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                              DIRNORM1_v, DIRNORM1_x, &
                              DIRNORM2_v, DIRNORM2_x)


   !-------------------------------------------------------------------------
   ! Purpose:   Evaluate norms for Pa and increment approximations
   !
   ! Called from da_compare_decomp_methods
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (xbx_type),intent(inout) :: xbx     ! Header & non-gridded vars.
   type (be_type), intent(in)    :: be      ! background error structure.
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: rank_rsvd_Pa, rankA, neign, ntrunc

   real*8,         intent(in) :: eignvec(neign, neign)
   real*8,         intent(in) :: eignval(neign)
   real*8,         intent(in) :: qhat(cv_size, 0:neign)

   real*8,         intent(in) :: eignvec_hess(rankA, rankA)
   real*8,         intent(in) :: eignval_hess(rankA)
   real*8,         intent(in) :: qhat_hess(cv_size, 0:rankA)

   real*8,         intent(in) :: xhat_hess(cv_size)
   real*8,         intent(in) :: ghat0(cv_size)

   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)

   !Check Variables
   real(wp), intent(in), allocatable        :: omega_Pav(:,:) ! Gaussian draw
   real(wp), intent(in), allocatable        :: yhat_Pav_LRU(:,:)  ! yhat = Pav * omega_Pa
   real(wp), intent(in), allocatable        :: yhat_Pav_LRA(:,:)  ! yhat = Pav * omega_Pa
   real(wp), intent(in), allocatable        :: yhat_Pa_LRU(:,:)  ! yhat = Pa * omega_Pa
   real(wp), intent(in), allocatable        :: yhat_Pa_LRA(:,:)  ! yhat = Pa * omega_Pa
   real(wp), intent(in), allocatable        :: La_reduced(:,:)  !  La^-1 = U^-1 Vm * (Im + LAM_m)^1/2 * Vm^T

   real(wp), intent(in), allocatable        :: dvnorm_hess(:), dxnorm_hess(:)

   !Metric outputs
   real*8, intent(out)        :: EUCNORM_x(3), FROBNORM_Pa(3), SPECNORM_Pa(3)
   real*8, intent(out)        :: EUCNORM_v(3), FROBNORM_Pav(3), SPECNORM_Pav(3)
   real*8, intent(out)        :: EUCNORM_xPa(3), FROBNORM_PaPa(3), SPECNORM_PaPa(3)
   real*8, intent(out)        :: DIRNORM1_v(3), DIRNORM1_x(3)
   real*8, intent(out)        :: DIRNORM2_v(3), DIRNORM2_x(3)

!!   real*8, intent(out)        :: EUCNORM_xB(3), FROBNORM_PaB(3), SPECNORM_PaB(3)



   real(wp)        :: UTomega(cv_size)
   real*8          :: xhat(cv_size)
   real*8          :: diffvec(cv_size)
   real*8          :: dvnorm_(cv_size)
   real*8          :: dxnorm_(cv_size)

   real(wp), allocatable  :: qhat_Pav(:,:)
   real(wp), allocatable  :: yhat_Pav_(:,:)

   real(wp), allocatable  :: qhat_Pa(:,:)
   real(wp), allocatable  :: yhat_Pa_(:,:)

   real(wp), allocatable  :: LaQ(:,:)
   real(wp), allocatable  :: LaVlS(:,:)
   real(wp), allocatable  :: LaVl(:,:)
   real(wp), allocatable  :: LaPaLa_reduced(:,:)


   real*8          :: temp
   real*8          :: eignvec_res(rank_rsvd_Pa, rank_rsvd_Pa)
   real*8          :: eignval_res(rank_rsvd_Pa), LAM(neign), LAM_hess(rankA)

   integer         :: iens, jens, i, j, k, ij, ic, itype, dummy, info, neign_trunc

   if (trace_use) call da_trace_entry("da_evaluate_decomp")

   write(unit=stdout,fmt=*) ' '
   write(unit=stdout,fmt='(A,I6,A,I6)') 'Evaluate Norms, k = ', ntrunc, ', k+p = ', neign


   EUCNORM_x=0.D0; FROBNORM_Pa=0.D0; SPECNORM_Pa=0.D0

   EUCNORM_v=0.D0; FROBNORM_Pav=0.D0; SPECNORM_Pav=0.D0

   EUCNORM_xPa=0.D0; FROBNORM_PaPa=0.D0; SPECNORM_PaPa=0.D0

   DIRNORM1_v=0.D0; DIRNORM1_x=0.D0
   DIRNORM2_v=0.D0; DIRNORM2_x=0.D0

   if (neign .lt. rankA) then
      neign_trunc = neign
   else
      neign_trunc = ntrunc
   end if

if (hess_comp_x) then
   !=====================================================================
   ! Compare xhat_approx to xhat_hess
   !  (1a)  v - increment in preconditioned v-space
   !  (1b)  x - increment in state x-space
   !  (1c)  x _ Pa^-1 - increment in state x-space normalized by Pa
   !=====================================================================

   !=====================================================================
   ! Compare xhat to xhat_hess and Uxhat to Uxhat_hess (vtox)
   !=====================================================================
   write(unit=stdout,fmt='(A)') 'Evaluating (xa_ - xa)'
   do itype = 1, 2
      xhat = 0.D0
      ! ||dva_ - dva||^2
      if (itype .eq. 1) then !LRA
         write(unit=stdout,fmt='(A)') 'LRA'
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
      end if

      if (itype .eq. 2) then !LRU
         write(unit=stdout,fmt='(A)') 'LRU'
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
      end if

      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             neign,neign,neign_trunc, &
                             LAM,eignvec,qhat,&
                             ghat0,xhat)

      if (itype .eq. 2) then !LRU
         xhat = xhat + ghat0
      end if

      !----------------------
      ! (1a) ||dva_ - dva||^2
      !----------------------
      diffvec = xhat - xhat_hess
      EUCNORM_v(itype) = da_dot_cv(cv_size, diffvec, diffvec, grid, mz, jp_start, jp_end)

      !----------------------
      ! (1b) ||dxa_ - dxa||^2
      !----------------------
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, &
                             grid%ep, diffvec(1:be%cv%size_jb), &
                             grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                             , grid%vch &
#endif
                             )

      temp = 0.D0
#if (WRF_CHEM == 1)
      !Only for chem so far
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
                  temp = temp + grid%xach%alpha_ant(i,j,k,ic)**2
               end do
            end do
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         do k = 1, num_bb_steps
            do j=jts,jte
               do i=its,ite
                  temp = temp + grid%xach%alpha_bb(i,j,k,ic)**2
               end do
            end do
         end do
      end do
#endif

#ifdef DM_PARALLEL
      EUCNORM_x(itype) = wrf_dm_sum_real(temp)
#else
      EUCNORM_x(itype) = temp
#endif

      !------------------------------
      ! (1c) ||dxa_ - dxa||^2_{Pa^-1}
      !------------------------------
      do i = 1, rankA
         LAM_hess(i) = eignval_hess(i) - 1.D0
      end do

      dvnorm_ = 0.D0
      !Multiply Pav^-1 * (dva_ - dva)
      ! L^T H^T R^-1 H L * diffvec
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             rankA,rankA,rankA, &
                             LAM_hess,eignvec_hess,qhat_hess,&
                             diffvec,dvnorm_,.true.)

      ! (I + L^T H^T R^-1 H L) * diffvec
      dvnorm_ = dvnorm_ + diffvec

      !Multiply (dva_ - dva)^T * Pav^-1 * (dva_ - dva) = ||dxa_ - dxa||^2_{Pa^-1}
      EUCNORM_xPa(itype) = da_dot_cv(cv_size, diffvec, dvnorm_, grid, mz, jp_start, jp_end)


   !=====================================================================
   ! Compare directions of xhat/xhat_hess and Uxhat/Uxhat_hess
   !  (1d)  v - increment in preconditioned v-space
   !  (1e)  x - increment in state x-space
   !=====================================================================
      !------------------------------------------------
      ! (1d) [ dv_^T/sqrt(dv_^T dv) dv/sqrt(dv^T dv) ]
      !------------------------------------------------
      if (allocated(dvnorm_hess)) then
         temp = sqrt(da_dot_cv(cv_size, xhat, xhat, grid, mz, jp_start, jp_end))
         dvnorm_ = xhat / temp
         DIRNORM1_v(itype) = da_dot_cv(cv_size, dvnorm_, dvnorm_hess, grid, mz, jp_start, jp_end)

         dvnorm_ = dvnorm_ - dvnorm_hess
         DIRNORM2_v(itype) = da_dot_cv(cv_size, dvnorm_, dvnorm_, grid, mz, jp_start, jp_end)
      end if

      !------------------------------------------------
      ! (1e) [ dx_^T/sqrt(dx_^T dx) dv/sqrt(dx^T dx) ]
      !------------------------------------------------
      if (allocated(dxnorm_hess)) then
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, &
                                grid%ep, xhat(1:be%cv%size_jb), &
                                grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                                , grid%vch &
#endif
                                )

         dxnorm_ = 0.D0
         ij = 0
#if (WRF_CHEM == 1)
         !Only for chem so far
         do ic = PARAM_FIRST_SCALAR, num_scaleant
            do k = 1, num_ant_steps
               do j=jts,jte
                  do i=its,ite
                     ij = ij + 1
                     if (ij .gt. cv_size) &
                        call da_error(__FILE__,__LINE__, &
                            (/"cv_size is incorrect"/))

                     dxnorm_(ij) = grid%xach%alpha_ant(i,j,k,ic)
                  end do
               end do
            end do
         end do

         do ic = PARAM_FIRST_SCALAR, num_scalebb
            do k = 1, num_bb_steps
               do j=jts,jte
                  do i=its,ite
                     ij = ij + 1
                     if (ij .gt. cv_size) &
                        call da_error(__FILE__,__LINE__, &
                            (/"cv_size is incorrect"/))

                     dxnorm_(ij) = grid%xach%alpha_bb(i,j,k,ic)
                  end do
               end do
            end do
         end do
#endif
         temp = sqrt(da_dot_cv(cv_size, dxnorm_, dxnorm_, grid, mz, jp_start, jp_end))
         dxnorm_ = dxnorm_ / temp
         DIRNORM1_x(itype) = da_dot_cv(cv_size, dxnorm_, dxnorm_hess, grid, mz, jp_start, jp_end)

         dxnorm_ = dxnorm_ - dxnorm_hess
         DIRNORM2_x(itype) = da_dot_cv(cv_size, dxnorm_, dxnorm_, grid, mz, jp_start, jp_end)
      end if

      write(unit=stdout,fmt='(A)') '    EUCNORM_v     EUCNORM_x   EUCNORM_xPa    DIRNORM1_v    DIRNORM1_x    DIRNORM2_v    DIRNORM2_x'
      write(unit=stdout,fmt='(7D14.4)') EUCNORM_v(itype), EUCNORM_x(itype), EUCNORM_xPa(itype), DIRNORM1_v(itype), DIRNORM1_x(itype), DIRNORM2_v(itype), DIRNORM2_x(itype)

   end do

   ! Adaptive increment
   if ( minval(eignval(1:ntrunc)) .lt. 2.D0 ) then
      EUCNORM_v(3) = EUCNORM_v(2)
      EUCNORM_x(3) = EUCNORM_x(2)
      EUCNORM_xPa(3) = EUCNORM_xPa(2)
      DIRNORM1_v(3) = DIRNORM1_v(2)
      DIRNORM1_x(3) = DIRNORM1_x(2)
      DIRNORM2_v(3) = DIRNORM2_v(2)
      DIRNORM2_x(3) = DIRNORM2_x(2)
   else
      EUCNORM_v(3) = EUCNORM_v(1)
      EUCNORM_x(3) = EUCNORM_x(1)
      EUCNORM_xPa(3) = EUCNORM_xPa(1)
      DIRNORM1_v(3) = DIRNORM1_v(1)
      DIRNORM1_x(3) = DIRNORM1_x(1)
      DIRNORM2_v(3) = DIRNORM2_v(1)
      DIRNORM2_x(3) = DIRNORM2_x(1)
   end if
end if

   !=====================================================================
   ! Compare Pa_approx to Pa_hess
   !  (2a)  Pav - posterior covariance in preconditioned v-space
   !  (2b)  Pa  - posterior covariance in state x-space
   !  (2c)  Pa _ Pa^-1 - posterior covariance in state x-space normalized by Pa
   !=====================================================================

   !---------------
   ! (2a) ||Pav_ - Pav||
   !---------------
if (.true.) then

   allocate(yhat_Pav_(cv_size,1:rank_rsvd_Pa))

   do itype = 1, 2
      !Generate yhat for approximate Hessian
      yhat_Pav_ = 0.D0

      if (itype .eq. 1) then !LRA
         if (.not.allocated(yhat_Pav_LRA)) cycle

         write(unit=stdout,fmt='(A)') 'LRA'
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
      end if

      if (itype .eq. 2) then !LRU
         if (.not.allocated(yhat_Pav_LRU)) cycle

         write(unit=stdout,fmt='(A)') 'LRU'
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
      end if
      write(unit=stdout,fmt='(A)') 'Perform RSVD on the error residual (Pav_ - Pav)'

      do iens = 1, rank_rsvd_Pa
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             neign,neign,neign_trunc, &
                             LAM,eignvec,qhat, &
                             omega_Pav(:,iens),yhat_Pav_(:,iens))
      end do

      !Subtract (Pav_ * OMEGA - Pav * OMEGA)
      do iens = 1, rank_rsvd_Pa
         if (itype .eq. 1) yhat_Pav_(:,iens) = yhat_Pav_(:,iens) - yhat_Pav_LRA(:,iens)
         if (itype .eq. 2) yhat_Pav_(:,iens) = yhat_Pav_(:,iens) - yhat_Pav_LRU(:,iens)
      end do

      allocate(qhat_Pav(cv_size,0:rank_rsvd_Pa))

      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            yhat_Pav_,qhat_Pav,cv_size,rank_rsvd_Pa)

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega_Pav, yhat_Pav_, qhat_Pav, cv_size, &
                     eignval_res, eignvec_res, rank_rsvd_Pa, dummy)

      deallocate(qhat_Pav)

      !Evaluate Spectral/Frobenius norms from eigenvalues
      SPECNORM_Pav(itype) = maxval(abs(eignval_res(1:dummy)))
      FROBNORM_Pav(itype) = sqrt(sum(eignval_res(1:dummy)*eignval_res(1:dummy)))

      write(unit=stdout,fmt='(A)') '    SPECNORM_Pav   FROBNORM_Pav'
      write(unit=stdout,fmt='(2D15.5)')SPECNORM_Pav(itype),FROBNORM_Pav(itype)

   end do

   deallocate(yhat_Pav_)

   ! Adaptive Pa
   if ( minval(eignval(1:ntrunc)) .lt. 2.D0 ) then
      FROBNORM_Pav(3) = FROBNORM_Pav(2)
      SPECNORM_Pav(3) = SPECNORM_Pav(2)
   else
      FROBNORM_Pav(3) = FROBNORM_Pav(1)
      SPECNORM_Pav(3) = SPECNORM_Pav(1)
   end if
end if

#if (WRF_CHEM == 1)
   !---------------
   ! (2b) ||Pa_ - Pa||
   !---------------
   !Only for chem so far
if (.true.) then
   allocate(yhat_Pa_(cv_size,1:rank_rsvd_Pa))

   do itype = 1, 2
      !Generate yhat for approximate Hessian
      yhat_Pa_ = 0.D0

      if (itype .eq. 1) then !LRA
         if (.not.allocated(yhat_Pa_LRA)) cycle

         write(unit=stdout,fmt='(A)') 'LRA'
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do

      end if

      if (itype .eq. 2) then !LRU
         if (.not.allocated(yhat_Pa_LRU)) cycle

         write(unit=stdout,fmt='(A)') 'LRU'
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
      end if
      write(unit=stdout,fmt='(A)') 'Perform RSVD on the error residual (Pa_ - Pa)'

      !Assumes 1-to-1 size between, e.g., scaleant+scalebb and xhat (SQUARE U)
      do iens = 1, rank_rsvd_Pa
         !Transfer random cv vector to vv and vch
         call da_cv_to_vv(cv_size, omega_Pav(1:be%cv%size_jb, iens), &
                          mz, grid%vv &
#if (WRF_CHEM == 1)
                          , grid%vch  &
#endif
                          )

         ! Transfer vch (and vv later) to xach (and xa later)
#if (WRF_CHEM == 1)
         !$OMP PARALLEL DO &
         !$OMP PRIVATE ( ij, ic, k, j, i)
         do ij = 1 , grid%num_tiles
            do ic = PARAM_FIRST_SCALAR, num_scaleant
               grid%xach%alpha_ant(its:ite,jts:jte,:,ic) = 0.D0
               do k = 1, num_ant_steps    !LEVELS
                  do j = jts, jte
                     do i = its, ite
                         grid%xach%alpha_ant (i,j,k,ic) = grid%vch%alpha_ant (i,j,k,ic)
                     end do
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO &
         !$OMP PRIVATE ( ij, ic, k, j, i)
         do ij = 1 , grid%num_tiles
            do ic = PARAM_FIRST_SCALAR, num_scalebb
               grid%xach%alpha_bb(its:ite,jts:jte,:,ic) = 0.D0
               do k = 1, num_bb_steps    !LEVELS
                  do j = jts, jte
                     do i = its, ite
                         grid%xach%alpha_bb (i,j,k,ic) = grid%vch%alpha_bb (i,j,k,ic)
                     end do
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO
#endif

         ! Calculate U^T * omega (inputs: xa and xach; output: UTomega)
         UTomega = 0.D0
         call da_transform_vtox_adj(grid, be%cv%size_jb, xbx, be, grid%ep, grid%vp, grid%vv, &
                                   UTomega(1:be%cv%size_jb) &
#if (WRF_CHEM == 1)
                                   , grid%vch &
#endif
                                   )

         yhat_Pa_(:,iens) = 0.D0
         !yhat_v = Pav * U^T * omega
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             neign,neign,neign_trunc, &
                             LAM,eignvec,qhat, &
                             UTomega,yhat_Pa_(:,iens))

         !Pa * omega = U * Pav * U^T * omega
         ! outputs: xa, xach
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                yhat_Pa_(1:be%cv%size_jb,iens), grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                                , grid%vch &
#endif
                                )

         yhat_Pa_(:,iens) = 0.D0
         ! Transfer vch (and vv later) from xach (and xa later)
#if (WRF_CHEM == 1)
         !$OMP PARALLEL DO &
         !$OMP PRIVATE ( ij, ic, k, j, i)
         do ij = 1 , grid%num_tiles
            do ic = PARAM_FIRST_SCALAR, num_scaleant
               grid%vch%alpha_ant(its:ite,jts:jte,:,ic) = 0.D0
               do k = 1, num_ant_steps    !LEVELS
                  do j = jts, jte
                     do i = its, ite
                         grid%vch%alpha_ant (i,j,k,ic) = grid%xach%alpha_ant (i,j,k,ic)
                     end do
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO &
         !$OMP PRIVATE ( ij, ic, k, j, i)
         do ij = 1 , grid%num_tiles
            do ic = PARAM_FIRST_SCALAR, num_scalebb
               grid%vch%alpha_bb(its:ite,jts:jte,:,ic) = 0.D0
               do k = 1, num_bb_steps    !LEVELS
                  do j = jts, jte
                     do i = its, ite
                         grid%vch%alpha_bb (i,j,k,ic) = grid%xach%alpha_bb (i,j,k,ic)
                     end do
                  end do
               end do
            end do
         end do
         !$OMP END PARALLEL DO
#endif
         !Transfer vch (and later vv) to yhat_Pa_
         call da_vv_to_cv(grid%vv, grid%xp, mz, cv_size, &
                          yhat_Pa_(1:be%cv%size_jb, iens) &
#if (WRF_CHEM == 1)
                          , grid%vch &
#endif
                           )

        !Subtract (Pa_ * OMEGA - Pa * OMEGA)
         if (itype .eq. 1) yhat_Pa_(:,iens) = yhat_Pa_(:,iens) - yhat_Pa_LRA(:,iens)
         if (itype .eq. 2) yhat_Pa_(:,iens) = yhat_Pa_(:,iens) - yhat_Pa_LRU(:,iens)
      end do

      allocate(qhat_Pa(cv_size,0:rank_rsvd_Pa))

      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            yhat_Pa_,qhat_Pa,cv_size,rank_rsvd_Pa)

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega_Pav, yhat_Pa_, qhat_Pa, cv_size, &
                     eignval_res, eignvec_res, rank_rsvd_Pa, dummy)

      deallocate(qhat_Pa)

      !Evaluate Spectral/Frobenius norms from eigenvalues
      SPECNORM_Pa(itype) = maxval(abs(eignval_res(1:dummy)))
      FROBNORM_Pa(itype) = sqrt(sum(eignval_res(1:dummy)*eignval_res(1:dummy)))

      write(unit=stdout,fmt='(A)') '    SPECNORM_Pa    FROBNORM_Pa'
      write(unit=stdout,fmt='(2D15.5)')SPECNORM_Pa(itype),FROBNORM_Pa(itype)

   end do

   deallocate(yhat_Pa_)

   ! Adaptive Pa
   if ( minval(eignval(1:ntrunc)) .lt. 2.D0 ) then
      FROBNORM_Pa(3) = FROBNORM_Pa(2)
      SPECNORM_Pa(3) = SPECNORM_Pa(2)
   else
      FROBNORM_Pa(3) = FROBNORM_Pa(1)
      SPECNORM_Pa(3) = SPECNORM_Pa(1)
   end if
end if
#endif


   !--------------------------------------------------
   ! (2c) ||Pa_ - Pa||_{Pa^-1} = ||La^-1 (Pa_ - Pa) La^-1||
   !--------------------------------------------------
if (.true. .and. &
    allocated(La_reduced) ) then

   !! Perform decomp in reduced-dimensional space (much cheaper, equal accuracy)
   ! LaQ = La_reduced^T * Q
   allocate(LaQ(rankA,neign_trunc))
   do iens = 1, rankA
      do jens = 1, neign_trunc
         LaQ(iens,jens) = da_dot_cv(cv_size, La_reduced(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
      end do
   end do

   !LaVl = LaQ * E = La * Vl
   allocate(LaVl(rankA,neign_trunc))
   do jens = 1, neign_trunc
      do iens = 1, rankA
         LaVl(iens,jens) = da_dot(neign_trunc, LaQ(iens,1:neign_trunc), eignvec(1:neign_trunc,jens))
      end do
   end do

!!   !This would need to be used if CG produces neign_trunc Ritz vectors of length = neign
!!   ! LaQ = La_reduced^T * Q
!!   allocate(LaQ(rankA,neign))
!!   do iens = 1, rankA
!!      do jens = 1, neign
!!         LaQ(iens,jens) = da_dot_cv(cv_size, La_reduced(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
!!      end do
!!   end do
!!
!!   !LaVl = LaQ * E = La * Vl
!!   allocate(LaVl(rankA,neign_trunc))
!!   do jens = 1, neign_trunc
!!      do iens = 1, rankA
!!         LaVl(iens,jens) = da_dot(neign, LaQ(iens,1:neign), eignvec(1:neign,jens))
!!      end do
!!   end do


   write(unit=stdout,fmt='(A)') 'Perform RSVD on the error residual Pa^-1/2 * (Pa_ - Pa) * Pa^-1/2'
   do itype = 1, 2
      if (itype .eq. 1) then !LRA
         if (.not.hess_comp_PaPa_LRA) cycle
         write(unit=stdout,fmt='(A)') 'LRA'
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
      end if

      if (itype .eq. 2) then !LRU
         if (.not.hess_comp_PaPa_LRU) cycle
         write(unit=stdout,fmt='(A)') 'LRU'
         do i = 1, neign
            LAM(i) = (eignval(i) - 1.D0) / eignval(i)
         end do
      end if


      !LaVlS = LaQ * E * LAM^1/2 = La * Vl * S
      allocate(LaVlS(rankA,neign_trunc))
      do jens = 1, neign_trunc
         LaVlS(:,jens) = LaVl(:,jens) * sqrt(LAM(jens))
      end do

      !LaPaLa_reduced = LaVlS * LaVlS^T
      allocate(LaPaLa_reduced(1:rankA,1:rankA))
      do iens = 1, rankA
         do jens = 1, rankA
            LaPaLa_reduced(iens,jens) = da_dot(neign_trunc, LaVlS(iens,1:neign_trunc), LaVlS(jens,1:neign_trunc))
         end do
      end do
      deallocate(LaVlS)

      !Subtract truth component before spectral decomp
      do iens = 1, rankA
         if (itype .eq. 1) then !LRA
            LaPaLa_reduced(iens,iens) = LaPaLa_reduced(iens,iens) - 1.D0
         end if
         if (itype .eq. 2) then !LRU
            LaPaLa_reduced(iens,iens) = LaPaLa_reduced(iens,iens) - (eignval_hess(iens) - 1.D0)
         end if
      end do

      ! Determine eigenvalues of the symmetric matrix
      !-----------------------------------------------------------------
      info = 0
      call SYEV(LaPaLa_reduced, eignval_res(1:rankA), JOBZ = 'N' , INFO = info) !F95

      !Evaluate Spectral/Frobenius norms from eigenvalues
      SPECNORM_PaPa(itype) = maxval(abs(eignval_res(1:rankA)))
      FROBNORM_PaPa(itype) = sqrt(sum(eignval_res(1:rankA)*eignval_res(1:rankA)))

      write(unit=stdout,fmt='(A)') '    SPECNORM_PaPa   FROBNORM_PaPa'
      write(unit=stdout,fmt='(2D15.5)')SPECNORM_PaPa(itype),FROBNORM_PaPa(itype)

      deallocate(LaPaLa_reduced)
   end do
   deallocate(LaQ)
   deallocate(LaVl)

   ! Adaptive Pa
   if ( minval(eignval(1:ntrunc)) .lt. 2.D0 ) then
      FROBNORM_PaPa(3) = FROBNORM_PaPa(2)
      SPECNORM_PaPa(3) = SPECNORM_PaPa(2)
   else
      FROBNORM_PaPa(3) = FROBNORM_PaPa(1)
      SPECNORM_PaPa(3) = SPECNORM_PaPa(1)
   end if
end if

   ! ||Pa_ - Pa||_{B^-1}
   ! Can not do these in reduced space
   !SPECNORM_PaB(itype) = ??maxval(abs(eignval_res))
   !FROBNORM_PaB(itype) = ??sqrt(sum(eignval_res*eignval_res))

   if (trace_use) call da_trace_exit("da_evaluate_decomp")

end subroutine da_evaluate_decomp
#endif
