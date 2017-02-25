#if defined(LAPACK)
subroutine da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                              rank_rsvd_Pa, rankA, neign, neign_trunc, &
                              xhat_hess, ghat0, &
                              omega_Pav, yhat_Pav, yhat_Pa, &
                              eignval, eignvec, qhat, &
                              eignval_hess, eignvec_hess, qhat_hess, &
                              EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)

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
   integer,        intent(in) :: rank_rsvd_Pa, rankA, neign, neign_trunc

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

   real*8, intent(out)        :: EUCNORM_x(2), FROBNORM_Pa(2), SPECNORM_Pa(2)
   real*8, intent(out)        :: EUCNORM_v(2), FROBNORM_Pav(2), SPECNORM_Pav(2)
   real*8, intent(out)        :: EUCNORM_xPa(2), FROBNORM_PaPa(2), SPECNORM_PaPa(2)
!!   real*8, intent(out)        :: EUCNORM_xB(2), FROBNORM_PaB(2), SPECNORM_PaB(2)

   !Check Variables
   real(wp), intent(in)        :: yhat_Pav(cv_size,rank_rsvd_Pa)  ! yhat = Pav * omega_Pa
   real(wp), intent(in)        :: omega_Pav(cv_size,rank_rsvd_Pa) ! Gaussian draw
   real(wp), intent(in)        :: yhat_Pa(cv_size,rank_rsvd_Pa)  ! yhat = Pa * omega_Pa

   real(wp)        :: UTomega(cv_size)
   real*8          :: xhat(cv_size)
   real*8          :: diffvec(cv_size)

   real(wp), allocatable  :: qhat_Pav(:,:)
   real(wp), allocatable  :: yhat_Pav_(:,:)
   real(wp), allocatable  :: qhat_Pa(:,:)
   real(wp), allocatable  :: yhat_Pa_(:,:)

   real*8          :: temp
   real*8          :: eignvec_res(rank_rsvd_Pa, rank_rsvd_Pa)
   real*8          :: eignval_res(rank_rsvd_Pa), LAM(neign)

   integer         :: iens, jens, i, j, k, ij, ic, itype, dummy

   if (trace_use) call da_trace_entry("da_evaluate_decomp")

   write(unit=stdout,fmt=*) ' '
   write(unit=stdout,fmt='(A,I6)') 'Evaluate Norms, k = ',neign


   EUCNORM_x=0.D0; FROBNORM_Pa=0.D0; SPECNORM_Pa=0.D0

   EUCNORM_v=0.D0; FROBNORM_Pav=0.D0; SPECNORM_Pav=0.D0

   EUCNORM_xPa=0.D0; FROBNORM_PaPa=0.D0; SPECNORM_PaPa=0.D0

   !=====================================================================
   ! Compare xhat to xhat_hess and Uxhat to Uxhat_hess (vtox)
   !=====================================================================
   write(unit=stdout,fmt='(A)') 'Evaluating (xa_ - xa)'
   do itype = 1, 2

      ! ||dva_ - dva||
      if (itype .eq. 1) then !LRA
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
         xhat = 0.D0
      end if

      if (itype .eq. 2) then !LRU
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
         xhat = ghat0
      end if

      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             neign,neign,neign_trunc, &
                             LAM,eignvec,qhat,&
                             ghat0,xhat)


      diffvec = xhat - xhat_hess
      EUCNORM_v(itype) = da_dot_cv(cv_size, diffvec, diffvec, grid, mz, jp_start, jp_end)

      ! ||dxa_ - dxa||
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, &
                             grid%ep, diffvec(1:be%cv%size_jb), &
                             grid%vv, grid%vp, grid%vch)

      temp = 0.D0
      !Only for chem so far
#if (WRF_CHEM == 1)
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
      !Should I really be comparing scaleant/scalebb or log(scaleant) and log(scalebb) across
      ! the entire grid?  Or should I only compare where emissions are non-zero? 
      ! Or where emisions are largest?  Above a threshold? Can I use Fischer information?
      ! Or simply posterior variance reduction?  I suppose that's where ||xa_-xa||_Pa^-1
      ! comes in. Although there can be variance reduction even in places with zero emissions,
      ! which is probably arbitrary.


      ! ||dxa_ - dxa||_{Pa^-1}
      xhat = 0.D0
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             rankA,rankA,rankA, &
                             eignval_hess,eignvec_hess,qhat_hess,&
                             diffvec,xhat)

      EUCNORM_xPa(itype) = da_dot_cv(cv_size, diffvec, xhat, grid, mz, jp_start, jp_end)

      !WHAT ABOUT dx-DIRECTION for DGN/LM nonlinear problems?
   end do

   !=====================================================================
   ! Compare Pa_approx to Pa_hess
   !    Pav - posterior covariance in preconditioned v-space
   !    Pa  - posterior covariance in state x-space
   !=====================================================================
   write(unit=stdout,fmt='(A)') 'Perform RSVD on the error residual (Pa_ - Pa)'

   allocate(qhat_Pav(cv_size,0:rank_rsvd_Pa))
   allocate(yhat_Pav_(cv_size,1:rank_rsvd_Pa))

   ! ||Pav_ - Pav||
   do itype = 1, 2
      !Generate yhat for approximate Hessian
      yhat_Pav_ = 0.D0

      if (itype .eq. 1) then !LRA
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
      end if

      if (itype .eq. 2) then !LRU
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
!         !Skipping this step, will cancel between yhat_Pav_ and yhat_Pav 
!         yhat_Pav_(:,iens) = omega_Pav(:,iens)
      end if

      do iens = 1, rank_rsvd_Pa
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                neign,neign,neign_trunc, &
                                LAM,eignvec,qhat, &
                                omega_Pav(:,iens),yhat_Pav_(:,iens))
      end do

      !Subtract (Pav_ * OMEGA - Pav * OMEGA)
      do iens = 1, rank_rsvd_Pa
         yhat_Pav_(:,iens) = yhat_Pav_(:,iens) - yhat_Pav(:,iens)
      end do

      qhat_Pav(:,0) = 0.D0
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            yhat_Pav_,qhat_Pav,cv_size,rank_rsvd_Pa)

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega_Pav, yhat_Pav_, qhat_Pav, cv_size, &
                     eignval_res, eignvec_res, rank_rsvd_Pa, dummy)

      !Evaluate Spectral/Frobenius norms from eigenvalues
      SPECNORM_Pav(itype) = maxval(eignval_res)
      FROBNORM_Pav(itype) = sqrt(sum(eignval_res*eignval_res))
   end do

   deallocate(qhat_Pav)
   deallocate(yhat_Pav_)

   allocate(qhat_Pa(cv_size,0:rank_rsvd_Pa))
   allocate(yhat_Pa_(cv_size,1:rank_rsvd_Pa))

   ! ||Pa_ - Pa||
   do itype = 1, 2
      !Generate yhat for approximate Hessian
      yhat_Pa_ = 0.D0

      if (itype .eq. 1) then !LRA
         do i = 1, neign
            LAM(i) = 1.D0 / eignval(i)
         end do
      end if

      if (itype .eq. 2) then !LRU
         do i = 1, neign
            LAM(i) = - (eignval(i) - 1.D0) / eignval(i)
         end do
!         !Skipping this step, will cancel between yhat_Pa_ and yhat_Pa
!         yhat_Pa_(:,iens) = omega_Pa(:,iens)
      end if


      !Assumes 1-to-1 size between, e.g., scaleant+scalebb and xhat (SQUARE U)
      do iens = 1, rank_rsvd_Pa
         !Transfer random cv vector to vv and vch
         call da_cv_to_vv(cv_size, omega_Pav(1:be%cv%size_jb, iens), &
                          mz, grid%vv, &
#if (WRF_CHEM == 1)
                          grid%vch  &
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
                                   UTomega(1:be%cv%size_jb), &
#if (WRF_CHEM == 1)
                                   grid%vch &
#endif
                                   )

         !yhat_v = Pav * U^T * omega
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                neign,neign,neign_trunc, &
                                LAM,eignvec,qhat, &
                                UTomega,yhat_Pa_(:,iens))

         !Pa * omega = U * Pav * U^T * omega
         ! outputs: xa, xach
         yhat_Pa_ = 0.D0
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                yhat_Pa_(1:be%cv%size_jb,iens), grid%vv, grid%vp, grid%vch)

         yhat_Pa_ = 0.D0
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
         yhat_Pa_(:, iens) = yhat_Pa_(:, iens) &
                              - yhat_Pa(:, iens)
      end do

      qhat_Pa(:,0) = 0.D0
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            yhat_Pa_,qhat_Pa,cv_size,rank_rsvd_Pa)

      call da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega_Pav, yhat_Pa_, qhat_Pa, cv_size, &
                     eignval_res, eignvec_res, rank_rsvd_Pa, dummy)


      !Evaluate Spectral/Frobenius norms from eigenvalues
      SPECNORM_Pa(itype) = maxval(eignval_res)
      FROBNORM_Pa(itype) = sqrt(sum(eignval_res*eignval_res))
   end do

   deallocate(qhat_Pa)
   deallocate(yhat_Pa_)


   ! ||Pa_ - Pa||_{Pa^-1}
   ! Should be able to do these in reduced space
   !SPECNORM_PaPa(itype) = ??maxval(eignval_res)
   !FROBNORM_PaPa(itype) = ??sqrt(sum(eignval_res*eignval_res))

   ! ||Pa_ - Pa||_{B^-1}
   ! Can not do these in reduced space
   !SPECNORM_PaB(itype) = ??maxval(eignval_res)
   !FROBNORM_PaB(itype) = ??sqrt(sum(eignval_res*eignval_res))

   if (trace_use) call da_trace_exit("da_evaluate_decomp")

end subroutine da_evaluate_decomp
#endif
