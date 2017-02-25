#if defined(LAPACK)
subroutine da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA,  &
                                   omega_Pav, yhat_Pav, yhat_Pa, &
                                   eignval_hess, eignvec_hess, qhat_hess)

   !-------------------------------------------------------------------------
   ! Purpose:   Initialize decomposition comparison for full Hessian
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
   integer,        intent(in) :: rank_rsvd_Pa, rankA

   real*8,         intent(in) :: eignvec_hess(rankA, rankA)
   real*8,         intent(in) :: eignval_hess(rankA)
   real*8,         intent(in) :: qhat_hess(cv_size, 0:rankA)

   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)

   real(wp), allocatable, intent(out)        :: yhat_Pav(:,:)  ! yhat = Pav * omega_Pa
   real(wp), allocatable, intent(out)        :: omega_Pav(:,:) ! Gaussian draw
   real(wp), allocatable, intent(out)        :: yhat_Pa(:,:)  ! yhat = Pa * omega_Pa

   real(wp)        :: UTomega(1:cv_size)
   real(wp)        :: LAM_hess(1:rankA)
   integer         :: iens, jens, nens, neign, i, j, k, ij, ic

   if (trace_use) call da_trace_entry("da_evaluate_decomp_init")
   write(unit=stdout,fmt='(A)') 'Initialize Comparison'

   allocate(omega_Pav(1:cv_size,1:rank_rsvd_Pa))
   do iens = 1, rank_rsvd_Pa
      call da_set_randomcv(cv_size, omega_Pav(:,iens))
   end do

   !!!Generate yhat_v and yhat_x for true inverse Hessian (use LRU)
   !Same LAM_hess for both
   do i = 1, rankA
      LAM_hess(i) = - (eignval_hess(i) - 1.D0) / eignval_hess(i)
   end do

   !! Initialize Pav comparison
   allocate(yhat_Pav(1:cv_size,1:rank_rsvd_Pa))
   yhat_Pav = 0.D0
   do iens = 1, rank_rsvd_Pa
!      !Skipping this step, will cancel between yhat_Pav_ and yhat_Pav
!      yhat_Pav(:,iens) = omega_Pav(:,iens)
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             rankA,rankA,rankA, &
                             LAM_hess,eignvec_hess,qhat_hess, &
                             omega_Pav(:,iens),yhat_Pav(:,iens))
   end do


   !! Initialize Pa comparison
   !Assumes 1-to-1 size between, e.g., scaleant+scalebb and xhat (SQUARE U)
   allocate(yhat_Pa(1:cv_size,1:rank_rsvd_Pa))
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
                             rankA,rankA,rankA, &
                             LAM_hess,eignvec_hess,qhat_hess, &
                             UTomega,yhat_Pa(:,iens))

      !Pa * omega = U * Pav * U^T * omega
      ! outputs: xa, xach
      yhat_Pa = 0.D0
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                             yhat_Pa(1:be%cv%size_jb,iens), grid%vv, grid%vp, grid%vch)

      yhat_Pa = 0.D0
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

      !Transfer vch (and later vv) to yhat_Pa
      call da_vv_to_cv(grid%vv, grid%xp, mz, cv_size, &
                       yhat_Pa(1:be%cv%size_jb, iens) &
#if (WRF_CHEM == 1)
                       , grid%vch &
#endif
                        )

   end do

   if (trace_use) call da_trace_exit("da_evaluate_decomp_init")

end subroutine da_evaluate_decomp_init
#endif
