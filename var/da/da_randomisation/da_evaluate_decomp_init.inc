#if defined(LAPACK)
subroutine da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA,  &
                                   omega_Pav, &
                                   yhat_Pav_LRU, yhat_Pav_LRA, &
                                   yhat_Pa_LRU, yhat_Pa_LRA, &
                                   La_reduced, &
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

   real(wp), allocatable, intent(out)        :: omega_Pav(:,:) ! Gaussian draw

   real(wp), allocatable, intent(out)        :: yhat_Pav_LRU(:,:)  ! yhat = Pav * omega_Pa
   real(wp), allocatable, intent(out)        :: yhat_Pav_LRA(:,:)  ! yhat = Pav * omega_Pa
   real(wp), allocatable, intent(out)        :: yhat_Pa_LRU(:,:)  ! yhat = Pa * omega_Pa
   real(wp), allocatable, intent(out)        :: yhat_Pa_LRA(:,:)  ! yhat = Pa * omega_Pa

   real(wp), allocatable, intent(out)        :: La_reduced(:,:)  !  La^-1 = U^-1 Vm * (Im + LAM_m)^1/2 * Vm^T


   real(wp), allocatable        :: omega_g(:,:) ! Global Gaussian draw
   real(wp), allocatable        :: omega1d(:) ! Gaussian draw
   real(wp), allocatable        :: dumcv(:)
   real(wp), allocatable        :: dumcv_g(:)
   real(wp), allocatable        :: yhat_Pav(:,:)
   real(wp), allocatable        :: yhat_Pa(:,:)

   real(wp)        :: UTomega(1:cv_size)
   real(wp)        :: LAM_hess(1:rankA)
   integer         :: iens, jens, nens, neign, i, j, k, ij, ic, cv_s, cv_f, itype

   if (trace_use) call da_trace_entry("da_evaluate_decomp_init")
   write(unit=stdout,fmt='(A)') 'Initialize Comparison'

   write(unit=stdout,fmt='(A)') 'Initialize OMEGA_Pa'
   call da_gen_omega(omega_Pav, cv_size, rank_rsvd_Pa, grid, mz)

   !-------------------------------------
   !  Initialize Pav comparison
   !-------------------------------------
   do itype = 1, 2
      if (itype .eq. 1) then
         if (.not.hess_comp_Pav_LRA) cycle

         allocate(yhat_Pav_LRA(1:cv_size,1:rank_rsvd_Pa))
         do i = 1, rankA
            LAM_hess(i) = 1.D0 / eignval_hess(i)
         end do
      end if

      if (itype .eq. 2) then
         if (.not.hess_comp_Pav_LRU) cycle

         allocate(yhat_Pav_LRU(1:cv_size,1:rank_rsvd_Pa))
         do i = 1, rankA
            LAM_hess(i) = - (eignval_hess(i) - 1.D0) / eignval_hess(i)
         end do
      end if

      allocate(yhat_Pav(1:cv_size,1:rank_rsvd_Pa))

      yhat_Pav = 0.D0

      do iens = 1, rank_rsvd_Pa
         write(unit=stdout,fmt='(A,I5,A,I5)') 'Initialize Pav*OMEGA_i, i = ', iens ,' of ',rank_rsvd_Pa

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                omega_Pav(:,iens),yhat_Pav(:,iens),.true.)

         if (itype .eq. 1) yhat_Pav_LRA(:,iens) = yhat_Pav(:,iens)
         if (itype .eq. 2) yhat_Pav_LRU(:,iens) = yhat_Pav(:,iens)
      end do
      deallocate( yhat_Pav )
   end do



#if (WRF_CHEM == 1)
   !-------------------------------------
   !  Initialize Pa comparison
   !-------------------------------------
   !Only for chem so far
   !Assumes 1-to-1 size between, e.g., scaleant+scalebb and xhat (SQUARE U)
   do itype = 1, 2
      if (itype .eq. 1) then
         if (.not.hess_comp_Pa_LRA) cycle

         allocate(yhat_Pa_LRA(1:cv_size,1:rank_rsvd_Pa))
         do i = 1, rankA
            LAM_hess(i) = 1.D0 / eignval_hess(i)
         end do
      end if

      if (itype .eq. 2) then
         if (.not.hess_comp_Pa_LRU) cycle

         allocate(yhat_Pa_LRU(1:cv_size,1:rank_rsvd_Pa))
         do i = 1, rankA
            LAM_hess(i) = - (eignval_hess(i) - 1.D0) / eignval_hess(i)
         end do
      end if

      allocate(yhat_Pa(1:cv_size,1:rank_rsvd_Pa))

      yhat_Pa = 0.D0

      do iens = 1, rank_rsvd_Pa
         write(unit=stdout,fmt='(A,I5,A,I5)') 'Initialize Pa*OMEGA_i, i = ', iens ,' of ',rank_rsvd_Pa

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

         yhat_Pa(:,iens) = 0.D0
         !yhat_v = Pav * U^T * omega
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                UTomega,yhat_Pa(:,iens),.true.)

         !Pa * omega = U * Pav * U^T * omega
         ! outputs: xa, xach
         call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                yhat_Pa(1:be%cv%size_jb,iens), grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                                , grid%vch &
#endif
                                )

         yhat_Pa(:,iens) = 0.D0
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

         if (itype .eq. 1) yhat_Pa_LRA(:,iens) = yhat_Pa(:,iens)
         if (itype .eq. 2) yhat_Pa_LRU(:,iens) = yhat_Pa(:,iens)
      end do

      deallocate( yhat_Pa )
   end do
#endif


   !-------------------------------------
   !  Initialize Pa - Pa-norm comparison
   !-------------------------------------
   if (hess_comp_PaPa_LRA .or. hess_comp_PaPa_LRU) then
      allocate(La_reduced(cv_size,1:rankA))

      ! La_reduced = Vm * (Im + LAM_m)^1/2
      do iens = 1, rankA
         do i = 1, cv_size
            La_reduced(i,iens) = da_dot(rankA, qhat_hess(i,1:rankA), eignvec_hess(1:rankA,iens))
         end do
         La_reduced(:,iens) = La_reduced(:,iens) * sqrt(eignval_hess(iens))
      end do
   end if
   if (trace_use) call da_trace_exit("da_evaluate_decomp_init")

end subroutine da_evaluate_decomp_init
#endif
