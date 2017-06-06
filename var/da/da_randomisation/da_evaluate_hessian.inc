#if defined(LAPACK)
subroutine da_evaluate_hessian(grid, xbx, be, mz, cv_size, &
                              rankA, neign, neign_trunc, &
                              eignval_all, eignvec_all, qhat)

   !-------------------------------------------------------------------------
   ! Purpose:   Evaluate inverse Hessian convergence (CHEM only)
   !
   ! Called from da_main
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (xbx_type),intent(inout) :: xbx     ! Header & non-gridded vars.
   type (be_type), intent(in)    :: be      ! background error structure.
   integer,        intent(in) :: cv_size

   integer, intent(in)        :: mz(:)

   integer, intent(in) :: rankA, neign, neign_trunc
   real*8,  intent(in) :: eignvec_all(neign, neign)
   real*8,  intent(in) :: eignval_all(neign)
   real*8,  intent(in) :: qhat(cv_size, 0:neign)

   !Check Variables
   real(wp)        :: UTomega(cv_size)
   real(wp), allocatable  :: omega_Pa(:,:)
   real(wp), allocatable  :: yhat_Pa(:,:)
   real(wp), allocatable  :: yhat_Pa_prev(:,:)
   real(wp), allocatable  :: yhat_Pa_curr(:,:)
   real(wp), allocatable  :: qhat_Pa(:,:)


   real*8          :: eignvec_(neign,neign)
   real*8          :: eignval_(neign)
   real*8          :: temp_l, temp_e(neign)
   logical         :: swapped
   real*8          :: temp
   real*8          :: eignvec_res(rankA,rankA)
   real*8          :: eignval_res(rankA)

   real*8          :: LAM(neign)

   integer         :: iens, i, j, k, ij, ic, itype, dummy, iter
   real*8          :: FROBNORM_Pa, SPECNORM_Pa

   integer         :: jp_start, jp_end

   if (trace_use) call da_trace_entry("da_evaluate_hessian")

   jp_start   = be % cv % size_jb + be % cv % size_je + 1
   jp_end     = be % cv % size_jb + be % cv % size_je + be % cv % size_jp

   FROBNORM_Pa=0.D0; SPECNORM_Pa=0.D0

   !=====================================================================
   ! Compare Pa_approx to Pa_prev
   !    Pa  - posterior covariance in state x-space
   !=====================================================================

   !---------------
   ! ||Pa_ - Pa||
   !---------------
   write(unit=stdout,fmt='(A)') 'Initialize OMEGA_Pa'

   call da_gen_omega(omega_Pa, cv_size, rankA, grid, mz)

   allocate(yhat_Pa(cv_size,1:rankA))
   allocate(yhat_Pa_curr(cv_size,1:rankA))
   allocate(yhat_Pa_prev(cv_size,1:rankA))

   !First sort eigenvalues and eigenvectors in descending order
   eignval_ = eignval_all
   eignvec_ = eignvec_all
   do j = neign-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (eignval_(i) .lt. eignval_(i+1)) then
         temp_l = eignval_(i)
         eignval_(i) = eignval_(i+1)
         eignval_(i+1) = temp_l

         temp_e = eignvec_(:,i)
         eignvec_(:,i) = eignvec_(:,i+1)
         eignvec_(:,i+1) = temp_e

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do

   write(unit=stdout,fmt='(A)') 'Perform RSVD on the error residual (Pa_j - Pa_j-1)'
!   do itype = 1, 2
   do itype = 2, 2
      if (itype .eq. 1) then !LRA
         write(unit=stdout,fmt='(A)') 'LRA'
         do i = 1, neign
            LAM(i) = 1.D0 / eignval_(i)
         end do

      end if

      if (itype .eq. 2) then !LRU
         write(unit=stdout,fmt='(A)') 'LRU'
         do i = 1, neign
            LAM(i) = - (eignval_(i) - 1.D0) / eignval_(i)
         end do
      end if

      write(unit=stdout,fmt='(A)') 'HESS:  iter    SPECNORM_Pa    FROBNORM_Pa'

      do iter = 1, neign_trunc
         if (iter .gt. 1) then
            do iens = 1, rankA
               yhat_Pa_prev(:,iens) = yhat_Pa_curr(:,iens)
            end do
         else
            yhat_Pa_prev = 0.D0
         end if

         !Generate yhat for current rank (iter)
         yhat_Pa_curr = 0.D0

         !Assumes 1-to-1 size between, e.g., scaleant+scalebb and xhat (SQUARE U)
         do iens = 1, rankA
            !Transfer random cv vector to vv and vch
            call da_cv_to_vv(cv_size, omega_Pa(1:be%cv%size_jb, iens), &
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

            yhat_Pa_curr(:,iens) = 0.D0
            !yhat_v = Pav * U^T * omega
            call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                   neign,neign,iter, &
                                   LAM,eignvec_,qhat, &
                                   UTomega,yhat_Pa_curr(:,iens))

            !Pa * omega = U * Pav * U^T * omega
            ! outputs: xa, xach
            call da_transform_vtox(grid, be%cv%size_jb, xbx, be, grid%ep, &
                                   yhat_Pa_curr(1:be%cv%size_jb,iens), grid%vv, grid%vp &
#if (WRF_CHEM == 1)
                                   , grid%vch &
#endif
                                   )

            yhat_Pa_curr(:,iens) = 0.D0
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

            !Transfer vch (and later vv) to yhat_Pa_curr
            call da_vv_to_cv(grid%vv, grid%xp, mz, cv_size, &
                             yhat_Pa_curr(1:be%cv%size_jb, iens) &
#if (WRF_CHEM == 1)
                             , grid%vch &
#endif
                              )

         end do

         do iens = 1, rankA
           !Subtract (Pa_i * OMEGA - Pa_i-1 * OMEGA)
           yhat_Pa(:,iens) = yhat_Pa_curr(:,iens) - yhat_Pa_prev(:,iens)
         end do

         allocate(qhat_Pa(cv_size,0:rankA))

         call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                               yhat_Pa,qhat_Pa,cv_size,rankA)

         call da_rsvd56(grid, mz, jp_start, jp_end, &
                        omega_Pa, yhat_Pa, qhat_Pa, cv_size, &
                        eignval_res, eignvec_res, rankA, dummy, .false. )

         deallocate(qhat_Pa)

         !Evaluate Spectral/Frobenius norms from eigenvalues
         SPECNORM_Pa = maxval(abs(eignval_res(1:dummy)))
         FROBNORM_Pa = sqrt(sum(eignval_res(1:dummy)*eignval_res(1:dummy)))

         write(unit=stdout,fmt='(A,I5,2D15.5)') 'HESS: ', iter, SPECNORM_Pa, FROBNORM_Pa
      end do
   end do

   deallocate(omega_Pa)
   deallocate(yhat_Pa)
   deallocate(yhat_Pa_prev)
   deallocate(yhat_Pa_curr)

   if (trace_use) call da_trace_exit("da_evaluate_hessian")

end subroutine da_evaluate_hessian
#endif
