#if defined(LAPACK)
subroutine da_compare_decomp_methods(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                      rankA, eignval_hess, eignvec_hess, qhat_hess, &
                                      ghat0,xhat_hess)

   !-------------------------------------------------------------------------
   ! Purpose:         Compare RSVD and Lanczos increments and Pa to exact Hessian
   !
   ! Called from da_calculate_hessian in final stage
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   ! References: Halko et al. (2011); Bousserez and Henze (2016); Golub and van Loan(1996)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (xbx_type),intent(inout) :: xbx     ! Header & non-gridded vars.
   type (be_type), intent(in)    :: be      ! background error structure.
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: rankA
   real*8,         intent(in) :: eignvec_hess(rankA, rankA)
   real*8,         intent(in) :: eignval_hess(rankA)
   real,           intent(in) :: qhat_hess(cv_size, 0:rankA) ! Ritz vectors
   real,           intent(in) :: ghat0(cv_size)               ! Input vector to multiply by A
   real,           intent(in) :: xhat_hess(cv_size)           ! Input vector: xhat from FULL HESS
   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)


   ! RSVD eval vars
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
   integer, parameter           :: max_oversamp=15
   integer, parameter           :: noversamp=4
   integer, parameter           :: oversamps(noversamp) = (/0, 5, 10, max_oversamp/)
   integer                      :: rank_rsvd_Pa

   ! Lanczos eval vars
   real(wp)                     :: ghat(1:cv_size)
   real(wp)                     :: fhat(1:cv_size)
   real*8                       :: alpha(nens_compare), beta(0:nens_compare)
   real*8                       :: subdiag(nens_compare)
   real*8                       :: sstwrk(2*nens_compare-2)
   integer                      :: nsstwrk, info

   !Shared eval vars
   real(wp), allocatable        :: yhat_Pav(:,:)  ! yhat = Pav * omega
   real(wp), allocatable        :: omega_Pav(:,:) ! Gaussian draw
   real(wp), allocatable        :: yhat_Pa(:,:)  ! yhat = Pa * omega

   real*8, allocatable          :: eignvec(:, :)
   real*8, allocatable          :: eignval(:)
   real(wp)                     :: qhat(1:cv_size,0:nens_compare+max_oversamp)
   real(wp)                     :: LAM(1:nens_compare+max_oversamp)
   real(wp)                     :: LAM_hess(1:rankA)
   integer                      :: iens, jens, nens, neign, neign_p, i, j, p, dummy, iter, isamp
   integer                      :: ounit, iostat

   real                         :: gdot

   real*8        :: EUCNORM_x(2), FROBNORM_Pa(2), SPECNORM_Pa(2)
   real*8        :: EUCNORM_v(2), FROBNORM_Pav(2), SPECNORM_Pav(2)
   real*8        :: EUCNORM_xPa(2), FROBNORM_PaPa(2), SPECNORM_PaPa(2)

   character(len=filename_len)  :: filename
   character*3  :: stat
   character*2  :: pchar
   character*10 :: suffix
 
   if (trace_use) call da_trace_entry("da_compare_decomp_methods")

   write(unit=stdout,fmt='(A)') 'Compare Hessian decompositions'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------

   rank_rsvd_Pa=2*rankA

   call da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA, &
                                   omega_Pav, yhat_Pav, yhat_Pa, &
                                   eignval_hess, eignvec_hess, qhat_hess)

   !-------------------------------------------------------------------------
   ! [2.0] TSVD:
   !-------------------------------------------------------------------------

   nens = nens_compare

   ! Compare to TSVD for different ranks of the approximate Hessian
   suffix = "tsvd"
   do neign = 1, nens_compare
      call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                              rank_rsvd_Pa, rankA, rankA, neign, &
                              xhat_hess, ghat0, &
                              omega_Pav, yhat_Pav, yhat_Pa, &
                              eignval_hess, eignvec_hess, qhat_hess(:,0:neign), &
                              eignval_hess, eignvec_hess, qhat_hess, &
                              EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)

      call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                              neign, suffix)
   end do



   !-------------------------------------------------------------------------
   ! [3.0] RSVD5.6:
   !-------------------------------------------------------------------------

   nens = nens_compare + max_oversamp

   !Initialize yhat, qhat, omega for RSVD

   allocate(yhat(1:cv_size,1:nens))
   allocate(omega(1:cv_size,1:nens))
   do iens = 1, nens
      call da_set_randomcv(cv_size, omega(:,iens))
   end do

   !Multiply obs Hessian by omega (Y=A_obs * omega)
   LAM_hess = eignval_hess - 1.D0

   yhat = 0.D0
   do iens = 1, nens
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             rankA,rankA,rankA, &
                             LAM_hess,eignvec_hess,qhat_hess, &
                             omega(:,iens),yhat(:,iens))

   end do

   qhat(:,0) = 0.D0

   call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         yhat,qhat,cv_size,nens)

   ! Perform RSVD for different ranks of the approximate Hessian
   do neign_p = 1, nens_compare

      !Loop over p (oversample values) (NOT VERY EFFICIENT IF NEED TIME SAVINGS)
      do isamp = 1, noversamp
         p = oversamps(isamp)
   
         neign = neign_p + p
   
         allocate(eignvec(1:neign,1:neign))
         allocate(eignval(1:neign))
   
         ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
         call da_rsvd56(grid, mz, jp_start, jp_end, &
                        omega(:,1:neign), yhat(:,1:neign), qhat(:,0:neign), cv_size, &
                        eignval, eignvec, neign, dummy)
   
         call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                 rank_rsvd_Pa, rankA, neign, min(neign_p, dummy), &
                                 xhat_hess, ghat0, &
                                 omega_Pav, yhat_Pav, yhat_Pa, &
                                 eignval, eignvec, qhat(:,0:neign), &
                                 eignval_hess, eignvec_hess, qhat_hess, &
                                 EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)

         deallocate(eignval)
         deallocate(eignvec)

         write(pchar,fmt='(I2.2)') isamp
         suffix = "rsvd_p="//pchar

         call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                 neign, suffix)


      end do
   end do

   deallocate(yhat)
   deallocate(omega)


   !-------------------------------------------------------------------------
   ! [4.0] CG-Lanczos:
   !-------------------------------------------------------------------------

   nens = nens_compare

   ghat = ghat0
   beta(0)   = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
   qhat(:,:) = 0.D0

   LAM_hess = eignval_hess

   do iter=1, nens
      qhat(:,iter) = ghat / beta(iter-1)

      !!!fhat = qhat + AMAT * AMAT^T * qhat (AMAT is sqrt obs HESSIAN)

      fhat = 0.D0
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                             rankA,rankA,rankA, &
                             LAM_hess,eignvec_hess,qhat_hess, &
                             qhat(:,iter),fhat)



      alpha(iter) = da_dot_cv(cv_size, qhat(:,iter), fhat, grid, mz, jp_start, jp_end)

      ghat        = fhat - alpha(iter)*qhat(:,iter) - beta(iter-1)*qhat(:,iter-1)

      do j = 1,2 !Second iteration improves orthogonalization for numerical precision (Bjorck, 1994)
         do i = 1, iter
            gdot = da_dot_cv(cv_size, ghat, qhat(:,i), grid, mz, jp_start, jp_end)
            ghat = ghat - gdot * qhat(:,i)
         end do
      end do

      beta(iter)  = SQRT(da_dot_cv (cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
   end do

   suffix = "cg-lanczos"
   do neign=1, nens
    ! Determine eigenvalues and eigenvectors of the Lanczos tri-diagonal matrix
    !--------------------------------------------------------------------------
      allocate(eignvec(1:neign,1:neign))
      allocate(eignval(1:neign))

      eignval(1:neign)   = alpha(1:neign)
      subdiag(1:neign-1) = beta(1:neign-1)
      nsstwrk           = MAX(2*neign-2,1)
      info              = 0
      call DSTEQR('I',neign,eignval(1:neign),subdiag(1:neign-1),eignvec,neign,&
                   sstwrk(1:nsstwrk),info)
      if (info /=0) write(stdout,*) 'Error in Lanczos minimization: SSTEQR returned ',info

      !Both LRU and LRA for CG? Or just LRA?
      call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                              rank_rsvd_Pa, rankA, neign, neign, &
                              xhat_hess, ghat0, &
                              omega_Pav, yhat_Pav, yhat_Pa, &
                              eignval, eignvec, qhat(:,0:neign), &
                              eignval_hess, eignvec_hess, qhat_hess, &
                              EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)
      deallocate(eignval)
      deallocate(eignvec)

      call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                              neign, suffix)
   end do

   if (allocated(yhat_Pav)) deallocate(yhat_Pav)
   if (allocated(omega_Pav)) deallocate(omega_Pav)
   if (allocated(yhat_Pa)) deallocate(yhat_Pa)

   if (trace_use) call da_trace_exit("da_compare_decomp_methods")

end subroutine da_compare_decomp_methods
#endif
