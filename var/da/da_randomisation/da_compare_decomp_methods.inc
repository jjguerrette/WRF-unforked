#if defined(LAPACK)
subroutine da_compare_decomp_methods(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                      rankA, eignval_hess, eignvec_hess, qhat_hess, &
                                      ghat0,xhat_hess,AMAT)

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
   real*8,         intent(in) :: qhat_hess(cv_size, 0:rankA) ! Ritz vectors
   real*8,         intent(in) :: ghat0(cv_size)               ! Input vector to multiply by A
   real*8,         intent(in) :: xhat_hess(cv_size)           ! Input vector: xhat from FULL HESS
   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)
   real*8,         intent(in) :: AMAT(cv_size, 1:rankA) ! Half Hessian


   ! RSVD eval vars
   real(wp), allocatable        :: yhat(:,:)  ! yhat = A * omega
   real(wp), allocatable        :: omega(:,:) ! Gaussian draw
   real(wp), allocatable        :: omega1d(:) ! Gaussian draw
   real(wp), allocatable        :: omega_g(:,:)
   real(wp), allocatable        :: dumcv(:)
   real(wp), allocatable        :: dumcv_g(:)

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
   real(wp), allocatable        :: omega_Pav(:,:) ! Gaussian draw
   real(wp), allocatable        :: yhat_Pav_LRU(:,:)  ! yhat = Pav * omega
   real(wp), allocatable        :: yhat_Pav_LRA(:,:)  ! yhat = Pav * omega
   real(wp), allocatable        :: yhat_Pa_LRU(:,:)  ! yhat = Pa * omega
   real(wp), allocatable        :: yhat_Pa_LRA(:,:)  ! yhat = Pa * omega
   real(wp), allocatable        :: La_reduced(:,:)  !  La^-1 = U^-1 Vm * (Im + LAM_m)^1/2 * Vm^T
   real(wp), allocatable        :: dvnorm_hess(:), dxnorm_hess(:)


   real*8, allocatable          :: eignvec(:, :)
   real*8, allocatable          :: eignval(:)
   real(wp), allocatable        :: qhat(:,:)
   real(wp)                     :: LAM_hess(1:rankA)
   integer                      :: iens, jens, nens, neign, ntrunc, &
                                   i, j, p, dummy, iter, isamp, cv_s, cv_f
   integer                      :: ounit, iostat

   real(wp)                     :: gdot

   real*8        :: EUCNORM_x(3), FROBNORM_Pa(3), SPECNORM_Pa(3)
   real*8        :: EUCNORM_v(3), FROBNORM_Pav(3), SPECNORM_Pav(3)
   real*8        :: EUCNORM_xPa(3), FROBNORM_PaPa(3), SPECNORM_PaPa(3)
   real*8        :: DIRNORM1_v(3), DIRNORM1_x(3)
   real*8        :: DIRNORM2_v(3), DIRNORM2_x(3)



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

!   rank_rsvd_Pa = 2 * rankA

   rank_rsvd_Pa=floor(scale_hcheck*real(rankA,8))

   write(unit=stdout,fmt='(A)') 'Initialize OMEGA_Pa'
   call da_gen_omega(omega_Pav, cv_size, rank_rsvd_Pa, grid, mz)

   call da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA, &
                                   eignval_hess, eignvec_hess, qhat_hess, &
                                   xhat_hess, &
                                   omega_Pav, &
                                   yhat_Pav_LRU, yhat_Pav_LRA, &
                                   yhat_Pa_LRU, yhat_Pa_LRA, &
                                   La_reduced, &
                                   dvnorm_hess, dxnorm_hess)

   !-------------------------------------------------------------------------
   ! [2.0] TSVD:
   !-------------------------------------------------------------------------
   if (tsvd_compare) then
      nens = nens_compare

      ! Compare to TSVD for different ranks of the approximate Hessian
      suffix = "tsvd"
      write(unit=stdout,fmt='(3A)') 'Evaluating ',suffix ,' dv, dx, Pav, and Pa'

      do neign = 1, nens
         call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                 rank_rsvd_Pa, rankA, rankA, neign, &
                                 xhat_hess, ghat0, &
                                 omega_Pav, &
                                 yhat_Pav_LRU, yhat_Pav_LRA, &
                                 yhat_Pa_LRU, yhat_Pa_LRA, &
                                 La_reduced, &
                                 eignval_hess, eignvec_hess, qhat_hess(:,0:rankA), &
                                 eignval_hess, eignvec_hess, qhat_hess, &
                                 dvnorm_hess, dxnorm_hess, &
                                 EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                 DIRNORM1_v, DIRNORM1_x, &
                                 DIRNORM2_v, DIRNORM2_x )


         call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                 DIRNORM1_v, DIRNORM1_x , &
                                 DIRNORM2_v, DIRNORM2_x , &
                                 neign, suffix)
      end do
   end if

   !-------------------------------------------------------------------------
   ! [3.0] CG-Lanczos:
   !-------------------------------------------------------------------------

   if (cglz_compare) then
      write(unit=stdout,fmt='(A)') 'Evaluating CG-Lanczos dv, dx, Pav, and Pa'

      nens = nens_compare

      allocate(qhat(1:cv_size,0:nens))

      ghat = ghat0
      beta(0)   = SQRT(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))
      qhat(:,:) = 0.D0

      LAM_hess = eignval_hess - 1.D0

      do iter=1, nens
         qhat(:,iter) = ghat / beta(iter-1)

         !fhat = qhat + AMAT * AMAT^T * qhat (AMAT is sqrt obs HESSIAN)
         fhat = qhat(:,iter)

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                qhat(:,iter),fhat,.true.)


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
                                 omega_Pav, &
                                 yhat_Pav_LRU, yhat_Pav_LRA, &
                                 yhat_Pa_LRU, yhat_Pa_LRA, &
                                 La_reduced, &
                                 eignval, eignvec, qhat(:,0:neign), &
                                 eignval_hess, eignvec_hess, qhat_hess, &
                                 dvnorm_hess, dxnorm_hess, &
                                 EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                 DIRNORM1_v, DIRNORM1_x, &
                                 DIRNORM2_v, DIRNORM2_x )

         deallocate(eignval)
         deallocate(eignvec)

         call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                 DIRNORM1_v, DIRNORM1_x , &
                                 DIRNORM2_v, DIRNORM2_x , &
                                 neign, suffix)
      end do

      deallocate(qhat)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] RSVD5.6:
   !-------------------------------------------------------------------------

   if (rsvd_compare) then

      nens = nens_compare + max_oversamp
!      nens = nens_compare

      write(unit=stdout,fmt='(A)') 'Evaluating RSVD5.6 dv, dx, Pav, and Pa'

      !Initialize yhat, qhat, omega for RSVD
      allocate(qhat(1:cv_size,0:nens))
      allocate(yhat(1:cv_size,1:nens))

      allocate(omega(cv_size,nens))
      call da_gen_omega(omega, cv_size, nens, grid, mz)

      !Multiply obs Hessian by omega (Y=A_obs * omega)
      LAM_hess = eignval_hess - 1.D0

      yhat = 0.D0
      do iens = 1, nens
         write(unit=stdout,fmt='(A,I5,A,I5)') 'Initialize Pav*OMEGA_i, i = ', iens ,' of ',nens
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                omega(:,iens),yhat(:,iens),.true.)

      end do

      qhat(:,1:nens) = yhat(:,1:nens)
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nens)

      do neign = 1, nens
         allocate(eignvec(1:neign,1:neign))
         allocate(eignval(1:neign))

      
         ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
         call da_rsvd56(grid, mz, jp_start, jp_end, &
                        omega(:,1:neign), yhat(:,1:neign), qhat(:,0:neign), cv_size, &
                        eignval, eignvec, neign, dummy)

         eignval = eignval + 1.D0

         do isamp = 1, noversamp

            p = oversamps(isamp)
            ntrunc = neign - p
            if  (ntrunc.gt.0 .and. ntrunc.le.nens_compare ) then
 
               call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                       rank_rsvd_Pa, rankA, neign, ntrunc, &
                                       xhat_hess, ghat0, &
                                       omega_Pav, &
                                       yhat_Pav_LRU, yhat_Pav_LRA, &
                                       yhat_Pa_LRU, yhat_Pa_LRA, &
                                       La_reduced, &
                                       eignval, eignvec, qhat(:,0:neign), &
                                       eignval_hess, eignvec_hess, qhat_hess, &
                                       dvnorm_hess, dxnorm_hess, &
                                       EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                       EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                       EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                       DIRNORM1_v, DIRNORM1_x, &
                                       DIRNORM2_v, DIRNORM2_x )

               write(pchar,fmt='(I2.2)') p
               suffix = "rsvd_p="//pchar

               call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                       EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                       EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                                       DIRNORM1_v, DIRNORM1_x , &
                                       DIRNORM2_v, DIRNORM2_x , &
                                       neign, suffix)
            end if

         end do

         deallocate(eignval)
         deallocate(eignvec)

      end do

      deallocate(qhat)
      deallocate(yhat)
      deallocate(omega)

   end if

   if (allocated(omega_Pav)) deallocate(omega_Pav)

   if (allocated(yhat_Pav_LRU)) deallocate(yhat_Pav_LRU)
   if (allocated(yhat_Pav_LRA)) deallocate(yhat_Pav_LRA)

   if (allocated(yhat_Pa_LRU)) deallocate(yhat_Pa_LRU)
   if (allocated(yhat_Pa_LRA)) deallocate(yhat_Pa_LRA)

   if (allocated(La_reduced)) deallocate(La_reduced)

   if (allocated(dvnorm_hess)) deallocate(dvnorm_hess)
   if (allocated(dxnorm_hess)) deallocate(dxnorm_hess)

   if (trace_use) call da_trace_exit("da_compare_decomp_methods")

end subroutine da_compare_decomp_methods
#endif
