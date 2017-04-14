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
   real(wp), allocatable        :: yhat_LaPaLa_LRU(:,:)  ! yhat = La^-1/2 * Pa * La^-1/2 * omega_Pa
   real(wp), allocatable        :: yhat_LaPaLa_LRA(:,:)  ! yhat = La^-1/2 * Pa * La^-1/2 * omega_Pa
   real(wp), allocatable        :: yhat_LaPaLa_half(:,:)  ! yhat = La^-1/2 * omega_Pa

   real*8, allocatable          :: eignvec(:, :)
   real*8, allocatable          :: eignval(:)
   real(wp), allocatable        :: qhat(:,:)
   real(wp)                     :: LAM(1:nens_compare+max_oversamp)
   real(wp)                     :: LAM_hess(1:rankA)
   integer                      :: iens, jens, nens, neign, neign_p, &
                                   i, j, p, dummy, iter, isamp, cv_s, cv_f
   integer                      :: ounit, iostat

   real(wp)                     :: gdot

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

!   rank_rsvd_Pa = 2 * rankA

   rank_rsvd_Pa=floor(scale_hcheck*real(rankA,8))

   call da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA, &
                                   omega_Pav, &
                                   yhat_Pav_LRU, yhat_Pav_LRA, &
                                   yhat_Pa_LRU, yhat_Pa_LRA, &
                                   yhat_LaPaLa_LRU, yhat_LaPaLa_LRA, yhat_LaPaLa_half, &
                                   eignval_hess, eignvec_hess, qhat_hess)

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
                                 yhat_LaPaLa_LRU, yhat_LaPaLa_LRA, yhat_LaPaLa_half, &
                                 eignval_hess, eignvec_hess, qhat_hess(:,0:rankA), &
                                 eignval_hess, eignvec_hess, qhat_hess, &
                                 EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)

         call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                 EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                 EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
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

!      LAM_hess = eignval_hess

!Should be much more accurate for high iter
      LAM_hess = eignval_hess - 1.D0


      do iter=1, nens
         qhat(:,iter) = ghat / beta(iter-1)

         !!!fhat = qhat + AMAT * AMAT^T * qhat (AMAT is sqrt obs HESSIAN, using full eigenval_hess above)

!         fhat = 0.D0
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
                                 yhat_LaPaLa_LRU, yhat_LaPaLa_LRA, yhat_LaPaLa_half, &
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

      deallocate(qhat)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] RSVD5.6:
   !-------------------------------------------------------------------------

   if (rsvd_compare) then

      nens = nens_compare + max_oversamp

      write(unit=stdout,fmt='(A)') 'Evaluating RSVD5.6 dv, dx, Pav, and Pa'

      !Initialize yhat, qhat, omega for RSVD
      allocate(qhat(1:cv_size,0:nens))
      allocate(yhat(1:cv_size,1:nens))

      call da_gen_omega(omega, cv_size, nens, grid, mz)

!      allocate(omega(1:cv_size,1:nens))
!
!!da_set_randomcv sometimes uses the same seed for two different draw, giving identical cv's on a particular processor
!
!!Either the seed needs to be forced to change between multiple draws or it needs to be the same and only a single draw can be taken for all cv values.
!
!!May want to use a similar approach to global below in the first call to da_wrfvar.exe for RIOT.  This may improve LRU result significantly.
!
!!! Cheap enough, but numbers are less random due to changing of seed between nens*nproc draws
!!      do iens = 1, nens
!!         write(unit=stdout,fmt='(A,I5,A,I5)') 'Initialize OMEGA_i, i = ', iens ,' of ',nens
!!         do
!!            call da_set_randomcv(cv_size, omega(:,iens))
!!
!!            !Ensures independent draws
!!            if (iens .eq. 1) then
!!               exit
!!            else
!!               if ( sum(abs(omega(:,iens) - omega(:,iens-1))) .gt. 1.0D-2 ) exit
!!            end if
!!         end do
!!      end do
!
!! Should make this into a subroutine --> inputs: cv_size, nens; outputs: omega
!#ifdef DM_PARALLEL
!      ! Generate global random omega matrix
!      ! Sequence from single global seed (non-repeating)
!      cv_size_domain = wrf_dm_sum_integer(cv_size)
!      if (rootproc) then
!         allocate(omega_g(1:cv_size_domain,nens))
!         allocate(omega1d(1:cv_size_domain*nens))
!         call da_set_randomcv(cv_size_domain*nens, omega1d)
!         cv_s = 1
!         do iens = 1,nens
!            cv_f = cv_s + cv_size_domain-1
!            omega_g(:,iens) = omega1d(cv_s:cv_f)
!            cv_s = cv_f+1
!         end do
!         deallocate(omega1d)
!      else
!         allocate(omega_g(1,1:nens))
!      end if
!
!      call mpi_barrier(comm, ierr)
!
!      !Transfer to local omega
!      allocate(dumcv(1:cv_size))
!      allocate(dumcv_g(1:cv_size_domain))
!      do iens = 1, nens
!         if(rootproc) then
!            dumcv_g = omega_g(:,iens)
!         else
!            dumcv_g = 0.D0
!         end if
!         dumcv = 0.D0
!
!         call da_global_to_cv(cv_size, cv_size_domain, dumcv, grid, mz, dumcv_g)
!         omega(:,iens) = dumcv
!
!         call mpi_barrier(comm, ierr)
!      end do
!      deallocate(dumcv)
!      deallocate(dumcv_g)
!      deallocate(omega_g)
!#else
!      ! Generate local=global omega
!      ! Single seed
!      allocate(omega1d(1:cv_size*nens))
!      call da_set_randomcv(cv_size*nens, omega1d)
!      cv_s = 1
!      do iens = 1,nens
!         cv_f = cv_s + cv_size-1
!         omega(:,iens) = omega1d(cv_s:cv_f)
!         cv_s = cv_f+1
!      end do
!      deallocate(omega1d)
!#endif

      ! Check if each draw is independent or sometimes repeated
      do iens = 1, nens
         write(unit=stdout,fmt='(A,I3,A,5E25.15)') 'omega(1:5,',iens,') = ', omega(1:5,iens)
      end do


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

      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            yhat,qhat,cv_size,nens)

!   !Loop over p (oversample values) (NOT VERY EFFICIENT IF NEED TIME SAVINGS)
!   do isamp = 1, noversamp
!      p = oversamps(isamp)
!
!      ! Perform RSVD for different ranks of the approximate Hessian
!      do neign_p = 1, nens_compare

      do neign_p = 1, nens
!         neign = neign_p + p
!         allocate(eignvec(1:neign,1:neign))
!         allocate(eignval(1:neign))

         allocate(eignvec(1:neign_p,1:neign_p))
         allocate(eignval(1:neign_p))

      
         ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
         call da_rsvd56(grid, mz, jp_start, jp_end, &
                        omega(:,1:neign_p), yhat(:,1:neign_p), qhat(:,0:neign_p), cv_size, &
                        eignval, eignvec, neign_p, dummy)


         do isamp = 1, noversamp

            p = oversamps(isamp)
            neign = neign_p - p
            if  (neign.gt.0 .and. neign.le.nens_compare ) then
 
               call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                       rank_rsvd_Pa, rankA, neign_p, neign, &
                                       xhat_hess, ghat0, &
                                       omega_Pav, &
                                       yhat_Pav_LRU, yhat_Pav_LRA, &
                                       yhat_Pa_LRU, yhat_Pa_LRA, &
                                       yhat_LaPaLa_LRU, yhat_LaPaLa_LRA, yhat_LaPaLa_half, &
                                       eignval, eignvec, qhat(:,0:neign_p), &
                                       eignval_hess, eignvec_hess, qhat_hess, &
                                       EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                       EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                       EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa)


               write(pchar,fmt='(I2.2)') p
               suffix = "rsvd_p="//pchar

               call da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                                       EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                                       EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
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

   if (allocated(yhat_LaPaLa_LRU)) deallocate(yhat_LaPaLa_LRU)
   if (allocated(yhat_LaPaLa_LRA)) deallocate(yhat_LaPaLa_LRA)
   if (allocated(yhat_LaPaLa_half)) deallocate(yhat_LaPaLa_half)

   if (trace_use) call da_trace_exit("da_compare_decomp_methods")

end subroutine da_compare_decomp_methods
#endif
