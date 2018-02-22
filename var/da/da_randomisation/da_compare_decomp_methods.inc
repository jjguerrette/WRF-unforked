#if defined(LAPACK)
subroutine da_compare_decomp_methods(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                      rankA, eignval_hess, eignvec_hess, qhat_hess, &
                                      ghat0, xhat_hess, AMAT)

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
   real(wp), allocatable        :: dumcv(:)
   real(wp), allocatable        :: dumcv_g(:)

   integer, parameter           :: max_oversamp=15
   integer, parameter           :: noversamp=4
   integer, parameter           :: oversamps(noversamp) = (/0, 5, 10, max_oversamp/)
   integer                      :: rank_rsvd_Pa

   ! Lanczos eval vars
   real(wp)                     :: ghat(1:cv_size)
   real(wp)                     :: ghat_CG(1:cv_size)
   real(wp)                     :: grad(1:cv_size)
   real(wp)                     :: cv_dummy(1:cv_size)
   real(wp)                     :: xhat_REF(1:cv_size,2)
   real*8                       :: alpha(nsamp_compare), beta(0:nsamp_compare)
   real*8                       :: subdiag(nsamp_compare)
   real*8                       :: sstwrk(2*nsamp_compare-2)
   integer                      :: nsstwrk, info

   !Shared eval vars
   real(wp), allocatable        :: omega_Pav(:,:) ! Gaussian draw
   real(wp), allocatable        :: yhat_Pav_LRU(:,:)  ! yhat = Pav * omega
   real(wp), allocatable        :: yhat_Pav_LRA(:,:)  ! yhat = Pav * omega
   real(wp), allocatable        :: yhat_Pa_LRU(:,:)  ! yhat = Pa * omega
   real(wp), allocatable        :: yhat_Pa_LRA(:,:)  ! yhat = Pa * omega
   real(wp), allocatable        :: La_reduced(:,:)  !  La^-1 = U^-1 Vm * (Im + LAM_m)^1/2 * Vm^T
   real(wp), allocatable        :: dvnorm_hess(:), dxnorm_hess(:)


   real*8, allocatable          :: eignvec(:,:)
   real*8, allocatable          :: eignval(:)
   real(wp), allocatable        :: qhat(:,:)
   real(wp)                     :: LAM_hess(1:rankA)
   real*8, allocatable          :: Vk_exact(:,:), Vk_approx(:,:)
   real*8                       :: EIGDOT(rankA)

   integer                      :: iens, jens, nsamp, neign, ntrunc, &
                                   i, j, p, dummy, iter, isamp, cv_s, cv_f
   integer                      :: itype, l
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
   character*5  :: nchar

   character*10 :: suffix

   real, allocatable                 :: temp_e(:)
   real                              :: temp_l
   logical                           :: swapped

   if (trace_use) call da_trace_entry("da_compare_decomp_methods")

   write(unit=stdout,fmt='(A)') 'Compare Hessian decompositions'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------

!   rank_rsvd_Pa = 2 * rankA

   rank_rsvd_Pa=floor(scale_hcheck*real(rankA,8))

   write(unit=stdout,fmt='(A)') 'Initialize OMEGA_Pa'
   if ( hess_comp_Pa_LRA  .OR. hess_comp_Pa_LRU &
   .OR. hess_comp_Pav_LRA .OR. hess_comp_Pav_LRU ) then
      allocate(omega_Pav(cv_size,rank_rsvd_Pa))
      call da_gen_omega(omega_Pav, cv_size, rank_rsvd_Pa, grid, mz)
   end if

   call da_evaluate_decomp_init(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                   rank_rsvd_Pa, rankA, &
                                   eignval_hess, eignvec_hess, qhat_hess, &
                                   xhat_hess, &
                                   omega_Pav, &
                                   yhat_Pav_LRU, yhat_Pav_LRA, &
                                   yhat_Pa_LRU, yhat_Pa_LRA, &
                                   La_reduced, &
                                   dvnorm_hess, dxnorm_hess)

   ghat = ghat0

   if (check_precon) then
      neign = rankA

      l = max(min(neign-nsamp_compare-1,spectral_trunc(1)),0)
      spectral_trunc(1) = l

      write(nchar,fmt='(I0)')rankA

      allocate(qhat(cv_size,0:neign))
      allocate(eignvec(1:neign,1:neign))
      allocate(eignval(1:neign))

      eignvec = eignvec_hess
      eignval = eignval_hess
      qhat = qhat_hess

      call da_hessian_io ('w',grid, cv_size,1,neign,neign,&
                           eignvec,eignval, &
                           qhat,mz,ghat)

!      eignvec = 0.D0
!      eignvec(:,1:neign-l) = eignvec_hess(:,1+l:neign)
!
!      eignval = 0.D0
!      eignval(1:neign-l) = eignval_hess(1+l:neign)
!
!      qhat = 0.D0
!      qhat(:,1:neign-l) = qhat_hess(:,1+l:neign)
!
!      !Limit reference increment to spectra of Hessian outside preconditioned zone
!      do itype = 1, 2
!         if (itype .eq. 1) then !LRA
!            do i = 1, neign-l
!               LAM_hess(i) = -(eignval(i) - 1.D0) / eignval(i)
!            end do
!         end if
!
!         if (itype .eq. 2) then !LRU
!            do i = 1, neign-l
!               LAM_hess(i) = 1.D0 / eignval(i)
!            end do
!         end if
!
!         cv_dummy = 0.D0
!         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
!                                neign,neign,neign-l,&
!                                LAM_hess,eignvec,qhat,&
!                                ghat,cv_dummy,.true.)
!
!         if (itype .eq. 2) then !LRU
!            cv_dummy = cv_dummy + ghat
!         end if
!
!         xhat_REF(:,itype) = cv_dummy
!
!!         xhat_REF(:,itype) = xhat_hess - cv_dummy
!!
!!         xhat_REF(:,itype) = xhat_hess
!      end do

      deallocate(qhat)
      deallocate(eignvec)
      deallocate(eignval)

      xhat_REF(:,1) = xhat_hess
      xhat_REF(:,2) = xhat_hess

      !Precondition gradient
      cv_dummy = ghat
      call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                           1, cv_dummy, ghat, 'AD' )

      ! Vk_exact = Q * E
      allocate(Vk_exact(1:cv_size,neign))
      Vk_exact = 0.D0
      do iens = 1, neign
         Vk_exact(:,iens) = qhat_hess(:,iens)
      end do

   else

      xhat_REF(:,1) = xhat_hess
      xhat_REF(:,2) = xhat_hess

      !-------------------------------------------------------------------------
      ! [2.0] TSVD:
      !-------------------------------------------------------------------------
      if (tsvd_compare) then
         nsamp = nsamp_compare

         ! Compare to TSVD for different ranks of the approximate Hessian
         suffix = "tsvd"
         write(unit=stdout,fmt='(3A)') 'Evaluating ',suffix ,' dv, dx, Pav, and Pa'

         do neign = 1, nsamp
            call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                    rank_rsvd_Pa, rankA, rankA, neign, &
                                    xhat_REF, ghat, &
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
                                    neign, 0, suffix)
         end do
      end if
   end if

   !-------------------------------------------------------------------------
   ! [3.0] CG-Lanczos:
   !-------------------------------------------------------------------------

   if (cglz_compare) then
      write(unit=stdout,fmt='(A)') 'Evaluating CG-Lanczos dv, dx, Pav, and Pa'

      nsamp = nsamp_compare

      allocate(qhat(1:cv_size,0:nsamp))

      ghat_CG = ghat
      beta(0)   = SQRT(da_dot_cv(cv_size, ghat_CG, ghat_CG, grid, mz, jp_start, jp_end))
      qhat(:,:) = 0.D0

      LAM_hess = eignval_hess - 1.D0

      do iter=1, nsamp
         qhat(:,iter) = ghat_CG / beta(iter-1)

         !grad = PMAT^T*PMAT * qhat + PMAT^T * AMAT * AMAT^T * PMAT *qhat 
         ! - AMAT is sqrt obs HESSIAN
         ! - PMAT is spectral preconditioner

         if (check_precon) then
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 1, qhat(:,iter), cv_dummy, 'TL' )
         else
            cv_dummy = qhat(:,iter)
         end if

         grad = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                cv_dummy,grad,.true.)

         cv_dummy = cv_dummy + grad

         if (check_precon) then
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 1, cv_dummy, grad, 'AD' )
         else
            grad = cv_dummy
         end if

         alpha(iter) = da_dot_cv(cv_size, qhat(:,iter), grad, grid, mz, jp_start, jp_end)

         ghat_CG        = grad - alpha(iter)*qhat(:,iter) - beta(iter-1)*qhat(:,iter-1)

         do j = 1,2 !Second iteration improves orthogonalization for numerical precision (Bjorck, 1994)
            do i = 1, iter
               gdot = da_dot_cv(cv_size, ghat_CG, qhat(:,i), grid, mz, jp_start, jp_end)
               ghat_CG = ghat_CG - gdot * qhat(:,i)
            end do
         end do

         beta(iter)  = SQRT(da_dot_cv (cv_size, ghat_CG, ghat_CG, grid, mz, jp_start, jp_end))
      end do

      suffix = "cg-lanczos"
      do neign=1, nsamp
       ! Determine eigenvalues and eigenvectors of the Lanczos tri-diagonal matrix
       !--------------------------------------------------------------------------
         allocate(eignvec(1:neign,1:neign))
         allocate(eignval(1:neign))

         eignval(1:neign)   = alpha(1:neign)
         subdiag(1:neign-1) = beta(1:neign-1)
         nsstwrk           = MAX(2*neign-2,1)
         info              = 0
         call DSTEQR('I',neign,eignval,subdiag(1:neign-1),eignvec,neign,&
                      sstwrk(1:nsstwrk),info)
         if (info /=0) write(stdout,*) 'Error in Lanczos minimization: SSTEQR returned ',info

         !Sort eigenvalues and eigenvectors in descending order
         allocate(temp_e(1:neign))
         do j = neign-1, 1, -1
           swapped = .false.
           do i = 1, j
             if (eignval(i) .lt. eignval(i+1)) then
               temp_l = eignval(i)
               eignval(i) = eignval(i+1)
               eignval(i+1) = temp_l

               temp_e = eignvec(1:neign,i)
               eignvec(1:neign,i) = eignvec(1:neign,i+1)
               eignvec(1:neign,i+1) = temp_e

               swapped = .true.
             end if
           end do
           if (.NOT. swapped) exit
         end do
         deallocate(temp_e)

         if (check_precon .and. mod(neign,5).eq.0) then
            if(rootproc) then
               ! =========
               ! open file
               ! =========
               call da_get_unit(ounit)

               write(filename,fmt='("EIGDOT_CG_n=",I0)') neign
               open(unit   = ounit,     &
                  FILE   = trim(filename), &
                  FORM   = 'FORMATTED',  &
                  ACCESS = 'SEQUENTIAL', &
                  iostat =  iostat,     &
                  STATUS = 'REPLACE')
            end if

            write(unit=stdout,fmt='("Check Eigenvector Alignment of Preconditioned Hessian (CG):")') 
            allocate(Vk_approx(1:cv_size,1))
            do iens = 1, neign

               ! Vk_approx = PMAT * Q * E
               Vk_approx = 0.D0
               do jens = 1, neign
                  Vk_approx(:,1) = Vk_approx(:,1) + qhat(:,jens) * eignvec(jens,iens)
               end do
               call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                    1, Vk_approx(:,1), Vk_approx(:,1), 'IA' )

               do jens = 1, rankA
                  EIGDOT(jens) = da_dot_cv(cv_size, Vk_approx(:,1), Vk_exact(:,jens), grid, mz, jp_start, jp_end)
               end do

               if (rootproc) write(unit=ounit,fmt='(I8,'//trim(nchar)//'E16.6)') iens, EIGDOT
            end do
            deallocate(Vk_approx)

            if(rootproc) then
               close(ounit)
               call da_free_unit(ounit)
            end if
         end if

         !Both LRU and LRA for CG? Or just LRA?
         call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                 rank_rsvd_Pa, rankA, neign, neign, &
                                 xhat_REF, ghat, &
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
                                 neign, 0, suffix)
      end do

      deallocate(qhat)

   end if

   !-------------------------------------------------------------------------
   ! [4.0] RSVD5.6:
   !-------------------------------------------------------------------------

   if (rsvd_compare) then

      nsamp = nsamp_compare + max_oversamp
!      nsamp = nsamp_compare

      write(unit=stdout,fmt='(A)') 'Evaluating RSVD5.6 dv, dx, Pav, and Pa'

      !Initialize yhat, qhat, omega for RSVD
      allocate(qhat(1:cv_size,0:nsamp))
      allocate(yhat(1:cv_size,1:nsamp))

      allocate(omega(cv_size,nsamp))
      call da_gen_omega(omega, cv_size, nsamp, grid, mz)

      !Multiply obs Hessian by omega (Y=A_obs * omega)
      LAM_hess = eignval_hess - 1.D0

      yhat = 0.D0
      do iens = 1, nsamp
         write(unit=stdout,fmt='(A,I5,A,I5)') 'Initialize Pav*OMEGA_i, i = ', iens ,' of ',nsamp

         !grad = PMAT^T*PMAT * qhat + PMAT^T * AMAT*AMAT^T * PMAT * qhat - I * qhat
         ! - AMAT is sqrt obs HESSIAN
         ! - PMAT is spectral preconditioner

         if (check_precon) then
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 1, omega(:,iens), cv_dummy, 'TL' )
         else
            cv_dummy = omega(:,iens)
         end if

         grad = 0.D0

         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size, &
                                rankA,rankA,rankA, &
                                LAM_hess,eignvec_hess,qhat_hess, &
                                cv_dummy,grad,.true.)

         if (check_precon) then
            cv_dummy = cv_dummy + grad
            call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                 1, cv_dummy, grad, 'AD' )

            grad = grad - omega(:,iens)
         end if

         yhat(:,iens) = grad
      end do

      qhat(:,1:nsamp) = yhat(:,1:nsamp)
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nsamp)

      do neign = 1, nsamp
         allocate(eignvec(1:neign,1:neign))
         allocate(eignval(1:neign))

         ! Eigenvalue decomposition of A (Algorithm 5.6 of Halko)
         call da_rsvd56(grid, mz, jp_start, jp_end, &
                        omega(:,1:neign), yhat(:,1:neign), qhat(:,0:neign), cv_size, &
                        eignval, eignvec, neign)

         eignval = eignval + 1.D0

         if (check_precon .and. mod(neign,5).eq.0) then
            if(rootproc) then
               ! =========
               ! open file
               ! =========
               call da_get_unit(ounit)

               write(filename,fmt='("EIGDOT_RSVD_n=",I0)') neign
               open(unit   = ounit,     &
                  FILE   = trim(filename), &
                  FORM   = 'FORMATTED',  &
                  ACCESS = 'SEQUENTIAL', &
                  iostat =  iostat,     &
                  STATUS = 'REPLACE')
            end if

            write(unit=stdout,fmt='("Check Eigenvector Alignment of Preconditioned Hessian (RSVD):")') 
            allocate(Vk_approx(1:cv_size,1))
            do iens = 1, neign

               ! Vk_approx = PMAT * Q * E
               Vk_approx = 0.D0
               do jens = 1, neign
                  Vk_approx(:,1) = Vk_approx(:,1) + qhat(:,jens) * eignvec(jens,iens)
               end do
               call da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, &
                                    1, Vk_approx(:,1), Vk_approx(:,1), 'IA' )

               do jens = 1, rankA
                  EIGDOT(jens) = da_dot_cv(cv_size, Vk_approx(:,1), Vk_exact(:,jens), grid, mz, jp_start, jp_end)
               end do

               if (rootproc) write(unit=ounit,fmt='(I8,'//trim(nchar)//'E16.6)') iens, EIGDOT
            end do
            deallocate(Vk_approx)

            if(rootproc) then
               close(ounit)
               call da_free_unit(ounit)
            end if
         end if

         do isamp = 1, noversamp

            p = oversamps(isamp)
            ntrunc = neign - p
            if  (ntrunc.gt.0 .and. ntrunc.le.nsamp_compare ) then
 
               call da_evaluate_decomp(grid, xbx, be, mz, jp_start, jp_end, cv_size, &
                                       rank_rsvd_Pa, rankA, neign, ntrunc, &
                                       xhat_REF, ghat, &
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
                                       neign, p, suffix)
            end if

         end do

         deallocate(eignval)
         deallocate(eignvec)

      end do

      deallocate(qhat)
      deallocate(yhat)
      deallocate(omega)

   end if

   if (allocated(Vk_exact)) deallocate(Vk_exact)


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
