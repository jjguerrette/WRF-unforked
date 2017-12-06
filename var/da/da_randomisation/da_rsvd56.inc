#if defined(LAPACK)
subroutine da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega, yhat, qhat, cv_size, &
                     eignval, eignvec, nens, neign, print_messages)

   !-------------------------------------------------------------------------
   ! Purpose:         Second stage of RSVD 5.6
   !
   ! Called from da_randomise_svd and da_compare_decomp_methods
   !
   ! Solve nearly-symmetric low-dimensional system A * x = b, &
   !  where KMAT * F1^T = F2^T, F1^T = Q^T * OMEGA, F2^T = Q^T * Y
   !  or    F1 * KMAT^T = F1 * KMAT = F2
   !  D = F1
   !  b = F2(:,iens)
   !  x = KMAT(iens,:)
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   ! References: Halko et al. (2011)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(in) :: grid
   integer, intent(in)      :: jp_start, jp_end
   integer, intent(in)      :: mz(:)

   integer, intent(in)      :: cv_size        ! Total cv size

   integer, intent(in)     :: nens                     ! Number of total eigenpairs
   real, intent(in)        :: yhat(1:cv_size, 1:nens)  ! yhat = A * omega
   real, intent(in)        :: omega(1:cv_size, 1:nens) ! Gaussian draw
   real, intent(in)        :: qhat(1:cv_size, 0:nens)  ! orthonormal basis set 

   real, intent(out)       :: eignvec(nens, nens)
   real, intent(out)       :: eignval(nens)

   integer, intent(out)    :: neign                    ! Number of non-truncated eigenpairs
   logical, optional, intent(in) :: print_messages

   real(wp)                     :: LAM(nens)
   integer                      :: iens, jens, i, j, p1
   real(wp)                     :: F1(nens, nens) ! intermediate matrix
   real(wp)                     :: F2(nens, nens) ! intermediate matrix
   real(wp)                     :: U(nens, nens)
   real(wp)                     :: VT(nens, nens)
   real(wp)                     :: SIG(nens)
   real(wp)                     :: UTF2(nens,nens), VSIGinvT(nens,nens)
   real(wp)                     :: KMAT(nens, nens)  ! KMAT = Q^T * A * Q

   integer                      :: piv(nens)
   integer                      :: info

   real(wp)                     :: temp_l
   real(wp)                     :: temp_e(nens)
   character(len=3) :: cnens=""
   logical          :: swapped
  
   if (trace_use) call da_trace_entry("da_rsvd56")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

   if(.not.present(print_messages)) write(unit=stdout,fmt='(A,I5,A)') 'Performing Stage 2 of RSVD5.6 for ', nens ,' ensembles'
   if(.not.present(print_messages)) write(unit=stdout,fmt=*) ' '

   do iens = 1, nens
      do jens = 1, nens
         !F1 = OMEGA^T * Q
         F1(iens,jens) = da_dot_cv(cv_size, omega(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)

         !F2 = Y^T * Q
         F2(iens,jens) = da_dot_cv(cv_size, yhat(:,iens), qhat(:,jens), grid, mz, jp_start, jp_end)
      end do
   end do

   write(cnens,fmt='(I0)') nens

   info = 0

   !! Solve F1*KMAT = F2

   ! DECOMPOSE F1 = U * SIG * VT 
   call GESVD(F1(1:nens,1:nens), &
              SIG(1:nens), &
              U = U(1:nens,1:nens), &
              VT = VT(1:nens,1:nens), &
              INFO = info)

   call mpi_barrier(comm, ierr)

   if (info .ne. 0) then
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F1: ',F1(iens,1:nens)
      end do

      write(stdout,*) 'Error in randomised SVD: GESV returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if

   ! KMAT = VT^T * SIG^-1 * U^T * F2 
   do iens = 1, nens
      VSIGinvT(iens,1:nens) = VT(iens,1:nens) / SIG(iens)
      do jens = 1, nens
         UTF2(iens,jens) = da_dot(nens,U(1:nens,iens),F2(1:nens,jens))
      end do
   end do

   do iens = 1, nens
      do jens = 1, nens
         KMAT(jens,iens) = da_dot(nens,VSIGinvT(:,iens),UTF2(:,jens))
      end do
   end do


   ! Force symmetry in KMAT
   KMAT = (KMAT + TRANSPOSE(KMAT)) / 2.D0

 ! Determine eigenvalues and eigenvectors of the symmetric KMAT matrix
 !-----------------------------------------------------------------
   info = 0
   call SYEV(KMAT, LAM, JOBZ = 'V' , INFO = info) !F95

   call mpi_barrier(comm, ierr)

   if (info .ne. 0) then
      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F2: ',F2(iens,1:nens)
      end do

      do iens = 1, nens
         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'KMAT: ',KMAT(iens,1:nens)
      end do

      write(stdout,*) 'Error in randomised SVD: SYEV returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if

   do iens = 1, nens
      do jens = 1, nens 
         eignvec(iens,jens) = KMAT(iens,jens)
      end do
   end do

   temp_e = 1.D0
   if (all(LAM .le. 0.D0)) then
      if(.not.present(print_messages)) write(unit=stdout,fmt='(A)') &
         'WARNING: Preconditioned Hessian: all(eigs <= 0), setting them positive.'
      LAM = - LAM
   else if ( sum(sign(temp_e,LAM)) .lt. 0.D0 ) then
      if(.not.present(print_messages)) write(unit=stdout,fmt='(A,F10.1,A)') &
         'WARNING: Preconditioned Hessian: More negative than positive eigenvalues,', sum(sign(temp_e,LAM)), ', setting LAM=-LAM'
      LAM = - LAM
   end if
  
   !Sort eigenvalues and eigenvectors in descending order
   do j = nens-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (LAM(i) .lt. LAM(i+1)) then
         temp_l = LAM(i)
         LAM(i) = LAM(i+1)
         LAM(i+1) = temp_l

         temp_e = eignvec(:,i)
         eignvec(:,i) = eignvec(:,i+1)
         eignvec(:,i+1) = temp_e

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do

   !Preconditioned observation Hessian should be positive semi-definite
   neign = nens
   do iens = 1, nens
      LAM(iens) = max(LAM(iens),0.D0)
   end do
   eignval = LAM

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: rsvd56')
#endif

   if (trace_use) call da_trace_exit("da_rsvd56")

end subroutine da_rsvd56
#endif
