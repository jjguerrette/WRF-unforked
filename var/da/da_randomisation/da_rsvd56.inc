#if defined(LAPACK)
subroutine da_rsvd56(grid, mz, jp_start, jp_end, &
                     omega, yhat, qhat, cv_size, &
                     eignval, eignvec, nens, neign, print_messages)

   !-------------------------------------------------------------------------
   ! Purpose:         Second stage of RSVD 5.6
   !
   ! Called from da_randomise_svd
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   ! References: Halko et al. (2011)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)       :: grid
   integer, intent(in)     :: jp_start, jp_end
   integer, intent(in)     :: mz(:)

   integer, intent(in)               :: cv_size        ! Total cv size

   integer, intent(in)     :: nens                     ! Number of total eigenpairs
   real(wp), intent(in)    :: yhat(1:cv_size, 1:nens)  ! yhat = A * omega
   real(wp), intent(in)    :: omega(1:cv_size, 1:nens) ! Gaussian draw
   real(wp), intent(in)    :: qhat(1:cv_size, 0:nens)  ! orthonormal basis set 

   real*8, intent(out)     :: eignvec(nens, nens)
   real*8, intent(out)     :: eignval(nens)
   integer, intent(out)    :: neign                    ! Number of non-truncated eigenpairs
   logical, optional, intent(in) :: print_messages

   real(wp), dimension(1:nens) :: LAM
 
   integer                      :: iens, jens, i, j, p

   real(wp)                     :: F1(nens, nens) ! intermediate matrix
   real(wp)                     :: F2(nens, nens) ! intermediate matrix
   real(wp)                     :: U(nens, nens)
   real(wp)                     :: VT(nens, nens)
   real(wp)                     :: SIG(nens)
   real(wp), allocatable        :: TEMP1(:,:), TEMP2(:,:)


   real(wp)                     :: KMAT(nens, nens)  ! KMAT = Q^T * A * Q
   real(wp)                     :: KMAT_(nens, nens)  ! KMAT = Q^T * A * Q
   integer                      :: piv(nens)
   integer                      :: info, cyc

   real                         :: gdot, beta, temp_l
   real                         :: temp_e(nens)
   character(len=3) :: cnens=""
   logical          :: swapped
  
   if (trace_use) call da_trace_entry("da_rsvd56")

   if(.not.present(print_messages)) write(unit=stdout,fmt='(A)') 'Performing Stage 2 of RSVD5.6'
   if(.not.present(print_messages)) write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] ADD SECTION TITLES:
   !-------------------------------------------------------------------------

   do iens = 1, nens
      do jens = 1, nens
         F1(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), omega(:,jens), grid, mz, jp_start, jp_end)
         F2(jens,iens) = da_dot_cv(cv_size, qhat(:,iens), yhat(:,jens), grid, mz, jp_start, jp_end)
      end do
   end do

!REMOVE EVENTUALLY or add option
!   write(cnens,fmt='(I0)') nens
!   do iens = 1, nens
!      write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F1: ',F1(iens,1:nens)
!   end do
!   do iens = 1, nens
!      write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'F2: ',F2(iens,1:nens)
!   end do
!
!REMOVE EVENTUALLY

 ! Solve system D * x = b, &
 !  where KMAT * F1^T = F2^T, F1^T = Q^T * OMEGA, F2^T = Q^T * Y
 !  D = F1
 !  b = F2(:,iens)
 !  x = KMAT(iens,:)
 !-----------------------------------------------------------
   write(cnens,fmt='(I0)') nens

   info = 0
   p = svd_p
   if (p .eq. 0) then
      call GESV(F1, F2, piv, info) !F95

      ! GESV returns KMAT^T in F2
      do iens = 1, nens
         do jens = 1, nens
            KMAT(iens,jens) = F2(jens,iens)
         end do
      end do
   else !NEED TO ADJUST OVERSAMPLING METHOD!!!
      if (p .lt. 0) p = 0
      call GESVD(F1(1:nens,1:nens-p), &
                 SIG(1:nens-p), &
                 U = U(1:nens,1:nens-p), &
                 VT = VT(1:nens-p,1:nens-p), &
                 INFO = info)

      allocate(TEMP1(1:nens-p,1:nens-p))
      allocate(TEMP2(1:nens-p,1:nens-p))
      do iens = 1, nens-p
         TEMP2(iens,:) = (1.D0 / SIG(iens)) * VT(iens,:)
         do jens = 1, nens-p
            TEMP1(jens,iens) = sum(F2(1:nens,iens) * U(1:nens,jens))
         end do
      end do


      do iens = 1, nens-p
         do jens = 1, nens-p
               KMAT(iens,jens) = sum(TEMP1(:,iens) * TEMP2(:,jens))
         end do
      end do
      deallocate(TEMP1)
      deallocate(TEMP2)
   end if
   if (info /=0) then
      write(stdout,*) 'Error in randomised SVD: GESV returned ',info
      CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
   end if
   ! Force symmetry in KMAT
   do iens = 1, nens
      do jens = min(iens+1,nens), nens
         KMAT(iens,jens) = (KMAT(iens,jens) + KMAT(jens,iens)) / 2.D0
         KMAT(jens,iens) = KMAT(iens,jens)
      end do
   end do

!REMOVE EVENTUALLY or add option
!   write(cnens,fmt='(I0)') nens
!   do iens = 1, nens
!      write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'KMAT: ',KMAT(iens,1:nens)
!   end do
!REMOVE EVENTUALLY

 ! Determine eigenvalues and eigenvectors of the symmetric KMAT matrix
 !-----------------------------------------------------------------
   info = 0
   KMAT_ = KMAT
   call SYEV(KMAT, LAM, JOBZ = 'V' , INFO = info) !F95

   if (info /=0) then
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

! This doesn't really work for Hessian with poor symmetry.  Really want to know how many eigenvalues are positive and how many are negative
!   else if ( abs(maxval(LAM)) .lt. abs(minval(LAM)) ) then
!      if(.not.present(print_messages)) write(unit=stdout,fmt='(A)') &
!         'WARNING: Preconditioned Hessian: ABS(MAX(eigs)) < ABS(MIN(eigs)), setting LAM=-LAM'
!      LAM = - LAM

   else if ( sum(sign(temp_e,LAM)) .lt. 0.D0 ) then
      if(.not.present(print_messages)) write(unit=stdout,fmt='(A,F10.1,A)') &
         'WARNING: Preconditioned Hessian: More negative than positive eigenvalues,', sum(sign(temp_e,LAM)), ', setting LAM=-LAM'
      LAM = - LAM
   end if
  
   !First sort eigenvalues and eigenvectors in descending order
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

   neign = nens
!   cyc = 0

   eignval = LAM
   do iens = 1, nens
!!      eignval(iens) = LAM(iens) !!!+ 1.D0
!      if (cyc .eq. 1) cycle
      if (neign .eq. nens .and. LAM(iens) .lt. 0.D0) then
         if (iens .gt. 1) then
            neign = iens - 1
!         else
!!NEED TO ADJUST THIS LANGUAGE
!            write(unit=stdout,fmt='(A)') &
!               'WARNING: There need to be some positive eigenvalues of KMAT'
!            cyc = 1
!
!!!            CALL wrf_error_fatal('Canceling RANDOMIZED SVD')
         end if
      end if
   end do


!REMOVE EVENTUALLY or add option
!   if (cyc .eq. 1) then
!      write(cnens,fmt='(I0)') nens
!      do iens = 1, nens
!         write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'KMAT: ',KMAT_(iens,1:nens)
!      end do
!   end if
!   write(unit=stdout,fmt='(A,'//trim(cnens)//'E27.16E3)') 'Eigenvalues: ',eignval(1:nens)
!REMOVE EVENTUALLY


   if (trace_use) call da_trace_exit("da_rsvd56")

end subroutine da_rsvd56
#endif
