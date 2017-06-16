subroutine da_amat_mul_trunc(grid, mz, jp_start, jp_end, cv_size, &
                             nmodes, neign, ntrunc, &
                             eignval, eignvec, qhat, &
                             shat, xhat, nuflag)

   !-------------------------------------------------------------------------
   ! Purpose:  Multiply a control vector by the Analysis Error Cov Matrix A 
   !
   ! Called from da_solve
   !
   ! History: 09/01/2016  Creation (JJ Guerrette)
   !
   ! Modified from da_amat_mul
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(in) :: grid
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: nmodes
   integer,        intent(in) :: neign
   integer,        intent(in) :: ntrunc
   real*8,         intent(in) :: eignvec(nmodes, nmodes)  
   real*8,         intent(in) :: eignval(nmodes)
   real*8,         intent(in) :: qhat(cv_size, 0:nmodes) ! Ritz vectors
   real*8,         intent(in) :: shat(cv_size)          ! Input vector to multiply by A
   real*8,      intent(inout) :: xhat(cv_size)          ! Output vector: xhat = A.shat
   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)
   logical, intent(in), optional :: nuflag

   real*8                     :: nuhat(1:cv_size,1:ntrunc)
!   real*8                     :: QTs(1:neign)
!   real*8                     :: ETQTs(1:neign_trunc)
   real*8                     :: gam(ntrunc), gam_(ntrunc), temp_l, temp_e(cv_size)
   integer                    :: i, j, start, fin, step, order(ntrunc), temp_i
   real*8                     :: dot_cv
   logical                    :: swapped

   if (trace_use) call da_trace_entry("da_amat_mul_trunc")

   if (present(nuflag)) then
      nuhat = qhat(1:cv_size,1:ntrunc)
   else
      ! Is this slow due to non-parallel da_dot across cv_size?
      do j = 1, ntrunc
         do i = 1, cv_size
            nuhat(i,j) = da_dot(neign, qhat(i,1:neign), eignvec(1:neign,j))
         end do
      end do

!      do i = 1, neign
!         QTs(i) = da_dot_cv(cv_size, qhat(:,i), shat, grid, mz, jp_start, jp_end)
!      end do
!      do i = 1, ntrunc
!         ETQTs(i) = da_dot(neign,QTs(1:neign),eignvec(1:neign,i)
!         gam(i) = eignval(i) * ETQTs(i)
!         order(i) = i
!      end do
   end if

   !Sort gamma values in descending order to improve mat-vec precision
   do i = 1, ntrunc 
      dot_cv = da_dot_cv(cv_size, nuhat(:,i), shat, grid, mz, jp_start, jp_end)
      gam(i) = eignval(i) * dot_cv
      order(i) = i
   end do
   gam_ = gam

! This ordering significantly improves numerical precision
   do j = ntrunc-1, 1, -1
     swapped = .false.
     do i = 1, j
       if (abs(gam_(i)) .lt. abs(gam_(i+1))) then
         temp_l = gam_(i)
         gam_(i) = gam_(i+1)
         gam_(i+1) = temp_l

         temp_i = order(i)
         order(i) = order(i+1)
         order(i+1) = temp_i

!! Many more reassignments
!         temp_e = nuhat(:,i)
!         nuhat(:,i) = nuhat(:,i+1)
!         nuhat(:,i+1) = temp_e

         swapped = .true.
       end if
     end do
     if (.NOT. swapped) exit
   end do

   !!Initial xhat value defined in calling subroutine
   do i = ntrunc, 1, -1 !From smallest to largest gamma
!      write(unit=stdout,fmt='(A,I0,A,E25.15)') 'gamma_', i, ' = ', gam(i)
!      write(unit=stdout,fmt='(A,I0,A,2E25.15)') 'gamma_', i, ' = ', gam(i), gam(i)*minval(nuhat(:,i))
!      write(unit=stdout,fmt='(A,I0,A,2E25.15,I4)') 'gamma_', i, ' = ', gam(i), gam(i)*minval(nuhat(:,i)), order(i)

      xhat = xhat + gam(order(i)) * nuhat(:,order(i))
   end do

!   !!Initial xhat value defined in calling subroutine
!   if (eignval(1) .lt. eignval(ntrunc)) then
!      start = 1
!      fin = ntrunc
!      step = 1
!   else
!      start = ntrunc
!      fin = 1
!      step = -1
!   end if
!
!!   do i = start, fin, step
!!!   do i = 1, ntrunc !What about when lam_k < lam_k+1
!   do i = ntrunc, 1, -1 !This ordering eliminates some precision issues when lam_k>lam_k+1
!
!      dot_cv = da_dot_cv(cv_size, nuhat(:,i), shat, grid, mz, jp_start, jp_end)
!      xhat = xhat + eignval(i) * nuhat(:,i) * dot_cv
!
!   end do

   if (trace_use) call da_trace_exit ("da_amat_mul_trunc")

end subroutine da_amat_mul_trunc
