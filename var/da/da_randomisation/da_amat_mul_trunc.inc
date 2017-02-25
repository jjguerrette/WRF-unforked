subroutine da_amat_mul_trunc(grid, mz, jp_start, jp_end, cv_size, &
                             nmodes, neign, ntrunc, &
                             eignval, eignvec, qhat, &
                             shat, xhat)

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
   real,           intent(in) :: qhat(cv_size, 0:nmodes) ! Ritz vectors
   real,           intent(in) :: shat(cv_size)          ! Input vector to multiply by A
   real,           intent(inout):: xhat(cv_size)          ! Output vector: xhat = A.shat
   integer, intent(in)        :: jp_start, jp_end
   integer, intent(in)        :: mz(:)
   real                       :: nuhat(1:cv_size,1:ntrunc)

   integer                    :: i, j
   real                       :: dot_cv

   if (trace_use) call da_trace_entry("da_amat_mul_trunc")

   do j = 1, ntrunc
      do i = 1, cv_size
         nuhat(i,j) = da_dot(neign, qhat(i,1:neign), eignvec(1:neign,j))
      end do
   end do

   !!Initial xhat value defined in calling subroutine
   do i = 1, ntrunc
      dot_cv = da_dot_cv(cv_size, nuhat(:,i), shat, grid, mz, jp_start, jp_end)
      xhat = xhat + eignval(i) * nuhat(:,i) * dot_cv
   end do

   if (trace_use) call da_trace_exit ("da_amat_mul_trunc")

end subroutine da_amat_mul_trunc
