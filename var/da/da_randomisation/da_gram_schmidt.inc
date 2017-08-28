subroutine da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nvec)

   !-------------------------------------------------------------------------
   ! Purpose:        Perform Gram-Schmidt Orthogonalization 
   ! Input: qhat
   ! Output: qhat
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)     :: grid
   integer, intent(in)             :: jp_start, jp_end
   integer, intent(in)             :: mz(:)
   integer, intent(in)             :: cv_size, nvec
   real(wp), intent(inout)         :: qhat(1:cv_size,0:nvec)               ! control variable
 
   integer                         :: ivec, jvec, i
   real(wp)                        :: gdot, beta
 
   if (trace_use) call da_trace_entry("da_gram_schmidt")

      qhat(:,0) = 0.D0
      ! Generate qhat basis vectors by Modified Gram-Schmidt orthonormalization
      do ivec = 1, nvec
         if (ivec .gt. 1) then
            !Second iteration improves numerical precision (Bjorck, 1994) 
            do i = 1,2
            do jvec = 1, ivec-1
               gdot = da_dot_cv(cv_size, qhat(:,ivec), qhat(:,jvec), &
                                grid, mz, jp_start, jp_end)
               qhat(:,ivec) = qhat(:,ivec) - gdot * qhat(:,jvec)
            end do
            end do
         end if

         beta  = SQRT(da_dot_cv(cv_size, qhat(:,ivec), qhat(:,ivec), &
                                 grid, mz, jp_start, jp_end))
         qhat(:,ivec) = qhat(:,ivec) / beta
      end do

end subroutine da_gram_schmidt


