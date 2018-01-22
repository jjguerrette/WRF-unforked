subroutine da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat,cv_size,nvec,qstart)

   !-------------------------------------------------------------------------
   ! Purpose:        Perform Gram-Schmidt Orthogonalization 
   ! Input: qhat(:,1:nvec) = yhat
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
   real*8, intent(inout)           :: qhat(1:cv_size,0:nvec)  ! orthogonal basis
   integer, optional, intent(in)   :: qstart

   integer                         :: ivec, jvec, i, i0
   real*8                          :: gdot, beta
 
   if (trace_use) call da_trace_entry("da_gram_schmidt")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      qhat(:,0) = 0.D0
      ! Generate qhat basis vectors by Modified Gram-Schmidt orthonormalization
      i0 = 1
      if (present(qstart)) then
         if (qstart .gt. 0) i0 = qstart
      end if

      do ivec = i0, nvec
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

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: gram_schmidt')
#endif

   if (trace_use) call da_trace_exit("da_gram_schmidt")

end subroutine da_gram_schmidt


