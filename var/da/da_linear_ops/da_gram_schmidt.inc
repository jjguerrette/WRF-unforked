subroutine da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            Q,cv_size,nvec,qstart,R)

   !-------------------------------------------------------------------------
   ! Purpose:        Perform Gram-Schmidt Orthogonalization 
   ! Input: Q(:,1:nvec) = yhat
   ! Output: Q
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)     :: grid
   integer, intent(in)             :: jp_start, jp_end
   integer, intent(in)             :: mz(:)
   integer, intent(in)             :: cv_size, nvec
   real*8, intent(inout)           :: Q(1:cv_size,0:nvec)  ! orthogonal basis
   integer, optional, intent(in)   :: qstart
   real*8, optional, intent(out)   :: R(1:nvec,1:nvec)     ! upper triangular matrix

   integer      :: ivec, jvec, i0, icount
   real*8       :: gdot, beta, tt

   if (trace_use) call da_trace_entry("da_gram_schmidt")

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL start_timing
#endif

      Q(:,0) = 0.D0
      ! Generate Q basis vectors by Modified Gram-Schmidt orthonormalization
      i0 = 1
      if (present(qstart)) then
         if (qstart .gt. 0) i0 = qstart
      end if
      if (present(R)) R = 0.D0

      do ivec = i0, nvec
         tt  = SQRT(da_dot_cv(cv_size, Q(:,ivec), Q(:,ivec), &
                                 grid, mz, jp_start, jp_end))
         icount = 0                                                   !REORTH2

         if (ivec .gt. 1) then
            !!-------------------------
            !! Reorthogonalization Loop
            !!-------------------------

            ! Second iteration improves numerical precision (Bjorck, 1994) 
            ! (depracated)
!!            do icount = 1,2 

            !REORTH2
            ! The do loop is made more robust by adding convergence criteria,
            ! which reduces QR wall time and still improves precision.
            ! [Gander, 1980, Algorithms for the QR-Decomposition, REPORT NO. 80-02]

            do                                                        !REORTH2
               icount = icount + 1                                    !REORTH2
               do jvec = 1, ivec-1
                  gdot = da_dot_cv(cv_size, Q(:,ivec), Q(:,jvec), &
                                   grid, mz, jp_start, jp_end)

!!                  if (icount.eq.1 &
                  if (i0.eq.1 .and. present(R)) then
                     if (icount.eq.1) then
                        R(jvec,ivec) = gdot
                     else                                             !REORTH2
                        R(jvec,ivec) = R(jvec,ivec) + gdot            !REORTH2
                     end if
                  end if

                  Q(:,ivec) = Q(:,ivec) - gdot * Q(:,jvec)
               end do

               beta  = SQRT(da_dot_cv(cv_size, Q(:,ivec), Q(:,ivec), &
                                       grid, mz, jp_start, jp_end))

               !Tolerance for continued reorthogonalization           !REORTH2
               if (beta .ge. tt*0.1D0) then                           !REORTH2
                  exit                                                !REORTH2
               else                                                   !REORTH2
                  tt=beta                                             !REORTH2
               end if                                                 !REORTH2
            end do
         else
            beta = tt
         end if
         if (i0.eq.1 .and. present(R)) R(ivec,ivec) = beta

         Q(:,ivec) = Q(:,ivec) / beta
      end do

#ifdef DM_PARALLEL
   IF ( rootproc ) CALL end_timing('da_end_timing: gram_schmidt')
#endif

   if (trace_use) call da_trace_exit("da_gram_schmidt")

end subroutine da_gram_schmidt


