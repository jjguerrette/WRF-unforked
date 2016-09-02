subroutine da_amat_mul_trunc(be, grid, cv_size, nmodes, neign, ntrunc, eignval, eignvec, qhat, shat, xhat)

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

   type (be_type), intent(in) :: be                     ! Background error structure.
   type (domain),  intent(in) :: grid
   integer,        intent(in) :: cv_size
   integer,        intent(in) :: nmodes
   integer,        intent(in) :: neign
   integer,        intent(in) :: ntrunc
   real*8,         intent(in) :: eignvec(nmodes, nmodes)  
   real*8,         intent(in) :: eignval(nmodes)
   real,           intent(in) :: qhat(cv_size, 0:nmodes) ! Ritz vectors
   real,           intent(in) :: shat(cv_size)          ! Input vector to multiply by A
   real,           intent(out):: xhat(cv_size)          ! Output vector: xhat = A.shat

   real                       :: nuhat(1:cv_size,1:ntrunc)

#ifdef CLOUD_CV
   integer                    :: mz(13)                 ! mz for each variable
#else
   integer                    :: mz(7)
#endif
   integer                    :: jp_start, jp_end       ! Start/end indices of Jp.
   integer                    :: i, j
   real                       :: dot_cv

   if (trace_use) call da_trace_entry("da_amat_mul_trunc")

#ifdef CLOUD_CV
   mz       = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%v6%mz, be%v7%mz, be%v8%mz, be%v9%mz,  be%v10%mz, be%v11%mz, be%alpha%mz, be % ne /)
#else
   mz       = (/ be%v1%mz, be%v2%mz, be%v3%mz, be%v4%mz, be%v5%mz, be%alpha%mz, be % ne /)
#endif
   jp_start = be % cv % size_jb + be % cv % size_je + 1
   jp_end   = be % cv % size_jb + be % cv % size_je + be % cv % size_jp


   do j = 1, ntrunc
      do i = 1, cv_size
         nuhat(i,j) = da_dot(neign, qhat(i,1:neign), eignvec(1:neign,j))
      end do
   end do

   xhat     = 0.0
   do i = 1, ntrunc
      dot_cv = da_dot_cv(cv_size, nuhat(:,i), shat, grid, mz, jp_start, jp_end)
      xhat = xhat + eignval(i) * nuhat(:,i) * dot_cv
   end do

   if (trace_use) call da_trace_exit ("da_amat_mul_trunc")

end subroutine da_amat_mul_trunc
