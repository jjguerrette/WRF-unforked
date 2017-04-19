#if defined(LAPACK)
subroutine da_evaluate_increment(grid, xbx, be, mz, cv_size, &
                                 xhat_prev, xhat_curr)

   !-------------------------------------------------------------------------
   ! Purpose:   Evaluate increment convergence
   !
   ! Called from da_main
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   type (domain),  intent(inout) :: grid
   type (xbx_type),intent(inout) :: xbx     ! Header & non-gridded vars.
   type (be_type), intent(in)    :: be      ! background error structure.
   integer,        intent(in) :: cv_size

   integer, intent(in)        :: mz(:)

   real*8, intent(in) :: xhat_prev(cv_size)
   real*8, intent(in) :: xhat_curr(cv_size)

   !Check Variables
   real*8          :: diffvec(cv_size)
   real*8          :: temp

   integer         :: iens, i, j, k, ij, ic, itype
   real*8          :: EUCNORM_x(2)

   if (trace_use) call da_trace_entry("da_evaluate_increment")

   EUCNORM_x=0.D0

   !=====================================================================
   ! Compare Uxhat to Uxhat_prev (vtox)
   !=====================================================================
   do itype = 1,2
      if (itype .eq. 1) then
         write(unit=stdout,fmt='(A)') 'Evaluating (dxa_i - dxa_i-1)'
         diffvec = xhat_curr - xhat_prev
      end if

      if (itype .eq. 2) then
         write(unit=stdout,fmt='(A)') 'Evaluating (xa_i - xa_i-1)'
         diffvec = xhat_curr
      end if

      ! ||dxa_i - dxa_i-1||^2
      call da_transform_vtox(grid, be%cv%size_jb, xbx, be, &
                             grid%ep, diffvec(1:be%cv%size_jb), &
                             grid%vv, grid%vp, grid%vch)

      temp = 0.D0
      !Only for chem so far
#if (WRF_CHEM == 1)
      do ic = PARAM_FIRST_SCALAR, num_scaleant
         do k = 1, num_ant_steps
            do j=jts,jte
               do i=its,ite
                  temp = temp + grid%xach%alpha_ant(i,j,k,ic)**2
               end do
            end do
         end do
      end do

      do ic = PARAM_FIRST_SCALAR, num_scalebb
         do k = 1, num_bb_steps
            do j=jts,jte
               do i=its,ite
                  temp = temp + grid%xach%alpha_bb(i,j,k,ic)**2
               end do
            end do
         end do
      end do
#endif
#ifdef DM_PARALLEL
      EUCNORM_x(itype) = wrf_dm_sum_real(temp)
#else
      EUCNORM_x(itype) = temp
#endif
   end do

   write(unit=stdout,fmt='(A,2D20.10)') 'EUCNORM_x = ', EUCNORM_x

   if (trace_use) call da_trace_exit("da_evaluate_increment")

end subroutine da_evaluate_increment
#endif
