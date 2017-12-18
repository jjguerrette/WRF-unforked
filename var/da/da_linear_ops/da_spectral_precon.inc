subroutine da_spectral_precon(grid, cv_size, mz, jp_start, jp_end, nouter_prev, cv_v, cv_z, mode )

   !-------------------------------------------------------------------------
   ! Purpose:   Apply cummulative square-root spectral preconditioner
   !
   ! Called from ???
   !
   ! History: 10/31/2017  Creation (JJ Guerrette)
   !
   ! References: Tshimanga et al. (2008) and Gratton et al. (2011)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)  :: grid
   integer, intent(in)          :: cv_size                ! Total cv size
   integer, intent(in)          :: mz(:)
   integer, intent(in)          :: jp_start, jp_end
   integer, intent(in)          :: nouter_prev
   real*8, intent(in)           :: cv_v(cv_size)    ! Background square-root preconditioned CV
   real*8, intent(out)          :: cv_z(cv_size)    ! Spectral square-root preconditioned CV
   character*2, intent(in)      :: mode

   real*8                       :: cv_z0(cv_size)
   integer                      :: it, it_inner, l
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   real*8, allocatable   :: eignvec_prev(:,:)
   real*8, allocatable   :: eignval_prev(:)
   real*8, allocatable   :: qhat_prev(:,:)

   integer  :: nmode_it, neign_it, outer_start, outer_end, increment
   real     :: expon

   if (trace_use) call da_trace_entry("da_spectral_precon")

   write(unit=stdout,fmt='(2A)') 'Apply Spectral preconditioner in mode: ', trim(mode)
   write(unit=stdout,fmt=*) ' '

   if (trim(mode) == 'TL') then
      outer_start = nouter_prev
      outer_end = 1
      increment = -1
      expon = -0.5D0
   else if (trim(mode) == 'AD') then
      outer_start = 1
      outer_end = nouter_prev
      increment = 1
      expon = -0.5D0
   else if (trim(mode) == 'IT') then
      outer_start = 1
      outer_end = nouter_prev
      increment = 1
      expon = 0.5D0
   else if (trim(mode) == 'IA') then
      outer_start = nouter_prev
      outer_end = 1
      increment = -1
      expon = 0.5D0
   else
      call da_error(__FILE__,__LINE__, &
       (/"Spectral preconditioning must use TL, AD, IT, or IA mode."/))
   end if

   cv_z = cv_v

   do it = outer_start, outer_end, increment
      if (.not.any((/use_lanczos,use_randomblock/))) &
         call da_error(__FILE__,__LINE__, &
          (/"Spectral preconditioning can only be used with Lanczos or RSVD."/))

      call da_hessian_io ('r',grid,cv_size,it,nmode_it,neign_it)

      allocate(eignval_prev(nmode_it))
      allocate(eignvec_prev(nmode_it,nmode_it))
      allocate(qhat_prev(cv_size,0:nmode_it))

      call da_hessian_io ('r',grid,cv_size,it,nmode_it,neign_it,&
                         eignvec_prev,eignval_prev,qhat_prev,mz)

      !Need to add capability for user to set l <= neign_it for each set of eigenpairs
      l = neign_it
!      l = min(neign_it,l_user_input(it))

      do it_inner = 1, l
         eignval_prev(it_inner) = eignval_prev(it_inner)**expon - 1.D0
      end do

      ! cv_z^it = (I_n + V (LAM^-0.5 - I_l) V^T ) * cv_z^it-1
      cv_z0 = 0.D0
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             nmode_it,nmode_it,l, &
                             eignval_prev,eignvec_prev,qhat_prev,&
                             cv_z,cv_z0)
      cv_z = cv_z + cv_z0

      deallocate(eignval_prev)
      deallocate(eignvec_prev)
      deallocate(qhat_prev)
   end do

   if (trace_use) call da_trace_exit("da_spectral_precon")

end subroutine da_spectral_precon
