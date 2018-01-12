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
   integer                      :: it, imode, l
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   real*8, allocatable   :: eignvec(:,:)
   real*8, allocatable   :: eignval(:)
   real*8, allocatable   :: qhat(:,:)

   integer  :: nmodes, neign, outer_start, outer_end, increment
   real     :: expon

   if (trace_use) call da_trace_entry("da_spectral_precon")

   write(unit=stdout,fmt='(2A)') 'Apply Spectral preconditioner in mode: ', trim(mode)
   write(unit=stdout,fmt=*) ' '

   if (trim(mode) == 'TL') then !Tangent preconditioner
      outer_start = nouter_prev
      outer_end = 1
      increment = -1
      expon = -0.5D0
   else if (trim(mode) == 'AD') then !Adjoint or Transpose of preconditioner
      outer_start = 1
      outer_end = nouter_prev
      increment = 1
      expon = -0.5D0
   else if (trim(mode) == 'IT') then !Inverse of tangent preconditioner
      outer_start = 1
      outer_end = nouter_prev
      increment = 1
      expon = 0.5D0
   else if (trim(mode) == 'IA') then !Inverse of adjoint preconditioner
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

      call da_hessian_io ('r',grid,cv_size,it,nmodes,neign)

      allocate(eignval(nmodes))
      allocate(eignvec(nmodes,nmodes))
      allocate(qhat(cv_size,0:nmodes))

      call da_hessian_io ('r',grid,cv_size,it,nmodes,neign,&
                         eignvec,eignval,qhat,mz)

      ! Precondition with the l leading eigen modes from iteration it
      !  (requires eignval/eignvec to be sorted by decreasing value)
      l = max(min(neign,spectral_trunc(it)),0)
      if (l .gt. 0) then
         do imode = 1, l
            eignval(imode) = eignval(imode)**expon - 1.D0
         end do

         ! cv_z^it = (I_n + V (LAM^-0.5 - I_l) V^T ) * cv_z^it-1
         cv_z0 = 0.D0
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                                nmodes,neign,l, &
                                eignval,eignvec,qhat,&
                                cv_z,cv_z0)
         cv_z = cv_z + cv_z0
      end if

      deallocate(eignval)
      deallocate(eignvec)
      deallocate(qhat)
   end do

   if (trace_use) call da_trace_exit("da_spectral_precon")

end subroutine da_spectral_precon
