subroutine da_rotate_omega(grid, cv_size, mz, jp_start, jp_end, nouter_prev, cv_z, cv_v )

   !-------------------------------------------------------------------------
   ! Purpose:   Rotate random vectors into the space complementary to earlier eigenvectors
   !
   ! Called from da_randomise_svd.inc
   !
   ! History: 10/31/2017  Creation (JJ Guerrette)
   !
   ! References: Bousserez, Guerrette, and Henze; unpublished
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)  :: grid
   integer, intent(in)          :: cv_size                ! Total cv size
   integer, intent(in)          :: mz(:)
   integer, intent(in)          :: jp_start, jp_end
   integer, intent(in)          :: nouter_prev
   real*8, intent(in)           :: cv_z(cv_size)    ! Spectral square-root preconditioned CV
   real*8, intent(out)          :: cv_v(cv_size)    ! Background square-root preconditioned CV

   real*8                       :: cv_v0(cv_size)
   integer                      :: it, it_inner, l
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   real*8, allocatable   :: eignvec_prev(:,:)
   real*8, allocatable   :: eignval_prev(:)
   real*8, allocatable   :: qhat_prev(:,:)

   real*8, allocatable   :: nuhat_all(:,:)
   real*8, allocatable   :: lam_all(:)
   real*8, allocatable   :: nuhat(:,:)

   integer  :: i, nmode_it, neign_it, &
               outer_start, outer_end, increment, &
               eig_count, eig_total

   if (trace_use) call da_trace_entry("da_rotate_omega")

   write(unit=stdout,fmt='(A)') 'Prerotate Gaussian random vector with previous Hessian eigenvectors for RSVD. '
   write(unit=stdout,fmt=*) ' '

   eig_total = 0
   do it = 1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmode_it,neign_it)
      eig_total = eig_total + neign_it
   end do
   allocate(nuhat_all(cv_size,0:eig_total))
   allocate(lam_all(eig_total))

   eig_count = 0

  cv_v = cv_z

   do it = 1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmode_it,neign_it)

      allocate(eignvec_prev(nmode_it,nmode_it))
      allocate(eignval_prev(nmode_it))
      allocate(qhat_prev(cv_size,0:nmode_it))

      call da_hessian_io ('r',grid, cv_size,it,nmode_it,neign_it,&
                         eignvec_prev,eignval_prev,qhat_prev,mz)
      do it_inner = 1, neign_it
         eig_count = eig_count + 1

         lam_all(eig_count) = eignval_prev(it_inner)
         do i = 1, cv_size
            nuhat_all(i,eig_count) = da_dot(nmode_it, qhat_prev(i,1:nmode_it), eignvec_prev(1:nmode_it,it_inner))
         end do
      end do
      deallocate(eignvec_prev)
      deallocate(eignval_prev)
      deallocate(qhat_prev)
   end do

   allocate(eignval_prev(eig_total))
   allocate(eignvec_prev(eig_total,eig_total))
   eignval_prev = 1.D0
  
   call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         nuhat_all,cv_size,eig_total)

   ! cv_v = (I_n - V V^T ) * cv_v
   cv_v0 = 0.D0
   call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                          eig_total,eig_total,eig_total, &
                          eignval_prev,eignvec_prev,nuhat_all,&
                          cv_v,cv_v0, .true.)
   cv_v = cv_v - cv_v0
   deallocate(eignval_prev)
   deallocate(eignvec_prev)


   !! Multiply by inverse preconditioner to recover OMEGA in v-space
   eig_count = 0
   do it=1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmode_it,neign_it)

      allocate(eignval_prev(neign_it))
      allocate(eignvec_prev(neign_it,neign_it))
      allocate(nuhat(cv_size,0:neign_it))

      !Need to add capability for user to set l <= neign_it for each set of eigenpairs
      l = neign_it
!      l = min(neign_it,l_user_input(it))

      eignval_prev = 0.D0
      do it_inner = 1, neign_it
         eig_count = eig_count + 1
         if (it_inner .le. l) then
            eignval_prev(it_inner) = lam_all(eig_count)**0.5D0 - 1.D0
            nuhat(:,it_inner) = nuhat_all(:,eig_count)
         end if
      end do

      ! cv_v^it = (I_n + V (LAM^0.5 - I_l) V^T ) * cv_v^it-1
      cv_v0 = 0.D0
      call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                             neign_it,neign_it,l, &
                             eignval_prev,eignvec_prev,nuhat,&
                             cv_v,cv_v0,.true.)
      cv_v = cv_v + cv_v0

      deallocate(eignval_prev)
      deallocate(eignvec_prev)
      deallocate(nuhat)
   end do

   deallocate(nuhat_all)
   deallocate(lam_all)

   if (trace_use) call da_trace_exit("da_rotate_omega")

end subroutine da_rotate_omega
