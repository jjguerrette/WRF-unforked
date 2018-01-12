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
   integer                      :: it, imode, l
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   real*8, allocatable   :: eignvec(:,:)
   real*8, allocatable   :: eignval(:)
   real*8, allocatable   :: qhat(:,:)

   real*8, allocatable   :: nuhat_all(:,:)
   real*8, allocatable   :: lam_all(:)
   real*8, allocatable   :: nuhat(:,:)

   integer  :: i, nmodes, neign, &
               outer_start, outer_end, increment, &
               eig_count, eig_total

   if (trace_use) call da_trace_entry("da_rotate_omega")

   write(unit=stdout,fmt='(A)') 'Prerotate Gaussian random vector with previous Hessian eigenvectors for RSVD. '
   write(unit=stdout,fmt=*) ' '

   eig_total = 0
   do it = 1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmodes,neign)
      eig_total = eig_total + neign
   end do
   allocate(nuhat_all(cv_size,0:eig_total))
   allocate(lam_all(eig_total))

   nuhat_all = 0.D0
   eig_count = 0

   cv_v = cv_z

   do it = 1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmodes,neign)

      allocate(eignvec(nmodes,nmodes))
      allocate(eignval(nmodes))
      allocate(qhat(cv_size,0:nmodes))

      call da_hessian_io ('r',grid, cv_size,it,nmodes,neign,&
                         eignvec,eignval,qhat,mz)

      ! Precondition with the l leading eigen modes from iteration it
      !  (requires eignval/eignvec to be sorted by decreasing value)
      l = max(min(neign,spectral_trunc(it)),0)
      if (l .gt. 0) then
         do imode = 1, l
            eig_count = eig_count + 1
            lam_all(eig_count) = eignval(imode)
            do i = 1, nmodes
               nuhat_all(:,eig_count) = nuhat_all(:,eig_count) + qhat(:,i) * eignvec(i,imode)
            end do
         end do
      end if

      deallocate(eignvec)
      deallocate(eignval)
      deallocate(qhat)
   end do

   eig_total = eig_count

   allocate(eignval(eig_total))
   allocate(eignvec(eig_total,eig_total))
   eignval = 1.D0
  
   call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                         nuhat_all(:,0:eig_total),cv_size,eig_total)

   ! cv_v = (I_n - V V^T ) * cv_v
   cv_v0 = 0.D0
   call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                          eig_total,eig_total,eig_total, &
                          eignval,eignvec,nuhat_all(:,0:eig_total),&
                          cv_v,cv_v0, .true.)
   cv_v = cv_v - cv_v0
   deallocate(eignval)
   deallocate(eignvec)


   !! Multiply by inverse preconditioner to recover OMEGA in v-space
   eig_count = 0
   do it=1, nouter_prev
      call da_hessian_io ('r',grid, cv_size,it,nmodes,neign)

      allocate(eignval(neign))
      allocate(eignvec(neign,neign))
      allocate(nuhat(cv_size,0:neign))

      ! Precondition with the l leading eigen modes from iteration it
      !  (requires eignval/eignvec to be sorted by decreasing value)
      l = max(min(neign,spectral_trunc(it)),0)
      if (l .gt. 0) then
         eignval = 0.D0
         do imode = 1, l
            eig_count = eig_count + 1
            eignval(imode) = lam_all(eig_count)**0.5D0 - 1.D0
            nuhat(:,imode) = nuhat_all(:,eig_count)
         end do

         ! cv_v^it = (I_n + V (LAM^0.5 - I_l) V^T ) * cv_v^it-1
         cv_v0 = 0.D0
         call da_amat_mul_trunc(grid,mz,jp_start,jp_end,cv_size,&
                                neign,neign,l, &
                                eignval,eignvec,nuhat,&
                                cv_v,cv_v0,.true.)
         cv_v = cv_v + cv_v0
      end if

      deallocate(eignval)
      deallocate(eignvec)
      deallocate(nuhat)
   end do

   deallocate(nuhat_all)
   deallocate(lam_all)

   if (trace_use) call da_trace_exit("da_rotate_omega")

end subroutine da_rotate_omega
