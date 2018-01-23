subroutine da_precon_omega(grid, cv_size, mz, jp_start, jp_end, nsamp, ghat )

   !-------------------------------------------------------------------------
   ! Purpose:   Select preconditioned omega vectors for next RIOT outer loop
   !
   ! Called from da_randomise_svd and da_randomise_svd_51
   !
   ! History: 04/24/2017  Creation (JJ Guerrette)
   !
   ! References: Bousserez, Guerrette, and Henze (2017) [unpublished]
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)  :: grid
   integer, intent(in)          :: cv_size                ! Total cv size
   integer, intent(in)          :: mz(:)
   integer, intent(in)          :: jp_start, jp_end
   real*8, intent(in)           :: ghat(1:cv_size)        ! control variable (local).
   integer, intent(in)          :: nsamp                   ! Total number of eigenpairs

   real*8                       :: omega(cv_size,nsamp) ! Selected Hessian eigenvectors
   real*8                       :: nuhat(cv_size,1) ! Eigenvectors of A
   real*8                       :: ghat_prev(cv_size)

   integer                      :: iens, jens, i, j
   character*10                 :: cproc, cens
   character(len=filename_len)  :: filename

   real*8, allocatable   :: eignvec_prev(:,:)
   real*8, allocatable   :: eignval_prev(:)
   real*8, allocatable   :: qhat_prev(:,:)

   real*8, allocatable   :: crit(:)
   real*8, allocatable   :: crit_prev(:)
   real*8, allocatable   :: lam_all(:)

   real*8   :: gdot, gdot_prev
   real*8   :: crit_save(nsamp)
   integer  :: precon_list(nsamp)
   integer    :: it, eig_count, eig_total, crit_imin(1), &
                 nsamp_it, neign_it, num_omega_it, nouter_seed

   if (trace_use) call da_trace_entry("da_precon_omega")

   write(unit=stdout,fmt='(A)') 'Minimize cost function using randomized SVD method'
   write(unit=stdout,fmt=*) ' '

   !-------------------------------------------------------------------------
   ! [1.0] Initialization:
   !-------------------------------------------------------------------------

   eig_count = 0

   gdot = sqrt(da_dot_cv(cv_size, ghat, ghat, grid, mz, jp_start, jp_end))

   !-------------------------------------------------------------------------
   ! [2.0] Select preconditioned omega vectors
   !-------------------------------------------------------------------------

   nouter_seed = rand_outer-1
   if (any(riot_precon .eq. (/12,13,14,15/))) nouter_seed = rand_outer

   if ( any(riot_precon .eq. (/14,15/)) ) then
      eig_total = 0
      do it=nouter_seed, 1, -1
         call da_hessian_io ('r',grid, cv_size,it,nsamp_it,neign_it)
         if(riot_precon .eq. 14) eig_total = eig_total + neign_it
         if(riot_precon .eq. 15) eig_total = eig_total + nsamp_it
      end do
   end if

   allocate(crit(sum(ntmax(1:nouter_seed))))
   allocate(lam_all(sum(ntmax(1:nouter_seed))))

   do it=nouter_seed, 1, -1
      call da_hessian_io ('r',grid, cv_size,it,nsamp_it,neign_it)

      allocate(eignvec_prev(nsamp_it,nsamp_it))
      allocate(eignval_prev(nsamp_it))
      allocate(qhat_prev(cv_size,0:nsamp_it))

      call da_hessian_io ('r',grid, cv_size,it,nsamp_it,neign_it,&
                         eignvec_prev,eignval_prev,qhat_prev,mz)

      if (any(riot_precon .eq. (/3,4,15/))) then
         num_omega_it = nsamp_it
      else
         num_omega_it = neign_it
      end if
      do iens = 1, num_omega_it
         if (any(riot_precon .eq. (/3,4,15/))) then
            nuhat(:,1) = qhat_prev(:,iens)
         else
            do i = 1, cv_size
               nuhat(i,1) = da_dot(nsamp_it, qhat_prev(i,1:nsamp_it), eignvec_prev(1:nsamp_it,iens))
            end do
         end if
         eig_count = eig_count + 1

         if (any(riot_precon .eq. (/1,2,3,4,12/))) &
            crit(eig_count) = ( da_dot_cv(cv_size, nuhat(:,1), ghat, &
                                           grid, mz, jp_start, jp_end) &
                                            / gdot ) **2

         if (riot_precon .eq. 13) &
            crit(eig_count) = eignval_prev(iens)

         if (any(riot_precon .eq. (/14,15/))) &
            crit(eig_count) = real(eig_total - eig_count + 1,8)

         lam_all(eig_count) = eignval_prev(iens)

         if (eig_count .le. nsamp) then
            crit_save(eig_count) = crit(eig_count)
            precon_list(eig_count) = eig_count
            omega(:,eig_count) = nuhat(:,1)
         else
            crit_imin(1:1) = minloc(crit_save)

            if(crit(eig_count) .ge. crit_save(crit_imin(1))) then
               crit_save(crit_imin(1)) = crit(eig_count)
               precon_list(crit_imin(1)) = eig_count
               omega(:,crit_imin(1)) = nuhat(:,1)
            end if
         end if
      end do

      deallocate(eignvec_prev)
      deallocate(eignval_prev)
      deallocate(qhat_prev)
   end do

   if (eig_count .lt. nsamp) then
      crit_save(eig_count+1:nsamp) = 0.D0

      nouter_seed = nsamp-eig_count

      allocate(qhat_prev(cv_size,nouter_seed))
      call da_gen_omega(qhat_prev, cv_size, nouter_seed, grid, mz)

      omega(:,eig_count+1:nsamp) = qhat_prev

      deallocate(qhat_prev)
   end if

   write(unit=stdout,fmt='("    crit                      lam")')
   do iens = 1, eig_count
      write(unit=stdout,fmt='(2D26.16)') crit(iens), lam_all(iens)
   end do

   do iens = 1, nsamp
      write(unit=stdout,fmt='("iens = ",i4,", ieig = ",i4,", crit(i) = ",D26.16)') iens, precon_list(iens), crit_save(iens)
   end do

   deallocate(crit)
   deallocate(lam_all)

   !Prepend with gradient vector
   if (any(riot_precon .eq. (/2,4/))) then
      omega(:,2:nsamp) = omega(:,1:nsamp-1)
      omega(:,1) = ghat / gdot
   end if

   !Orthonormalization - should do this every time
!   if (any(riot_precon .eq. (/1,2,3,4,15/))) then
      allocate(qhat_prev(cv_size,0:nsamp))
      qhat_prev(:,1:nsamp) = omega(:,1:nsamp)
      call da_gram_schmidt(grid,mz,jp_start,jp_end, &
                            qhat_prev,cv_size,nsamp)

      omega(:,1:nsamp) = qhat_prev(:,1:nsamp)
      deallocate(qhat_prev)
!   end if

   !-------------------------------------------------------------------------
   ! [3.0] Write preconditioned omega vectors
   !-------------------------------------------------------------------------

   ! Write preconditioned omega vectors from all ensembles on each processor
   write(cproc,fmt='(i4.4)') myproc
   do iens = 1, nsamp
      write(unit=stdout,fmt='("Writing preconditioned omega ensemble; iens = ",i4)') iens
      write(cens,fmt='(i4.4)') iens

      filename = 'omega_precon.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
      call da_cv_io ('w', filename, cv_size, omega(1:cv_size,iens), grid, mz)
   end do


   if (trace_use) call da_trace_exit("da_precon_omega")


end subroutine da_precon_omega


