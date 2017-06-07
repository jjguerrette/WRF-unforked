subroutine da_precon_omega(grid, cv_size, mz, jp_start, jp_end, nens, ghat )

   !-------------------------------------------------------------------------
   ! Purpose:   Select preconditioned omega vectors for next RIOT outer loop
   !
   ! Called from da_randomise_svd and da_randomise_svd_51
   !
   ! History: 04/24/2017  Creation (JJ Guerrette)
   !
   ! References: Bousserez et al. (2017)
   !
   !-------------------------------------------------------------------------

   implicit none

   type(domain), intent(inout)  :: grid
   integer, intent(in)          :: cv_size                ! Total cv size
   integer, intent(in)          :: mz(:)
   integer, intent(in)          :: jp_start, jp_end
   real*8, intent(in)           :: ghat(1:cv_size)        ! control variable (local).
   integer, intent(in)          :: nens                   ! Total number of eigenpairs

   real*8                       :: omega(cv_size,nens) ! Selected Hessian eigenvectors
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
   real*8   :: crit_save(nens)
   integer  :: precon_list(nens)
   integer    :: it, eig_count, eig_total, crit_imin(1), nens_it, neign_it, nouter_seed

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

!   if (riot_precon .eq. 1) then
!
!      allocate(crit(ensdim_svd*(svd_outer-1)))
!      allocate(crit_prev(ensdim_svd*(svd_outer-1)))
!
!      do it=svd_outer-1, 1, -1
!
!         nens_it = ensdim_svd
!         if (it.gt.1) then
!            nens_it = min(ntmax(it),ensdim_svd)
!         end if
!
!         allocate(eignvec_prev(nens_it,nens_it))
!         allocate(eignval_prev(nens_it))
!         allocate(qhat_prev(cv_size,0:nens_it))
!
!         if(eig_count .eq. 0) then !!.and. riot_precon.eq.2) then
!            call da_hessian_io ('r',grid, cv_size,it,nens_it,neign_it,&
!                               eignvec_prev,eignval_prev,qhat_prev,mz,ghat_prev)
!            gdot_prev = sqrt(da_dot_cv(cv_size, ghat_prev, ghat_prev, grid, mz, jp_start, jp_end))
!         else
!            call da_hessian_io ('r',grid, cv_size,it,nens_it,neign_it,&
!                               eignvec_prev,eignval_prev,qhat_prev,mz)
!         end if
!
!         do iens = 1, neign_it
!            do i = 1, cv_size
!               nuhat(i,1) = da_dot(nens_it, qhat_prev(i,1:nens_it), eignvec_prev(1:nens_it,iens))
!            end do
!
!            eig_count = eig_count + 1
!
!!Criteria 1
!            crit(eig_count) = &
!               (da_dot_cv(cv_size, nuhat(:,1), ghat, grid, &
!                          mz, jp_start, jp_end) / gdot ) **2
!            crit_prev(eig_count) = &
!               (da_dot_cv(cv_size, nuhat(:,1), ghat_prev, grid, &
!                          mz, jp_start, jp_end) / gdot_prev ) **2
!
!!Criteria 2
!!!            crit(eig_count) = &
!!!               (da_dot_cv(cv_size, nuhat(:,1), ghat, grid, &
!!!                        mz, jp_start, jp_end) / gdot)**2 / eignval_prev(iens)
!
!!Criteria 3
!!!            crit(eig_count) = &
!!!               (da_dot_cv(cv_size, nuhat(:,1), ghat, grid, &
!!!                        mz, jp_start, jp_end) / gdot / eignval_prev(iens))**2
!
!
!            if (eig_count .le. nens) then
!               crit_save(eig_count) = crit(eig_count)
!               precon_list(eig_count) = eig_count
!               omega(:,eig_count) = nuhat(:,1)
!            else
!               crit_imin(1:1) = minloc(crit_save)
!
!               if(crit(eig_count) .ge. crit_save(crit_imin(1))) then
!                  crit_save(crit_imin(1)) = crit(eig_count)
!                  precon_list(crit_imin(1)) = eig_count
!                  omega(:,crit_imin(1)) = nuhat(:,1)
!               end if
!            end if
!
!         end do
!
!         deallocate(eignvec_prev)
!         deallocate(eignval_prev)
!         deallocate(qhat_prev)
!      end do
!      write(cens,fmt='(I0)') eig_count
!      write(unit=stdout,fmt='("crit_prev=",'//trim(cens)//'D15.5)') crit_prev(1:eig_count)
!      write(unit=stdout,fmt='("crit=",'//trim(cens)//'D15.5)') crit(1:eig_count)
!
!      do iens = 1, nens
!         write(unit=stdout,fmt='("iens = ",i3,", ieig = ",i3,", crit(i) = ",D15.5)') iens, precon_list(iens), crit_save(iens)
!      end do
!
!      deallocate(crit)
!      deallocate(crit_prev)
!
!   end if

   if (riot_precon .eq. 1) nouter_seed = svd_outer-1
   if (any(riot_precon .eq. (/2,3,4/))) nouter_seed = svd_outer

   if(riot_precon .eq. 4) then
      eig_total = 0
      do it=nouter_seed, 1, -1
         nens_it = ensdim_svd
         if (it.gt.1) then
            nens_it = min(ntmax(it),ensdim_svd)
         end if

         call da_hessian_io ('r',grid, cv_size,it,nens_it,neign_it)
         eig_total = eig_total + neign_it
      end do
   end if

   allocate(crit(ensdim_svd*nouter_seed))
   allocate(lam_all(ensdim_svd*nouter_seed))

   do it=nouter_seed, 1, -1
      nens_it = ensdim_svd
      if (it.gt.1) then
         nens_it = min(ntmax(it),ensdim_svd)
      end if

      allocate(eignvec_prev(nens_it,nens_it))
      allocate(eignval_prev(nens_it))
      allocate(qhat_prev(cv_size,0:nens_it))

      call da_hessian_io ('r',grid, cv_size,it,nens_it,neign_it,&
                         eignvec_prev,eignval_prev,qhat_prev,mz)

      do iens = 1, neign_it
         do i = 1, cv_size
            nuhat(i,1) = da_dot(nens_it, qhat_prev(i,1:nens_it), eignvec_prev(1:nens_it,iens))
         end do

         eig_count = eig_count + 1

         if (riot_precon .eq. 1 .or. riot_precon .eq. 2) &
            crit(eig_count) = ( da_dot_cv(cv_size, nuhat(:,1), ghat, &
                                           grid, mz, jp_start, jp_end) &
                                            / gdot ) **2

         if (riot_precon .eq. 3) &
            crit(eig_count) = eignval_prev(iens)

         if (riot_precon .eq. 4) &
            crit(eig_count) = real(eig_total - eig_count + 1,8)

         lam_all(eig_count) = eignval_prev(iens)

         if (eig_count .le. nens) then
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

   if (eig_count .lt. nens) then
      crit_save(eig_count+1:nens) = 0.D0

      nouter_seed = nens-eig_count

      allocate(qhat_prev(cv_size,nouter_seed))

      call da_gen_omega(qhat_prev, cv_size, nouter_seed, grid, mz)

      omega(:,eig_count+1:nens) = qhat_prev

      deallocate(qhat_prev)
   end if

!   write(cens,fmt='(I0)') eig_count
   write(unit=stdout,fmt='("crit                    lam")')
   do iens = 1, eig_count
      write(unit=stdout,fmt='(2D26.16)') crit(iens), lam_all(iens)
   end do

   do iens = 1, nens
      write(unit=stdout,fmt='("iens = ",i3,", ieig = ",i3,", crit(i) = ",D26.16)') iens, precon_list(iens), crit_save(iens)
   end do

   deallocate(crit)
   deallocate(lam_all)


   !-------------------------------------------------------------------------
   ! [3.0] Write preconditioned omega vectors
   !-------------------------------------------------------------------------

   ! Write preconditioned omega vectors from all ensembles on each processor
   write(cproc,fmt='(i4.4)') myproc
   do iens = 1, nens
      write(unit=stdout,fmt='("Writing preconditioned omega ensemble; iens = ",i3)') iens
      write(cens,fmt='(i4.4)') iens

      filename = 'omega_precon.e'//trim(adjustl(cens))//'.p'//trim(adjustl(cproc))
      call da_cv_io ('w', filename, cv_size, omega(1:cv_size,iens), grid, mz)
   end do


   if (trace_use) call da_trace_exit("da_precon_omega")


end subroutine da_precon_omega


