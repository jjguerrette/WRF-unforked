subroutine da_gen_omega(omega_mat, vec_len, nvec, grid, mz )

   !-------------------------------------------------------------------------
   ! Purpose:   Generate Gaussian omega_mat with correct stats, single global seed
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !-------------------------------------------------------------------------

   implicit none
   real(wp), allocatable, intent(inout) :: omega_mat(:,:) ! Gaussian draw
   integer, intent(in)                :: vec_len, nvec
   type (domain),  intent(inout)      :: grid
   integer, intent(in)                :: mz(:)

   real(wp), allocatable        :: omega_g(:,:) ! Global Gaussian draw
   real(wp), allocatable        :: omega1d(:) ! Gaussian draw
   real(wp), allocatable        :: dumcv(:)
   real(wp), allocatable        :: dumcv_g(:)

   integer         :: ivec, i, j, k, cv_s, cv_f, vec_len_domain

   if (trace_use) call da_trace_entry("da_gen_omega")

   if(allocated(omega_mat)) deallocate(omega_mat)
   allocate(omega_mat(1:vec_len,1:nvec))

! Solves this problem: da_set_randomcv often uses the same seed for two different draws, giving identical cv's for multiple ensemble members on a particular processor.  Either the seed needs to be forced to change between multiple draws or it needs to be the same and only a single call to da_set_randomcv can be used for all cv values.

#ifdef DM_PARALLEL
   ! Generate global random omega matrix
   ! Sequence from single global seed (non-repeating)
   vec_len_domain = wrf_dm_sum_integer(vec_len)

   if (rootproc) then
      allocate(omega_g(1:vec_len_domain,nvec))
      allocate(omega1d(1:vec_len_domain*nvec))
      call da_set_randomcv(vec_len_domain*nvec, omega1d)
      cv_s = 1
      do ivec = 1,nvec
         cv_f = cv_s + vec_len_domain-1
         omega_g(:,ivec) = omega1d(cv_s:cv_f)
         cv_s = cv_f+1
      end do
      deallocate(omega1d)
   else
      allocate(omega_g(1,1:nvec))
   end if

   call mpi_barrier(comm, ierr)

   !Transfer to local omega
   allocate(dumcv(1:vec_len))
   allocate(dumcv_g(1:vec_len_domain))
   do ivec = 1, nvec
      if(rootproc) then
         dumcv_g = omega_g(:,ivec)
      else
         dumcv_g = 0.D0
      end if
      dumcv = 0.D0

      call da_global_to_cv(vec_len, vec_len_domain, dumcv, grid, mz, dumcv_g)
      omega_mat(:,ivec) = dumcv

      call mpi_barrier(comm, ierr)
   end do
   deallocate(dumcv)
   deallocate(dumcv_g)
   deallocate(omega_g)
#else
   ! Generate local=global omega
   ! Single seed
   allocate(omega1d(1:vec_len*nvec))
   call da_set_randomcv(vec_len*nvec, omega1d)
   cv_s = 1
   do ivec = 1,nvec
      cv_f = cv_s + vec_len-1
      omega_mat(:,ivec) = omega1d(cv_s:cv_f)
      cv_s = cv_f+1
   end do
   deallocate(omega1d)
#endif

   if (trace_use) call da_trace_exit("da_gen_omega")

end subroutine da_gen_omega

