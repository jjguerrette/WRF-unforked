subroutine da_gen_omega(omega, vec_len, nvec, grid, mz )

   !-------------------------------------------------------------------------
   ! Purpose:   Generate Gaussian omega with correct stats, single global seed
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !-------------------------------------------------------------------------

   implicit none
   integer, intent(in)             :: vec_len, nvec
   real(wp), intent(inout)         :: omega(1:vec_len,1:nvec) ! Local Gaussian draw
   type (domain),  intent(inout)   :: grid
   integer, intent(in)             :: mz(:)

   real(wp), allocatable        :: omega1d_g(:) ! 1-D Global Gaussian draw
   real(wp), allocatable        :: dumcv(:)
   real(wp), allocatable        :: dumcv_g(:)

   integer         :: ivec, i, j, k, cv_s, cv_f, vec_len_domain

   if (trace_use) call da_trace_entry("da_gen_omega")

! Solves this problem: da_set_randomcv often uses the same seed for two different draws, giving identical cv's for multiple ensemble members on a particular processor.  Either the seed needs to be forced to change between multiple draws or it needs to be the same and only a single call to da_set_randomcv can be used for all cv values.

#ifdef DM_PARALLEL
   ! Generate global random omega matrix
   ! Sequence from single global seed (non-repeating)
   vec_len_domain = wrf_dm_sum_integer(vec_len)

   if (rootproc) then
      allocate(omega1d_g(1:vec_len_domain*nvec))
      allocate(dumcv_g(1:vec_len_domain))
      call da_set_randomcv(vec_len_domain*nvec, omega1d_g)
   else
      allocate(omega1d_g(1))
      allocate(dumcv_g(1))
   end if

   !Transfer to local omega
   allocate(dumcv(1:vec_len))
   cv_f = 0
   do ivec = 1, nvec
      cv_s = cv_f + 1
      cv_f = cv_s + vec_len_domain-1
      if(rootproc) then
         dumcv_g = omega1d_g(cv_s:cv_f)
      else
         dumcv_g = 0.D0
      end if
      dumcv = 0.D0

      call mpi_barrier(comm, ierr)

      call da_global_to_cv(vec_len, vec_len_domain, dumcv, grid, mz, dumcv_g)
      omega(:,ivec) = dumcv
   end do
   deallocate(omega1d_g)
   deallocate(dumcv)
   deallocate(dumcv_g)
#else
   ! Generate local=global omega
   ! Single seed
   allocate(omega1d_g(1:vec_len*nvec))
   call da_set_randomcv(vec_len*nvec, omega1d_g)
   cv_s = 1
   do ivec = 1,nvec
      cv_f = cv_s + vec_len-1
      omega(:,ivec) = omega1d_g(cv_s:cv_f)
      cv_s = cv_f+1
   end do
   deallocate(omega1d_g)
#endif

   if (trace_use) call da_trace_exit("da_gen_omega")

end subroutine da_gen_omega

