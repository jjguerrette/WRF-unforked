subroutine da_cv_indices(xp, mzs, cv_size, icv, jcv, kcv &
#if (WRF_CHEM == 1)
                        , vch &
#endif
                        )

   !---------------------------------------------------------------------------
   ! Purpose: Store i,j,k locations of each cv index.
   !---------------------------------------------------------------------------

   implicit none

#if (WRF_CHEM == 1)
   type (xch_type), optional, intent(inout) :: vch          ! Grdipt/EOF rcv for chem
#endif
   type (xpose_type), intent(in) :: xp          ! Dimensions and xpose buffers.
   integer,        intent(in)    :: mzs(:)      ! Background error structure levels.
   integer,        intent(in)    :: cv_size     ! Length of CV array.
   integer,        intent(inout) :: icv(1:cv_size) ! Control variables loc i
   integer,        intent(inout) :: jcv(1:cv_size) ! Control variables loc j
   integer,        intent(inout) :: kcv(1:cv_size) ! Control variables loc k

   integer   :: is,ie       ! Local grid range in y coordinate.
   integer   :: js,je       ! Local grid range in x coordinate.
   integer   :: ix,jy       ! Local grid horizontal dimensions.
   integer   :: mz          ! Max vertical coordinate for v1 through v5 arrays.
   integer   :: ne          ! Ensemble size.
   integer   :: cv_s,cv_e   ! Starting and ending indices into CV array.
   integer   :: ijm         ! Size of interior of v1 through v5 arrays.
   integer   :: ijmn        ! Size of interior of alpha cv arrays.
#if (WRF_CHEM == 1)
   integer   :: ma, mb      ! Max anthro and bioburn time coordinates
   integer   :: ic          ! chem loop iterator
#endif

   integer   :: i,j,k,ijk,m,n

   if (trace_use) call da_trace_entry("da_cv_indices")

   if( use_rf )then
      is = xp % its
      ie = xp % ite
      js = xp % jts
      je = xp % jte
   else
      call da_error(__FILE__,__LINE__,(/"This subroutine should not be called for use_rf = .false."/))
   endif
   ix = ie-is+1
   jy = je-js+1
   cv_e = 0
   ijk=0

   ! Store v1
   mz = mzs(1)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
        end do
     end do
   end if
   ! Store v2
   mz = mzs(2)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v3
   mz = mzs(3)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v4
   mz = mzs(4)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v5
   mz = mzs(5)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
#ifdef CLOUD_CV
   ! Store v6
   mz = mzs(6)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v7
   mz = mzs(7)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v8
   mz = mzs(8)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v9
   mz = mzs(9)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v10
   mz = mzs(10)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
   ! Store v11
   mz = mzs(11)
   if (mz > 0) then
      ijm = ix * jy * mz
      cv_s = cv_e + 1
      cv_e = cv_s + ijm - 1
      do k=1,mz
         do j=js,je
            do i=is,ie
               ijk = ijk + 1
               icv(ijk) = i
               jcv(ijk) = j
               kcv(ijk) = k
            end do
         end do
      end do
   end if
#endif
#ifdef CLOUD_CV
   ! Store alpha:
   mz = mzs(12)
   ne = mzs(13)
#else
   ! Store alpha:
   mz = mzs(6)
   ne = mzs(7)
#endif
   if ( ne > 0 ) then
      ix = ite_int - its_int + 1
      jy = jte_int - jts_int + 1
      ijmn = ix * jy * mz * ne
      cv_s = cv_e + 1
      cv_e = cv_s + ijmn - 1
      do n = 1, ne
         do m = 1, mz
            do j = jts_int, jte_int !js, je
               do i = its_int, ite_int ! is, ie
                  ijk = ijk + 1
                  icv(ijk) = i
                  jcv(ijk) = j
                  kcv(ijk) = k
               end do
            end do
         end do
      end do
   end if

#if (WRF_CHEM == 1)
   if (present(vch) ) then
      ix = ie-is+1
      jy = je-js+1
#ifdef CLOUD_CV
      ma = mzs(14)
      mb = mzs(15)
#else
      ma = mzs(8)
      mb = mzs(9)
#endif
      if ( ma > 0 ) then
         ijm = ix * jy * ma
         do ic = PARAM_FIRST_SCALAR, num_scaleant
            cv_s = cv_e + 1
            cv_e = cv_s + ijm - 1
            do k=1,ma
               do j=js,je
                  do i=is,ie
                     ijk = ijk + 1
                     icv(ijk) = i
                     jcv(ijk) = j
                     kcv(ijk) = k
                  end do
               end do
            end do
         end do
      end if

      if ( mb > 0 ) then
         ijm = ix * jy * mb
         do ic = PARAM_FIRST_SCALAR, num_scalebb
            cv_s = cv_e + 1
            cv_e = cv_s + ijm - 1
            do k=1,mb
               do j=js,je
                  do i=is,ie
                     ijk = ijk + 1
                     icv(ijk) = i
                     jcv(ijk) = j
                     kcv(ijk) = k
                  end do
               end do
            end do
         end do
      end if
   end if

#endif

   if (trace_use) call da_trace_exit("da_cv_indices")

end subroutine da_cv_indices


