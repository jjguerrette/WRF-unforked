#if defined(LAPACK)
subroutine da_write_norms(    EUCNORM_x, FROBNORM_Pa, SPECNORM_Pa, &
                              EUCNORM_v, FROBNORM_Pav, SPECNORM_Pav, &
                              EUCNORM_xPa, FROBNORM_PaPa, SPECNORM_PaPa, &
                              neign, suffix)

   !-------------------------------------------------------------------------
   ! Purpose:   Write norms for Pa and increment approximations
   !
   ! Called from da_compare_decomp_methods
   !
   ! History: 02/22/2017  Creation (JJ Guerrette)
   !
   !-------------------------------------------------------------------------

   implicit none

   real*8, intent(in)        :: EUCNORM_x(3), FROBNORM_Pa(3), SPECNORM_Pa(3)
   real*8, intent(in)        :: EUCNORM_v(3), FROBNORM_Pav(3), SPECNORM_Pav(3)
   real*8, intent(in)        :: EUCNORM_xPa(3), FROBNORM_PaPa(3), SPECNORM_PaPa(3)
!!   real*8, intent(in)        :: EUCNORM_xB(3), FROBNORM_PaB(3), SPECNORM_PaB(3)

   integer, intent(in)       :: neign
   character(*), intent(in)  :: suffix

   real*8                       :: vals(3)
   character(len=filename_len)  :: filename
   integer                      :: itype
   integer                      :: ounit, iostat

   if (trace_use) call da_trace_entry("da_write_norms")

   do itype = 1, 9
      select case (itype)
         case (1)
            vals = EUCNORM_x
            filename = "EUCNORM_x_"//trim(suffix)
         case (2)
            vals = EUCNORM_v
            filename = "EUCNORM_v_"//trim(suffix)
         case (3)
            vals = EUCNORM_xPa
            filename = "EUCNORM_xPa_"//trim(suffix)
         case (4)
            vals = FROBNORM_Pa
            filename = "FROBNORM_Pa_"//trim(suffix)
         case (5)
            vals = FROBNORM_Pav
            filename = "FROBNORM_Pav_"//trim(suffix)
         case (6)
            vals = FROBNORM_PaPa
            filename = "FROBNORM_PaPa_"//trim(suffix)
         case (7)
            vals = SPECNORM_Pa
            filename = "SPECNORM_Pa_"//trim(suffix)
         case (8)
            vals = SPECNORM_Pav
            filename = "SPECNORM_Pav_"//trim(suffix)
         case (9)
            vals = SPECNORM_PaPa
            filename = "SPECNORM_PaPa_"//trim(suffix)
      end select
      if(rootproc) then
         ! =========
         ! open file
         ! =========
         call da_get_unit(ounit)
         if (neign .eq. 1) then
            open(unit   = ounit,     &
               FILE   = trim(filename), &
               FORM   = 'FORMATTED',  &
               ACCESS = 'SEQUENTIAL', &
               iostat =  iostat,     &
               STATUS = 'REPLACE')
         else
            open(unit   = ounit,     &
               FILE   = trim(filename), &
               FORM   = 'FORMATTED',  &
               POSITION = 'APPEND', &
               iostat =  iostat,     &
               STATUS = 'UNKNOWN')

         end if

         if (iostat /= 0) then
            write(unit=message(1),fmt='(A,I5, I5, A)') "Error message: ",iostat, itype, trim(filename)
            call da_warning(__FILE__,__LINE__,message(1:1))
            call da_wrf_error_fatal('Error opening NORM file')
         end if
         write(ounit,fmt='(I6,3E27.16E3)') neign, vals
         close(ounit)
         call da_free_unit(ounit)
      end if
#ifdef DM_PARALLEL
      call mpi_barrier(comm, ierr)
#endif
   end do

   if (trace_use) call da_trace_exit("da_write_norms")

end subroutine da_write_norms
#endif
