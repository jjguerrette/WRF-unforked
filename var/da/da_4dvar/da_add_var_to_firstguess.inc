subroutine da_add_var_to_firstguess(overwrite_name) 
!---------------------------------------------------------------------------
!  Purpose: Rewrite fg file including variables initialized in WRFDA 
!
!  Creator :   JJ Guerrette, Dec 10, 2017
!         
!  The following WRF fields are copied/created:
!  All variables in input stream 0
!---------------------------------------------------------------------------

   implicit none

   INTERFACE
     integer(c_int32_t) function copyfile(ifile, ofile) bind(c)
       import :: c_int32_t, C_CHAR
       CHARACTER(KIND=C_CHAR), DIMENSION(*), intent(in) :: ifile, ofile
     END function copyfile
   END INTERFACE

   INTERFACE
     integer(c_int32_t) function deletefile(dfile) bind(c)
       import :: c_int32_t, C_CHAR
       CHARACTER(KIND=C_CHAR), DIMENSION(*), intent(in) :: dfile
     END function deletefile
   END INTERFACE

   character(*), intent(in), optional    :: overwrite_name               ! I/O filename

   character(len=120) :: file_name 
   character(len=4) ::  fgname
   integer :: ierr

   if ( present(overwrite_name) ) then
      fgname = trim(overwrite_name)
   else
      fgname = 'fg'
   end if
   file_name = trim(fgname)//'_temp'

   config_flags%input_outname = trim(file_name)
   config_flags%write_input = .true.

   call med_filter_out (head_grid, config_flags)

   if (rootproc) then
      ierr = deletefile(trim(fgname)//C_NULL_CHAR)

      ierr = copyfile(trim(file_name)//C_NULL_CHAR, trim(fgname)//C_NULL_CHAR)

      ierr = ierr + deletefile(trim(file_name)//C_NULL_CHAR)

      if ( ierr /= 0 ) then
         write(unit=message(1),fmt='(a)') "Failed to create "//trim(file_name)//" from "//trim(fgname)
         call da_error(__FILE__,__LINE__,message(1:1))
      endif
   end if

   config_flags%input_outname = "wrfinput_d<domain>"
   config_flags%write_input = .false.

end subroutine da_add_var_to_firstguess 
