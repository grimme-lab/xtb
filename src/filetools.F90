! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
!
! xtb is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! xtb is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine print_filelist(iunit)
   use xtb_mctc_global, only : persistentEnv
   use xtb_readin
   implicit none
   integer,intent(in) :: iunit

   call persistentEnv%io%list(iunit)

end subroutine print_filelist

subroutine delete_file(name)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   character(len=*),intent(in)  :: name

   call persistentEnv%io%deleteFile(name)

end subroutine delete_file

subroutine touch_file(name)
   use xtb_mctc_global, only : persistentEnv
   use xtb_setparam
   implicit none
   character(len=*),intent(in)  :: name

   call persistentEnv%io%touchFile(name)

end subroutine touch_file

subroutine open_binary(iunit,name,status)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   integer, intent(out) :: iunit
   character(len=*),intent(in)  :: name
   character(len=1),intent(in)  :: status

   select case(status)
   case default
      iunit = -1
   case('r','R')
      call persistentEnv%io%readBinary(iunit, name)
   case('w','W')
      call persistentEnv%io%writeBinary(iunit, name)
   end select

end subroutine open_binary

subroutine open_file(iunit,name,status)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   character(len=*),intent(in)  :: name
   character(len=1),intent(in)  :: status
   integer,intent(out) :: iunit

   select case(status)
   case default
      iunit = -1
   case('r','R')
      call persistentEnv%io%readFile(iunit, name)
   case('w','W')
      call persistentEnv%io%writeFile(iunit, name)
   end select

end subroutine open_file

subroutine close_file(unit)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   integer,intent(in) :: unit

   call persistentEnv%io%closeFile(unit)

end subroutine close_file

subroutine remove_file(unit)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   integer,intent(in) :: unit

   call persistentEnv%io%closeFile(unit, remove=.true.)

end subroutine remove_file
