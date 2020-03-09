! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Manipulation of character variables
module xtb_mctc_chartools
   implicit none
   private

   public :: toLowercase


   !> ASCII offset between lowercase and uppercase letters
   integer, parameter :: offset = iachar('A') - iachar('a')


contains


!> Convert input string to lowercase
function toLowercase(str) result(lcStr)

   !> Input string
   character(len=*), intent(in) :: str

   !> Lowercase version of string
   character(len=len_trim(str)):: lcStr

   integer :: ilen, iquote, i, iav, iqc

   ilen = len_trim(str)
   iquote = 0
   lcstr = str

   do i = 1, ilen
      iav = iachar(str(i:i))
      if (iquote == 0 .and. (iav == 34 .or.iav == 39)) then
         iquote = 1
         iqc = iav
         cycle
      end if
      if (iquote == 1 .and. iav==iqc) then
         iquote=0
         cycle
      end if
      if (iquote == 1) cycle
      if (iav >= iachar('A') .and. iav <= iachar('Z')) then
         lcstr(i:i) = achar(iav - offset)
      else
         lcstr(i:i) = str(i:i)
      end if
   end do

end function toLowercase


end module xtb_mctc_chartools
