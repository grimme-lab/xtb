! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

subroutine fodenmak(uhf,nmo,eps,occ,efermi)
   use mctc_econv, only : autoev
   implicit none
   integer nmo
   logical uhf
   real*8 efermi                      
   real*8 occ(*)
   real*8 eps(*)

   integer i
   real*8 inte

   inte=2.0d0
   if(uhf) inte=1.0d0

   do i=1,nmo
      if(eps(i)*autoev.le.efermi) then
         occ(i)=inte-occ(i)
      endif
   enddo

end subroutine fodenmak
