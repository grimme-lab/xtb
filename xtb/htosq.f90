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

subroutine htosq(n,a,b)
   implicit none
! ---------------------------------------------------------------------
!      expand trigonal matrix b to full matrix a
!      a and b may refer to the same address
! ---------------------------------------------------------------------
   integer, intent(in)  :: n
   real(wp),intent(out) :: a(n,n)
   real(wp),intent(in)  :: b(n*(n+1)/2)
   integer :: i,j,ioff

   do i=n,1,-1
      ioff=i*(i-1)/2
      do j=i,1,-1
         a(j,i)=b(ioff+j)
      enddo
   enddo

   do i=1,n
      do j=1,i-1
         a(i,j)=a(j,i)
      enddo
   enddo

end subroutine htosq
