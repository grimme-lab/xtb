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

pure subroutine blowsy(ity,a,b,n)
   use iso_fortran_env, wp => real64
   implicit none
   ! blow up symmetric or antisymmetric matrix to full size
   integer, intent(in)  :: ity
   integer, intent(in)  :: n
   real(wp),intent(in)  :: a(*)
   real(wp),intent(out) :: b(n,n)
   integer :: ij,i,j
   ! determine if we have an antisymmetric integral

   ij=0
   if (ity.eq.-1) then
      do i=1,n
         do j=1,i-1
            ij=ij+1
            b(j,i)=-a(ij)
            b(i,j)=a(ij)
         enddo
         ij=ij+1
         b(i,i)=0.d0
      enddo
   else
      do i=1,n
         do j=1,i-1
            ij=ij+1
            b(j,i)=a(ij)
            b(i,j)=a(ij)
         enddo
         ij=ij+1
         b(i,i)=a(ij)
      enddo
   endif
end subroutine blowsy

