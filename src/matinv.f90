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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

subroutine smatinv (a,ldm,n,d)
   implicit none
   integer, intent(in) :: ldm, n
   real*8, intent(out)   :: d
   real*4, intent(inout) :: a(ldm,*)
   integer             :: i, j, k, l(n), m(n)
   real*8                :: biga, temp
   real*8, parameter     :: tol = 1.0d-12

   d = 1.0

   do k = 1, n
      l(k) = k
      m(k) = k
      biga = a(k,k)
      do j = k, n
         do i = k, n
            if ( abs(biga).lt.abs(a(j,i)) ) then
               biga = a(j,i)
               l(k) = i
               m(k) = j
            end if
         end do
      end do
      j = l(k)
      if ( j.gt.k ) then
         do i = 1, n
            temp = -a(i,k)
            a(i,k) = a(i,j)
            a(i,j) = temp
         end do
      end if
      i = m(k)
      if ( i.gt.k ) then
         do j = 1, n
            temp = -a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = temp
         end do
      end if
      if ( abs(biga).lt.tol ) then
         d = 0.0
         return
      end if
      do i = 1, n
         if ( i.ne.k ) a(k,i) = a(k,i)/(-biga)
      end do
      do i = 1, n
         do j = 1, n
            if ( i.ne.k ) then
               if ( j.ne.k ) a(j,i) = a(k,i)*a(j,k) + a(j,i)
            end if
         end do
      end do
      do j = 1, n
         if ( j.ne.k ) a(j,k) = a(j,k)/biga
      end do
      d = max(-1.0d25,min(1.0d25,d))
      d = d*biga
      a(k,k) = 1.0/biga
   end do
   !
   k = n
   do
      !
      k = k - 1
      if ( k.le.0 ) exit
      i = l(k)
      if ( i.gt.k ) then
         do j = 1, n
            temp = a(k,j)
            a(k,j) = -a(i,j)
            a(i,j) = temp
         end do
      end if
      j = m(k)
      if ( j.gt.k ) then
         do i = 1, n
            temp = a(i,k)
            a(i,k) = -a(i,j)
            a(i,j) = temp
         end do
      end if
   end do
   !
end subroutine

subroutine dmatinv (a,ldm,n,d)
   implicit none
   integer, intent(in) :: ldm, n
   real*8, intent(out)   :: d
   real*8, intent(inout) :: a(ldm,*)
   integer             :: i, j, k, l(n), m(n)
   real*8                :: biga, temp
   real*8, parameter     :: tol = 1.0d-12
   !
   d = 1.0
   !
   do k = 1, n
      l(k) = k
      m(k) = k
      biga = a(k,k)
      do j = k, n
         do i = k, n
            if ( abs(biga).lt.abs(a(j,i)) ) then
               biga = a(j,i)
               l(k) = i
               m(k) = j
            end if
         end do
      end do
      j = l(k)
      if ( j.gt.k ) then
         do i = 1, n
            temp = -a(i,k)
            a(i,k) = a(i,j)
            a(i,j) = temp
         end do
      end if
      i = m(k)
      if ( i.gt.k ) then
         do j = 1, n
            temp = -a(k,j)
            a(k,j) = a(i,j)
            a(i,j) = temp
         end do
      end if
      if ( abs(biga).lt.tol ) then
         d = 0.0
         return
      end if
      do i = 1, n
         if ( i.ne.k ) a(k,i) = a(k,i)/(-biga)
      end do
      do i = 1, n
         do j = 1, n
            if ( i.ne.k ) then
               if ( j.ne.k ) a(j,i) = a(k,i)*a(j,k) + a(j,i)
            end if
         end do
      end do
      do j = 1, n
         if ( j.ne.k ) a(j,k) = a(j,k)/biga
      end do
      d = max(-1.0d25,min(1.0d25,d))
      d = d*biga
      a(k,k) = 1.0/biga
   end do
   !
   k = n
   do
      !
      k = k - 1
      if ( k.le.0 ) exit
      i = l(k)
      if ( i.gt.k ) then
         do j = 1, n
            temp = a(k,j)
            a(k,j) = -a(i,j)
            a(i,j) = temp
         end do
      end if
      j = m(k)
      if ( j.gt.k ) then
         do i = 1, n
            temp = a(i,k)
            a(i,k) = -a(i,j)
            a(i,j) = temp
         end do
      end if
   end do
   !
end subroutine
