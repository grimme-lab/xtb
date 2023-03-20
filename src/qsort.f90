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

recursive subroutine qsort(a, first, last, ind)
  use xtb_mctc_accuracy, only : wp
  implicit none
  real(wp) :: a(*), x, t
  integer  :: ind(*)
  integer  :: first, last
  integer  :: i, j, ii

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     ii=ind(i); ind(i) = ind(j);  ind(j) = ii
     i=i+1
     j=j-1
  end do
  if (first < i-1) call qsort(a, first, i-1, ind)
  if (j+1 < last)  call qsort(a, j+1, last, ind)
end subroutine qsort

recursive subroutine qqsort(a, first, last)
   use xtb_mctc_accuracy, only : wp
  implicit none
  real(wp) :: a(*), x, t
  integer  :: first, last
  integer  :: i, j, ii

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     i=i+1
     j=j-1
  end do
  if (first < i-1) call qqsort(a, first, i-1)
  if (j+1 < last)  call qqsort(a, j+1, last)
end subroutine qqsort
