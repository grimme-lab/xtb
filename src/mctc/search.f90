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

!> Implements search algorithms
module xtb_mctc_search
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: bisectSearch


   interface bisectSearch
      module procedure :: bisectSearchReal
      module procedure :: bisectSearchInteger
   end interface bisectSearch


contains


!> Integer case for bisection search
pure subroutine bisectSearchInteger(j, xx, x)

   !> Located element such that xx(j) <= x < xx(j+1)
   integer, intent(out) :: j

   !> Array of values in monotonic order to search through
   integer, intent(in) :: xx(:)

   !> Value to locate j for
   integer, intent(in) :: x

   integer :: n
   integer :: jlower, jupper, jcurr

   n = size(xx)
   if (n == 0) then
      j = 0
      return
   end if

   if (x < xx(1)) then
      j = 0
   else if (x == xx(1)) then
      j = 1
   else if (x == xx(n)) then
      j = n - 1
   else if (x > xx(n)) then
      j = n
   else
      jlower=0
      jcurr=n+1
      do while ((jcurr-jlower) > 1)
         jupper=(jcurr+jlower)/2
         if((xx(n) >= xx(1)).eqv.(x >= xx(jupper)))then
            jlower = jupper
         else
            jcurr = jupper
         end if
      end do
      j = jlower
   end if

end subroutine bisectSearchInteger


!> Real case for bisection search
pure subroutine bisectSearchReal(j, xx, x, tol)

   !> Located element such that xx(j) <= x < xx(j+1)
   integer, intent(out) :: j

   !> Array of values in monotonic order to search through
   real(wp), intent(in) :: xx(:)

   !> Value to locate j for
   real(wp), intent(in) :: x

   !> Tolerance for equality comparision
   real(wp), intent(in), optional :: tol

   integer :: n
   integer :: jlower, jupper, jcurr
   real(wp) :: rTol
   logical :: ascending

   n = size(xx)
   if (n == 0) then
      j = 0
      return
   end if

   if (present(tol)) then
      rTol = tol
   else
      rTol = epsilon(0.0_wp)
   end if

   if (x < xx(1) - rTol) then
      j = 0
   else if (abs(x - xx(1)) <= rTol) then
      j = 1
   else if (abs(x - xx(n)) <= rTol) then
      j = n - 1
   else if (x > xx(n) + rTol) then
      j = n
   else
      ascending = (xx(n) >=  xx(1))
      jlower = 0
      jcurr = n+1
      do while ((jcurr-jlower) > 1)
         jupper = (jcurr+jlower)/2
         if (ascending .eqv. (x >= xx(jupper) + rTol)) then
            jlower = jupper
         else
            jcurr = jupper
         end if
      end do
      j = jlower
   end if

end subroutine bisectSearchReal


end module xtb_mctc_search
