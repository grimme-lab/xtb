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


contains


!> real case for bisection search
pure subroutine bisectSearch(j, xx, x, tol)

   !> located element such that xx(j) < x < xx(j+1)
   integer, intent(out) :: j

   !> array of values in monotonic order to search through
   real(wp), intent(in) :: xx(:)

   !> value to locate j for
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

end subroutine bisectSearch


end module xtb_mctc_search
