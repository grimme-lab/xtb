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

!> Implements sorting algorithms
module xtb_mctc_sort
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: indexHeapSort


contains


!> Real case heap sort returning an index.
!  Based on Numerical Recipes Software 1986-92
pure subroutine indexHeapSort(indx, array, tolerance)
   !> Indexing array on return
   integer, intent(out) :: indx(:)
   !> Array of values to be sorted
   real(wp), intent(in) :: array(:)
   !> Tolerance for equality of two elements
   real(wp), intent(in), optional :: tolerance

   integer :: n, ir, ij, il, ii
   integer :: indxTmp
   real(wp) :: arrayTmp, tol

   !:ASSERT(size(array)==size(indx))

   if (present(tolerance)) then
      tol = tolerance
   else
      tol = epsilon(0.0_wp)
   end if

   do ii = 1, size(indx)
      indx(ii) = ii
   end do
   n = size(array)
   if (n <= 1) return
   il = n/2 + 1
   ir = n
   do
      if (il > 1) then
         il = il - 1
         indxTmp = indx(il)
         arrayTmp = array(indxTmp)
      else
         indxTmp = indx(ir)
         arrayTmp = array(indxTmp)
         indx(ir) = indx(1)
         ir = ir - 1
         if (ir < 1) then
            indx(1) = indxTmp
            return
         end if
      end if
      ii = il
      ij = 2 * il
      do while (ij <= ir)
         if (ij < ir) then
            if (array(indx(ij)) < array(indx(ij+1)) - tol) then
               ij = ij + 1
            end if
         end if
         if(arrayTmp < array(indx(ij)) - tol) then
            indx(ii) = indx(ij)
            ii = ij
            ij = 2*ij
         else
            ij = ir + 1
         end if
      end do
      indx(ii)=indxTmp
   end do

end subroutine indexHeapSort


end module xtb_mctc_sort
