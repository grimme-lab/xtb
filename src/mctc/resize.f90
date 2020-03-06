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

!> Reallocation implementation for resizing arrays
module xtb_mctc_resize
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: resize


   !> Overloaded resize interface
   interface resize
      module procedure :: resizeChar
      module procedure :: resizeInt
      module procedure :: resizeReal
      module procedure :: resizeReal2d
   end interface resize


contains


!> Reallocate list of integers
pure subroutine resizeInt(var, n)

   !> TODO
   integer, allocatable, intent(inout) :: var(:)

   !> TODO
   integer, intent(in), optional :: n

   integer, allocatable :: tmp(:)
   integer :: length, currentLength

   currentLength = size(var)
   if (currentLength > 0) then
      if (present(n)) then
         if (n <= currentLength) return
         length = n
      else
         length = currentLength + currentLength/2 + 1
      endif
      allocate(tmp(length), source=0)
      tmp(:currentLength) = var(:currentLength)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length), source=0)
   endif

end subroutine resizeInt


!> Reallocate list of characters
pure subroutine resizeChar(var, n)

   !> TODO
   character(len=*), allocatable, intent(inout) :: var(:)

   !> TODO
   integer, intent(in), optional :: n

   character(len=:), allocatable :: tmp(:)
   integer :: length, currentLength

   currentLength = size(var)
   if (currentLength > 0) then
      if (present(n)) then
         if (n <= currentLength) return
         length = n
      else
         length = currentLength + currentLength/2 + 1
      endif
      allocate(tmp(length), mold=var)
      tmp(:currentLength) = var(:currentLength)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length), source=' ')
   endif

end subroutine resizeChar


!> Reallocate list of reals
pure subroutine resizeReal(var, n)

   !> TODO
   real(wp), allocatable, intent(inout) :: var(:)

   !> TODO
   integer, intent(in), optional :: n

   real(wp), allocatable :: tmp(:)
   integer :: length, order, currentLength

   currentLength = size(var)
   if (currentLength > 0) then
      if (present(n)) then
         if (n <= currentLength) return
         length = n
      else
         length = currentLength + currentLength/2 + 1
      endif
      allocate(tmp(length), source=0.0_wp)
      tmp(:currentLength) = var(:currentLength)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length), source=0.0_wp)
   endif

end subroutine resizeReal


!> Reallocate list of reals
pure subroutine resizeReal2d(var, n)

   !> TODO
   real(wp), allocatable, intent(inout) :: var(:,:)

   !> TODO
   integer, intent(in), optional :: n(2)

   real(wp), allocatable :: tmp(:,:)
   integer :: length, order, currentLength

   currentLength = size(var, 2)
   if (currentLength > 0) then
      if (present(n)) then
         if (n(2) <= currentLength) return
         order = n(1)
         length = n(2)
      else
         length = currentLength + currentLength/2 + 1
         order = size(var, 1)
      endif
      allocate(tmp(order, length), source=0.0_wp)
      tmp(:, :currentLength) = var(:, :currentLength)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         order = n(1)
         length = n(2)
      else
         order = 3
         length = 64
      endif
      allocate(var(order,length), source=0.0_wp)
   endif

end subroutine resizeReal2d


end module xtb_mctc_resize
