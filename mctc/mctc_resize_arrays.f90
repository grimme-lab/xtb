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

module mctc_resize_arrays
   use iso_fortran_env, only: wp => real64
   implicit none
   public :: resize
   private

   interface resize
      module procedure :: resize_char
      module procedure :: resize_int
      module procedure :: resize_real
      module procedure :: resize_real_2d
   end interface resize

contains

subroutine resize_int(var, n)
   integer, allocatable, intent(inout) :: var(:)
   integer, intent(in), optional :: n
   integer, allocatable :: tmp(:)
   integer :: length, current_length
   current_length = size(var)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(tmp(length), source=0)
      tmp(:current_length) = var(:current_length)
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
end subroutine resize_int

subroutine resize_char(var, n)
   character(len=*), allocatable, intent(inout) :: var(:)
   integer, intent(in), optional :: n
   character(len=:), allocatable :: tmp(:)
   integer :: length, current_length
   current_length = size(var)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(tmp(length), mold=var)
      tmp(:current_length) = var(:current_length)
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
end subroutine resize_char

subroutine resize_real(var, n)
   real(wp), allocatable, intent(inout) :: var(:)
   integer, intent(in), optional :: n
   real(wp), allocatable :: tmp(:)
   integer :: length, order, current_length
   current_length = size(var)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(tmp(length), source=0.0_wp)
      tmp(:current_length) = var(:current_length)
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
end subroutine resize_real

subroutine resize_real_2d(var, n)
   real(wp), allocatable, intent(inout) :: var(:,:)
   integer, intent(in), optional :: n(2)
   real(wp), allocatable :: tmp(:,:)
   integer :: length, order, current_length
   current_length = size(var, 2)
   if (current_length > 0) then
      if (present(n)) then
         if (n(2) <= current_length) return
         order = n(1)
         length = n(2)
      else
         length = current_length + current_length/2 + 1
         order = size(var, 1)
      endif
      allocate(tmp(order, length), source=0.0_wp)
      tmp(:, :current_length) = var(:, :current_length)
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
end subroutine resize_real_2d

end module mctc_resize_arrays
