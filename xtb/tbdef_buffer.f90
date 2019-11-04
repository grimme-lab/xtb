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

module tbdef_buffer
   use iso_fortran_env
   implicit none
   public :: tb_buffer
   public :: len, size
   private


   type :: tb_buffer
      private
      integer :: error_code = 0
      integer :: length = 0
      integer :: pos = 0
      character(len=:), allocatable :: raw
      integer, allocatable :: table(:,:)
   contains
      generic :: allocate => new_default
      procedure, private :: new_default => buffer_new_default
      procedure :: get_error => buffer_get_error
      procedure :: push_back => buffer_push_back
      procedure :: reset => buffer_reset
      procedure :: next => buffer_next
      procedure :: getline => buffer_getline
      procedure :: resize => buffer_resize
      procedure :: deallocate => buffer_destroy
   end type tb_buffer


   interface len
      module procedure :: buffer_length
   end interface

   interface size
      module procedure :: buffer_size
   end interface


contains


subroutine buffer_new_default(self)
   class(tb_buffer), intent(out) :: self
end subroutine buffer_new_default


integer pure elemental function buffer_size(self)
   class(tb_buffer), intent(in) :: self
   if (allocated(self%table)) then
      buffer_size = size(self%table, 2)
   else
      buffer_size = 0
   endif
end function buffer_size

integer pure elemental function buffer_length(self)
   class(tb_buffer), intent(in) :: self
   if (allocated(self%table)) then
      buffer_length = self%length
   else
      buffer_length = 0
   endif
end function buffer_length

pure subroutine buffer_resize(self, n)
   class(tb_buffer), intent(inout) :: self
   integer, intent(in), optional :: n
   integer, allocatable :: table(:,:)
   integer :: length, current_length
   current_length = size(self)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(table(2, length), source=0)
      table(:, :current_length) = self%table(:, :current_length)
      deallocate(self%table)
      call move_alloc(table, self%table)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(self%table(2, length), source=0)
   endif
end subroutine buffer_resize


subroutine buffer_push_back(self, string)
   class(tb_buffer), intent(inout) :: self
   character(len=*), intent(in) :: string
   integer :: pos, last, length

   pos = len(self%raw)
   length = len(string)
   last = len(self)
   if (last >= size(self)) call self%resize

   ! increment counter
   self%length = self%length+1
   ! save position of string in buffer
   self%table(:,last+1) = pos+[1,length]
   ! append string to buffer
   self%raw = self%raw // string
end subroutine buffer_push_back


subroutine buffer_reset(self)
   class(tb_buffer), intent(inout) :: self
   self%pos = 0
end subroutine buffer_reset

logical function buffer_next(self) result(next)
   class(tb_buffer), intent(inout) :: self
   next = self%pos < self%length
   if (next) self%pos = self%pos+1
end function buffer_next

subroutine buffer_getline(self, line)
   class(tb_buffer), intent(in) :: self
   character(len=:), allocatable, intent(out) :: line
   line = self%raw(self%table(1,self%pos):self%table(2,self%pos))
end subroutine buffer_getline


logical pure elemental function buffer_get_error(self) result(error)
   class(tb_buffer), intent(in) :: self
   error = self%error_code /= 0
end function buffer_get_error

subroutine buffer_destroy(self)
   class(tb_buffer), intent(inout) :: self
   if (allocated(self%raw)) deallocate(self%raw)
   if (allocated(self%table)) deallocate(self%table)
end subroutine buffer_destroy


end module tbdef_buffer
