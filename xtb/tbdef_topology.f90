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

module tbdef_topology
   implicit none
   public :: tb_topology
   public :: len, size, assignment(=)
   private


   type :: tb_topology
      private
      integer :: order = 2
      integer :: length = 0
      integer, allocatable :: list(:,:)
   contains
      generic :: allocate => new_default, new_from_list
      procedure, private :: new_default => top_new_default
      procedure, private :: new_from_list => top_new_from_list
      procedure, pass(self) :: to_list => list_assign_top
      procedure :: push_back => top_push_back
      procedure :: set_item => top_set_item
      procedure :: get_item => top_get_item
      procedure :: resize => top_resize
      procedure :: deallocate => top_destroy
   end type tb_topology


   interface len
      module procedure :: top_length
   end interface len

   interface size
      module procedure :: top_size
   end interface size

   interface assignment(=)
      module procedure :: top_new_from_list
      module procedure :: list_assign_top
   end interface assignment(=)


contains


subroutine top_new_default(self, order, size)
   class(tb_topology), intent(out) :: self
   integer, intent(in), optional :: order
   integer, intent(in), optional :: size
   if (present(order)) self%order = order
   if (present(size)) call self%resize(size)
end subroutine top_new_default

subroutine top_new_from_list(self, list)
   class(tb_topology), intent(out) :: self
   integer, intent(in) :: list(:,:)
   self%order = size(list, 1)
   call self%resize(size(list, 2))
   self%list = list
end subroutine top_new_from_list

subroutine list_assign_top(list, self)
   class(tb_topology), intent(in) :: self
   integer, allocatable, intent(out) :: list(:,:)
   list = self%list(:,:len(self))
end subroutine list_assign_top


integer pure elemental function top_length(self) result(length)
   class(tb_topology), intent(in) :: self
   if (allocated(self%list)) then
      length = self%length
   else
      length = 0
   endif
end function top_length

integer pure elemental function top_size(self) result(length)
   class(tb_topology), intent(in) :: self
   if (allocated(self%list)) then
      length = size(self%list, 2)
   else
      length = 0
   endif
end function top_size


subroutine top_get_item(self, at, item)
   class(tb_topology), intent(in) :: self
   integer, intent(in) :: at
   integer, intent(out) :: item(:)
   if (at > len(self)) item = 0
   item = self%list(:, at)
end subroutine top_get_item

subroutine top_set_item(self, at, item)
   class(tb_topology), intent(inout) :: self
   integer, intent(in) :: at
   integer, intent(in) :: item(:)
   if (at > len(self)) call self%resize(at)
   self%list(:, at) = item
end subroutine top_set_item


subroutine top_push_back(self, item)
   class(tb_topology), intent(inout) :: self
   integer, intent(in) :: item(:)
   integer :: pos, last, length

   last = len(self)
   if (last >= size(self)) call self%resize

   self%length = last+1
   self%list(:,last+1) = item
end subroutine top_push_back

subroutine top_resize(self, n)
   class(tb_topology), intent(inout) :: self
   integer, intent(in), optional :: n
   integer, allocatable :: list(:,:)
   integer :: length, current_length
   current_length = size(self)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(list(self%order, length), source=0)
      list(:, :current_length) = self%list(:, :current_length)
      deallocate(self%list)
      call move_alloc(list, self%list)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(self%list(self%order, length), source=0)
   endif
end subroutine top_resize

subroutine top_destroy(self)
   class(tb_topology), intent(inout) :: self
   if (allocated(self%list)) deallocate(self%list)
end subroutine top_destroy

end module tbdef_topology
