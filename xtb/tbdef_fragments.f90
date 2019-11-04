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

module tbdef_fragments
   implicit none
   public :: tb_fragments
   public :: len, size
   private


   type :: tb_fragments
      integer, allocatable :: list(:)
      integer :: n = 0
   contains
      generic :: allocate => new_default, new_from_list
      procedure, private :: new_default => frag_new_default
      procedure, private :: new_from_list => frag_new_from_list
      procedure :: get_list => frag_get_list
      procedure :: deallocate => frag_destroy
   end type tb_fragments


   interface len
      module procedure :: frag_length
   end interface len


contains


subroutine frag_new_default(self, number_of_atoms)
   class(tb_fragments), intent(out) :: self
   integer, intent(in) :: number_of_atoms
   allocate(self%list(number_of_atoms), source=0)
end subroutine frag_new_default

subroutine frag_new_from_list(self, list)
   class(tb_fragments), intent(out) :: self
   integer, intent(in) :: list(:)
   if (all(list > 0)) then
      self%list = list
      self%n = maxval(list)
   endif
end subroutine frag_new_from_list


integer pure elemental function frag_length(self) result(length)
   class(tb_fragments), intent(in) :: self
   if (allocated(self%list)) then
      length = maxval(self%list)
   else
      length = 0
   endif
end function frag_length


subroutine frag_get_list(self, fragment, list)
   class(tb_fragments), intent(in) :: self
   integer, intent(in) :: fragment
   integer, allocatable, intent(out) :: list(:)
   integer :: i
   list = pack([(i, i=1, size(self%list))], mask=self%list.eq.fragment)
end subroutine frag_get_list


subroutine frag_destroy(self)
   class(tb_fragments), intent(inout) :: self
   if (allocated(self%list)) deallocate(self%list)
end subroutine frag_destroy


end module tbdef_fragments
