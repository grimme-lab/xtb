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

!> implemenation of a list of atoms
module xtb_type_atomlist
   use xtb_mctc_accuracy, only : wp
   implicit none
   public :: TAtomList
   public :: size, len, assignment(=)
   private

   character, parameter :: p_delimiter = ','
   character, parameter :: p_skip = '-'

   type :: TAtomList
      logical, allocatable :: list(:)
      logical :: default = .false.
      character :: delimiter = p_delimiter
      character :: skip = p_skip
      logical :: error = .false.
   contains
      generic :: new => from_integers, from_logicals, from_string, from_defaults
      ! Should be private, but must be declared as public due to PGI bug #28452
      procedure, non_overridable :: from_defaults
      procedure, private :: from_integers => atomlist_assign_integers
      procedure, private :: from_logicals => atomlist_assign_logicals
      procedure, private :: from_string => atomlist_assign_string
      procedure, pass(self) :: to_string => string_assign_atomlist
      procedure, pass(self) :: to_list => list_assign_atomlist
      procedure :: switch_truth => atomlist_switch_truth
      procedure :: set_truth => atomlist_set_truth
      procedure :: get_truth => atomlist_get_truth
      procedure :: get_error => atomlist_get_error
      generic :: remove => remove_integer, remove_integers, remove_logicals
      procedure, private :: remove_integer => atomlist_remove_integer
      procedure, private :: remove_integers => atomlist_remove_integers
      procedure, private :: remove_logicals => atomlist_remove_logicals
      generic :: add => add_integer, add_integers, add_logicals, add_string
      procedure, private :: add_integer => atomlist_add_integer
      procedure, private :: add_integers => atomlist_add_integers
      procedure, private :: add_logicals => atomlist_add_logicals
      procedure, private :: add_string => atomlist_assign_string
      procedure, private :: parse => atomlist_parse_string
      procedure :: resize => atomlist_resize
      generic :: gather => gather_int, gather_real, gather_real_2d
      procedure, private :: gather_int => atomlist_gather_int
      procedure, private :: gather_real => atomlist_gather_real
      procedure, private :: gather_real_2d => atomlist_gather_real_2d
      generic :: scatter => scatter_real_2d
      procedure, private :: scatter_real_2d => atomlist_scatter_real_2d
      procedure :: destroy => atomlist_destroy
      final :: atomlist_finalizer
   end type TAtomList

   interface TAtomList
      module procedure :: atomlist_from_logicals
      module procedure :: atomlist_from_integers
      module procedure :: atomlist_from_string
   end interface TAtomList

   interface len
      module procedure :: atomlist_length
   end interface len

   interface size
      module procedure :: atomlist_size
   end interface size

   interface assignment(=)
      module procedure :: atomlist_assign_logicals
      module procedure :: atomlist_assign_integers
      module procedure :: atomlist_assign_string
      module procedure :: string_assign_atomlist
      module procedure :: list_assign_atomlist
   end interface assignment(=)

contains

pure function atomlist_from_logicals(list, truth, delimiter, skip) result(self)
   logical, intent(in) :: list(:)
   logical, intent(in), optional :: truth
   character, intent(in), optional :: delimiter, skip
   type(TAtomList) :: self
   if (present(truth)) self%default = .not.truth
   if (present(delimiter)) self%delimiter = delimiter
   if (present(skip)) self%skip = skip
   call self%new(list)
end function atomlist_from_logicals

pure function atomlist_from_integers(list, truth, delimiter, skip) result(self)
   integer, intent(in) :: list(:)
   logical, intent(in), optional :: truth
   character, intent(in), optional :: delimiter, skip
   type(TAtomList) :: self
   if (present(truth)) self%default = .not.truth
   if (present(delimiter)) self%delimiter = delimiter
   if (present(skip)) self%skip = skip
   call self%new(list)
end function atomlist_from_integers

pure function atomlist_from_string(list, truth, delimiter, skip) result(self)
   character(len=*), intent(in) :: list
   logical, intent(in), optional :: truth
   character, intent(in), optional :: delimiter, skip
   type(TAtomList) :: self
   if (present(truth)) self%default = .not.truth
   if (present(delimiter)) self%delimiter = delimiter
   if (present(skip)) self%skip = skip
   call self%new(list)
end function atomlist_from_string

subroutine from_defaults(self)
   class(TAtomList), intent(out) :: self
end subroutine from_defaults

pure elemental subroutine atomlist_switch_truth(self)
   class(TAtomList), intent(inout) :: self
   self%default = .not.self%default
end subroutine atomlist_switch_truth

pure elemental subroutine atomlist_set_truth(self, truth)
   class(TAtomList), intent(inout) :: self
   logical, intent(in) :: truth
   self%default = .not.truth
end subroutine atomlist_set_truth

pure elemental function atomlist_get_truth(self) result(truth)
   class(TAtomList), intent(in) :: self
   logical :: truth
   truth = .not.self%default
end function atomlist_get_truth

pure elemental function atomlist_get_error(self) result(error)
   class(TAtomList), intent(in) :: self
   logical :: error
   error = self%error
end function atomlist_get_error

subroutine atomlist_write_formatted(self, unit, iotype, v_list, iostat, iomsg)
   class(TAtomList), intent(in) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in) :: iotype
   integer, intent(in) :: v_list(:)
   integer, intent(out) :: iostat
   character(len=*), intent(inout) :: iomsg
   character(len=:), allocatable :: buffer
   call self%to_string(buffer)
   write(unit, '(a)', iostat=iostat, iomsg=iomsg) buffer
end subroutine atomlist_write_formatted

subroutine atomlist_read_formatted(self, unit, iotype, v_list, iostat, iomsg)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: unit
   character(len=*), intent(in) :: iotype
   integer, intent(in) :: v_list(:)
   integer, intent(out) :: iostat
   character(len=*), intent(inout) :: iomsg
   character(len=:), allocatable :: string

   call get_line(unit, string, iostat)
   call self%new(string)

contains

subroutine get_line(unit, line, iostat)
   use, intrinsic :: iso_fortran_env, only : iostat_eor
   integer, intent(in) :: unit
   character(len=:), allocatable, intent(out) :: line
   integer, intent(out), optional :: iostat

   integer, parameter :: buffersize = 256
   character(len=buffersize) :: buffer
   integer :: size
   integer :: err

   line = ''
   do
      read(unit,'(a)',advance='no',iostat=err,size=size)  &
      &    buffer
      if (err.gt.0) then
         if (present(iostat)) iostat=err
         return ! an error occurred
      endif
      line = line // buffer(:size)
      if (err.lt.0) then
         if (err.eq.iostat_eor) err = 0
         if (present(iostat)) iostat=err
         return
      endif
   enddo

end subroutine get_line
end subroutine atomlist_read_formatted

integer pure elemental function atomlist_size(self)
   class(TAtomList), intent(in) :: self
   if (allocated(self%list)) then
      atomlist_size = size(self%list)
   else
      atomlist_size = 0
   endif
end function atomlist_size

integer pure elemental function atomlist_length(self)
   class(TAtomList), intent(in) :: self
   if (allocated(self%list)) then
      atomlist_length = count(self%list.neqv.self%default)
   else
      atomlist_length = 0
   endif
end function atomlist_length

pure subroutine atomlist_add_integer(self, item)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: item
   integer :: i
   call self%resize(item)
   self%list(item) = .not.self%default
end subroutine atomlist_add_integer

pure subroutine atomlist_add_integers(self, list)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: list(:)
   integer :: i
   call self%resize(maxval(list))
   do i = 1, size(list)
      self%list(list(i)) = .not.self%default
   enddo
end subroutine atomlist_add_integers

pure subroutine atomlist_add_logicals(self, list)
   class(TAtomList), intent(inout) :: self
   logical, intent(in) :: list(:)
   integer :: i
   call self%resize(size(list))
   self%list(:size(list)) = list .or. self%list(:size(list))
end subroutine atomlist_add_logicals

pure subroutine atomlist_remove_integer(self, item)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: item
   integer :: i
   call self%resize(item)
   self%list(item) = self%default
end subroutine atomlist_remove_integer

pure subroutine atomlist_remove_integers(self, list)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: list(:)
   integer :: i
   call self%resize(maxval(list))
   do i = 1, size(list)
      self%list(list(i)) = self%default
   enddo
end subroutine atomlist_remove_integers

pure subroutine atomlist_remove_logicals(self, list)
   class(TAtomList), intent(inout) :: self
   logical, intent(in) :: list(:)
   integer :: i
   call self%resize(size(list))
   self%list(:size(list)) = list .and. self%list(:size(list))
end subroutine atomlist_remove_logicals

pure subroutine list_assign_atomlist(list, self)
   integer, allocatable, intent(out) :: list(:)
   class(TAtomList), intent(in) :: self
   integer :: i, j
   allocate(list(len(self)), source=0)
   j = 0
   do i = 1, size(self)
      if (self%list(i).neqv.self%default) then
         j = j+1
         list(j) = i
      endif
   enddo
end subroutine list_assign_atomlist

pure subroutine string_assign_atomlist(string, self)
   character(len=:), allocatable, intent(out) :: string
   class(TAtomList), intent(in) :: self
   character(len=10) :: buffer
   integer :: i, last
   logical :: state, first
   last = -1
   first = .true.
   state = .not.self%default
   do i = 1, size(self)
      if (state.eqv.self%list(i)) then
         state = .not.state
         if (state.eqv.self%default) then
            last = i
            write(buffer,'(i0)') i
            if (first) then
               first = .false.
               string = trim(buffer)
            else
               string = string // self%delimiter // trim(buffer)
            endif
         else
            if (i-1 .ne. last) then
               write(buffer,'(i0)') i-1
               string = string // self%skip // trim(buffer)
            endif
         endif
      endif
   enddo
   if (state.eqv.self%default .and.last.ne.size(self)) then
      write(buffer,'(i0)') size(self)
      string = string // self%skip // trim(buffer)
   endif
end subroutine string_assign_atomlist

pure subroutine atomlist_assign_logicals(self, list)
   class(TAtomList), intent(inout) :: self
   logical, intent(in) :: list(:)
   call self%resize(size(list))
   self%list = list
end subroutine atomlist_assign_logicals

pure subroutine atomlist_assign_integers(self, list)
   class(TAtomList), intent(inout) :: self
   integer, intent(in) :: list(:)
   integer :: i
   call self%resize(maxval(list))
   do i = 1, size(list)
      self%list(list(i)) = .not.self%default
   enddo
end subroutine atomlist_assign_integers

pure subroutine atomlist_assign_string(self, string)
   class(TAtomList), intent(inout) :: self
   character(len=*), intent(in) :: string
   character(len=:), allocatable :: buffer
   integer, allocatable :: list(:)
   integer :: pos, last, n

   last = 0
   do
      pos = index(string(last+1:), self%delimiter)
      if (pos > 0) then
         call self%parse(string(last+1:last+pos-1))
         if (self%error) exit
         last = last+pos
      else
         call self%parse(string(last+1:))
         exit
      endif
   enddo

end subroutine atomlist_assign_string

pure subroutine atomlist_parse_string(self, string)
   class(TAtomList), intent(inout) :: self
   character(len=*),intent(in) :: string
   integer :: pos, item, begin, last, err

   pos = index(string,self%skip)
   if (pos.eq.0) then
      read(string,*,iostat=err) item
      if (err.ne.0 .or.item.eq.0) then
         self%error = .true.
         return
      endif
      call self%add(item)
   else
      read(string(:pos-1),*,iostat=err) begin
      if (err.ne.0 .or. begin.eq.0) then
         self%error = .true.
         return
      endif
      read(string(pos+1:),*,iostat=err) last
      if (err.ne.0 .or. last.eq.0) then
         self%error = .true.
         return
      endif
      if (last.lt.begin) then
         self%error = .true.
         return
      endif
      do item = begin, last
         call self%add(item)
      enddo
   endif

end subroutine atomlist_parse_string

pure subroutine atomlist_resize(self, n)
   class(TAtomList), intent(inout) :: self
   integer, intent(in), optional :: n
   logical, allocatable :: list(:)
   integer :: length, current_length
   current_length = size(self)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(list(length), source=self%default)
      list(:current_length) = self%list(:current_length)
      deallocate(self%list)
      call move_alloc(list, self%list)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(self%list(length), source=self%default)
   endif
end subroutine atomlist_resize

subroutine atomlist_gather_int(self, array, selected)
   class(TAtomList), intent(in) :: self
   integer, intent(in) :: array(:)
   integer, allocatable, intent(out) :: selected(:)
   selected = pack(array, mask=self%list.neqv.self%default)
end subroutine atomlist_gather_int

subroutine atomlist_gather_real(self, array, selected)
   class(TAtomList), intent(in) :: self
   real(wp), intent(in) :: array(:)
   real(wp), allocatable, intent(out) :: selected(:)
   selected = pack(array, mask=self%list.neqv.self%default)
end subroutine atomlist_gather_real

subroutine atomlist_gather_real_2d(self, array, selected, dim)
   class(TAtomList), intent(in) :: self
   real(wp), intent(in) :: array(:,:)
   integer, intent(in), optional :: dim
   real(wp), allocatable, intent(out) :: selected(:,:)
   integer :: selected_dim, i, j
   if (present(dim)) then
      selected_dim = dim
   else
      selected_dim = 2
   endif
   if (selected_dim == 1) then
      allocate(selected(len(self), size(array, 2)), source=0.0_wp)
      j = 0
      do i = 1, min(size(array, 1), size(self))
         if (self%list(i).neqv.self%default) then
            j = j+1
            selected(j,:) = array(i,:)
         endif
      enddo
   else
      allocate(selected(size(array, 1), len(self)), source=0.0_wp)
      j = 0
      do i = 1, min(size(array, 2), size(self))
         if (self%list(i).neqv.self%default) then
            j = j+1
            selected(:,j) = array(:,i)
         endif
      enddo
   endif
end subroutine atomlist_gather_real_2d

subroutine atomlist_scatter_real_2d(self, reduced, array, dim, scale)
   class(TAtomList), intent(in) :: self
   real(wp), intent(inout) :: array(:,:)
   integer, intent(in), optional :: dim
   real(wp), intent(in), optional :: scale
   real(wp), intent(in) :: reduced(:,:)
   integer :: selected_dim, i, j
   real(wp) :: alpha
   if (present(dim)) then
      selected_dim = dim
   else
      selected_dim = 2
   endif
   if (present(scale)) then
      alpha = scale
   else
      alpha = 0.0_wp
   endif
   if (selected_dim == 1) then
      j = 0
      do i = 1, min(size(array, 1), size(self))
         if (self%list(i).neqv.self%default) then
            j = j+1
            array(i,:) = reduced(j,:) + alpha * array(i,:)
         endif
      enddo
   else
      j = 0
      do i = 1, min(size(array, 2), size(self))
         if (self%list(i).neqv.self%default) then
            j = j+1
            array(:,i) = reduced(:,j) + alpha * array(:,i)
         endif
      enddo
   endif
end subroutine atomlist_scatter_real_2d

pure elemental subroutine atomlist_destroy(self)
   class(TAtomList), intent(inout) :: self
   if (allocated(self%list)) deallocate(self%list)
end subroutine atomlist_destroy

pure elemental subroutine atomlist_finalizer(self)
   type(TAtomList), intent(inout) :: self
   call self%destroy
end subroutine atomlist_finalizer

end module xtb_type_atomlist
