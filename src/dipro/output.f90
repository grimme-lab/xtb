! This file is part of the dipro module in xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

!> Helper routines to for output generation
module xtb_dipro_output
   implicit none
   private

   public :: format_list, to_string

contains

!> Format a list of logicals as integer ranges
function format_list(list, default, delimiter, skip) result(string)
   !> Generated string representation of the list
   character(len=:), allocatable :: string
   !> Actual list to format as string
   logical, intent(in) :: list(:)
   !> Default value of the list
   logical, intent(in), optional :: default
   !> Delimiter between two ranges
   character(len=*), intent(in), optional :: delimiter
   !> Symbol to represent skipping integers in a range
   character(len=*), intent(in), optional :: skip

   integer :: i, last
   logical :: state, first
   logical :: default_
   character(len=:), allocatable :: delimiter_
   character(len=:), allocatable :: skip_

   if (present(default)) then
      default_ = default
   else
      default_ = .false.
   end if
   if (present(delimiter)) then
      delimiter_ = delimiter
   else
      delimiter_ = ", "
   end if
   if (present(skip)) then
      skip_ = skip
   else
      skip_ = "-"
   end if

   last = -1
   first = .true.
   state = .not.default_
   do i = 1, size(list)
      if (state.eqv.list(i)) then
         state = .not.state
         if (state.eqv.default_) then
            last = i
            if (first) then
               first = .false.
               string = to_string(i)
            else
               string = string // delimiter_ // to_string(i)
            endif
         else
            if (i-1 .ne. last) then
               string = string // skip_ // to_string(i-1)
            endif
         endif
      endif
   enddo
   if (state.eqv.default_ .and.last.ne.size(list)) then
      string = string // skip_ // to_string(size(list))
   endif
end function format_list

!> Format an integer as string without using internal IO
pure function to_string(val) result(string)
   !> Integer value to format
   integer, intent(in) :: val
   !> Resulting string representation
   character(len=:), allocatable :: string

   character(len=1), parameter :: numbers(0:9) = &
      ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
   integer, parameter :: buffer_len = range(val)+2
   character(len=buffer_len) :: buffer
   integer :: pos
   integer :: n

   if (val == 0) then
      string = numbers(0)
      return
   end if

   n = abs(val)
   buffer = ""

   pos = buffer_len + 1
   do while (n > 0)
      pos = pos - 1
      buffer(pos:pos) = numbers(mod(n, 10))
      n = n/10
   end do
   if (val < 0) then
      pos = pos - 1
      buffer(pos:pos) = '-'
   end if

   string = buffer(pos:)
end function to_string

end module xtb_dipro_output
