! This file is part of dipro.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> Helper routines to for output generation
module dipro_output
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

end module dipro_output
