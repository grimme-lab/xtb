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

!> Version information for this project
module xtb_dipro_version
   implicit none
   private

   public :: dipro_version_string, dipro_version_compact
   public :: get_dipro_version


   !> String representation of the dipro version
   character(len=*), parameter :: dipro_version_string = "0.1.0"

   !> Numeric representation of the dipro version
   integer, parameter :: dipro_version_compact(3) = [0, 1, 0]

contains

!> Getter function to retrieve dipro version
subroutine get_dipro_version(major, minor, patch, string)
   !> Major version number of the dipro version
   integer, intent(out), optional :: major
   !> Minor version number of the dipro version
   integer, intent(out), optional :: minor
   !> Patch version number of the dipro version
   integer, intent(out), optional :: patch
   !> String representation of the dipro version
   character(len=:), allocatable, intent(out), optional :: string

   if (present(major)) then
      major = dipro_version_compact(1)
   end if
   if (present(minor)) then
      minor = dipro_version_compact(2)
   end if
   if (present(patch)) then
      patch = dipro_version_compact(3)
   end if
   if (present(string)) then
      string = dipro_version_string
   end if

end subroutine get_dipro_version

end module xtb_dipro_version
