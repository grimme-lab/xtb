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

!> Version information for this project
module dipro_version
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

end module dipro_version
