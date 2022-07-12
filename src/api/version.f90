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

!> Defines the API version for xtb
module xtb_api_version
   use, intrinsic :: iso_c_binding
   implicit none
   private

   public :: getAPIVersion_api

   integer(c_int), parameter :: apiMajor = 1
   integer(c_int), parameter :: apiMinor = 0
   integer(c_int), parameter :: apiPatch = 0


contains


function getAPIVersion_api() result(version) &
      & bind(C, name="xtb_getAPIVersion")
   !DEC$ ATTRIBUTES DLLEXPORT :: getAPIVersion_api
   integer(c_int) :: version

   version = 10000_c_int * apiMajor + 100_c_int * apiMinor + apiPatch

end function getAPIVersion_api


end module xtb_api_version
