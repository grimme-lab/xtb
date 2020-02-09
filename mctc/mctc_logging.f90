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

!> Error handling module.
module mctc_logging
   use, intrinsic :: iso_fortran_env
   implicit none
   private

   public :: mctc_environment
   public :: mctc_error

   type :: mctc_environment
      integer :: unit = output_unit
   end type mctc_environment

   type :: mctc_error
      character(len=:), allocatable :: msg
      integer :: code = 1
      logical :: fatal = .true.
   end type mctc_error

contains

end module mctc_logging
