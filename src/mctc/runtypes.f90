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

!> Supported run types in this program
module xtb_mctc_runtypes
   implicit none
   private

   public :: runType


   !> Possible run types
   type :: TRunTypeEnum

      !> Single point calculation
      integer :: sp = 1

      !> Geometry optimization
      integer :: opt = 2

      !> Frequency calculation
      integer :: freq = 3

      !> Molecular dynamics
      integer :: md = 4

   end type TRunTypeEnum

   !> Run type enumerator
   type(TRunTypeEnum), parameter :: runType = TRunTypeEnum()


contains


end module xtb_mctc_runtypes
