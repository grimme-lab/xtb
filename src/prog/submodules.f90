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

!> Available submodule in this program
module xtb_prog_submodules
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: xtbSubmodule, getSubmodule


   !> Possible submodules for xtb
   type :: TSubmoduleEnum

      !> Unknown run mode
      integer :: invalid = 0

      !> Main run mode of xtb
      integer :: main = 1

      !> Info run mode of xtb
      integer :: info = 2

      !> Thermodynamic function evaluator
      integer :: thermo = 3

      !> Force field topology generator
      integer :: topo = 4

      !> Fragment docking
      integer :: dock = 5

      !> IR spectra from DFTB+ output
      integer :: ir = 6

   end type TSubmoduleEnum

   !> Actual enumerator for the submodules
   type(TSubmoduleEnum), parameter :: xtbSubmodule = TSubmoduleEnum()


contains


!> Translate a string in a submodule enumerator
function getSubmodule(argument) result(submod)

   !> String identifying the submodule name
   character(len=*), intent(in) :: argument

   !> Identifier of the submodule
   integer :: submod

   select case(argument)
   case default
      submod = xtbSubmodule%invalid

   case('info')
      submod = xtbSubmodule%info

   case('thermo')
      submod = xtbSubmodule%thermo

   case('topo')
      submod = xtbSubmodule%topo

   case('dock')
      submod = xtbSubmodule%dock

   case('ir')
      submod = xtbSubmodule%ir

   end select

end function getSubmodule


end module xtb_prog_submodules
