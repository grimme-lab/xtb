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

!> Available boundary conditions in this program
module xtb_mctc_boundaryconditions
   implicit none
   private

   public :: boundaryCondition


   !> Possible boundary conditions
   type :: TBoundaryCondEnum

      !> Finite molecular structure
      integer :: cluster = 0

      !> Three dimensional infinite periodic boundary conditions
      integer :: pbc3d = 3

   end type 

   !> Boundary condition enumerator
   type(TBoundaryCondEnum), parameter :: boundaryCondition = TBoundaryCondEnum()


contains


end module xtb_mctc_boundaryconditions
