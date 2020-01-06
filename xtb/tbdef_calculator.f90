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

!> abstract calculator that hides implementation details
!  from calling codes
module tbdef_calculator
   use iso_fortran_env, only: wp => real64
   use tbdef_basisset
   use tbdef_param
   implicit none

   public :: tb_calculator
   private

   type :: tb_calculator
      type(tb_basisset), allocatable :: basis
      type(scc_parameter), allocatable :: param
   end type tb_calculator

end module tbdef_calculator
