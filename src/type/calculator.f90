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
module xtb_type_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_type_basisset
   use xtb_type_param
   use xtb_xtb_data, only : TxTBData
   implicit none

   public :: tb_calculator
   private

   type :: tb_calculator
      type(TBasisset), allocatable :: basis
      type(scc_parameter), allocatable :: param
      type(TxTBData), allocatable :: xtb
   end type tb_calculator

end module xtb_type_calculator
