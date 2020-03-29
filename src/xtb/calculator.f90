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

!> Implementation of the calculator interface of the xTB methods
module xtb_xtb_calculator
   use xtb_mctc_accuracy, only : wp
   use xtb_type_identitymap, only : TIdentityMap
   use xtb_type_molecule, only : TMolecule
   use xtb_type_calculator, only : TCalculator
   use xtb_xtb_data
   implicit none
   private

   public :: TxTBCalculator, init


   !> Container for performing calculations with xTB methods
   type, extends(TCalculator) :: TxTBCalculator
   end type TxTBCalculator


   !> Initialize the xTB calculator instance from a database
   interface init
      module procedure :: initCalculator
   end interface


contains


!> Initialize the xTB calculator instance from a database
subroutine initCalculator(self, input)

   !> Instance of the calculator
   type(TxTBCalculator), intent(out) :: self

   !> Parameter data base
   type(TxTBData), intent(in) :: input

end subroutine initCalculator


end module xtb_xtb_calculator
