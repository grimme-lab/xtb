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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Helper routines to handle calculations with the tblite library
module xtb_dipro_xtb
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
#if WITH_TBLITE
   use tblite_xtb_calculator, only : xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator
   use tblite_xtb_gfn1, only : new_gfn1_calculator
   use tblite_xtb_ipea1, only : new_ipea1_calculator
   implicit none
   private

   public :: get_calculator

contains

!> Create a new calculator for a given method string
subroutine get_calculator(xcalc, mol, method, error)
   !> Instance of the new calculator
   type(xtb_calculator), intent(out) :: xcalc
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Name of the requested method
   character(len=*), intent(in) :: method
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select case(method)
   case default
      call fatal_error(error, "Unknown method '"//method//"' requested")
!      error stop
   case("gfn2")
      call new_gfn2_calculator(xcalc, mol)
   case("gfn1")
      call new_gfn1_calculator(xcalc, mol)
   case("ipea1")
      call new_ipea1_calculator(xcalc, mol)
   end select
end subroutine get_calculator
#endif

end module xtb_dipro_xtb
