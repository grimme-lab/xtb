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

!> Helper routines to handle calculations with the tblite library
module dipro_xtb
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
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
      error stop
   case("gfn2")
      call new_gfn2_calculator(xcalc, mol)
   case("gfn1")
      call new_gfn1_calculator(xcalc, mol)
   case("ipea1")
      call new_ipea1_calculator(xcalc, mol)
   end select
end subroutine get_calculator

end module dipro_xtb
