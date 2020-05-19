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
module xtb_calculators
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment
   implicit none
   private

   public :: gfn0_calculation, gfn1_calculation, gfn2_calculation

   interface
      module subroutine gfn0_calculation &
            (iunit,env,opt,mol,hl_gap,energy,gradient,stress,lattice_gradient)
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(peeq_options),   intent(in)    :: opt
         type(TEnvironment), intent(inout)    :: env
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
         real(wp), intent(out) :: stress(3,3)
         real(wp), intent(out) :: lattice_gradient(3,3)
      end subroutine gfn0_calculation
      module subroutine gfn1_calculation &
            (iunit,env,opt,mol,pcem,wfn,hl_gap,energy,gradient)
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         use xtb_type_pcem
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(scc_options),    intent(in)    :: opt
         type(TEnvironment), intent(inout)    :: env
         type(tb_pcem),        intent(inout) :: pcem
         type(TWavefunction),intent(inout) :: wfn
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
      end subroutine gfn1_calculation
      module subroutine gfn2_calculation &
            (iunit,env,opt,mol,pcem,wfn,hl_gap,energy,gradient)
         use xtb_type_options
         use xtb_type_molecule
         use xtb_type_wavefunction
         use xtb_type_basisset
         use xtb_type_param
         use xtb_type_data
         use xtb_type_pcem
         integer, intent(in) :: iunit
         type(TMolecule),    intent(inout) :: mol
         type(TWavefunction),intent(inout) :: wfn
         type(scc_options),    intent(in)    :: opt
         type(TEnvironment), intent(inout)    :: env
         type(tb_pcem),        intent(inout) :: pcem
         real(wp), intent(out) :: energy
         real(wp), intent(out) :: hl_gap
         real(wp), intent(out) :: gradient(3,mol%n)
      end subroutine gfn2_calculation
   end interface
end module xtb_calculators
