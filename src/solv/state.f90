! This file is part of xtb.
!
! Copyright (C) 2020 Sebastian Ehlert
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

!> Define reference states of solution
module xtb_solv_state
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : boltzmann => kB
   implicit none
   private

   public :: solutionState, getStateShift


   !> Possible reference states for the solution
   type :: TSolutionStateEnum

      !> 1 l of ideal gas and 1 l of liquid solution
      integer :: gsolv = 1

      !> 1 bar of ideal gas and 1 mol/L of liquid solution at infinite dilution
      integer :: reference = 2

      !> 1 bar of ideal gas and 1 mol/L of liquid solution
      integer :: mol1bar = 3

   end type TSolutionStateEnum

   !> Actual solvation state enumerator
   type(TSolutionStateEnum), parameter :: solutionState = TSolutionStateEnum()


   real(wp), parameter :: refDensity = 1.0e-3_wp ! kg__au/(1.0e10_wp*AA__Bohr)**3
   real(wp), parameter :: refMolecularMass = 1.0_wp ! amu__au
   real(wp), parameter :: idealGasMolVolume = 24.79_wp
   real(wp), parameter :: ambientTemperature = 298.15_wp ! * Boltzman


contains


pure function getStateShift(state, temperature, density, molecularMass) &
      & result(freeEnergyShift)

   !> Reference state
   integer, intent(in) :: state

   !> Temperature of the solution
   real(wp), intent(in) :: temperature

   !> Mass density of the solvent
   real(wp), intent(in) :: density

   !> Molecular mass of the solvent
   real(wp), intent(in) :: molecularMass

   !> Resulting shift to the solvation free energy
   real(wp) :: freeEnergyShift

   select case(state)
   case default
      freeEnergyShift = 0.0_wp
   case(solutionState%reference)
      freeEnergyShift = temperature * boltzmann &
         & * (log(idealGasMolVolume * temperature / ambientTemperature) &
         & + log(density/refDensity * refMolecularMass/molecularMass))
   case(solutionState%mol1bar)
      freeEnergyShift = temperature * boltzmann &
         & * log(idealGasMolVolume * temperature / ambientTemperature)
   end select

end function getStateShift


end module xtb_solv_state
