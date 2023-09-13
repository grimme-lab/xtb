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

!> Gathers the input data needed to create solvation models
module xtb_solv_input
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_state, only : solutionState
   implicit none
   private

   public :: TSolvInput


   !> Input data for the solvation model
   type :: TSolvInput

      !> Selected solvent
      character(len=:), allocatable :: solvent

      !> Use the analytical linearized Poisson Boltzmann model
      logical :: alpb = .true.

      !> Generalized Born interaction kernel
      integer :: kernel = gbKernel%p16

      !> Use COSMO solvation model instead of Born based one
      logical :: cosmo = .false.
      
      !> Use TM convention for COSMO solvation model
      logical :: tmcosmo = .false.

      !> Additional CPCM-X mode
      character(len=:),allocatable :: cpxsolvent

      !> Reference state
      integer :: state = solutionState%gsolv

      !> Temperature
      real(wp) :: temperature = 298.15_wp

      !> Number of grid points
      integer :: nAng = 230

      !> Ion strength
      real(wp) :: ionStrength = 0.0_wp

      !> Ion radius
      real(wp) :: ionRad = 0.0_wp

      !> Onufriev-Bashford-Case correction
      real(wp) :: obcPar(3) = [1.0_wp, 0.8_wp, 4.85_wp]

      !> Dielectric smoothing parameter
      real(wp) :: smoothPar = 0.3_wp * aatoau

      !> Born radii integration cutoff
      real(wp) :: bornCutoff = 35.0_wp * aatoau

      !> SASA integration offset
      real(wp) :: surfaceOffset = 2.0_wp * aatoau

      !> Threshold for surface elements
      real(wp) :: surfaceThr = 1.0e-6_wp

   end type TSolvInput


#ifdef __PGIC__
   !> Create a custom constructor due to PGI 20.7 failing for derived type
   !  constructors using deferred-length character variables
   interface TSolvInput
      module procedure :: pgiWrapper
   end interface TSolvInput
#endif


contains


function pgiWrapper(solvent, alpb, kernel, state, temperature, nAng, &
      & ionStrength, ionRad) result(self)

   !> Selected solvent
   character(len=*), intent(in) :: solvent

   !> Use the analytical linearized Poisson Boltzmann model
   logical, intent(in), optional :: alpb

   !> Generalized Born interaction kernel
   integer, intent(in), optional :: kernel

   !> Reference state
   integer, intent(in), optional :: state

   !> Temperature
   real(wp), intent(in), optional :: temperature

   !> Number of grid points
   integer, intent(in), optional :: nAng

   !> Ion strength
   real(wp), intent(in), optional :: ionStrength

   !> Ion radius
   real(wp), intent(in), optional :: ionRad

   !> Instance of the solvation input data
   type(TSolvInput) :: self

   self%solvent = solvent
   if (present(alpb)) self%alpb = alpb
   if (present(kernel)) self%kernel = kernel
   if (present(state)) self%state = state
   if (present(temperature)) self%temperature = temperature
   if (present(nAng)) self%nAng = nAng
   if (present(ionStrength)) self%ionStrength = ionStrength
   if (present(ionRad)) self%ionRad = ionRad

end function pgiWrapper


end module xtb_solv_input
