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

!> Parametrisation data for xTB methods
module xtb_xtb_data
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   use xtb_type_param, only : dftd_parameter
   implicit none
   private

   public :: TxTBData, init
   public :: TRepulsionData, TCoulombData, THamiltonianData, TDispersionData
   public :: THalogenData, TMultipoleData, TShortRangeData
   public :: generateValenceShellData, angToShellData


   !> Data for the dispersion contribution
   type :: TDispersionData

      !> Damping parameters
      type(dftd_parameter) :: dpar

      !> Weighting factor for Gaussian interpolation
      real(wp) :: wf

      !> Charge steepness
      real(wp) :: g_a

      !> Charge height
      real(wp) :: g_c

   end type TDispersionData


   !> Data for the repulsion contribution
   type :: TRepulsionData

      !> Repulsion exponent for heavy elements
      real(wp) :: kExp

      !> Repulsion exponent for light elements
      real(wp) :: kExpLight

      !> Repulsion exponent
      real(wp) :: rExp

      !> Electronegativity scaling of repulsion
      real(wp) :: enScale

      !> Exponents of repulsion term
      real(wp), allocatable :: alpha(:)

      !> Effective nuclear charge
      real(wp), allocatable :: zeff(:)

      !> Electronegativitity for scaling of repulsion
      real(wp), allocatable :: electronegativity(:)

      !> FIXME: real space cutoff should not be part of data
      real(wp) :: cutoff

   end type TRepulsionData


   !> Data for the evaluation of the xTB core Hamiltonian
   type :: THamiltonianData

      !> Scaling factors for different interacting shells
      real(wp) :: kScale(0:3, 0:3)

      !> Scaling factor for diffuse or polarisation function
      real(wp) :: kDiff

      !> Shell dependence of the EN polynom
      real(wp) :: enScale(0:3, 0:3)

      !> Quartic contribution to EN polynom
      real(wp) :: enscale4

      !> Exponent for shell exponent weighting
      real(wp) :: wExp

      !> Principal quantum number of each shell
      integer, allocatable :: principalQuantumNumber(:, :)

      !> Angular momentum of each shell
      integer, allocatable :: angShell(:, :)

      !> Valence character of each shell
      integer, allocatable :: valenceShell(:, :)

      !> Number of primitives for expansion of Slater functions
      integer, allocatable :: numberOfPrimitives(:, :)

      !> Exponent of the Slater function
      real(wp), allocatable :: slaterExponent(:, :)

      !> Atomic level information
      real(wp), allocatable :: selfEnergy(:, :)

      !> Reference occupation of the atom
      real(wp), allocatable :: referenceOcc(:, :)

      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: kCN(:, :)

      !> Electronegativity used in the shell polynomials
      real(wp), allocatable :: electronegativity(:)

      !> Atomic radii used in the shell polynomials
      real(wp), allocatable :: atomicRad(:)

      !> Shell polynomials to scale Hamiltonian elements
      real(wp), allocatable :: shellPoly(:, :)

      !> Pair parameters to scale Hamiltonian elements
      real(wp), allocatable :: pairParam(:, :)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQShell(:, :)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQAtom(:)

   end type THamiltonianData


   !> Data for the evalutation of the Coulomb interactions
   type :: TCoulombData

      !> Use electronegativity equlibration for reference density
      logical :: enEquilibration

      !> Include second order electrostatics
      logical :: secondOrder

      !> Include third order electrostatics
      logical :: thirdOrder

      !> Third order electrostatics is shell resolved
      logical :: shellResolved

      !> Exponent of the generalized gamma function
      real(wp) :: gExp

      !> Atomic hardnesses used in second order electrostatics
      real(wp), allocatable :: chemicalHardness(:)

      !> Scaling factors for shell electrostatics
      real(wp), allocatable :: shellHardness(:, :)

      !> Third order Hubbard derivatives
      real(wp), allocatable :: thirdOrderAtom(:)

      !> Shell resolved third order Hubbard derivatives
      real(wp), allocatable :: thirdOrderShell(:, :)

      !> Charge widths for EEQ model
      real(wp), allocatable :: chargeWidth(:)

      !> Electronegativity for EEQ model
      real(wp), allocatable :: electronegativity(:)

      !> Coordination number dependence of the EN
      real(wp), allocatable :: kCN(:)

   end type TCoulombData


   !> Data for the evaluation of the multipole electrostatics
   type :: TMultipoleData

      !> Coordination number shift
      real(wp) :: cnShift

      !> Coordination number exponent for radii
      real(wp) :: cnExp

      !> Maximum radius
      real(wp) :: cnRMax

      !> Damping parameter for charge-dipole interactions
      real(wp) :: dipDamp

      !> Damping parameter for dipole-dipole, charge-quadrupole interactions
      real(wp) :: quadDamp

      !> Valence coordination number for radii
      real(wp), allocatable :: valenceCN(:)

      !> Cutoff radii for multipole electrostatics
      real(wp), allocatable :: multiRad(:)

      !> Dipole exchange-correlation kernel
      real(wp), allocatable :: dipKernel(:)

      !> Quadrupole exchange-correlation kernel
      real(wp), allocatable :: quadKernel(:)

   end type TMultipoleData


   !> Data for halogen bond correction
   type :: THalogenData

      !> Scaling factor of the atomic radii
      real(wp) :: radScale

      !> Damping parameter for the halogen bond interactions
      real(wp) :: dampingPar

      !> Strength of the halogen bond
      real(wp), allocatable :: bondStrength(:)

      !> Atomic radii
      real(wp), allocatable :: atomicRad(:)

   end type THalogenData


   !> Short range basis correction
   type TShortRangeData

      !> Additional offset for the reference bond lengths
      real(wp) :: shift

      !> Scaling factor for the energy contribution
      real(wp) :: prefactor

      !> Steepness of the EN dependence
      real(wp) :: steepness

      !> Scaling factor for electronegativity differences
      real(wp) :: enScale

   end type TShortRangeData


   !> Parametrisation data for the xTB method
   type :: TxTBData

      !> Internal version number
      integer :: level

      !> Number of shells
      integer, allocatable :: nShell(:)

      !> Parametrisation data for repulsive interactions
      type(TRepulsionData) :: repulsion

      !> Parametrisation data for core Hamiltonian
      type(THamiltonianData) :: hamiltonian

      !> Parametrisation data for Coulombic interactions
      type(TCoulombData) :: coulomb

      !> Parametrisation data for dispersion interactions
      type(TDispersionData) :: dispersion

      !> Parametrisation data for multipole electrostatics (optional)
      type(TMultipoleData), allocatable :: multipole

      !> Parametrisation data for halogen bond correction (optional)
      type(THalogenData), allocatable :: halogen

      !> Parametrisation data for the short range basis correction (optional)
      type(TShortRangeData), allocatable :: srb

      !> Shift for IP/EA calculations
      real(wp) :: ipeashift

   end type TxTBData


   !> Default constructor for the data types
   interface init
      module procedure :: initRepulsion
      module procedure :: initHalogen
      module procedure :: initMultipole
      module procedure :: initCoulomb
   end interface init


   ! ========================================================================
   ! MULTIPOLE DATA
   !> Valence coordination number for radii
   real(wp), parameter :: valenceCN(1:86) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, &
      & 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp]

   !> Cutoff radii for multipole electrostatics
   real(wp), parameter :: multiRad(1:86) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 4.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp]


contains


!> Generator for valence shell data from the angular momenta of the shells
subroutine generateValenceShellData(valenceShell, nShell, angShell)

   !> Valency character of each shell
   integer, intent(out) :: valenceShell(:, :)

   !> Number of shells for each atom
   integer, intent(in) :: nShell(:)

   !> Angular momenta of each shell
   integer, intent(in) :: angShell(:, :)

   integer :: lAng, iZp, iSh
   logical :: valShell(0:3)

   valenceShell(:, :) = 0
   do iZp = 1, size(nShell, dim=1)
      valShell(:) = .true.
      do iSh = 1, nShell(iZp)
         lAng = angShell(iSh, iZp)
         if (valShell(lAng)) then
            valShell(lAng) = .false.
            valenceShell(iSh, iZp) = 1
         end if
      end do
   end do

end subroutine generateValenceShellData


!> Initialize halogen bond data
subroutine initHalogen(self, radScale, dampingPar, halogenBond)

   !> Data instance
   type(THalogenData), intent(out) :: self

   !> Scaling parameter for the atomic radii
   real(wp), intent(in) :: radScale

   !> Damping parameter
   real(wp), intent(in) :: dampingPar

   !> Halogen bond strength
   real(wp), intent(in) :: halogenBond(:)

   integer :: maxElem

   maxElem = size(halogenBond)

   self%radScale = radScale
   self%dampingPar = dampingPar
   self%atomicRad = atomicRad(:maxElem)
   self%bondStrength = halogenBond(:maxElem)

end subroutine initHalogen


subroutine initMultipole(self, cnShift, cnExp, cnRMax, dipDamp, quadDamp, &
      & dipKernel, quadKernel)

   !> Data instance
   type(TMultipoleData), intent(out) :: self

   !>
   real(wp), intent(in) :: cnShift

   !>
   real(wp), intent(in) :: cnExp

   !>
   real(wp), intent(in) :: cnRMax

   !>
   real(wp), intent(in) :: dipDamp

   !>
   real(wp), intent(in) :: quadDamp

   !>
   real(wp), intent(in) :: dipKernel(:)

   !>
   real(wp), intent(in) :: quadKernel(:)

   integer :: maxElem

   maxElem = min(size(dipKernel), size(quadKernel))

   self%cnShift = cnShift
   self%cnExp = cnExp
   self%cnRMax = cnRMax
   self%dipDamp = dipDamp
   self%quadDamp = quadDamp
   self%dipKernel = dipKernel(:maxElem)
   self%quadKernel = quadKernel(:maxElem)
   self%valenceCN = valenceCN(:maxElem)
   self%multiRad = multiRad(:maxElem)

end subroutine initMultipole


subroutine initRepulsion(self, kExp, kExpLight, rExp, enScale, alpha, zeff, &
      & electronegativity)

   !> Data instance
   type(TRepulsionData), intent(out) :: self

   !>
   real(wp), intent(in) :: kExp

   !>
   real(wp), intent(in) :: kExpLight

   !>
   real(wp), intent(in) :: rExp

   !>
   real(wp), intent(in) :: enScale

   !>
   real(wp), intent(in) :: alpha(:)

   !>
   real(wp), intent(in) :: zeff(:)

   !>
   real(wp), intent(in), optional :: electronegativity(:)

   integer :: maxElem

   maxElem = min(size(alpha), size(zeff))
   if (present(electronegativity)) then
      maxElem = min(maxElem, size(electronegativity))
   end if

   self%cutoff = 40.0_wp
   self%kExp = kExp
   self%kExpLight = kExpLight
   self%rExp = rExp
   self%enScale = enScale
   self%alpha = alpha(:maxElem)
   self%zeff = zeff(:maxElem)
   if (present(electronegativity)) then
      self%electronegativity = electronegativity(:maxElem)
   else
      self%electronegativity = paulingEN(:maxElem)
   end if

end subroutine initRepulsion


subroutine initCoulomb(self, nShell, chemicalHardness, shellHardness, &
      & thirdOrderAtom, electronegativity, kCN, chargeWidth)

   !> Data instance
   type(TCoulombData), intent(out) :: self

   !>
   integer, intent(in) :: nShell(:)

   !>
   real(wp), intent(in) :: chemicalHardness(:)

   !>
   real(wp), intent(in), optional :: shellHardness(:, :)

   !>
   real(wp), intent(in), optional :: thirdOrderAtom(:)

   !>
   real(wp), intent(in), optional :: electronegativity(:)

   !>
   real(wp), intent(in), optional :: kCN(:)

   !>
   real(wp), intent(in), optional :: chargeWidth(:)

   integer :: maxElem

   maxElem = size(chemicalHardness)
   if (present(shellHardness)) then
      maxElem = min(maxElem, size(shellHardness, dim=2))
   end if
   if (present(thirdOrderAtom)) then
      maxElem = min(maxElem, size(thirdOrderAtom))
   end if
   if (present(electronegativity).and.present(kCN).and.present(chargeWidth)) then
      maxElem = min(maxElem, size(electronegativity), size(kCN), size(chargeWidth))
   end if

   self%chemicalHardness = chemicalHardness(:maxElem)
   if (present(shellHardness)) then
      self%shellHardness = shellHardness(:, :maxElem)
   end if
   if (present(thirdOrderAtom)) then
      self%thirdOrderAtom = thirdOrderAtom(:maxElem)
   end if
   if (present(electronegativity).and.present(kCN).and.present(chargeWidth)) then
      self%electronegativity = electronegativity(:maxElem)
      self%kCN = kCN(:maxElem)
      self%chargeWidth = chargeWidth(:maxElem)
   end if

end subroutine initCoulomb


!> Transform a data array from angular momenta to shell number references
subroutine angToShellData(kDat, nShell, angShell, angDat)

   !> Data in terms of shell number of each species
   real(wp), intent(out) :: kDat(:, :)

   !> Number of shells for each species
   integer, intent(in) :: nShell(:)

   !> Angular momenta of each shell
   integer, intent(in) :: angShell(:, :)

   !> Data in terms of angular momenta of each shell
   real(wp), intent(in) :: angDat(0:, :)

   integer :: nElem, iZp, iSh, lAng, iKind

   nElem = min(size(kDat, dim=2), size(nShell), size(angShell, dim=2), &
      & size(angDat, dim=2))

   kDat(:, :) = 0.0_wp
   do iZp = 1, nElem
      do iSh = 1, nShell(iZp)
         lAng = angShell(iSh, iZp)
         kDat(iSh, iZp) = angDat(lAng, iZp)
      end do
   end do

end subroutine angToShellData


end module xtb_xtb_data
