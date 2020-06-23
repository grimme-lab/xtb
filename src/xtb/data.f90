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
   use xtb_type_dispersionmodel, only : TDispersionModel
   implicit none
   private

   public :: TxTBData, init
   public :: TRepulsionData, TCoulombData, THamiltonianData, TDispersionData
   public :: THalogenData, TMultipoleData, TShortRangeData
   public :: newData, getData
   public :: generateValenceShellData, angToShellData


   interface newData
      module procedure :: newAtomicData
      module procedure :: newShellData
   end interface newData


   interface getData
      module procedure :: getAtomicData
      module procedure :: getShellData
   end interface getData


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

      !> Reference data for the dispersion
      type(TDispersionModel) :: dispm

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

      !> Name of the parametrisation
      character(len=:), allocatable :: name

      !> Reference to the publication
      character(len=:), allocatable :: doi

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

   contains

      !> Write informative printout for the parametrisation data
      procedure :: writeInfo

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
   real(wp), parameter :: valenceCN(1:118) = [&
      & 1.0_wp, 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, &
      & 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 6.0_wp, 4.0_wp, 4.0_wp, 2.0_wp, 3.0_wp, &
      & 3.0_wp, 3.0_wp, 3.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp]

   !> Cutoff radii for multipole electrostatics
   real(wp), parameter :: multiRad(1:118) = [&
      & 1.4_wp, 3.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.0_wp, 1.9_wp, 1.8_wp, 2.4_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 2.1_wp, 3.1_wp, 2.5_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 3.9_wp, 4.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, 5.0_wp, &
      & 5.0_wp]


contains


subroutine writeInfo(self, unit, num)

   !> Instance of the parametrisation data
   class(TxTBData), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   !> Atomic numbers
   integer, intent(in), optional :: num(:)

   character(len=*), parameter :: rnum = '(f12.6)'
   character(len=*), parameter :: offset = '(8x)'
   character(len=*), parameter :: head = '(6x,"*",1x,a,":")'
   character(len=*), parameter :: rfmt = '('//offset//',a,t36,'//rnum//')'
   character(len=*), parameter :: afmt = '('//offset//',a,t36,4x,a)'
   character(len=:), allocatable :: name
   integer :: ii, jj

   write(unit, '(a)')
   if (allocated(self%name)) then
      allocate(character(len=2*len(self%name)-1) :: name)
      name = repeat(' ', len(name))
      do ii = 1, len(self%name)
         jj = 2*ii-1
         name(jj:jj) = self%name(ii:ii)
      end do
   else
      name = repeat(' ', 10)
      write(name, '(i0)') self%level
      name = 'xTB level '//trim(name)
   end if
   call generic_header(unit, name, 49, 10)

   write(unit, '(a)')
   if (allocated(self%doi)) then
      write(unit, afmt) "Reference", self%doi
   end if

   write(unit, head) "Hamiltonian"
   write(unit, rfmt, advance='no') "H0-scaling (s, p, d)"
   do ii = 0, 2
      write(unit, rnum, advance='no') self%hamiltonian%kScale(ii, ii)
   end do
   write(unit, '(a)')
   write(unit, rfmt) "zeta-weighting", self%hamiltonian%wExp

   write(unit, head) "Dispersion"
   write(unit, rfmt) "s8", self%dispersion%dpar%s8
   write(unit, rfmt) "a1", self%dispersion%dpar%a1
   write(unit, rfmt) "a2", self%dispersion%dpar%a2
   write(unit, rfmt) "s9", self%dispersion%dpar%s9

   write(unit, head) "Repulsion"
   write(unit, rfmt, advance='no') "kExp", self%repulsion%kExp
   if (self%repulsion%kExpLight /= self%repulsion%kExp) then
      write(unit, rnum, advance='no') self%repulsion%kExpLight
   end if
   write(unit, '(a)')
   write(unit, rfmt) "rExp", self%repulsion%rExp

   write(unit, head) "Coulomb"
   write(unit, rfmt) "alpha", self%coulomb%gExp
   if (allocated(self%coulomb%thirdOrderShell)) then
      write(unit, afmt) "third order", "shell-resolved"
   else if (allocated(self%coulomb%thirdOrderAtom)) then
      write(unit, afmt) "third order", "atomic"
   else
      write(unit, afmt) "third order", "false"
   end if

   if (allocated(self%multipole)) then
      write(unit, afmt) "anisotropic", "true"
      write(unit, rfmt) "a3", self%multipole%dipDamp
      write(unit, rfmt) "a5", self%multipole%quadDamp
      write(unit, rfmt) "cn-shift", self%multipole%cnShift
      write(unit, rfmt) "cn-exp", self%multipole%cnExp
      write(unit, rfmt) "max-rad", self%multipole%cnRMax
   else
      write(unit, afmt) "anisotropic", "false"
   end if

   if (allocated(self%halogen)) then
      write(unit, head) "Halogen bond correction"
      write(unit, rfmt) "rad-scale", self%halogen%radScale
      write(unit, rfmt) "damping", self%halogen%dampingPar
   end if

   if (allocated(self%srb)) then
      write(unit, head) "Polar bond correction"
      write(unit, rfmt) "rad-shift", self%srb%shift
      write(unit, rfmt) "strength", self%srb%prefactor
      write(unit, rfmt) "en-exp", self%srb%steepness
      write(unit, rfmt) "en-scale", self%srb%enScale
   end if
   write(unit, '(a)')

end subroutine writeInfo


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


subroutine newAtomicData(vec, num, data)

   real(wp), allocatable, intent(out) :: vec(:)

   integer, intent(in) :: num(:)

   real(wp), intent(in) :: data(:)

   allocate(vec(size(num)))
   call getAtomicData(vec, num, data)

end subroutine newAtomicData


subroutine getAtomicData(vec, num, data)

   real(wp), intent(out) :: vec(:)

   integer, intent(in) :: num(:)

   real(wp), intent(in) :: data(:)

   integer :: ii, izp

   do ii = 1, size(vec, dim=1)
      izp = num(ii)
      vec(ii) = data(izp)
   end do

end subroutine getAtomicData


subroutine newShellData(vec, num, nshell, data)

   real(wp), allocatable, intent(out) :: vec(:, :)

   integer, intent(in) :: num(:)

   integer, intent(in) :: nshell(:)

   real(wp), intent(in) :: data(:, :)

   allocate(vec(maxval(nshell), size(num)))
   call getShellData(vec, num, nshell, data)

end subroutine newShellData


subroutine getShellData(vec, num, nshell, data)

   real(wp), intent(out) :: vec(:, :)

   integer, intent(in) :: num(:)

   integer, intent(in) :: nshell(:)

   real(wp), intent(in) :: data(:, :)

   integer :: ii, ish, izp

   vec(:, :) = 0.0_wp
   do ii = 1, size(vec, dim=2)
      izp = num(ii)
      do ish = 1, nshell(izp)
         vec(ish, ii) = data(ish, izp)
      end do
   end do

end subroutine getShellData


end module xtb_xtb_data
