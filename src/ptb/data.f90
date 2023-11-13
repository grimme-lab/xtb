! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
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

!> Parametrisation data for the PTB method
module xtb_ptb_data
   use xtb_mctc_accuracy, only: wp
   use xtb_param_atomicrad, only: atomicRad
   use xtb_param_paulingen, only: paulingEN
   use xtb_type_param, only: dftd_parameter
   use xtb_type_dispersionmodel, only: TDispersionModel
   implicit none
   private

   public :: TPTBData, init
   public :: TRepulsionData, TCoulombData, THamiltonianData, TDispersionData
   public :: TShortRangeData, TCorePotentialData, TEEQData
   public :: newData, getData
   public :: angToShellData

   interface newData
      module procedure :: newAtomicData
      module procedure :: newShellData
      module procedure :: newAngShellData
   end interface newData

   interface getData
      module procedure :: getAtomicData
      module procedure :: getShellData
      module procedure :: getAngShellData
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

   !> Data for the effective core potential
   type :: TCorePotentialData
      !> Maximum number of core shells
      integer :: max_shell

      !> Maximum number of core primitives
      integer :: max_prim
      
      !> Number of core shells for each atom
      integer, allocatable :: nshell(:)

      !> Principal quantum number of each shell
      integer, allocatable :: pqn(:,:)

      !> Angular momentum of each shell
      integer, allocatable :: angshell(:,:)
      
      !> HF level of each shell
      real(wp), allocatable :: hflev(:,:)

      !> Slater exponents of each shell
      real(wp), allocatable :: sl_exp(:,:)

      !> Effective core potential scaling factors
      real(wp), allocatable :: kecpepsilon(:)

   end type TCorePotentialData

   !> EEQ parameters
   type :: TEEQData

      !> Alpha parameter for EEQ model
      real(wp), allocatable :: alp(:)

      !> EN parameter for EEQ model
      real(wp), allocatable :: chi(:)

      !> CN dependency parameter of EEQ model
      real(wp), allocatable :: cnf(:)

      !> Chemical hardness parameters for EEQ model
      real(wp), allocatable :: gam(:)

   end type TEEQData

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

   !> Data for the evaluation of the PTB core Hamiltonian
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

      !> Number of primitives for expansion of Slater functions
      integer, allocatable :: numberOfPrimitives(:, :)

      !> Reference occupation of the atom
      real(wp), allocatable :: refocc(:, :)

      !> Electronegativity used in the shell polynomials
      real(wp), allocatable :: electronegativity(:)

      !> Atomic radii used in the shell polynomials
      real(wp), allocatable :: atomicRad(:)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQShell(:, :)

      !> Charge dependence of the atomic levels
      real(wp), allocatable :: kQAtom(:)

      !----------------
      !> H0 information
      !----------------
      !> Atomic level information
      real(wp), allocatable :: selfEnergy(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), allocatable :: klh(:, :)
      !> Coordination number-star dependence of the atomic levels
      real(wp), allocatable :: kcnstar(:)
      !> Shift of atomic levels depending only on coordination number-star
      real(wp), allocatable :: kshift(:)
      !> Wolfsberg parameter for EHT
      real(wp), allocatable :: kla(:,:)
      !> Atomic radii for H0
      real(wp), allocatable :: kr(:)
      !> One-center off-diagonal scaling of S_H0 (squared)
      real(wp), allocatable :: kocod(:)
      !> OCOD scaling fo s-s', p-p', d-d'...
      real(wp), allocatable :: ksla(:,:)

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

   !> Parametrisation data for the PTB method
   type :: TPTBData

      !> Name of the parametrisation
      character(len=:), allocatable :: name

      !> Reference to the publication
      character(len=:), allocatable :: doi

      !> Internal version number
      integer :: level

      !> Number of shells
      integer, allocatable :: nshell(:)

      !> Parametrisation data for repulsive interactions
      type(TRepulsionData) :: repulsion

      !> Parametrisation data for core Hamiltonian
      type(THamiltonianData) :: hamiltonian

      !> Parametrisation data for Coulombic interactions
      type(TCoulombData) :: coulomb

      !> Parametrisation data for dispersion interactions
      type(TDispersionData) :: dispersion

      !> Parametrisation data for the effective core potential
      type(TCorePotentialData) :: corepotential

      !> Parametrisation data for the EEQ model
      type(TEEQData) :: eeq

      !> Parametrisation data for the short range basis correction (optional)
      type(TShortRangeData), allocatable :: srb


      !> Shift for IP/EA calculations
      real(wp) :: ipeashift

   contains

      !> Write informative printout for the parametrisation data
      procedure :: writeInfo

   end type TPTBData

   !> Default constructor for the data types
   interface init
      module procedure :: initRepulsion
      module procedure :: initCoulomb
      module procedure :: initCorepotential
      module procedure :: initEEQ
      module procedure :: initHamiltonian
   end interface init

contains

   subroutine writeInfo(self, unit, num)

      !> Instance of the parametrisation data
      class(TPTBData), intent(in) :: self

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

      write (unit, '(a)')
      if (allocated(self%name)) then
         allocate (character(len=2*len(self%name) - 1) :: name)
         name = repeat(' ', len(name))
         do ii = 1, len(self%name)
            jj = 2*ii - 1
            name(jj:jj) = self%name(ii:ii)
         end do
      else
         name = repeat(' ', 10)
         write (name, '(i0)') self%level
         name = 'PTB level '//trim(name)
      end if
      call generic_header(unit, name, 49, 10)

      write (unit, '(a)')
      if (allocated(self%doi)) then
         write (unit, afmt) "Reference", self%doi
      end if

      write (unit, head) "Hamiltonian"
      write (unit, rfmt, advance='no') "H0-scaling (s, p, d)"
      do ii = 0, 2
         write (unit, rnum, advance='no') self%hamiltonian%kScale(ii, ii)
      end do
      write (unit, '(a)')
      write (unit, rfmt) "zeta-weighting", self%hamiltonian%wExp

      write (unit, head) "Dispersion"
      write (unit, rfmt) "s8", self%dispersion%dpar%s8
      write (unit, rfmt) "a1", self%dispersion%dpar%a1
      write (unit, rfmt) "a2", self%dispersion%dpar%a2
      write (unit, rfmt) "s9", self%dispersion%dpar%s9

      write (unit, head) "Repulsion"
      write (unit, rfmt, advance='no') "kExp", self%repulsion%kExp
      if (self%repulsion%kExpLight /= self%repulsion%kExp) then
         write (unit, rnum, advance='no') self%repulsion%kExpLight
      end if
      write (unit, '(a)')
      write (unit, rfmt) "rExp", self%repulsion%rExp

      write (unit, head) "Coulomb"
      write (unit, rfmt) "alpha", self%coulomb%gExp
      if (allocated(self%coulomb%thirdOrderShell)) then
         write (unit, afmt) "third order", "shell-resolved"
      else if (allocated(self%coulomb%thirdOrderAtom)) then
         write (unit, afmt) "third order", "atomic"
      else
         write (unit, afmt) "third order", "false"
      end if

      if (allocated(self%srb)) then
         write (unit, head) "Polar bond correction"
         write (unit, rfmt) "rad-shift", self%srb%shift
         write (unit, rfmt) "strength", self%srb%prefactor
         write (unit, rfmt) "en-exp", self%srb%steepness
         write (unit, rfmt) "en-scale", self%srb%enScale
      end if
      write (unit, '(a)')

   end subroutine writeInfo

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

   subroutine initCorepotential(self, max_core_prim, max_core_shell, & !> Array sizes
      & cbas_sl_exp, cbas_nshell, & !> Slater exponents, number of shells
      & cbas_pqn, cbas_angshell, cbas_hflev, & !> Principal quantum number, angular momentum, HF level
      & kecpepsilon) !> Effective core potential scaling factors

      integer, intent(in) :: max_core_prim, max_core_shell
      integer, intent(in) :: cbas_nshell(:), cbas_pqn(:,:), cbas_angshell(:,:)
      real(wp), intent(in) :: cbas_sl_exp(:,:), cbas_hflev(:, :), kecpepsilon(:)

      !> Data instance
      type(TCorePotentialData), intent(out) :: self

      self%max_prim = max_core_prim
      self%max_shell = max_core_shell
      self%nshell = cbas_nshell(:)
      self%pqn = cbas_pqn(:, :)
      self%angshell = cbas_angshell(:, :)
      self%sl_exp = cbas_sl_exp(:, :)
      self%hflev = cbas_hflev(:, :)
      self%kecpepsilon = kecpepsilon(:)

   end subroutine initCorepotential

   subroutine initEEQ(self, num, alp, chi, cnf, gam)

      !> Data instance
      type(TEEQData), intent(out) :: self

      !> Atomic numbers for unique species
      integer, intent(in) :: num(:)

      !> Alpha parameter for EEQ model
      real(wp), intent(in) :: alp(:)

      !> EN parameter for EEQ model
      real(wp), intent(in) :: chi(:)

      !> CN dependency parameter of EEQ model
      real(wp), intent(in) :: cnf(:)

      !> Chemical hardness parameters for EEQ model
      real(wp), intent(in) :: gam(:)

      call newData(self%alp, num, alp)
      call newData(self%chi, num, chi)
      call newData(self%cnf, num, cnf)
      call newData(self%gam, num, gam)

   end subroutine initEEQ

   subroutine initHamiltonian(self, num, nshell, ang_shell, h0_levels, &
      & cn_dependency, cnstar_dependency, cnstar_shift, wolfsberg_par, &
      & atom_radii_h0, onecenteroffdiagonal, ocod_l)

      !> Data instance
      type(THamiltonianData), intent(out) :: self
      !> Atomic numbers for unique species
      integer, intent(in) :: num(:)
      !> Number of shells for each atom
      integer, intent(in) :: nshell(:)
      !> Angular momenta of each shell
      integer, intent(in) :: ang_shell(:, :)
      !> Atomic level information
      real(wp), intent(in) :: h0_levels(:, :)
      !> Coordination number dependence of the atomic levels
      real(wp), intent(in) :: cn_dependency(:,:)
      !> Coordination number-star dependence of the atomic levels
      real(wp), intent(in) :: cnstar_dependency(:)
      !> Shift of atomic levels depending only on coordination number-star
      real(wp), intent(in) :: cnstar_shift(:)
      !> Wolfsberg parameter for EHT
      real(wp), intent(in) :: wolfsberg_par(:, :)
      !> Atomic radii for H0
      real(wp), intent(in) :: atom_radii_h0(:)
      !> One-center off-diagonal scaling of S_H0 (squared)
      real(wp), intent(in) :: onecenteroffdiagonal(:)
      !> Ang-mom-dependent OCOD scaling for s-s', p-p', d-d'...
      real(wp), intent(in) :: ocod_l(:,:)

      call newData(self%selfEnergy, num, nshell, h0_levels)
      call newData(self%klh, num, nshell, cn_dependency)
      call newData(self%kcnstar, num, cnstar_dependency)
      call newData(self%kshift, num, cnstar_shift)
      call newData(self%kla, num, nshell, ang_shell, wolfsberg_par)
      call newData(self%kr, num, atom_radii_h0)
      call newData(self%kocod, num, onecenteroffdiagonal)
      call newData(self%ksla, num, nshell, ang_shell, ocod_l)

   end subroutine initHamiltonian

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
      if (present(electronegativity) .and. present(kCN) .and. present(chargeWidth)) then
         maxElem = min(maxElem, size(electronegativity), size(kCN), size(chargeWidth))
      end if

      self%chemicalHardness = chemicalHardness(:maxElem)
      if (present(shellHardness)) then
         self%shellHardness = shellHardness(:, :maxElem)
      end if
      if (present(thirdOrderAtom)) then
         self%thirdOrderAtom = thirdOrderAtom(:maxElem)
      end if
      if (present(electronegativity) .and. present(kCN) .and. present(chargeWidth)) then
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

      integer :: nElem, iZp, iSh, lAng

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

      allocate (vec(size(num)), source=0.0_wp)
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

      allocate (vec(maxval(nshell), size(num)))
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

   subroutine newAngShellData(vec, num, nshell, angmompershell, data)

      real(wp), allocatable, intent(out) :: vec(:, :)

      integer, intent(in) :: num(:)

      integer, intent(in) :: nshell(:)

      integer, intent(in) :: angmompershell(:,:)

      real(wp), intent(in) :: data(:, :)

      allocate (vec(maxval(nshell), size(num)))
      call getAngShellData(vec, num, nshell, angmompershell, data)

   end subroutine newAngShellData

   subroutine getAngShellData(vec, num, nshell, angmompershell, data)

      real(wp), intent(out) :: vec(:, :)

      integer, intent(in) :: num(:)

      integer, intent(in) :: nshell(:)

      integer, intent(in) :: angmompershell(:,:)

      real(wp), intent(in) :: data(:, :)

      integer :: ii, ish, izp, angmom

      vec(:, :) = 0.0_wp
      do ii = 1, size(vec, dim=2)
         izp = num(ii)
         do ish = 1, nshell(izp)
            angmom = angmompershell(ish, izp)
            vec(ish, ii) = data(angmom+1, izp)
         end do
      end do

   end subroutine getAngShellData

end module xtb_ptb_data
