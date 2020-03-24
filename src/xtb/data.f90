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

!> TODO
module xtb_xtb_data
   use xtb_mctc_accuracy, only : wp
   use xtb_param_atomicrad, only : atomicRad
   use xtb_param_paulingen, only : paulingEN
   implicit none
   private

   public :: TxTBData, init
   public :: TRepulsionData, TCoulombData, THamiltonianData
   public :: THalogenData, TMultipoleData
   public :: generateValenceShellData


   !>
   type :: TRepulsionData

      !>
      real(wp) :: kExp

      !>
      real(wp) :: kExpLight

      !>
      real(wp) :: rExp

      !>
      real(wp) :: enScale

      !>
      real(wp), allocatable :: alpha(:)

      !>
      real(wp), allocatable :: zeff(:)

      !>
      real(wp), allocatable :: electronegativity(:)

      !>
      real(wp) :: cutoff

   end type TRepulsionData


   !>
   type :: THamiltonianData

      !>
      integer, allocatable :: principalQuantumNumber(:, :)

      !>
      integer, allocatable :: angShell(:, :)

      !>
      integer, allocatable :: valenceShell(:, :)

      !>
      integer, allocatable :: numberOfPrimitives(:, :)

      !>
      real(wp), allocatable :: slaterExponent(:, :)

      !>
      real(wp), allocatable :: selfEnergy(:, :)

      !>
      real(wp), allocatable :: referenceOcc(:, :)

      !>
      real(wp), allocatable :: kCN(:, :)

      !>
      real(wp), allocatable :: electronegativity(:)

      !>
      real(wp), allocatable :: atomicRad(:)

      !>
      real(wp), allocatable :: shellPoly(:, :)

      !>
      real(wp), allocatable :: pairParam(:, :)

      !>
      real(wp), allocatable :: kQShell(:, :)

      !>
      real(wp), allocatable :: kQAtom(:)

   end type THamiltonianData


   !>
   type :: TCoulombData

      !>
      logical :: enEquilibration

      !>
      logical :: secondOrder

      !>
      logical :: thirdOrder

      !>
      logical :: shellResolved

      !>
      real(wp), allocatable :: chemicalHardness(:)

      !>
      real(wp), allocatable :: shellHardness(:, :)

      !>
      real(wp), allocatable :: thirdOrderAtom(:)

      !>
      real(wp), allocatable :: chargeWidth(:)

      !>
      real(wp), allocatable :: electronegativity(:)

      !>
      real(wp), allocatable :: kCN(:)

   end type TCoulombData


   !>
   type :: TMultipoleData

      !>
      real(wp) :: cnShift

      !>
      real(wp) :: cnExp

      !>
      real(wp) :: cnRMax

      !>
      real(wp) :: dipDamp

      !>
      real(wp) :: quadDamp

      !>
      real(wp), allocatable :: valenceCN(:)

      !>
      real(wp), allocatable :: multiRad(:)

      !>
      real(wp), allocatable :: dipKernel(:)

      !>
      real(wp), allocatable :: quadKernel(:)

   end type TMultipoleData


   !>
   type :: THalogenData

      !>
      real(wp) :: radScale

      !>
      real(wp) :: dampingPar

      !>
      real(wp), allocatable :: bondStrength(:)

      !>
      real(wp), allocatable :: atomicRad(:)

   end type THalogenData


   !>
   type :: TxTBData

      !>
      integer, allocatable :: nShell(:)

      !>
      type(TRepulsionData) :: repulsion

      !>
      type(THamiltonianData) :: hamiltonian

      !>
      type(TCoulombData) :: coulomb

      !>
      type(TMultipoleData), allocatable :: multipole

      !>
      type(THalogenData), allocatable :: halogen

   end type TxTBData


   interface init
      module procedure :: initRepulsion
      module procedure :: initHalogen
      module procedure :: initMultipole
      module procedure :: initCoulomb
   end interface init


   ! ========================================================================
   ! MULTIPOLE DATA
   !>
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

   !>
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


subroutine generateValenceShellData(valenceShell, nShell, angShell)

   integer, intent(out) :: valenceShell(:, :)

   integer, intent(in) :: nShell(:)

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


subroutine initHalogen(self, radScale, dampingPar, halogenBond)

   !>
   type(THalogenData), intent(out) :: self

   !>
   real(wp), intent(in) :: radScale

   !>
   real(wp), intent(in) :: dampingPar

   !>
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

   !>
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

   !>
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

   !>
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


end module xtb_xtb_data
