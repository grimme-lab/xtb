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
   implicit none
   private

   public :: TxTBData
   public :: TRepulsionData, TCoulombData, THamiltonianData
   public :: THalogenData, TMultipoleData


   !>
   type :: TRepulsionData

      !>
      real(wp) :: kExp

      !>
      real(wp) :: kExpLight

      !>
      real(wp) :: rExp

      !>
      real(wp), allocatable :: alpha(:)

      !>
      real(wp), allocatable :: zeff(:)

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
      real(wp), allocatable :: kCN(:, :)

      !>
      real(wp), allocatable :: electronegativity(:)

      !>
      real(wp), allocatable :: atomicRad(:)

      !>
      real(wp), allocatable :: shellPoly(:, :)

      !>
      real(wp), allocatable :: pairParam(:, :)

   end type THamiltonianData


   !>
   type :: TCoulombData

      !>
      real(wp), allocatable :: chemicalHardness(:)

      !>
      real(wp), allocatable :: shellHardness(:, :)

      !>
      real(wp), allocatable :: thirdOrderAtom(:)

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


contains


end module xtb_xtb_data
