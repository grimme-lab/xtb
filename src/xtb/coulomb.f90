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

!> Implementation of an isotropic electrostatics container
module xtb_xtb_coulomb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : blas_dot, blas_symv
   use xtb_xtb_data, only : TCoulombData
   use xtb_xtb_thirdorder, only : TThirdOrder, init
   implicit none
   private

   public :: TxTBCoulomb, init


   !> Isotropic electrostatics
   type :: TxTBCoulomb

      !> Third order electrostatics
      type(TThirdOrder) :: thirdOrder

      !> Coulomb matrix
      real(wp), allocatable :: jmat(:, :)

      !> Scratch memory for storing the shifts
      real(wp), allocatable :: shift(:)

   contains

      !> Add shifts from isotropic electrostatics
      procedure :: addShift

      !> Get energy from isotropic electrostatics
      procedure :: getEnergy

   end type TxTBCoulomb


   !> Initialize isotropic electrostatics
   interface init
      module procedure :: initCoulomb
   end interface init


contains


!> Initialize isotropic electrostatics from parametrisation data
subroutine initCoulomb(self, input, nshell, num)

   !> Instance of the isotropic electrostatics
   type(TxTBCoulomb), intent(out) :: self

   !> Parametrisation data for coulombic interactions
   type(TCoulombData), intent(in) :: input

   !> Number of shells for each species
   integer, intent(in) :: nshell(:)

   !> Atomic numbers of each element
   integer, intent(in) :: num(:)

   integer :: nsh

   call init(self%thirdOrder, input, nshell, num)

   nsh = sum(nshell(num))
   allocate(self%jmat(nsh, nsh))
   allocate(self%shift(nsh))

end subroutine initCoulomb


!> Add shifts from isotropic electrostatics
subroutine addShift(self, qat, qsh, atomicShift, shellShift)

   !> Instance of the isotropic electrostatics
   class(TxTBCoulomb), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Atomic potential shift
   real(wp), intent(inout) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   integer :: nsh

   nsh = size(shellShift)

   call self%thirdOrder%addShift(qat, qsh, atomicShift, shellShift)

   call blas_symv('l', nsh, 1.0_wp, self%jmat, nsh, qsh, 1, 1.0_wp, shellShift, 1)

end subroutine addShift


!> Get energy from isotropic electrostatics
pure subroutine getEnergy(self, qat, qsh, energy)

   !> Instance of the isotropic electrostatics
   class(TxTBCoulomb), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Third order contribution to the energy
   real(wp), intent(out) :: energy

   integer :: nsh
   real(wp) :: eThird

   nsh = size(qsh)

   call self%thirdOrder%getEnergy(qat, qsh, eThird)

   call blas_symv('l', nsh, 1.0_wp, self%jmat, nsh, qsh, 1, 0.0_wp, self%shift, 1)
   energy = 0.5_wp * blas_dot(nsh, self%shift, 1, qsh, 1) + eThird

end subroutine getEnergy


end module xtb_xtb_coulomb
