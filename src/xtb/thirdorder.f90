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

!> Implementation of the third order electrostatics
module xtb_xtb_thirdorder
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data, only : TCoulombData
   implicit none
   private

   public :: TThirdOrder, init


   !> Third order electrostatics
   type :: TThirdOrder

      !> Atomic Hubbard derivatives
      real(wp), allocatable :: atomicGam(:)

      !> Shell-resolved Hubbard derivatives
      real(wp), allocatable :: shellGam(:)

   contains

      !> Add shifts from third order contribution
      procedure :: addShift

      !> Get energy from third order contribution
      procedure :: getEnergy

   end type TThirdOrder


   !> Initialize third order electrostatics
   interface init
      module procedure :: initThirdOrder
   end interface init


contains


!> Initialize third order electrostatics
subroutine initThirdOrder(self, input, nshell, num)

   !> Instance of the third order electrostatics
   type(TThirdOrder), intent(out) :: self

   !> Parametrisation data for coulombic interactions
   type(TCoulombData), intent(in) :: input

   !> Number of shells for each species
   integer, intent(in) :: nshell(:)

   !> Atomic numbers of each element
   integer, intent(in) :: num(:)

   integer :: nat, nsh
   integer :: ii, iat, izp, ish

   nat = size(num, dim=1)
   nsh = sum(nshell(num))

   ! set 3rd order shell gammas
   if (allocated(input%thirdOrderShell)) then
      allocate(self%shellGam(nsh))
      ii = 0
      do iat = 1, nat
         izp = num(iat)
         do ish = 1, nShell(izp)
            self%shellGam(ii+ish) = input%thirdOrderShell(ish, izp)
         end do
         ii = ii + nShell(izp)
      end do
   else if (allocated(input%thirdOrderAtom)) then
      allocate(self%atomicGam(nat))
      do iat = 1, nat
         izp = num(iat)
         self%atomicGam(iat) = input%thirdOrderAtom(izp)
      end do
   end if

end subroutine initThirdOrder


pure subroutine addShift(self, qat, qsh, atomicShift, shellShift)

   !> Instance of the third order electrostatics
   class(TThirdOrder), intent(inout) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Atomic potential shift
   real(wp), intent(inout) :: atomicShift(:)

   !> Shell-resolved potential shift
   real(wp), intent(inout) :: shellShift(:)

   integer :: ii

   if (allocated(self%atomicGam)) then
      do ii = 1, size(atomicShift, dim=1)
         atomicShift(ii) = atomicShift(ii) + qat(ii)**2 * self%atomicGam(ii)
      end do
   end if

   if (allocated(self%shellGam)) then
      do ii = 1, size(shellShift, dim=1)
         shellShift(ii) = shellShift(ii) + qsh(ii)**2 * self%shellGam(ii)
      end do
   end if

end subroutine addShift


pure subroutine getEnergy(self, qat, qsh, energy)

   !> Instance of the third order electrostatics
   class(TThirdOrder), intent(in) :: self

   !> Atomic partial charges
   real(wp), intent(in) :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(in) :: qsh(:)

   !> Third order contribution to the energy
   real(wp), intent(out) :: energy

   integer :: ii

   energy = 0.0_wp

   if (allocated(self%atomicGam)) then
      do ii = 1, size(qat, dim=1)
         energy = energy + qat(ii)**3 * self%atomicGam(ii) / 3.0_wp
      end do
   end if

   if (allocated(self%shellGam)) then
      do ii = 1, size(qsh, dim=1)
         energy = energy + qsh(ii)**3 * self%shellGam(ii) / 3.0_wp
      end do
   end if

end subroutine getEnergy


end module xtb_xtb_thirdorder
