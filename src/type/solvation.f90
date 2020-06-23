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

!> Abstract solvation model
module xtb_type_solvation
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   implicit none
   private

   public :: TSolvation


   type, abstract :: TSolvation
   contains

      !> Update coordinates and internal state
      procedure(update), deferred :: update

      !> Add potential shift
      procedure(addShift), deferred :: addShift

      !> Calculate solvation energy
      procedure(getEnergy), deferred :: getEnergy

      !> Calculate derivatives of solvation energy
      procedure(addGradient), deferred :: addGradient

   end type TSolvation


   abstract interface
      !> Update coordinates and internal state
      subroutine update(self, env, num, xyz)
         import :: TSolvation, TEnvironment, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Computation environment
         type(TEnvironment), intent(inout) :: env

         !> Atomic numbers
         integer, intent(in) :: num(:)

         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)

      end subroutine update

      !> Add potential shift
      subroutine addShift(self, env, qat, qsh, atomicShift, shellShift)
         import :: TSolvation, TEnvironment, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Computation environment
         type(TEnvironment), intent(inout) :: env

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Atomic potential shift
         real(wp), intent(inout) :: atomicShift(:)

         !> Shell-resolved potential shift
         real(wp), intent(inout) :: shellShift(:)

      end subroutine addShift


      !> Calculate solvation energy
      subroutine getEnergy(self, env, qat, qsh, energy)
         import :: TSolvation, TEnvironment, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Computation environment
         type(TEnvironment), intent(inout) :: env

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Total solvation energy
         real(wp), intent(out) :: energy

      end subroutine getEnergy


      !> Calculate derivatives of solvation energy
      subroutine addGradient(self, env, num, xyz, qat, qsh, gradient)
         import :: TSolvation, TEnvironment, wp

         !> Instance of the solvation model
         class(TSolvation), intent(inout) :: self

         !> Computation environment
         type(TEnvironment), intent(inout) :: env

         !> Atomic numbers
         integer, intent(in) :: num(:)

         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(:, :)

         !> Atomic partial charges
         real(wp), intent(in) :: qat(:)

         !> Shell-resolved partial charges
         real(wp), intent(in) :: qsh(:)

         !> Molecular gradient
         real(wp), intent(inout) :: gradient(:, :)

      end subroutine addGradient

   end interface


end module xtb_type_solvation
