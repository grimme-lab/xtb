! This file is part of xtb.
! SPDX-Identifier: LGPL-3.0-or-later
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

module test_wsc
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_wsc

contains

!> Collect all exported unit tests
subroutine collect_wsc(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("0d", test_wigner_seitz_0d), &
      new_unittest("3d", test_wigner_seitz_3D) &
      ]

end subroutine collect_wsc


subroutine test_wigner_seitz_0d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   use xtb_type_wsc
   type(error_type), allocatable, intent(out) :: error
   type(TMolecule) :: mol
   type(tb_wsc) :: wsc

   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [ 0.00000000000000_wp,  0.00000000000000_wp, -0.73578586109551_wp, &
      & 1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp, &
      &-1.44183152868459_wp,  0.00000000000000_wp,  0.36789293054775_wp  &
      & ],shape(xyz))

   call init(mol, at, xyz)

   call generate_wsc(mol,wsc)

   call check(error, allocated(wsc%lattr))
   call check(error, allocated(wsc%itbl))

   call check(error, all(wsc%itbl <= 1))
   call check(error, all(wsc%lattr == 0))

end subroutine test_wigner_seitz_0D

subroutine test_wigner_seitz_3D(error)
   use xtb_mctc_accuracy, only : wp

   use xtb_type_molecule
   use xtb_type_wsc

   use xtb_pbc_tools

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-10_wp
   ! CaF2
   integer, parameter :: nat = 3
   integer, parameter :: at(nat) = [9,9,20]
   real(wp),parameter :: abc(3,nat) = reshape(&
      &[0.25_wp, 0.25_wp, 0.25_wp, &
      & 0.75_wp, 0.75_wp, 0.75_wp, &
      & 0.00_wp, 0.00_wp, 0.00_wp], shape(abc))
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))

   type(TMolecule) :: mol
   type(tb_wsc) :: wsc
   real(wp), allocatable :: xyz(:, :)

   allocate(xyz(3, nat))
   call coord_trafo(nat,lattice,abc,xyz)
   call init(mol, at, xyz, lattice=lattice)

   call generate_wsc(mol,wsc)

   call check(error, allocated(wsc%lattr))
   call check(error, allocated(wsc%itbl))

   call check(error, all(abs(wsc%lattr) <= 1))
   call check(error, wsc%itbl(1,1), 12)
   call check(error, wsc%itbl(2,2), 12)
   call check(error, wsc%itbl(3,3), 12)
   call check(error, wsc%itbl(1,2), 6)
   call check(error, wsc%itbl(1,3), 4)
   call check(error, wsc%itbl(2,3), 4)
   call check(error, wsc%itbl(1,2), wsc%itbl(2,1))
   call check(error, wsc%itbl(1,3), wsc%itbl(3,1))
   call check(error, wsc%itbl(2,3), wsc%itbl(3,2))

end subroutine test_wigner_seitz_3D

end module test_wsc
