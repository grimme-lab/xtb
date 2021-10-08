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

module test_latticepoint
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_latticepoint

contains

!> Collect all exported unit tests
subroutine collect_latticepoint(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("pbc3d", test_latticepoint_pbc3d) &
      ]

end subroutine collect_latticepoint


subroutine test_latticepoint_pbc3d(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_latticepoint, only : TLatticePoint, init
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lattice_SiO2(3,3) = reshape(&
      &[ 8.7413053236641_wp,  0.0000000000000_wp,  0.0000000000000_wp,  &
      &  0.0000000000000_wp,  8.7413053236641_wp,  0.0000000000000_wp,  &
      &  0.0000000000000_wp,  0.0000000000000_wp,  8.7413053236641_wp], &
      & shape(lattice_SiO2))

   real(wp), parameter :: lattice_CaF2(3,3) = reshape(&
      &[ 5.9598811567890_wp,  2.1071361905157_wp,  3.6496669404404_wp,  &
      &  0.0000000000000_wp,  6.3214085715472_wp,  3.6496669404404_wp,  &
      &  0.0000000000000_wp,  0.0000000000000_wp,  7.2993338808807_wp], &
      & shape(lattice_CaF2))

   real(wp), parameter :: lattice_ammonia(3,3) = reshape(&
      &[ 6.4411018522600_wp,  0.0492571261505_wp,  0.2192046129910_wp,  &
      &  0.0462076831749_wp,  6.6435057067500_wp,  0.1670513770770_wp,  &
      &  0.2262248220170_wp, -0.9573234940220_wp,  6.7608039126200_wp], &
      & shape(lattice_ammonia)) * aatoau

   type(TEnvironment) :: env
   type(TLatticePoint) :: latp
   logical :: fail
   real(wp), allocatable :: latticePoint(:, :)

   call init(env)

   call init(latp, env, lattice_SiO2, boundaryCondition%pbc3d, 40.0_wp)
   call env%check(fail)
   call check_(error, .not.fail)
   call check_(error, allocated(latp%trans))
   call check_(error, allocated(latp%dist2))
   call check_(error, latp%nTrans, 389)

   call latp%getLatticePoints(latticePoint)
   call check_(error, size(latticePoint, dim=2), 389)

   !> Simply switch lattice without reinitalization
   call latp%update(env, lattice_CaF2)
   call env%check(fail)
   call check_(error, .not.fail)
   call check_(error, allocated(latp%trans))
   call check_(error, allocated(latp%dist2))
   call check_(error, latp%nTrans, 959)

   call latp%getLatticePoints(latticePoint, 30.0_wp)
   call check_(error, size(latticePoint, dim=2), 381)

   call latp%getLatticePoints(latticePoint, 50.0_wp)
   call check_(error, size(latticePoint, dim=2), 959)

   !> Reinitialize with new lattice and new cutoff
   call init(latp, env, lattice_ammonia, boundaryCondition%pbc3d, 60.0_wp)
   call env%check(fail)
   call check_(error, .not.fail)
   call check_(error, allocated(latp%trans))
   call check_(error, allocated(latp%dist2))
   call check_(error, latp%nTrans, 451)

   call latp%getLatticePoints(latticePoint)
   call check_(error, size(latticePoint, dim=2), 451)

   call latp%getLatticePoints(latticePoint, 40.0_wp)
   call check_(error, size(latticePoint, dim=2), 143)

   !> Reinitialize generator to exclude inversion symmetry
   call init(latp, env, lattice_ammonia, boundaryCondition%pbc3d, 50.0_wp, &
      & excludeInversion=.true.)
   call env%check(fail)
   call check_(error, .not.fail)
   call check_(error, allocated(latp%trans))
   call check_(error, allocated(latp%dist2))
   call check_(error, latp%nTrans, 130)

   call latp%getLatticePoints(latticePoint, 40.0_wp)
   call check_(error, size(latticePoint, dim=2), 72) ! 143 / 2 + 1

   call latp%getLatticePoints(latticePoint, 33.0_wp)
   call check_(error, size(latticePoint, dim=2), 41)

end subroutine test_latticepoint_pbc3d

end module test_latticepoint
