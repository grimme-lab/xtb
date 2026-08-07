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

module test_detrotra
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : sp, wp
   use xtb_detrotra, only : detrotra4, detrotra8
   implicit none
   private

   public :: collect_detrotra

   interface build_diatomic_modes
      module procedure build_diatomic_modes_sp
      module procedure build_diatomic_modes_wp
   end interface build_diatomic_modes

contains

!> Collect all exported unit tests
subroutine collect_detrotra(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("single precision linear diatomic", test_detrotra4_diatomic), &
      new_unittest("double precision linear diatomic", test_detrotra8_diatomic), &
      new_unittest("preserve sufficient low-mode candidates", test_low_mode_candidates), &
      new_unittest("adaptively expand mode candidates", test_adaptive_mode_candidates) &
      ]

end subroutine collect_detrotra

subroutine test_detrotra4_diatomic(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 2
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      0.0_wp, 0.0_wp, 0.0_wp, &
      1.0_wp, 0.0_wp, 0.0_wp], shape(xyz))
   real(sp) :: hess(3*nat, 3*nat), eig(3*nat)

   call build_diatomic_modes(hess)
   eig = [1.0e-14_sp, 2.0e-14_sp, 3.0e-14_sp, 3.4_sp, 3.4_sp, 40.6_sp]

   call detrotra4(.true., nat, xyz, hess, eig)

   call check(error, count(eig == 0.0_sp), 5)
   call check(error, eig(6), 40.6_sp, thr=1.0e-5_sp)
end subroutine test_detrotra4_diatomic

subroutine test_detrotra8_diatomic(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 2
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      0.0_wp, 0.0_wp, 0.0_wp, &
      1.0_wp, 0.0_wp, 0.0_wp], shape(xyz))
   real(wp) :: hess(3*nat, 3*nat), eig(3*nat)

   call build_diatomic_modes(hess)
   eig = [1.0e-14_wp, 2.0e-14_wp, 3.0e-14_wp, 3.4_wp, 3.4_wp, 40.6_wp]

   call detrotra8(.true., nat, xyz, hess, eig)

   call check(error, count(eig == 0.0_wp), 5)
   call check(error, eig(6), 40.6_wp, thr=1.0e-12_wp)
end subroutine test_detrotra8_diatomic

subroutine test_low_mode_candidates(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      0.0_wp, 0.0_wp, 0.0_wp, &
      1.0_wp, 0.0_wp, 0.0_wp, &
      0.0_wp, 1.0_wp, 0.0_wp], shape(xyz))
   real(wp) :: hess(3*nat, 3*nat), eig(3*nat)
   integer :: i

   hess = 0.0_wp
   do i = 1, 6
      hess(i, i) = 1.0_wp
   end do
   hess([1, 4, 7], 7) = 1.0_wp
   hess(8, 8) = 1.0_wp
   hess(9, 9) = 1.0_wp
   eig = [0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, &
      3.0_wp, 4.0_wp, 5.0_wp]

   call detrotra8(.false., nat, xyz, hess, eig)

   call check(error, count(eig == 0.0_wp), 6)
   call check(error, eig(7), 3.0_wp, thr=1.0e-12_wp)
end subroutine test_low_mode_candidates

subroutine test_adaptive_mode_candidates(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      0.0_wp, 0.0_wp, 0.0_wp, &
      1.0_wp, 0.0_wp, 0.0_wp, &
      0.0_wp, 1.0_wp, 0.0_wp], shape(xyz))
   real(wp) :: hess(3*nat, 3*nat), eig(3*nat)
   integer :: i

   hess = 0.0_wp
   do i = 1, 6
      hess(i, i) = 1.0_wp
   end do
   eig = [0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.01_wp, 0.10_wp, &
      0.20_wp, 0.30_wp, 0.40_wp]

   call detrotra8(.false., nat, xyz, hess, eig)

   call check(error, count(eig == 0.0_wp), 6)
   call check(error, eig(6), 0.0_wp)
   call check(error, eig(7), 0.20_wp, thr=1.0e-12_wp)
end subroutine test_adaptive_mode_candidates

subroutine build_diatomic_modes_sp(hess)
   real(sp), intent(out) :: hess(6, 6)

   hess = 0.0_sp
   hess([1, 4], 1) = 1.0_sp
   hess([2, 5], 2) = 1.0_sp
   hess([3, 6], 3) = 1.0_sp
   hess([4, 5], 4) = [-1.0_sp, 1.0_sp]
   hess([4, 6], 5) = [-1.0_sp, 1.0_sp]
   hess(4, 6) = 1.0_sp
end subroutine build_diatomic_modes_sp

subroutine build_diatomic_modes_wp(hess)
   real(wp), intent(out) :: hess(6, 6)

   hess = 0.0_wp
   hess([1, 4], 1) = 1.0_wp
   hess([2, 5], 2) = 1.0_wp
   hess([3, 6], 3) = 1.0_wp
   hess([4, 5], 4) = [-1.0_wp, 1.0_wp]
   hess([4, 6], 5) = [-1.0_wp, 1.0_wp]
   hess(4, 6) = 1.0_wp
end subroutine build_diatomic_modes_wp

end module test_detrotra
