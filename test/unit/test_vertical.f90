! This file is part of xtb.
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

module test_vertical
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_vertical, only : vfukui

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction
   implicit none
   private

   public :: collect_fukui


contains

!> Collect all exported unit tests
subroutine collect_fukui(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1", test_gfn1_fukui), &
      new_unittest("gfn2", test_gfn2_fukui) &
      ]

end subroutine collect_fukui

subroutine test_gfn1_fukui(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 4
   real(wp),parameter :: thr = 1.0e-2_wp
   character(len=*), parameter :: sym(nat) = ["B", "F", "F", "F"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -4.41826178485383_wp, 2.18219865869875_wp, -0.29163266946828_wp, &
      & -3.14927320876938_wp, 1.86405526998014_wp,  1.78831473195343_wp, &
      & -4.31209733618080_wp, 4.33072818478984_wp, -1.47918413760343_wp, &
      & -5.79274623699377_wp, 0.35153057534898_wp, -1.18447939588312_wp], &
      & shape(xyz))
   real(wp), parameter :: fukui_ref(3, nat) = reshape([ &
      &  0.471_wp, 0.150_wp, 0.310_wp, &
      &  0.176_wp, 0.283_wp, 0.230_wp, &
      &  0.176_wp, 0.283_wp, 0.230_wp, &
      &  0.176_wp, 0.283_wp, 0.230_wp], &
     &  shape(fukui_ref))
   !real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i,j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:), fukui(:,:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3,mol%n), fukui(3,mol%n))

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 1, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call vfukui(env, mol, chk, calc, fukui)

   do i = 1, size(fukui_ref, 2)
      do j = 1, size(fukui_ref, 1)
         call check(error, fukui(j, i), fukui_ref(j, i), thr=thr)
      end do
   end do

end subroutine test_gfn1_fukui

subroutine test_gfn2_fukui(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 4
   real(wp),parameter :: thr = 1.0e-2_wp
   character(len=*), parameter :: sym(nat) = ["B", "F", "F", "F"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & -4.41826178485383_wp, 2.18219865869875_wp, -0.29163266946828_wp, &
      & -3.14927320876938_wp, 1.86405526998014_wp,  1.78831473195343_wp, &
      & -4.31209733618080_wp, 4.33072818478984_wp, -1.47918413760343_wp, &
      & -5.79274623699377_wp, 0.35153057534898_wp, -1.18447939588312_wp], &
      & shape(xyz))
   real(wp), parameter :: fukui_ref(3, nat) = reshape([ &
      &  0.300_wp, -0.005_wp, 0.148_wp, &
      &  0.233_wp, 0.335_wp, 0.284_wp, &
      &  0.233_wp, 0.335_wp, 0.284_wp, &
      &  0.233_wp, 0.335_wp, 0.284_wp], &
     &  shape(fukui_ref))

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i,j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:), fukui(:,:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3,mol%n), fukui(3,mol%n))

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 1, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call vfukui(env, mol, chk, calc, fukui)

   do i = 1, size(fukui_ref, 2)
      do j = 1, size(fukui_ref, 1)
         call check(error, fukui(j, i), fukui_ref(j, i), thr=thr)
      end do
   end do

end subroutine test_gfn2_fukui

end module test_vertical
