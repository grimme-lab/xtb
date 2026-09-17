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

module test_optimizer
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart, only : TRestart
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, &
      & newWavefunction
   use xtb_setparam, only : p_olev_normal
   use xtb_optimizer, only : ancopt
   use xtb_relaxation_engine, only : fire, l_ancopt
   implicit none
   private

   public :: collect_optimizer

contains

!> Collect all exported unit tests
subroutine collect_optimizer(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("single-atom-ancopt", test_single_atom_ancopt), &
      new_unittest("single-atom-fire", test_single_atom_fire), &
      new_unittest("single-atom-l-ancopt", test_single_atom_l_ancopt) &
      ]

end subroutine collect_optimizer

subroutine test_single_atom_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 1), sigma(3, 3)
   logical :: fail
   integer :: iter_needed

   call setup_single_atom(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & fail, iter_needed)
   call ancopt(env, stdout, mol, chk, calc, egap, 300.0_wp, 50, 50, &
      & energy, gradient, sigma, p_olev_normal, .false., fail, iter_needed)

   call check_single_atom_result(error, fail, iter_needed)
end subroutine test_single_atom_ancopt

subroutine test_single_atom_fire(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 1), sigma(3, 3)
   logical :: fail
   integer :: iter_needed

   call setup_single_atom(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & fail, iter_needed)
   call fire(env, stdout, mol, chk, calc, p_olev_normal, 50, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)

   call check_single_atom_result(error, fail, iter_needed)
end subroutine test_single_atom_fire

subroutine test_single_atom_l_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 1), sigma(3, 3)
   logical :: fail
   integer :: iter_needed

   call setup_single_atom(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & fail, iter_needed)
   call l_ancopt(env, stdout, mol, chk, calc, p_olev_normal, 50, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)

   call check_single_atom_result(error, fail, iter_needed)
end subroutine test_single_atom_l_ancopt

subroutine setup_single_atom(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & fail, iter_needed)
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TxTBCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(3, 1), sigma(3, 3)
   logical, intent(out) :: fail
   integer, intent(out) :: iter_needed
   character(len=*), parameter :: sym(1) = ["Li"]
   real(wp), parameter :: xyz(3, 1) = 0.0_wp

   call init(env)
   call init(mol, sym, xyz, chrg=1.0_wp)
   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   energy = 0.0_wp
   egap = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp
   fail = .true.
   iter_needed = -1
end subroutine setup_single_atom

subroutine check_single_atom_result(error, fail, iter_needed)
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in) :: fail
   integer, intent(in) :: iter_needed

   call check(error, .not.fail)
   if (allocated(error)) return
   call check(error, iter_needed, 0)
end subroutine check_single_atom_result

end module test_optimizer
