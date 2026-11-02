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
   use testdrive, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau
   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart, only : TRestart
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, &
      & newWavefunction
   use xtb_setparam, only : p_olev_normal
   use xtb_optimizer, only : ancopt
   use xtb_relaxation_engine, only : fire, l_ancopt
   use xtb_fixparam, only : fixset
   use xtb_type_data, only : scc_results
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
      new_unittest("single-atom-l-ancopt", test_single_atom_l_ancopt), &
      new_unittest("fixed-atoms-ancopt", test_fixed_atoms_ancopt), &
      new_unittest("fixed-atoms-l-ancopt", test_fixed_atoms_l_ancopt), &
      new_unittest("all-fixed-ancopt", test_all_fixed_ancopt), &
      new_unittest("all-fixed-l-ancopt", test_all_fixed_l_ancopt), &
      new_unittest("all-fixed-fire", test_all_fixed_fire) &
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

!> Water dimer with all atoms but one hydrogen fixed. The free atom has all
!> three Cartesian degrees of freedom, projecting rotations of the whole
!> molecule would restrict it to a line and stop with a large gradient.
subroutine test_fixed_atoms_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 6), sigma(3, 3)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1, 2, 3, 4, 6])
   open(newunit=ilog, status='scratch')
   call ancopt(env, ilog, mol, chk, calc, egap, 300.0_wp, 250, 200, &
      & energy, gradient, sigma, p_olev_normal, .false., fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_fixed_atoms_result(error, fail, gradient)
end subroutine test_fixed_atoms_ancopt

subroutine test_fixed_atoms_l_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 6), sigma(3, 3)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1, 2, 3, 4, 6])
   open(newunit=ilog, status='scratch')
   call l_ancopt(env, ilog, mol, chk, calc, p_olev_normal, 200, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_fixed_atoms_result(error, fail, gradient)
end subroutine test_fixed_atoms_l_ancopt

!> With all atoms fixed there are no degrees of freedom left, the optimizers
!> must return immediately without moving any atom.
subroutine test_all_fixed_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 6), sigma(3, 3)
   real(wp), allocatable :: xyz0(:, :)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1, 2, 3, 4, 5, 6])
   xyz0 = mol%xyz
   open(newunit=ilog, status='scratch')
   call ancopt(env, ilog, mol, chk, calc, egap, 300.0_wp, 250, 200, &
      & energy, gradient, sigma, p_olev_normal, .false., fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_all_fixed_result(error, env, fail, iter_needed, mol%xyz, xyz0)
end subroutine test_all_fixed_ancopt

subroutine test_all_fixed_l_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 6), sigma(3, 3)
   real(wp), allocatable :: xyz0(:, :)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1, 2, 3, 4, 5, 6])
   xyz0 = mol%xyz
   open(newunit=ilog, status='scratch')
   call l_ancopt(env, ilog, mol, chk, calc, p_olev_normal, 200, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_all_fixed_result(error, env, fail, iter_needed, mol%xyz, xyz0)
end subroutine test_all_fixed_l_ancopt

subroutine test_all_fixed_fire(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 6), sigma(3, 3)
   real(wp), allocatable :: xyz0(:, :)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1, 2, 3, 4, 5, 6])
   xyz0 = mol%xyz
   open(newunit=ilog, status='scratch')
   call fire(env, ilog, mol, chk, calc, p_olev_normal, 200, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_all_fixed_result(error, env, fail, iter_needed, mol%xyz, xyz0)
end subroutine test_all_fixed_fire

subroutine setup_fixed_water_dimer(env, mol, chk, calc, energy, egap, &
      & gradient, sigma, fixed)
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TxTBCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(3, 6), sigma(3, 3)
   integer, intent(in) :: fixed(:)
   integer, parameter :: nat = 6
   character(len=*), parameter :: sym(nat) = ["O", "H", "H", "O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      &  0.000_wp,  0.000_wp,  0.000_wp, &
      &  0.957_wp,  0.000_wp,  0.000_wp, &
      & -0.240_wp,  0.927_wp,  0.000_wp, &
      &  2.900_wp,  0.000_wp,  0.000_wp, &
      &  3.300_wp,  0.800_wp,  0.300_wp, &
      &  3.300_wp, -0.400_wp, -0.750_wp], shape(xyz))
   type(scc_results) :: res

   call init(env)
   call init(mol, sym, xyz*aatoau)
   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)

   call fixset%allocate(size(fixed))
   fixset%n = size(fixed)
   fixset%atoms(:) = fixed

   call calc%singlepoint(env, mol, chk, 0, .false., energy, gradient, sigma, &
      & egap, res)
end subroutine setup_fixed_water_dimer

subroutine check_fixed_atoms_result(error, fail, gradient)
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in) :: fail
   real(wp), intent(in) :: gradient(:, :)

   call check(error, .not.fail)
   if (allocated(error)) return
   ! gradient threshold of the normal optimization level
   if (norm2(gradient(:, 5)) > 1.0e-3_wp) then
      call test_failed(error, "Gradient of the free atom is not converged")
   end if
end subroutine check_fixed_atoms_result

subroutine check_all_fixed_result(error, env, fail, iter_needed, xyz, xyz0)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment), intent(inout) :: env
   logical, intent(in) :: fail
   integer, intent(in) :: iter_needed
   real(wp), intent(in) :: xyz(:, :), xyz0(:, :)
   logical :: exitRun

   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (allocated(error)) return
   call check(error, .not.fail)
   if (allocated(error)) return
   call check(error, iter_needed, 0)
   if (allocated(error)) return
   call check(error, maxval(abs(xyz - xyz0)), 0.0_wp)
end subroutine check_all_fixed_result

end module test_optimizer
