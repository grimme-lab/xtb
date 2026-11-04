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
      new_unittest("all-fixed-fire", test_all_fixed_fire), &
      new_unittest("fixed-atom-periodic-ancopt", test_fixed_atom_periodic_ancopt), &
      new_unittest("fixed-atom-periodic-l-ancopt", &
      & test_fixed_atom_periodic_l_ancopt), &
      new_unittest("all-fixed-periodic-fire", test_all_fixed_periodic_fire) &
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

!> Water in a periodic box with the oxygen fixed. Rotations are no symmetry
!> under periodic boundary conditions, therefore only the three coordinates of
!> the fixed atom are removed and the ANC generation must not drop the three
!> rotational modes on top.
subroutine test_fixed_atom_periodic_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 3), sigma(3, 3)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_periodic_water(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1])
   open(newunit=ilog, status='scratch')
   call ancopt(env, ilog, mol, chk, calc, egap, 300.0_wp, 250, 200, &
      & energy, gradient, sigma, p_olev_normal, .false., fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_fixed_atom_periodic_result(error, env, fail, gradient)
end subroutine test_fixed_atom_periodic_ancopt

subroutine test_fixed_atom_periodic_l_ancopt(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, egap, gradient(3, 3), sigma(3, 3)
   logical :: fail
   integer :: iter_needed, ilog

   call setup_periodic_water(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, [1])
   open(newunit=ilog, status='scratch')
   call l_ancopt(env, ilog, mol, chk, calc, p_olev_normal, 200, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call check_fixed_atom_periodic_result(error, env, fail, gradient)
end subroutine test_fixed_atom_periodic_l_ancopt

!> Strained CaF2 with all atoms fixed. FIRE still has to relax the lattice,
!> while the Cartesian positions of the atoms stay untouched.
subroutine test_all_fixed_periodic_fire(error)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, energy0, egap, gradient(3, 3), sigma(3, 3)
   real(wp), allocatable :: xyz0(:, :), lattice0(:, :)
   logical :: fail, exitRun
   integer :: iter_needed, ilog

   call setup_strained_caf2(env, mol, chk, calc, energy, egap, gradient, sigma)
   xyz0 = mol%xyz
   lattice0 = mol%lattice
   energy0 = energy
   open(newunit=ilog, status='scratch')
   call fire(env, ilog, mol, chk, calc, p_olev_normal, 200, energy, egap, &
      & gradient, sigma, 0, fail, iter_needed)
   close(ilog)
   call fixset%deallocate

   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (allocated(error)) return
   call check(error, .not.fail)
   if (allocated(error)) return
   if (iter_needed < 1) then
      call test_failed(error, "Lattice was not optimized")
      return
   end if
   if (maxval(abs(mol%lattice - lattice0)) < 1.0e-4_wp) then
      call test_failed(error, "Lattice did not change")
      return
   end if
   if (energy >= energy0) then
      call test_failed(error, "Energy was not lowered by the lattice relaxation")
      return
   end if
   if (maxval(abs(mol%xyz - xyz0)) > 1.0e-10_wp) then
      call test_failed(error, "Fixed atoms were moved")
   end if
end subroutine test_all_fixed_periodic_fire

subroutine setup_periodic_water(env, mol, chk, calc, energy, egap, &
      & gradient, sigma, fixed)
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TxTBCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(3, 3), sigma(3, 3)
   integer, intent(in) :: fixed(:)
   integer, parameter :: nat = 3
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      &  4.000_wp,  4.000_wp,  4.000_wp, &
      &  5.050_wp,  4.000_wp,  4.100_wp, &
      &  3.700_wp,  5.000_wp,  3.900_wp], shape(xyz))
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      &  8.0_wp, 0.0_wp, 0.0_wp, &
      &  0.0_wp, 8.0_wp, 0.0_wp, &
      &  0.0_wp, 0.0_wp, 8.0_wp], shape(lattice))

   call setup_periodic(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & sym, xyz*aatoau, lattice*aatoau, fixed)
end subroutine setup_periodic_water

subroutine setup_strained_caf2(env, mol, chk, calc, energy, egap, &
      & gradient, sigma)
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TxTBCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(3, 3), sigma(3, 3)
   integer, parameter :: nat = 3
   character(len=*), parameter :: sym(nat) = ["F ", "F ", "Ca"]
   ! fractional coordinates shifted away from the cell boundary
   real(wp), parameter :: abc(3, nat) = reshape([ &
      & 0.35_wp, 0.35_wp, 0.35_wp, &
      & 0.85_wp, 0.85_wp, 0.85_wp, &
      & 0.10_wp, 0.10_wp, 0.10_wp], shape(abc))
   ! equilibrium lattice expanded by 3 %
   real(wp), parameter :: lattice(3, 3) = 1.03_wp * reshape([ &
      & 5.9598811567890_wp, 2.1071361905157_wp, 3.6496669404404_wp, &
      & 0.0000000000000_wp, 6.3214085715472_wp, 3.6496669404404_wp, &
      & 0.0000000000000_wp, 0.0000000000000_wp, 7.2993338808807_wp], &
      & shape(lattice))

   call setup_periodic(env, mol, chk, calc, energy, egap, gradient, sigma, &
      & sym, matmul(lattice, abc), lattice, [1, 2, 3])
end subroutine setup_strained_caf2

subroutine setup_periodic(env, mol, chk, calc, energy, egap, gradient, &
      & sigma, sym, xyz, lattice, fixed)
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TxTBCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(:, :), sigma(3, 3)
   character(len=*), intent(in) :: sym(:)
   real(wp), intent(in) :: xyz(:, :), lattice(3, 3)
   integer, intent(in) :: fixed(:)
   type(scc_results) :: res

   call init(env)
   call init(mol, sym, xyz, lattice=lattice)
   call newXTBCalculator(env, mol, calc, method=0)
   call newWavefunction(env, mol, calc, chk)

   call fixset%allocate(size(fixed))
   fixset%n = size(fixed)
   fixset%atoms(:) = fixed

   call calc%singlepoint(env, mol, chk, 0, .false., energy, gradient, sigma, &
      & egap, res)
end subroutine setup_periodic

subroutine check_fixed_atom_periodic_result(error, env, fail, gradient)
   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment), intent(inout) :: env
   logical, intent(in) :: fail
   real(wp), intent(in) :: gradient(:, :)
   logical :: exitRun

   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (allocated(error)) return
   call check(error, .not.fail)
   if (allocated(error)) return
   ! gradient threshold of the normal optimization level
   if (norm2(gradient(:, 2:)) > 1.0e-3_wp) then
      call test_failed(error, "Gradient of the free atoms is not converged")
   end if
end subroutine check_fixed_atom_periodic_result

end module test_optimizer
