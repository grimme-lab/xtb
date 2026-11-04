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
      new_unittest("all-fixed-periodic-fire", test_all_fixed_periodic_fire), &
      new_unittest("default-periodic-optimizer", test_default_periodic_optimizer), &
      new_unittest("fixed-cell-periodic-optimizer", test_fixed_cell_periodic_optimizer), &
      new_unittest("nonperiodic-default-optimizer", test_nonperiodic_default_optimizer), &
      new_unittest("explicit-periodic-optimizer", test_explicit_periodic_optimizer), &
      new_unittest("periodic-optlevel-and-failures", test_periodic_optlevel_and_failures) &
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


subroutine setup_periodic_nh3_tblite(env, mol, chk, calc, energy, egap, gradient, sigma)
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, &
      & newTBLiteCalculator, newTBLiteWavefunction
   type(TEnvironment), intent(out) :: env
   type(TMolecule), intent(out) :: mol
   type(TRestart), intent(out) :: chk
   type(TTBLiteCalculator), intent(out) :: calc
   real(wp), intent(out) :: energy, egap, gradient(3, 16), sigma(3, 3)
   integer, parameter :: nat = 16
   character(len=*), parameter :: sym(nat) = [ &
      & "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", "H", &
      & "N", "N", "N", "N"]
   real(wp), parameter :: xyz_aa(3, nat) = reshape([ &
      & 2.2300000000_wp, 1.7800007600_wp, 0.9468238481_wp, &
      & 1.7800007600_wp, 0.9468238481_wp, 2.2596045809_wp, &
      & 0.9468238481_wp, 2.2596045809_wp, 1.7800007600_wp, &
      & 4.7671705534_wp, 1.6326980357_wp, 4.9661124380_wp, &
      & 4.2819692113_wp, 2.4606754877_wp, 3.6526280234_wp, &
      & 3.4527982039_wp, 1.1481588918_wp, 4.1381359401_wp, &
      & 4.1381359401_wp, 3.4527982039_wp, 1.1481588918_wp, &
      & 4.9661124380_wp, 4.7671705534_wp, 1.6326980357_wp, &
      & 3.6526280234_wp, 4.2819692113_wp, 2.4606754877_wp, &
      & 2.4606754877_wp, 3.6526280234_wp, 4.2819692113_wp, &
      & 1.1481588918_wp, 4.1381359401_wp, 3.4527982039_wp, &
      & 1.6326980357_wp, 4.9661124380_wp, 4.7671705534_wp, &
      & 1.4358572817_wp, 1.4358572817_wp, 1.4358572817_wp, &
      & 3.9411899768_wp, 1.9736199441_wp, 4.4787989040_wp, &
      & 4.4787989040_wp, 3.9411899768_wp, 1.9736199441_wp, &
      & 1.9736199441_wp, 4.4787989040_wp, 3.9411899768_wp], [3, nat])
   real(wp), parameter :: lattice_aa(3, 3) = reshape([ &
      & 5.01336_wp, 0.00000_wp, 0.00000_wp, &
      & 0.00000_wp, 5.01336_wp, 0.00000_wp, &
      & 0.00000_wp, 0.00000_wp, 5.01336_wp], [3, 3])
   type(scc_results) :: res

   call init(env)
   call init(mol, sym, xyz_aa * aatoau, lattice=lattice_aa * aatoau)
   call newTBLiteCalculator(env, mol, calc, TTBLiteInput(method="gfn2"))
   call newTBLiteWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 0, .false., energy, gradient, sigma, &
      & egap, res)
end subroutine setup_periodic_nh3_tblite


subroutine cleanup_opt_files()
   integer :: iunit, stat
   character(len=*), parameter :: files(5) = [character(len=20) :: &
      & "xtbopt.log", "xtboptlog.cif", "xtbopt_emin.coord", ".xtboptok", "NOT_CONVERGED"]
   integer :: i
   logical :: exist

   do i = 1, size(files)
      inquire(file=trim(files(i)), exist=exist)
      if (exist) then
         open(newunit=iunit, file=trim(files(i)), status='old', iostat=stat)
         if (stat == 0) close(iunit, status='delete')
      end if
   end do
end subroutine cleanup_opt_files


subroutine test_default_periodic_optimizer(error)
   use xtb_features, only : get_xtb_feature
   use testdrive, only : skip_test
   use xtb_tblite_calculator, only : TTBLiteCalculator
   use xtb_main_defaults, only : initDefaults
   use xtb_geoopt, only : geometry_optimization
   use xtb_setparam, only : set, p_engine_pbc_lbfgs
   use xtb_pbc_tools, only : dlat_to_cell, dlat_to_rlat, dlat_to_dvol, xyz_to_abc
   use xtb_io_writer, only : writeMolecule
   use xtb_io_reader, only : readMolecule
   use xtb_mctc_filetypes, only : fileType
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_constants, only : pi

   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol, mol_written
   type(TRestart) :: chk
   type(TTBLiteCalculator) :: calc
   real(wp) :: energy, energy0, egap, gradient(3, 16), sigma(3, 3)
   real(wp) :: xyz0(3, 16), lattice0(3, 3), cellpar0(6), vol0
   real(wp) :: cellpar_ref(6), rec_lat_ref(3, 3), abc_ref(3, 16), vol_ref
   real(wp) :: cellpar_w(6), vol_w
   logical :: fail, exitRun, optcell_save
   integer :: iter_needed, iunit, isum, stat
   integer, allocatable :: engine_save
   character(len=256) :: line, vol_str, cell_str

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   if (allocated(set%opt_engine)) engine_save = set%opt_engine
   if (allocated(set%opt_engine)) deallocate(set%opt_engine)
   optcell_save = set%optcell
   set%optcell = .true.

   call setup_periodic_nh3_tblite(env, mol, chk, calc, energy, egap, gradient, sigma)
   xyz0 = mol%xyz
   lattice0 = mol%lattice
   cellpar0 = mol%cellpar
   vol0 = mol%volume
   energy0 = energy

   ! Check that initDefaults automatically selects pbc_lbfgs for periodic systems
   call initDefaults(env, calc, mol, 0)
   call check(error, allocated(set%opt_engine))
   if (.not.allocated(error)) then
      call check(error, set%opt_engine, p_engine_pbc_lbfgs)
   end if
   if (allocated(error)) then
      call cleanup_opt_files()
      if (allocated(engine_save)) then
         set%opt_engine = engine_save
      else if (allocated(set%opt_engine)) then
         deallocate(set%opt_engine)
      end if
      set%optcell = optcell_save
      return
   end if

   ! Also verify geometry_optimization selects pbc_lbfgs even when set%opt_engine is unallocated
   deallocate(set%opt_engine)
   call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 0, &
      & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
   call cleanup_opt_files()

   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (.not.allocated(error)) call check(error, .not.fail)
   if (.not.allocated(error)) then
      if (maxval(abs(mol%xyz - xyz0)) < 1.0e-3_wp) then
         call test_failed(error, "Atomic positions were not optimized")
      end if
   end if
   if (.not.allocated(error)) then
      if (maxval(abs(mol%lattice - lattice0)) < 1.0e-3_wp) then
         call test_failed(error, "Lattice vectors were not optimized")
      end if
   end if
   if (.not.allocated(error)) then
      if (energy >= energy0) then
         call test_failed(error, "Energy was not lowered by variable-cell optimization")
      end if
   end if

   ! Test H: Verify mol derived periodic properties match optimized lattice and written structure
   if (.not.allocated(error)) then
      call dlat_to_cell(mol%lattice, cellpar_ref)
      call dlat_to_rlat(mol%lattice, rec_lat_ref)
      vol_ref = dlat_to_dvol(mol%lattice)
      call xyz_to_abc(mol%n, mol%lattice, mol%xyz, abc_ref, mol%pbc)

      call check(error, maxval(abs(mol%cellpar - cellpar_ref)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%rec_lat - rec_lat_ref)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, mol%volume, vol_ref, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%abc - abc_ref)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      if (abs(mol%volume - vol0) < 1.0e-4_wp) then
         call test_failed(error, "Cell volume remained at initial value after cell relaxation")
      end if
   end if

   ! Write optimized VASP structure, read it back, and compare against mol and summary printout
   if (.not.allocated(error)) then
      open(newunit=iunit, status='scratch')
      call writeMolecule(mol, iunit, format=fileType%vasp, energy=energy, gnorm=norm2(gradient))
      rewind(iunit)
      call readMolecule(env, mol_written, unit=iunit, ftype=fileType%vasp)
      close(iunit)

      call dlat_to_cell(mol_written%lattice, cellpar_w)
      vol_w = dlat_to_dvol(mol_written%lattice)
      call check(error, maxval(abs(cellpar_w - mol%cellpar)), 0.0_wp, thr=1.0e-9_wp)
      if (.not.allocated(error)) then
         call check(error, vol_w, mol%volume, thr=1.0e-9_wp)
      end if
   end if

   ! Verify printed geometry summary (main_geometry / print_pbcsum) contains updated cell parameters and volume
   if (.not.allocated(error)) then
      open(newunit=isum, status='scratch')
      call main_geometry(isum, mol)
      rewind(isum)
      write(vol_str, '(f14.7)') mol%volume * autoaa**3
      write(cell_str, '(3f14.7,1x,3f9.3)') &
         & mol%cellpar(1:3) * autoaa, mol%cellpar(4:6) * 180.0_wp / pi
      vol_str = adjustl(vol_str)
      cell_str = adjustl(cell_str)
      exitRun = .false.
      fail = .false.
      do
         read(isum, '(a)', iostat=stat) line
         if (stat /= 0) exit
         if (index(line, trim(vol_str)) > 0) exitRun = .true.
         if (index(line, trim(cell_str)) > 0) fail = .true.
      end do
      close(isum)
      if (.not.exitRun) then
         call test_failed(error, "Printed geometry summary does not contain optimized cell volume")
      else if (.not.fail) then
         call test_failed(error, "Printed geometry summary does not contain optimized cell parameters")
      end if
   end if

   if (allocated(engine_save)) then
      set%opt_engine = engine_save
   else if (allocated(set%opt_engine)) then
      deallocate(set%opt_engine)
   end if
   set%optcell = optcell_save
end subroutine test_default_periodic_optimizer


subroutine test_fixed_cell_periodic_optimizer(error)
   use xtb_features, only : get_xtb_feature
   use testdrive, only : skip_test
   use xtb_tblite_calculator, only : TTBLiteCalculator
   use xtb_geoopt, only : geometry_optimization
   use xtb_setparam, only : set

   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TTBLiteCalculator) :: calc
   real(wp) :: energy, energy0, egap, gradient(3, 16), sigma(3, 3)
   real(wp) :: xyz0(3, 16), lattice0(3, 3), cellpar0(6), vol0
   logical :: fail, exitRun, optcell_save
   integer :: iter_needed
   integer, allocatable :: engine_save

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   if (allocated(set%opt_engine)) engine_save = set%opt_engine
   if (allocated(set%opt_engine)) deallocate(set%opt_engine)
   optcell_save = set%optcell
   set%optcell = .false.

   call setup_periodic_nh3_tblite(env, mol, chk, calc, energy, egap, gradient, sigma)
   xyz0 = mol%xyz
   lattice0 = mol%lattice
   cellpar0 = mol%cellpar
   vol0 = mol%volume
   energy0 = energy

   call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 0, &
      & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
   call cleanup_opt_files()

   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (.not.allocated(error)) call check(error, .not.fail)
   if (.not.allocated(error)) then
      if (maxval(abs(mol%xyz - xyz0)) < 1.0e-3_wp) then
         call test_failed(error, "Atomic positions were not optimized in fixed-cell mode")
      end if
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%lattice - lattice0)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%cellpar - cellpar0)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, mol%volume, vol0, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      if (energy >= energy0) then
         call test_failed(error, "Energy was not lowered in fixed-cell optimization")
      end if
   end if
   if (.not.allocated(error)) then
      if (norm2(gradient) > 1.0e-3_wp) then
         call test_failed(error, "Nuclear gradient did not converge in fixed-cell optimization")
      end if
   end if

   if (allocated(engine_save)) then
      set%opt_engine = engine_save
   else if (allocated(set%opt_engine)) then
      deallocate(set%opt_engine)
   end if
   set%optcell = optcell_save
end subroutine test_fixed_cell_periodic_optimizer


subroutine test_nonperiodic_default_optimizer(error)
   use xtb_main_defaults, only : initDefaults
   use xtb_geoopt, only : geometry_optimization
   use xtb_setparam, only : set, p_engine_rf

   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res
   integer, parameter :: nat = 3
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz_aa(3, nat) = reshape([ &
      &  0.00_wp, 0.00_wp, 0.00_wp, &
      &  1.02_wp, 0.00_wp, 0.00_wp, &
      & -0.28_wp, 0.98_wp, 0.00_wp], [3, nat])
   real(wp) :: energy, energy0, egap, gradient(3, nat), sigma(3, 3)
   logical :: fail, exitRun
   integer :: iter_needed
   integer, allocatable :: engine_save

   if (allocated(set%opt_engine)) engine_save = set%opt_engine
   if (allocated(set%opt_engine)) deallocate(set%opt_engine)

   call init(env)
   call init(mol, sym, xyz_aa * aatoau)
   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 0, .false., energy0, gradient, sigma, &
      & egap, res)
   energy = energy0

   call initDefaults(env, calc, mol, 0)
   call check(error, allocated(set%opt_engine))
   if (.not.allocated(error)) then
      call check(error, set%opt_engine, p_engine_rf)
   end if
   if (.not.allocated(error)) then
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, .not.exitRun)
   end if
   if (.not.allocated(error)) call check(error, .not.fail)
   if (.not.allocated(error)) then
      if (energy >= energy0) then
         call test_failed(error, "Nonperiodic optimization did not lower energy")
      end if
   end if

   if (allocated(engine_save)) then
      set%opt_engine = engine_save
   else if (allocated(set%opt_engine)) then
      deallocate(set%opt_engine)
   end if
end subroutine test_nonperiodic_default_optimizer


subroutine test_explicit_periodic_optimizer(error)
   use xtb_main_defaults, only : initDefaults
   use xtb_geoopt, only : geometry_optimization
   use xtb_setparam, only : set, p_engine_rf, p_engine_lbfgs, p_engine_inertial

   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   real(wp) :: energy, energy0, egap, gradient(3, 3), sigma(3, 3)
   real(wp) :: xyz0(3, 3), lattice0(3, 3)
   logical :: fail, exitRun, optcell_save
   integer :: iter_needed
   integer, allocatable :: engine_save
   character(len=:), allocatable :: err_msg

   if (allocated(set%opt_engine)) engine_save = set%opt_engine
   optcell_save = set%optcell
   set%optcell = .true.

   ! 1. Incompatible explicit optimizer (rf / ancopt) must be rejected for variable-cell PBC
   call setup_strained_caf2(env, mol, chk, calc, energy, egap, gradient, sigma)
   call fixset%deallocate
   xyz0 = mol%xyz
   lattice0 = mol%lattice

   set%opt_engine = p_engine_rf
   call initDefaults(env, calc, mol, 0)
   call check(error, set%opt_engine, p_engine_rf)
   if (.not.allocated(error)) then
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, fail, .true.)
   end if
   if (.not.allocated(error)) call check(error, exitRun, .true.)
   if (.not.allocated(error)) then
      call env%getLog(err_msg)
      if (index(err_msg, "does not support periodic geometry optimization") == 0) then
         call test_failed(error, "Missing error message for incompatible periodic optimizer")
      end if
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%xyz - xyz0)), 0.0_wp, thr=1.0e-12_wp)
   end if
   if (.not.allocated(error)) then
      call check(error, maxval(abs(mol%lattice - lattice0)), 0.0_wp, thr=1.0e-12_wp)
   end if

   ! 2. Incompatible explicit optimizer (lbfgs / l_ancopt) must also be rejected
   if (.not.allocated(error)) then
      call init(env)
      set%opt_engine = p_engine_lbfgs
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, fail, .true.)
      if (.not.allocated(error)) call check(error, exitRun, .true.)
      call env%getLog(err_msg)
   end if

   ! 3. Incompatible explicit optimizers (rf and lbfgs) must also be rejected for fixed-cell PBC (--nocellopt)
   if (.not.allocated(error)) then
      set%optcell = .false.
      call init(env)
      set%opt_engine = p_engine_rf
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, fail, .true.)
      if (.not.allocated(error)) call check(error, exitRun, .true.)
      call env%getLog(err_msg)
      if (.not.allocated(error)) then
         if (index(err_msg, "does not support periodic geometry optimization") == 0) then
            call test_failed(error, "rf was not rejected for fixed-cell periodic optimization")
         end if
      end if
   end if
   if (.not.allocated(error)) then
      call init(env)
      set%opt_engine = p_engine_lbfgs
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, fail, .true.)
      if (.not.allocated(error)) call check(error, exitRun, .true.)
      call env%getLog(err_msg)
      set%optcell = .true.
   end if

   ! 4. Supported explicit periodic optimizer (inertial / FIRE) must be preserved and succeed
   if (.not.allocated(error)) then
      call setup_strained_caf2(env, mol, chk, calc, energy, egap, gradient, sigma)
      call fixset%deallocate
      lattice0 = mol%lattice
      energy0 = energy
      set%opt_engine = p_engine_inertial
      call initDefaults(env, calc, mol, 0)
      call check(error, set%opt_engine, p_engine_inertial)
   end if
   if (.not.allocated(error)) then
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy, gradient, sigma, p_olev_normal, .false., .false., fail, iter_needed)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, .not.exitRun)
      if (.not.allocated(error)) call check(error, .not.fail)
      if (.not.allocated(error)) then
         if (maxval(abs(mol%lattice - lattice0)) < 1.0e-4_wp) then
            call test_failed(error, "Explicitly selected inertial optimizer did not relax lattice")
         end if
      end if
   end if

   if (allocated(engine_save)) then
      set%opt_engine = engine_save
   else if (allocated(set%opt_engine)) then
      deallocate(set%opt_engine)
   end if
   set%optcell = optcell_save
end subroutine test_explicit_periodic_optimizer


subroutine test_periodic_optlevel_and_failures(error)
   use xtb_geoopt, only : geometry_optimization
   use xtb_setparam, only : set, p_engine_pbc_lbfgs, p_olev_crude, p_olev_vtight
   use xtb_pbc_tools, only : dlat_to_cell, dlat_to_dvol

   type(error_type), allocatable, intent(out) :: error
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res
   real(wp) :: energy_crude, energy_vtight, energy_maxc, energy_ref
   real(wp) :: egap, egap_ref, gradient(3, 3), grad_ref(3, 3), sigma(3, 3), sigma_ref(3, 3)
   real(wp) :: gnorm_crude, gnorm_vtight, xyz0(3, 3), lattice0(3, 3), cellpar_ref(6), vol_ref
   logical :: fail, exitRun, optcell_save
   integer :: iter_crude, iter_vtight, iter_maxc, iter_fail
   integer :: no_fixed(0)
   integer, allocatable :: engine_save
   character(len=:), allocatable :: err_msg

   if (allocated(set%opt_engine)) engine_save = set%opt_engine
   optcell_save = set%optcell
   set%opt_engine = p_engine_pbc_lbfgs
   set%optcell = .false.

   ! 1. Verify optlevel propagation: vtight must converge to a tighter gradient norm than crude
   call setup_periodic_water(env, mol, chk, calc, energy_crude, egap, gradient, sigma, no_fixed)
   call fixset%deallocate
   call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
      & energy_crude, gradient, sigma, p_olev_crude, .false., .false., fail, iter_crude)
   call cleanup_opt_files()
   call env%check(exitRun)
   call check(error, .not.exitRun)
   if (.not.allocated(error)) call check(error, .not.fail)
   gnorm_crude = norm2(gradient)

   if (.not.allocated(error)) then
      call setup_periodic_water(env, mol, chk, calc, energy_vtight, egap, gradient, sigma, no_fixed)
      call fixset%deallocate
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 200, &
         & energy_vtight, gradient, sigma, p_olev_vtight, .false., .false., fail, iter_vtight)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, .not.exitRun)
      if (.not.allocated(error)) call check(error, .not.fail)
      gnorm_vtight = norm2(gradient)
   end if

   if (.not.allocated(error)) then
      if (iter_vtight <= iter_crude) then
         call test_failed(error, "vtight optimization did not require more cycles than crude")
      else if (gnorm_vtight > 2.0e-4_wp) then
         call test_failed(error, "vtight nuclear gradient norm exceeds vtight threshold")
      else if (gnorm_vtight >= gnorm_crude) then
         call test_failed(error, "vtight optimization did not achieve smaller gradient norm than crude")
      else if (energy_vtight >= energy_crude) then
         call test_failed(error, "vtight optimization did not achieve lower energy than crude")
      end if
   end if

   ! 2. Verify maximum-cycle termination (maxcycle_in = 2):
   !    fail must be .true., exitRun must be .false. (rescued), iter_needed == 2,
   !    and mol/energy/gradient/sigma must match the last evaluated structure.
   if (.not.allocated(error)) then
      set%optcell = .true.
      call setup_strained_caf2(env, mol, chk, calc, energy_maxc, egap, gradient, sigma)
      call fixset%deallocate
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 2, &
         & energy_maxc, gradient, sigma, p_olev_vtight, .false., .false., fail, iter_maxc)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, .not.exitRun)
      if (.not.allocated(error)) call check(error, fail, .true.)
      if (.not.allocated(error)) call check(error, iter_maxc, 2)
   end if
   if (.not.allocated(error)) then
      call dlat_to_cell(mol%lattice, cellpar_ref)
      vol_ref = dlat_to_dvol(mol%lattice)
      call check(error, maxval(abs(mol%cellpar - cellpar_ref)), 0.0_wp, thr=1.0e-12_wp)
      if (.not.allocated(error)) then
         call check(error, mol%volume, vol_ref, thr=1.0e-12_wp)
      end if
   end if
   if (.not.allocated(error)) then
      call calc%singlepoint(env, mol, chk, 0, .false., energy_ref, grad_ref, sigma_ref, egap_ref, res)
      call check(error, energy_maxc, energy_ref, thr=1.0e-10_wp)
      if (.not.allocated(error)) then
         call check(error, maxval(abs(gradient - grad_ref)), 0.0_wp, thr=1.0e-10_wp)
      end if
      if (.not.allocated(error)) then
         call check(error, maxval(abs(sigma - sigma_ref)), 0.0_wp, thr=1.0e-10_wp)
      end if
   end if

   ! 3. Verify failed single-point evaluation during periodic optimization:
   !    Force singlepoint failure by setting calc%xtbData%level = 2 on a periodic structure.
   !    Outputs must be cleanly zeroed (not left as uninitialized garbage), mol restored,
   !    fail == .true., exitRun == .true., and final_sp skipped.
   if (.not.allocated(error)) then
      call setup_strained_caf2(env, mol, chk, calc, energy_maxc, egap, gradient, sigma)
      call fixset%deallocate
      xyz0 = mol%xyz
      lattice0 = mol%lattice
      calc%xtbData%level = 2
      energy_maxc = -999.0_wp
      gradient = -999.0_wp
      sigma = -999.0_wp
      egap = -999.0_wp
      call geometry_optimization(env, mol, chk, calc, egap, 300.0_wp, 250, 10, &
         & energy_maxc, gradient, sigma, p_olev_normal, .false., .true., fail, iter_fail)
      call cleanup_opt_files()
      call env%check(exitRun)
      call check(error, fail, .true.)
      if (.not.allocated(error)) call check(error, .not.exitRun)
      if (.not.allocated(error)) call check(error, iter_fail, 0)
      if (.not.allocated(error)) call check(error, energy_maxc, 0.0_wp, thr=1.0e-12_wp)
      if (.not.allocated(error)) call check(error, maxval(abs(gradient)), 0.0_wp, thr=1.0e-12_wp)
      if (.not.allocated(error)) call check(error, maxval(abs(sigma)), 0.0_wp, thr=1.0e-12_wp)
      if (.not.allocated(error)) call check(error, maxval(abs(mol%xyz - xyz0)), 0.0_wp, thr=1.0e-12_wp)
      if (.not.allocated(error)) call check(error, maxval(abs(mol%lattice - lattice0)), 0.0_wp, thr=1.0e-12_wp)
      call env%getLog(err_msg)
   end if

   if (allocated(engine_save)) then
      set%opt_engine = engine_save
   else if (allocated(set%opt_engine)) then
      deallocate(set%opt_engine)
   end if
   set%optcell = optcell_save
end subroutine test_periodic_optlevel_and_failures

end module test_optimizer
