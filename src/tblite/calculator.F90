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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Proxy module to integrate the [tblite] tight-binding framework with xtb
!>
!> [tblite]: https://tblite.readthedocs.io
module xtb_tblite_calculator
   use mctc_env, only : error_type, fatal_error
   use mctc_io, only : structure_type, read_structure, filetype
#if WITH_TBLITE
   use tblite_basis_type, only : basis_type
   use tblite_container, only : container_type
   use tblite_context, only : context_type, context_terminal, escape
   use tblite_external_field, only : electric_field
   use tblite_param, only : param_record
   use tblite_results, only : results_type
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_solvation, only : new_solvation, solvation_type
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, sad_guess, eeq_guess
   use tblite_xtb_calculator, only : xtb_calculator, new_xtb_calculator
   use tblite_xtb_gfn2, only : new_gfn2_calculator, export_gfn2_param
   use tblite_xtb_gfn1, only : new_gfn1_calculator, export_gfn1_param
   use tblite_xtb_ipea1, only : new_ipea1_calculator, export_ipea1_param
   use tblite_xtb_singlepoint, only : xtb_singlepoint
#endif
   use xtb_mctc_accuracy, only : wp
   use xtb_type_calculator, only : TCalculator
   use xtb_type_data
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, assignment(=)
   use xtb_type_restart, only : TRestart
   use xtb_setparam
   use xtb_fixparam
   use xtb_scanparam
   use xtb_sphereparam
   use xtb_metadynamic
   use xtb_constrainpot
   implicit none
   private

   public :: TTBLiteCalculator, TTBLiteInput, newTBLiteCalculator, newTBLiteWavefunction

   !> Input for tblite library
   type :: TTBLiteInput
      !> Perform a spin-polarized calculation
      logical :: spin_polarized = .false.
      !> Electronic temperature in Kelvin
      real(wp) :: etemp = 300.0_wp
      !> Numerical accuracy
      real(wp) :: accuracy = 1.0_wp
      !> Parameter file to read
      character(len=:), allocatable :: param
      !> Method used for computation
      character(len=:), allocatable :: method
      !> Colorful output
      logical :: color = .false.
   end type TTBLiteInput

   !> Calculator interface for xTB based methods
   type, extends(TCalculator) :: TTBLiteCalculator
      !> Number of spin channels
      integer :: nspin
      !> Electronic temperature in Hartree
      real(wp) :: etemp
      !> Guess for initial wavefunction
      character(len=:), allocatable :: guess
      !> Colorful output
      logical :: color
#if WITH_TBLITE
      !> Instance of tblite calculator
      type(xtb_calculator) :: tblite
#endif
   contains
      !> Perform calculator for a singlepoint
      procedure :: singlepoint
      !> Write information instance
      procedure :: writeInfo
   end type TTBLiteCalculator


   !> Conversion factor from temperature to energy
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

contains


!> Create a new instance of an tblite based xTB calculator
subroutine newTBLiteCalculator(env, mol, calc, input)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'tblite_calculator_newTBLiteCalculator'
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Instance of the new calculator
   type(TTBLiteCalculator), intent(out) :: calc
   !> Input for calculator settings
   type(TTBLiteInput), intent(in) :: input

#if WITH_TBLITE
   character(len=:), allocatable :: method
   type(error_type), allocatable :: error
   type(structure_type) :: struc
   type(param_record) :: param

   struc = mol

   calc%nspin = merge(2, 1, input%spin_polarized)
   calc%etemp = input%etemp * kt
   calc%guess = "sad"
   calc%accuracy = input%accuracy
   calc%color = input%color

   if (allocated(input%param)) then
      call param%load(input%param, error)
      if (.not. allocated(error)) then
         call new_xtb_calculator(calc%tblite, struc, param, error)
      end if
   else
      method = "gfn2"
      if (allocated(input%method)) method = input%method
      select case(method)
      case default
         call fatal_error(error, "Unknown method '"//method//"' requested")
      case("gfn2")
         call new_gfn2_calculator(calc%tblite, struc)
      case("gfn1")
         call new_gfn1_calculator(calc%tblite, struc)
      case("ipea1")
         call new_ipea1_calculator(calc%tblite, struc)
      end select
   end if
   if (allocated(error)) then
      call env%error(error%message, source)
      return
   end if

!   if (allocated(config%efield)) then
!      block
!         class(container_type), allocatable :: cont
!         cont = electric_field(config%efield*vatoau)
!         call calc%tblite%push_back(cont)
!      end block
!   end if
!
!   if (config%spin_polarized) then
!      block
!         class(container_type), allocatable :: cont
!         type(spin_polarization), allocatable :: spin
!         real(wp), allocatable :: wll(:, :, :)
!         allocate(spin)
!         ! TODO: get_spin_constants(wll, struc, calc%tblite%bas)
!         call new_spin_polarization(spin, struc, wll, calc%tblite%bas%nsh_id)
!         call move_alloc(spin, cont)
!         call calc%tblite%push_back(cont)
!      end block
!   end if
!
!   if (allocated(config%solvation)) then
!      block
!         class(container_type), allocatable :: cont
!         class(solvation_type), allocatable :: solv
!         call new_solvation(solv, struc, config%solvation, error)
!         if (.not.allocated(error)) then
!            call move_alloc(solv, cont)
!            call calc%tblite%push_back(cont)
!         end if
!      end block
!   end if
!
!   if (allocated(error)) then
!      call env%error(error%message, source)
!      return
!   end if
#else
    call feature_not_implemented(env)
#endif
end subroutine newTBLiteCalculator


!> Create new wavefunction restart data for tblite library
subroutine newTBLiteWavefunction(env, mol, calc, chk)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'tblite_calculator_newTBLiteWavefunction'
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   !> Instance of the new calculator
   type(TTBLiteCalculator), intent(in) :: calc
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk

#if WITH_TBLITE
   type(structure_type) :: struc
   type(error_type), allocatable :: error

   struc = mol

   associate(wfn => chk%tblite)
      call new_wavefunction(wfn, struc%nat, calc%tblite%bas%nsh, calc%tblite%bas%nao, &
         & calc%nspin, calc%etemp)

      select case(calc%guess)
      case default
         call fatal_error(error, "Unknown starting guess '"//calc%guess//"' requested")
      case("sad")
         call sad_guess(struc, calc%tblite, wfn)
      case("eeq")
         call eeq_guess(struc, calc%tblite, wfn)
      end select
   end associate
   if (allocated(error)) then
      call env%error(error%message, source)
      return
   end if
#else
    call feature_not_implemented(env)
#endif
end subroutine newTBLiteWavefunction


!> Perform calculation for single point
subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
      & energy, gradient, sigma, hlgap, results)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'tblite_calculator_singlepoint'
   !> Calculator instance
   class(TTBLiteCalculator), intent(inout) :: self
   !> Computational environment
   type(TEnvironment), intent(inout) :: env
   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol
   !> Wavefunction data
   type(TRestart), intent(inout) :: chk
   !> Print level for IO
   integer, intent(in) :: printlevel
   !> Restart from previous results
   logical, intent(in) :: restart
   !> Total energy
   real(wp), intent(out) :: energy
   !> Molecular gradient
   real(wp), intent(out) :: gradient(:, :)
   !> Strain derivatives
   real(wp), intent(out) :: sigma(:, :)
   !> HOMO-LUMO gap
   real(wp), intent(out) :: hlgap
   !> Detailed results
   type(scc_results), intent(out) :: results

#if WITH_TBLITE
   type(structure_type) :: struc
   integer :: spin, charge, stat, unit, nspin
   logical :: exist
   type(error_type), allocatable :: error
   type(context_type) :: ctx

   struc = mol
   ctx%unit = env%unit
   ctx%terminal = context_terminal(self%color)

   ! Needed to update atomic charges after reading restart file
   call get_qat_from_qsh(self%tblite%bas, chk%tblite%qsh, chk%tblite%qat)

   call xtb_singlepoint(ctx, struc, self%tblite, chk%tblite, self%accuracy, &
      & energy, gradient, sigma, printlevel)
   if (ctx%failed()) then
      do while(ctx%failed())
         call ctx%get_error(error)
         call env%error(error%message, source)
      end do
      return
   end if

   results%e_total = energy
   results%converged = .true.
   results%dipole = sum(chk%tblite%dpat(:, :, 1), 2) + matmul(struc%xyz, chk%tblite%qat(: ,1))
   results%gnorm = norm2(gradient)

!   if (printlevel > 2) then
!      call ascii_levels(ctx%unit, printlevel, chk%tblite%homo, chk%tblite%emo, &
!         & chk%tblite%focc, 7)
!   end if
#else
    call feature_not_implemented(env)
#endif
end subroutine singlepoint


subroutine writeInfo(self, unit, mol)
   !> Calculator instance
   class(TTBLiteCalculator), intent(in) :: self
   !> Unit for I/O
   integer, intent(in) :: unit
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

#if WITH_TBLITE
#endif
end subroutine writeInfo


#if WITH_TBLITE
subroutine get_qat_from_qsh(bas, qsh, qat)
   type(basis_type), intent(in) :: bas
   real(wp), intent(in) :: qsh(:, :)
   real(wp), intent(out) :: qat(:, :)

   integer :: ish, ispin

   qat(:, :) = 0.0_wp
   !$omp parallel do schedule(runtime) collapse(2) default(none) &
   !$omp reduction(+:qat) shared(bas, qsh) private(ish)
   do ispin = 1, size(qsh, 2)
      do ish = 1, size(qsh, 1)
         qat(bas%sh2at(ish), ispin) = qat(bas%sh2at(ish), ispin) + qsh(ish, ispin)
      end do
   end do
end subroutine get_qat_from_qsh
#endif


#if ! WITH_TBLITE
subroutine feature_not_implemented(env)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   call env%error("Compiled without support for tblite library")
end subroutine feature_not_implemented
#endif

end module xtb_tblite_calculator
