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
   use tblite_ceh_ceh, only : new_ceh_calculator
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
   use tblite_data_spin, only : get_spin_constant 
   use xtb_tblite_mapping, only : convert_tblite_to_wfn
#endif
   use xtb_tblite_mapping, only : convert_tblite_to_results
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
   public :: get_ceh, num_grad_chrg

   type, private :: ceh
      !> numerical gradients
      logical :: grad
      !> step size for numerical gradients
      real(wp) :: step = 0.00001_wp
   endtype

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
      !> CEH charges
      type(ceh), allocatable :: ceh
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
      !> Perform a spin-polarized calculation
      logical :: spin_polarized = .false.
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


!> Create a new instance of tblite based xTB calculator
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
   calc%spin_polarized = input%spin_polarized

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
         call new_gfn2_calculator(calc%tblite, struc, error)
      case("gfn1")
         call new_gfn1_calculator(calc%tblite, struc, error)
      case("ipea1")
         call new_ipea1_calculator(calc%tblite, struc, error)
      case("ceh")
         calc%guess = method
         calc%nspin = 1
         calc%etemp = 4000.0_wp * kt
         call new_ceh_calculator(calc%tblite, struc, error)
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
   if (calc%spin_polarized) then
      block
         class(container_type), allocatable :: cont
         type(spin_polarization), allocatable :: spin
         real(wp), allocatable :: wll(:, :, :)
         allocate(spin)
         call get_spin_constants(wll, struc, calc%tblite%bas)
         call new_spin_polarization(spin, struc, wll, calc%tblite%bas%nsh_id)
         call move_alloc(spin, cont)
         call calc%tblite%push_back(cont)
      end block
   end if
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
   type(TTBLiteCalculator), intent(inout) :: calc
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
      case("ceh")
         block 
            use tblite_context, only : context_type, context_terminal
            use tblite_context_terminal, only : escape
            use tblite_ceh_singlepoint, only : ceh_singlepoint
            use tblite_lapack_solver, only : lapack_solver 
            use tblite_lapack_solver, only : lapack_algorithm
            type(context_type) :: ctx

            ctx%solver = lapack_solver(lapack_algorithm%gvd)

            ! temporary turn off colorful output
            ! ctx%terminal = context_terminal(calc%color)
            ! write (env%unit, '(0x,a)') escape(ctx%terminal%cyan) // "Calculation of CEH charges" // &
            !    & escape(ctx%terminal%reset)

            call ceh_singlepoint(ctx, calc%tblite, struc, wfn, calc%accuracy, 1)
         end block
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

   ! convert tblite results into xtb data !
   call convert_tblite_to_wfn(env, self%tblite%bas, mol, chk)
   call convert_tblite_to_results(results,mol,chk,energy,.true.,gradient=gradient)
   hlgap = results%hl_gap

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

#if WITH_TBLITE
subroutine get_spin_constants(wll, mol, bas)
   real(wp), allocatable, intent(out) :: wll(:, :, :)
   type(structure_type), intent(in) :: mol
   type(basis_type), intent(in) :: bas

   integer :: izp, ish, jsh, il, jl

   allocate(wll(bas%nsh, bas%nsh, mol%nid), source=0.0_wp)

   do izp = 1, mol%nid
      do ish = 1, bas%nsh_id(izp)
         il = bas%cgto(ish, izp)%ang
         do jsh = 1, bas%nsh_id(izp)
            jl = bas%cgto(jsh, izp)%ang
            wll(jsh, ish, izp) = get_spin_constant(jl, il, mol%num(izp))
         end do
      end do
   end do
end subroutine get_spin_constants
#endif

!> get CEH charges via tblite
subroutine get_ceh(env,mol,tblite, ceh_chrg)

   use xtb_propertyoutput, only : print_charges

   !> computational environment
   type(TEnvironment), intent(inout) :: env

   !> molecular structure data
   type(TMolecule), intent(in) :: mol

   !> tblite input
   type(TTBLiteInput), intent(in) :: tblite

   !> CEH charges, if present assumed the numerical gradients are requested
   real(wp), allocatable, optional, intent(out) :: ceh_chrg(:)

   !> initialize temporary calculator for CEH
   type(TTBLiteCalculator) :: calc_ceh
            
   !> initialize temporary input for CEH
   type(TTBLiteInput) :: tblite_ceh
   
   !> initialize temporary wfn for CEH
   type(TRestart) :: chk_ceh

   integer :: ich

   tblite_ceh = tblite        ! copy the tblite input
   tblite_ceh%method = "ceh" 
#if WITH_TBLITE

   if (.not. present(ceh_chrg)) &
      write(env%unit, '(1x, a, /, a)') "Calculation of CEH charges",repeat('-', 36)
   
   call newTBLiteCalculator(env, mol, calc_ceh, tblite_ceh)
   call newTBLiteWavefunction(env, mol, calc_ceh, chk_ceh)
   
   if (present(ceh_chrg)) then
      allocate(ceh_chrg(mol%n))
      ceh_chrg = chk_ceh%tblite%qat(:,1)  
   else
      ! create ceh.charges file !
      call open_file(ich, 'ceh.charges', 'w') 
      call print_charges(ich, mol%n, chk_ceh%tblite%qat(:,1))
      call close_file(ich)
      write(env%unit, '(1x, a)') "CEH charges written to ceh.charges"
   endif

#else
    call feature_not_implemented(env)
#endif

end subroutine get_ceh

!> get numerical gradients for charges
subroutine num_grad_chrg(env, mol, tblite)
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> Molecular structure data
   type(TMolecule), intent(inout) :: mol
   
   !> step size for numerical gradients
   type(TTBLiteInput), intent(in) :: tblite
   
   !> numerical gradients
   real(wp) :: numgrad(3, mol%n, mol%n)
   
   real(wp), allocatable, dimension(:) :: cehr, cehl
   !
   integer :: i,j,k, ich

   real(wp) :: step, step2 ! for numerical gradient
   
   numgrad=0.0_wp
   step = tblite%ceh%step
   step2 = 0.5_wp / step
   call get_ceh(env,mol,tblite)
   
   !$omp parallel do private(j,cehr,cehl) shared(env, numgrad, mol, tblite, step, step2) 
   do i = 1, mol%n
      do j = 1, 3
         mol%xyz(j,i) = mol%xyz(j,i) + step
         call get_ceh(env, mol, tblite, cehr)

         mol%xyz(j,i) = mol%xyz(j,i) - 2*step
         call get_ceh(env, mol, tblite, cehl)
         
         numgrad(j,i,:) = step2 * (cehr - cehl) ! numerical gradient
         mol%xyz(j,i) = mol%xyz(j,i) + step ! reset the coordinates
         
      enddo
   enddo
   !$omp end parallel do

   ! write the numerical gradient to the ceh.charges.numgrad file
   call open_file(ich, 'ceh.charges.numgrad', 'w')
   do i = 1, mol%n
      do j = 1, mol%n
         do k = 1, 3
            write(ich, '(3f12.6)') numgrad(k,j,i)
         enddo
      enddo
   enddo
   call close_file(ich)
   write(env%unit, '(1x, a)') "CEH gradients written to ceh.charges.numgrad"

 end subroutine num_grad_chrg
 


#if ! WITH_TBLITE
subroutine feature_not_implemented(env)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   call env%error("Compiled without support for tblite library")
end subroutine feature_not_implemented
#endif

end module xtb_tblite_calculator
