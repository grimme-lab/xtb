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
   use mctc_io_constants, only : codata
   use mctc_io_convert, only : aatoau, ctoau, autoev
#if WITH_TBLITE
   use tblite_basis_type, only : basis_type
   use tblite_ceh_ceh, only : new_ceh_calculator
   use tblite_container, only : container_type
   use tblite_context, only : context_type, context_terminal, escape
   use tblite_external_field, only : electric_field
   use tblite_param, only : param_record
   use tblite_post_processing_list, only : add_post_processing, post_processing_list
   use tblite_results, only : results_type
   use tblite_spin, only : spin_polarization, new_spin_polarization
   use tblite_solvation, only : solvation_type, solvation_input, new_solvation, &
      & new_solvation_cds, new_solvation_shift, alpb_input, cds_input, &
      & shift_input, cpcm_input, born_kernel, solution_state, solvent_data, &
      & get_solvent_data
   use tblite_wavefunction, only : wavefunction_type, new_wavefunction, &
      & sad_guess, eeq_guess, eeqbc_guess
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

   public :: TTBLiteCalculator, TTBLiteInput, TTBLiteSolvationInput
   public :: newTBLiteCalculator, newTBLiteWavefunction
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
      !> Maximum number of cycles for SCC convergence
      integer :: max_iter = 250
      !> Parameter file to read
      character(len=:), allocatable :: param
      !> Method used for computation
      character(len=:), allocatable :: method
      !> Colorful output
      logical :: color = .false.
      !> CEH charges
      type(ceh), allocatable :: ceh
      !> Electric field in cartesian coordinates
      real(wp), allocatable :: efield(:)
      !> Solvation model input
      type(TTBLiteSolvationInput), allocatable :: solvation
   end type TTBLiteInput

   !> Input for the solvation model in the tblite library
   type :: TTBLiteSolvationInput
      !> Solvation model name
      character(len=:), allocatable :: solvation_model
      !> Solvent name
      character(len=:), allocatable :: solvent
      !> Solvation reference state
      character(len=:), allocatable :: reference_state
   end type TTBLiteSolvationInput

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
      !> Shift for IP/EA calculations
      real(wp) :: ipeashift
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

   !> Convert J to atomic units
   real(wp), parameter :: jtoau = 1.0_wp / (codata%me*codata%c**2*codata%alpha**2)
   !> Convert V/Å = J/(C·Å) to atomic units
   real(wp), parameter :: vatoau = jtoau / (ctoau * aatoau)

   !> Global IP/EA shift for all models in Hartree
   real(wp), parameter :: global_ipeashift = 0.17806900_wp

   character(len=*),private,parameter :: outfmt = &
      '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

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
      method = calc%tblite%method
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

   ! The maximum number of iterations is a property of the tblite calculator 
   ! itself and has to be set AFTER the calculator is created
   calc%tblite%max_iter = input%max_iter

   if (allocated(input%efield)) then
      block
         class(container_type), allocatable :: cont
         cont = electric_field(input%efield*vatoau)
         call calc%tblite%push_back(cont)
      end block
   end if

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

   if (allocated(input%solvation)) then 
      block 
         type(solvation_input), allocatable :: solv_input
         ! Construct a tblite solvation input from the provided data
         allocate(solv_input)
         call construct_solv_input(input%solvation, solv_input, error)
         if (allocated(error)) then
            call env%error(error%message, source)
            return
         end if

         block
            class(container_type), allocatable :: cont
            class(solvation_type), allocatable :: solv
            call new_solvation(solv, struc, solv_input, error, method)
            if (allocated(error)) then
               call env%error(error%message, source)
               return
            end if
            call move_alloc(solv, cont)
            call calc%tblite%push_back(cont)
         end block
         if (allocated(solv_input%cds)) then
            block
               class(container_type), allocatable :: cont
               class(solvation_type), allocatable :: cds
               call new_solvation_cds(cds, struc, solv_input, error, method)
               if (allocated(error)) then
                  call env%error(error%message, source)
                  return
               end if
               call move_alloc(cds, cont)
               call calc%tblite%push_back(cont)
            end block
         end if
         if (allocated(solv_input%shift)) then
            block
               class(container_type), allocatable :: cont
               class(solvation_type), allocatable :: shift
               call new_solvation_shift(shift, solv_input, error, method)
               if (allocated(error)) then
                  call env%error(error%message, source)
                  return
               end if
               call move_alloc(shift, cont)
               call calc%tblite%push_back(cont)
            end block
         end if
      end block
   end if

   ! Set the IP/EA shift globally for all models
   calc%ipeashift = global_ipeashift

#else
    call feature_not_implemented(env)
#endif
end subroutine newTBLiteCalculator


#if WITH_TBLITE
subroutine construct_solv_input(input, solv_input, error)
   !> Source of the generated errors
   character(len=*), parameter :: source = 'tblite_calculator_constructsolvinput'
   !> Source of the solvation input from the xtb input
   type(TTBLiteSolvationInput), intent(in) :: input
   !> Constructed solvation input for tblite
   type(solvation_input), intent(out) :: solv_input
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   logical :: parametrized_solvation, alpb
   type(solvent_data), allocatable :: solvent
   integer :: sol_state, kernel, iostat
   real(wp) :: val

   ! Collect solvent data if a solvent name is provided
   if (allocated(input%solvent)) then
      parametrized_solvation = .true.
      allocate(solvent)
      solvent = get_solvent_data(input%solvent)
      ! If the solvent name is not recognized, we try interpreting it as a dielectric constant
      if (solvent%eps <= 0.0_wp) then
         parametrized_solvation = .false.
         read(input%solvent,*,iostat=iostat) solvent%eps
         if (iostat .ne. 0) then
            call fatal_error(error, "Invalid solvent/dielectric constant '"//input%solvent//"'")
            return
         end if
      end if
   else
      call fatal_error(error, "No solvent/dielectric constant specified")
      return
   end if

   ! Check if a solvation state is provided
   if (allocated(input%reference_state)) then
      select case(input%reference_state)
         case default
            call fatal_error(error, "Unknown solution state '"//input%reference_state//"' requested")
            return
         case("gsolv")
            sol_state = solution_state%gsolv
         case("bar1mol", "bar1M")
            sol_state = solution_state%bar1mol
         case("reference")
            sol_state = solution_state%reference
      end select
   else
      sol_state = solution_state%gsolv
   end if

   ! Select the solvation model and construct the input
   if (allocated(input%solvation_model)) then
      select case(input%solvation_model)
      case default
         call fatal_error(error, "Unknown solvation model '"//input%solvation_model//"' requested")
         return
      case("gbsa", "alpb")
         ! Check if ALPB or GBSA is requested
         alpb = .false.
         if (input%solvation_model == "alpb") then
            alpb = .true.
         end if

         ! Select the default born kernel
         kernel = merge(born_kernel%still, born_kernel%p16, alpb)

         ! Construct parametrized solvation model input for tblite
         if (parametrized_solvation) then
            solv_input%alpb = alpb_input(solvent%eps, solvent=solvent%solvent, &
               & kernel=kernel, alpb=alpb)
            solv_input%cds = cds_input(alpb=alpb, solvent=solvent%solvent)
            solv_input%shift = shift_input(alpb=alpb, solvent=solvent%solvent, &
               & state=sol_state)
         else
            call fatal_error(error, "ALPB/GBSA solvation models require a named solvent.")
            return
         end if
      case("gbe", "gb")
         ! Check if GBE or GB is requested
         alpb = .false.
         if (input%solvation_model == "gbe") then
            alpb = .true.
         end if

         ! Select the default born kernel
         kernel = merge(born_kernel%still, born_kernel%p16, alpb)

         if (.not.parametrized_solvation .and. sol_state /= solution_state%gsolv) then
            call fatal_error(error, "Solution state shift is only supported for named solvents")
            return
         end if

         ! Construct purely electrostatic solvation model input for tblite
         solv_input%alpb = alpb_input(solvent%eps, kernel=kernel, alpb=alpb)
      case("cosmo")
         ! ddCOSMO solvation model (currently called CPCM in tblite)
         if (sol_state /= solution_state%gsolv) then 
            call fatal_error(error, "Solution state shift not supported for ddCOSMO")
            return
         end if
         solv_input%cpcm = cpcm_input(solvent%eps)
      end select
   else
      call fatal_error(error, "No solvation model specified")
   end if

end subroutine construct_solv_input
# endif


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
         call eeq_guess(struc, calc%tblite, wfn, error)
      case("eeqbc")
         call eeqbc_guess(struc, calc%tblite, wfn, error)
      case("ceh")
         block 
            use tblite_context, only : context_type, context_terminal
            use tblite_context_terminal, only : escape
            use tblite_ceh_singlepoint, only : ceh_singlepoint
            use tblite_lapack_solver, only : lapack_solver 
            use tblite_lapack_solver, only : lapack_algorithm
            type(context_type) :: ctx

            ctx%solver = lapack_solver(lapack_algorithm%gvd)

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
   integer :: spin, charge, stat, unit, nspin, i
   logical :: exist
   type(error_type), allocatable :: error
   type(context_type) :: ctx
   type(post_processing_list), allocatable :: post_proc
   real(wp) :: efix
   real(wp), allocatable :: dpmom(:), qpmom(:)
   character(len=:), allocatable :: wbo_label, molmom_label

   struc = mol
   ctx%unit = env%unit
   ctx%terminal = context_terminal(self%color)

   ! Setup the required post-processing
   allocate(post_proc)
   ! Wiberg-Mayer bond orders
   wbo_label = "bond-orders"
   call add_post_processing(post_proc, wbo_label, error)
   if (allocated(error)) then
      call env%error(error%message, source)
      return
   end if

   ! Molecular multipole moments
   molmom_label = "molmom"
   call add_post_processing(post_proc, molmom_label, error)
   if (allocated(error)) then
      call env%error(error%message, source)
      return
   end if

   ! Needed to update atomic charges after reading restart file
   call get_qat_from_qsh(self%tblite%bas, chk%tblite%qsh, chk%tblite%qat)

   call xtb_singlepoint(ctx, struc, self%tblite, chk%tblite, self%accuracy, &
      & energy, gradient, sigma, printlevel, results%tblite_results, post_proc)
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

   ! ------------------------------------------------------------------------
   !  various external potentials
   efix = 0.0_wp
   call constrain_pot(potset,mol%n,mol%at,mol%xyz,gradient,efix)
   call constrpot   (mol%n,mol%at,mol%xyz,gradient,efix)
   call cavity_egrad(mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (metaset,mol%n,mol%at,mol%xyz,efix,gradient)
   call metadynamic (rmsdset,mol%n,mol%at,mol%xyz,efix,gradient)

   ! ------------------------------------------------------------------------
   !  fixing of certain atoms
   energy = energy + efix
   results%e_total = energy
   results%gnorm = norm2(gradient)
   if (fixset%n.gt.0) then
      do i=1, fixset%n
         gradient(1:3,fixset%atoms(i))=0
      enddo
   endif

   if (printlevel.ge.2) then
      ! start with summary header
      if (.not.set%silent) then
         write(env%unit,'(9x,53(":"))')
         write(env%unit,'(9x,"::",21x,a,21x,"::")') "SUMMARY"
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,outfmt) "total energy      ", results%e_total,"Eh   "
      write(env%unit,outfmt) "gradient norm     ", results%gnorm,  "Eh/a0"
      write(env%unit,outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
      if (.not.set%silent) then
         if (set%verbose) then
            write(env%unit,'(9x,"::",49("."),"::")')
            write(env%unit,outfmt) "HOMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo),  "eV   "
            write(env%unit,outfmt) "LUMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo+1),"eV   "
         endif
         write(env%unit,'(9x,"::",49("."),"::")')
         select case(self%tblite%method)
         case default
            call fatal_error(error, "Unknown method '"//self%tblite%method//"'")
         case("gfn2")
            write(env%unit,outfmt) "GFN2-xTB energy   ", energy - efix,       "Eh   "
         case("gfn1")
            write(env%unit,outfmt) "GFN1-xTB energy   ", energy - efix,       "Eh   "
         case("ipea1")
            write(env%unit,outfmt) "IPEA1-xTB energy  ", energy - efix,       "Eh   "
         end select
         write(env%unit,outfmt) "add. restraining  ", efix,       "Eh   "
         write(env%unit,outfmt) "total charge      ", sum(chk%wfn%q), "e    "
      endif
      write(env%unit,'(9x,53(":"))')
      write(env%unit,'(a)')
   endif

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

subroutine print_charges(ifile, n, q)
   implicit none
   integer, intent(in) :: ifile
   integer, intent(in) :: n
   real(wp), intent(in) :: q(n)
   integer :: i
   if (ifile /= -1) then
      do i = 1, n
         write (ifile, '(f14.8)') q(i)
      end do
   end if
end subroutine print_charges

!> get numerical gradients for charges
subroutine num_grad_chrg(env, mol, tblite)
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> Molecular structure data
   type(TMolecule), intent(in) :: mol
   
   !> step size for numerical gradients
   type(TTBLiteInput), intent(in) :: tblite
   
   !> numerical gradients
   real(wp) :: numgrad(3, mol%n, mol%n)
   
   real(wp), allocatable, dimension(:) :: cehr, cehl
   type(TMolecule) :: local_mol
   integer :: i, j, k, ich

   real(wp) :: step, step2 ! for numerical gradient
   
   numgrad=0.0_wp
   step = tblite%ceh%step
   step2 = 0.5_wp / step
   call get_ceh(env,mol,tblite)
   
   !$omp parallel do private(j, cehr, cehl, local_mol) shared(env, numgrad, mol, tblite, step, step2) default(none)
   do i = 1, mol%n
      call local_mol%copy(mol) ! create a thread-private copy of mol
      do j = 1, 3
         local_mol%xyz(j,i) = local_mol%xyz(j,i) + step
         call get_ceh(env, local_mol, tblite, cehr)

         local_mol%xyz(j,i) = local_mol%xyz(j,i) - 2*step
         call get_ceh(env, local_mol, tblite, cehl)
         
         numgrad(j,i,:) = step2 * (cehr - cehl)  ! numerical gradient
         local_mol%xyz(j,i) = local_mol%xyz(j,i) + step  ! reset the coordinates
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
