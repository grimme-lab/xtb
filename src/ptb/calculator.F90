! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
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

!> Density matrix (P) tight binding (TB) calculator
module xtb_ptb_calculator
   use xtb_type_calculator, only: TCalculator
   use xtb_type_data
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule, assignment(=)
   use xtb_type_param, only: TPTBParameter
   use xtb_type_restart, only: TRestart
   use xtb_setparam
   use xtb_fixparam
   use xtb_mctc_systools, only: rdpath

   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: autoev

#if WITH_TBLITE
   use xtb_tblite_mapping, only : convert_tblite_to_wfn, convert_tblite_to_results

   use multicharge_model, only: new_mchrg_model, mchrg_model_type

   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_lapack_solver, only: lapack_solver
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction
   use tblite_integral_type, only: integral_type
   use tblite_adjlist, only: adjacency_list

   !> PTB modules that are dependent on tblite
   use xtb_ptb_param, only: initPTB, ptbGlobals
   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_scf, only: twostepscf
   use xtb_ptb_corebasis, only: add_core_basis
   use xtb_ptb_integral_types, only: aux_integral_type
   use xtb_ptb_response, only: numgrad_polarizability
#endif
   !> PTB data type that is not dependent on tblite
   use xtb_ptb_data, only: TPTBData
   implicit none

   private

   public :: TPTBCalculator, newPTBCalculator

   !> Calculator interface for PTB method
   type, extends(TCalculator) :: TPTBCalculator

#if WITH_TBLITE
      !> PTB vDZP basis set
      type(basis_type) :: bas, cbas

      !> EEQ Model
      type(mchrg_model_type) :: eeqmodel
#endif

      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData

      !> Electronic temperature
      real(wp) :: etemp

   contains

      !> Perform PTB single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TPTBCalculator

   character(len=*), private, parameter :: outfmt = &
                                           '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'
   character(len=*), parameter :: intfmt = &
                                  '(10x,":",2x,a,i18,      10x,":")'
   character(len=*), parameter :: chrfmt = &
                                  '(10x,":",2x,a,a18,      10x,":")'
   character(len=*), parameter :: scifmt = &
                                  '(10x,":",2x,a,e22.7,1x,a,1x,":")'
   character(len=*), parameter :: dblfmt = &
                                  '(10x,":",2x,a,f18.7,5x,a,1x,":")'
   !> Conversion factor from temperature to energy
   real(wp), parameter :: kt = 3.166808578545117e-06_wp
   !> Finite-field step for numerical dipole polarizability
   real(wp), parameter :: dF = 5.0E-5_wp

contains

   subroutine newPTBCalculator(env, struc, calc, accuracy)

      character(len=*), parameter :: source = 'xtb_ptb_calculator_newPTBCalculator'

      type(TEnvironment), intent(inout) :: env

      type(TMolecule), intent(in) :: struc

      type(TPTBCalculator), intent(out) :: calc

      real(wp), intent(in), optional :: accuracy

      character(len=:), allocatable :: filename
      integer :: ich
      logical :: exist
      logical :: exitRun

#if WITH_TBLITE
      !> mctc-io structure type
      type(structure_type) :: mol

      mol = struc

      call rdpath(env%xtbpath, 'param_ptb.txt', filename, exist)
      if (.not. exist) filename = 'param_ptb.txt'

      if (present(accuracy)) then
         calc%accuracy = accuracy
      else
         calc%accuracy = 1.0_wp
      end if

      calc%etemp = set%etemp

      !> Obtain the parameter file
      allocate (calc%ptbData)
      call open_file(ich, filename, 'r')
      exist = ich /= -1
      if (exist) &
         call env%warning('Parameter file '//filename//' not supported yet.!', source)
      
      ! read the hardcoded parameters !
      call initPTB(calc%ptbData, mol%num)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not load parameters", source)
         return
      end if

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(mol, calc%bas)

      !> Add the core basis set to 'cbas' basis set type
      call add_core_basis(mol, calc%ptbData%corepotential, calc%cbas)
      
      !> set up the EEQ model
      call new_mchrg_model(calc%eeqmodel, chi=calc%ptbData%eeq%chi, &
      & rad=calc%ptbData%eeq%alp, eta=calc%ptbData%eeq%gam, kcn=calc%ptbData%eeq%cnf)

      !> check for external point charge field
      !> not implemented yet
      ! if (allocated(set%pcem_file)) then
      !    call open_file(ich, set%pcem_file, 'r')
      !    if (ich /= -1) then
      !       call read_pcem(ich, env, calc%pcem, calc%xtbData%coulomb)
      !       call close_file(ich)
      !    end if
      ! end if
#else
      call feature_not_implemented(env)
#endif

   end subroutine newPTBCalculator

   subroutine singlepoint(self, env, mol, chk, printlevel, restart, &
         & energy, gradient, sigma, hlgap, results)

      !> Source of the generated errors
      character(len=*), parameter :: source = 'ptb_calculator_singlepoint'

      !> Calculator instance
      class(TPTBCalculator), intent(inout) :: self

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
      !#################################################
      !> PTB INDIVIDUAL
      !#################################################
      type(structure_type) :: mctcmol
      !> tblite calculation context
      type(context_type) :: ctx
      !> Error type
      type(error_type), allocatable :: error
      !> Wiberg bond order
      real(wp), allocatable :: wbo(:, :, :)
      !> Static homogenoues external electric field
      real(wp), allocatable :: efield(:)
      !> Approx. effective core potential
      real(wp), allocatable :: Vecp(:, :)
      !> Auxiliary integrals
      type(aux_integral_type) :: auxints
      !> Exact integrals
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: neighborlist
      !> Effective self-energies
      real(wp), allocatable :: selfenergies(:)
      !> Electrostatic potential in second iteration
      real(wp), allocatable :: v_ES_2nditer(:)
      !> Coordination number for +U contribution
      real(wp) :: CN_plusU(mol%n)

      logical :: exitRun
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2

      energy = 0.0_wp
      gradient = 0.0_wp
      sigma = 0.0_wp

      call mol%update
      mctcmol = mol

      ctx%solver = lapack_solver(gvd)
      ctx%unit = env%unit
      ctx%verbosity = printlevel

      !> Set new PTB wavefunction in tblite format
      call newPTBWavefunction(env, self, mctcmol, chk%tblite)

      !> Static Homogeneous External Electric Field
      if (sum(abs(set%efield)) > 1.0E-6_wp) then
         allocate (efield(3), source=0.0_wp)
         efield = set%efield
      end if

      allocate (wbo(mctcmol%nat, mctcmol%nat, chk%tblite%nspin))
      call twostepscf(ctx, chk%tblite, self%ptbData, mctcmol, self%bas, self%cbas, ints, auxints, self%eeqmodel, &
         & results%dipole, results%quadrupole, vecp, neighborlist, selfenergies, v_ES_2nditer, CN_plusU, wbo, efield)
      !> INFO ON RETURNED VARIABLES: On return, ints%hamiltonian contains the last Hamiltonian matrix that was solved
      !> including all potentials and contributions. I.e., it does NOT contain H0 as intended in the usual SCF procedure.

      if (ctx%failed()) then
         do while (ctx%failed())
            call ctx%get_error(error)
            call env%error(error%message, source)
         end do
         call env%error("PTB two-step SCF terminated", source)
         call env%check(exitRun)
         if (exitRun) then
            return
         end if
      end if
      call convert_tblite_to_wfn(env, self%bas, mol, chk, wbo=wbo)
      call convert_tblite_to_results(results,mol,chk,energy,.true.)
      hlgap = results%hl_gap

      !> polarizability by simple perturbative treatment
      !> this is only done in alpha,beta cases
      if (set%elprop == p_elprop_alpha) then
         call numgrad_polarizability(ctx, self%ptbData, mctcmol, self%bas, chk%tblite, &
            & ints, auxints, vecp, neighborlist, selfenergies, v_ES_2nditer, CN_plusU, dF, results%alpha)
         if (ctx%failed()) then
            do while (ctx%failed())
               call ctx%get_error(error)
               call env%error(error%message, source)
            end do
            call env%error("PTB response terminated", source)
            call env%check(exitRun)
            if (exitRun) then
               return
            end if
         end if
      end if
#else
      call feature_not_implemented(env)
#endif

   end subroutine singlepoint

   subroutine writeInfo(self, unit, mol)

      !> Calculator instance
      class(TPTBCalculator), intent(in) :: self

      !> Unit for I/O
      integer, intent(in) :: unit

      !> Molecular structure data
      type(TMolecule), intent(in) :: mol

      call self%ptbData%writeInfo(unit, mol%at)

   end subroutine writeInfo


#if WITH_TBLITE
!---------------------------------------------
! Initialize new wavefunction
!---------------------------------------------
!> Create new wavefunction restart data for tblite library
   subroutine newPTBWavefunction(env, calc, mol, wfn)
      !> Source of the generated errors
      character(len=*), parameter :: source = 'ptb_calculator_newPTBWavefunction'
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      !> Instance of the new calculator
      type(TPTBCalculator), intent(in) :: calc
      !> mctc-io structure type
      type(structure_type) :: mol
      !> Wavefunction data
      type(wavefunction_type) :: wfn
      !> mctc-env error type
      type(error_type), allocatable :: error

      call new_wavefunction(wfn, mol%nat, calc%bas%nsh, calc%bas%nao, &
         & nspin=1, kt=calc%etemp * kt)

      if (allocated(error)) then
         call env%error(error%message, source)
         return
      end if
   end subroutine newPTBWavefunction
#endif

#if ! WITH_TBLITE
   subroutine feature_not_implemented(env)
      !> Computational environment
      type(TEnvironment), intent(inout) :: env

      call env%error("PTB not available without 'tblite'. Compiled without support for 'tblite' library.")
      call env%error("Please recompile without '-Dtblite=disabled' option or change meson setup.")
   end subroutine feature_not_implemented
#endif
end module xtb_ptb_calculator
