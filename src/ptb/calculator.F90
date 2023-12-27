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

! #ifndef WITH_TBLITE
! #define WITH_TBLITE 0
! #endif

!> Density matrix (P) tight binding (TB) calculator
module xtb_ptb_calculator
   use xtb_type_calculator, only: TCalculator
   use xtb_type_data
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule, assignment(=)
   use xtb_type_param, only: TPTBParameter
   use xtb_type_restart, only: TRestart
   use xtb_ptb_data, only: TPTBData
   use xtb_setparam
   use xtb_fixparam
   use xtb_mctc_systools, only: rdpath
   use xtb_readin, only: bool2string, bool2int

   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new
   use mctc_io_convert, only: autoev

   use multicharge_model, only: new_mchrg_model, mchrg_model_type

   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_lapack_solver, only: lapack_solver
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction
   use tblite_integral_type, only: integral_type
   use tblite_adjlist, only: adjacency_list

   use xtb_ptb_param, only: initPTB, ptbGlobals
   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_scf, only: twostepscf
   use xtb_ptb_corebasis, only: add_core_basis
   use xtb_ptb_integral_types, only: aux_integral_type
   use xtb_ptb_response, only: numgrad_polarizability
   implicit none

   private

   public :: TPTBCalculator, newPTBCalculator

   !> Calculator interface for PTB method
   type, extends(TCalculator) :: TPTBCalculator

      !> Structure type
      type(structure_type), allocatable :: mol

      !> PTB vDZP basis set
      type(basis_type) :: bas, cbas

      !> EEQ Model
      type(mchrg_model_type) :: eeqmodel

      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData

      !> Electronic temperature
      real(wp) :: etemp

      !> Maximum number of cycles for SCC convergence
      integer :: maxiter

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
   real(wp), parameter :: dF = 1.0E-1_wp

contains

   subroutine newPTBCalculator(env, struc, calc, accuracy)

      character(len=*), parameter :: source = 'xtb_ptb_calculator_newPTBCalculator'

      type(TEnvironment), intent(inout) :: env

      type(TMolecule), intent(in) :: struc

      type(TPTBCalculator), intent(out) :: calc

      real(wp), intent(in), optional :: accuracy

      character(len=:), allocatable :: filename
      type(TPTBParameter) :: globpar
      integer :: ich
      logical :: exist
      logical :: exitRun

! #if WITH_TBLITE
      type(structure_type) :: mol

      mol = struc
      calc%mol = mol

      call rdpath(env%xtbpath, 'param_ptb.txt', filename, exist)
      if (.not. exist) filename = 'param_ptb.txt'

      if (present(accuracy)) then
         calc%accuracy = accuracy
      else
         calc%accuracy = 1.0_wp
      end if

      calc%etemp = 300.0_wp
      calc%maxiter = 250

      !> Obtain the parameter file
      allocate (calc%ptbData)
      call open_file(ich, filename, 'r')
      exist = ich /= -1
      if (exist) then
         error stop "Parameter file not supported yet."
         call close_file(ich)
      else ! no parameter file, check if we have one compiled into the code
         call initPTB(calc%ptbData, mol%num)
      end if

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not load parameters", source)
         return
      end if

      !> set up the basis set for the PTB-Hamiltonian
      call add_vDZP_basis(calc%mol, calc%bas)
      !> Add the core basis set to 'cbas' basis set type
      call add_core_basis(calc%mol, calc%ptbData%corepotential, calc%cbas)
      !> set up the EEQ model
      call new_mchrg_model(calc%eeqmodel, chi=calc%ptbData%eeq%chi, &
      & rad=calc%ptbData%eeq%alp, eta=calc%ptbData%eeq%gam, kcn=calc%ptbData%eeq%cnf)
      !##### DEV WRITE #####
      ! loop over all atoms and print the number of shells and primitives
      ! write (*, *) "Number of atoms: ", struc%nat
      ! write (*, *) "Number of shells: ", calc%bas%nsh
      ! write (*, *) calc%bas%nsh_id
      ! do i = 1, struc%nat
      !    write (*, *) "i: ", i
      !    write (*, *) "ID(i): ", struc%id(i)
      !    write(*,*) "num(ID(i)): ", struc%num(struc%id(i))
      !    do j = 1, calc%bas%nsh_id(struc%id(i))
      !       write (*, *) "  shell: ", j, "  prim: ", calc%bas%cgto(j, i)%nprim
      !    end do
      ! end do
      !#####################

      !> check for external point charge field
      ! if (allocated(set%pcem_file)) then
      !    call open_file(ich, set%pcem_file, 'r')
      !    if (ich /= -1) then
      !       call read_pcem(ich, env, calc%pcem, calc%xtbData%coulomb)
      !       call close_file(ich)
      !    end if
      ! end if
! #else
      ! call feature_not_implemented(env)
! #endif

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

      !#################################################
      !> PTB INDIVIDUAL
      !#################################################
      !> tblite calculation context
      type(context_type) :: ctx
      !> Wiberg bond order
      real(wp), allocatable :: wbo(:, :, :)
      !> Static homogenoues external electric field
      real(wp), allocatable :: efield(:)
      !> Approx. effective core potential
      real(wp), allocatable :: Vecp(:, :)
      !> Static dipole polarizability tensor
      real(wp) :: alpha(3, 3)
      !> Auxiliary integrals
      type(aux_integral_type) :: auxints
      !> Exact integrals
      type(integral_type) :: ints
      !> Adjacency list
      type(adjacency_list) :: neighborlist
      !> Effective self-energies
      real(wp), allocatable :: selfenergies(:)

      logical :: exitRun
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2

      energy = 0.0_wp
      gradient = 0.0_wp
      sigma = 0.0_wp

      call mol%update

      ctx%solver = lapack_solver(gvd)
      ctx%unit = env%unit
      ctx%verbosity = printlevel

      !> Set new PTB wavefunction in tblite format
      call newPTBWavefunction(env, self, chk%tblite)

      !> Static Homogeneous External Electric Field
      if (sum(abs(set%efield)) > 1.0E-6_wp) then
         allocate (efield(3), source=0.0_wp)
         efield = set%efield
      end if

      allocate (wbo(self%mol%nat, self%mol%nat, chk%tblite%nspin))
      call twostepscf(ctx, chk%tblite, self%ptbData, self%mol, self%bas, self%cbas, ints, auxints, self%eeqmodel, &
         & results%dipole, vecp, neighborlist, selfenergies, wbo, efield)
      !> INFO ON RETURNED VARIABLES: On return, ints%hamiltonian contains the last Hamiltonian matrix that was solved
      !> including all potentials and contributions. I.e., it does NOT contain H0 as intended in the usual SCF procedure.

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Electronic structure method terminated", source)
         return
      end if

      call chk%wfn%allocate(mol%n, self%bas%nsh, self%bas%nao)
      chk%wfn%n = self%mol%nat
      chk%wfn%nel = nint(chk%tblite%nocc)
      chk%wfn%nopen = self%mol%uhf
      chk%wfn%nshell = self%bas%nsh
      chk%wfn%nao = self%bas%nao
      chk%wfn%P = chk%tblite%density(:, :, 1)
      chk%wfn%q = chk%tblite%qat(:, 1)
      chk%wfn%qsh = chk%tblite%qsh(:, 1)
      chk%wfn%focca = chk%tblite%focc(:, 1)
      chk%wfn%foccb = 0.0_wp
      chk%wfn%focc(:) = chk%tblite%focc(:, 1)
      chk%wfn%emo = chk%tblite%emo(:, 1) * autoev
      chk%wfn%C = chk%tblite%coeff(:, :, 1)
      chk%wfn%ihomo = chk%tblite%homo(1)
      chk%wfn%ihomoa = chk%tblite%homo(1)
      chk%wfn%ihomob = chk%tblite%homo(2)
      chk%wfn%wbo = wbo(:, :, 1)

      results%hl_gap = (chk%tblite%emo(chk%tblite%homo(1) + 1, 1) - chk%tblite%emo(chk%tblite%homo(1), 1)) * autoev
      hlgap = results%hl_gap

      !> polarizability by simple perturbative treatment
      !> this is only done in alpha,beta cases
      if (set%runtyp == p_run_alpha) then
         call numgrad_polarizability(ctx, self%ptbData, self%mol, self%bas, chk%tblite, &
            & ints, auxints, vecp, neighborlist, selfenergies, dF, alpha)
      end if

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

!---------------------------------------------
! Initialize new wavefunction
!---------------------------------------------
!> Create new wavefunction restart data for tblite library
   subroutine newPTBWavefunction(env, calc, wfn)
      !> Source of the generated errors
      character(len=*), parameter :: source = 'ptb_calculator_newPTBWavefunction'
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      !> Instance of the new calculator
      type(TPTBCalculator), intent(in) :: calc

      !#if WITH_TBLITE
      !> Wavefunction data
      type(wavefunction_type) :: wfn
      !> mctc-env error type
      type(error_type), allocatable :: error

      call new_wavefunction(wfn, calc%mol%nat, calc%bas%nsh, calc%bas%nao, &
         & nspin=1, kt=calc%etemp * kt)

      if (allocated(error)) then
         call env%error(error%message, source)
         return
      end if
      !#else
      ! call feature_not_implemented(env)
      !#endif
   end subroutine newPTBWavefunction

! #if ! WITH_TBLITE
   subroutine feature_not_implemented(env)
      !> Computational environment
      type(TEnvironment), intent(inout) :: env

      call env%error("Compiled without support for tblite library")
   end subroutine feature_not_implemented
! #endif
end module xtb_ptb_calculator
