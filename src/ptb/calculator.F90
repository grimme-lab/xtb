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

   use mctc_env, only: wp, error_type
   use mctc_io, only: structure_type, new

   use multicharge_model, only: new_mchrg_model, mchrg_model_type

   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_lapack_solver, only: lapack_solver
   use tblite_wavefunction, only: wavefunction_type, new_wavefunction

   use xtb_ptb_param, only: initPTB, ptbGlobals
   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_scf, only: twostepscf
   use xtb_ptb_corebasis, only: add_core_basis
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
   !> Conversion factor from temperature to energy
   real(wp), parameter :: kt = 3.166808578545117e-06_wp

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

      type(context_type) :: ctx

      integer :: i
      real(wp) :: efix
      logical, parameter :: ccm = .true.
      logical :: exitRun
      !> Divide-and-conquer solver
      integer :: gvd = 1
      !> Relatively robust solver
      integer :: gvr = 2
      !> Wavefunction data
      type(wavefunction_type) :: wfn

      call mol%update

      ctx%solver = lapack_solver(gvd)

      energy = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      hlgap = 0.0_wp
      efix = 0.0_wp

      !> Set new PTB wavefunction in tblite format
      call newPTBWavefunction(env, self, wfn)

      call twostepscf(ctx, wfn, self%ptbData, self%mol, self%bas, self%cbas, self%eeqmodel)
      stop

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Electronic structure method terminated", source)
         return
      end if

      ! ------------------------------------------------------------------------
      !  post processing of gradient and energy

      ! point charge embedding gradient file
      ! if (allocated(set%pcem_grad) .and. self%pcem%n > 0) then
      !    call open_file(ich, set%pcem_grad, 'w')
      !    do i = 1, self%pcem%n
      !       write (ich, '(3f12.8)') self%pcem%grd(1:3, i)
      !    end do
      !    call close_file(ich)
      ! end if

      ! ------------------------------------------------------------------------
      !  fixing of certain atoms
      !  print*,abs(efix/etot)
      energy = energy + efix
      results%e_total = energy
      results%gnorm = norm2(gradient)
      if (fixset%n .gt. 0) then
         do i = 1, fixset%n
            !print*,i,fixset%atoms(i)
            gradient(1:3, fixset%atoms(i)) = 0
         end do
      end if

      ! save point charge gradients in results
      ! if (self%pcem%n > 0) then
      !    results%pcem = self%pcem
      ! end if

      if (printlevel .ge. 2) then
         ! start with summary header
         if (.not. set%silent) then
            write (env%unit, '(9x,53(":"))')
            write (env%unit, '(9x,"::",21x,a,21x,"::")') "SUMMARY"
         end if
         write (env%unit, '(9x,53(":"))')
         write (env%unit, outfmt) "total energy      ", results%e_total, "Eh   "
         if (.not. set%silent .and. allocated(self%solvation)) then
            write (env%unit, outfmt) "total w/o Gsasa/hb", &
               &  results%e_total - results%g_sasa - results%g_hb - results%g_shift, "Eh   "
         end if
         write (env%unit, outfmt) "gradient norm     ", results%gnorm, "Eh/a0"
         write (env%unit, outfmt) "HOMO-LUMO gap     ", results%hl_gap, "eV   "
         if (.not. set%silent) then
            if (set%verbose) then
               write (env%unit, '(9x,"::",49("."),"::")')
               write (env%unit, outfmt) "HOMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo), "eV   "
               write (env%unit, outfmt) "LUMO orbital eigv.", chk%wfn%emo(chk%wfn%ihomo + 1), "eV   "
            end if
            write (env%unit, '(9x,"::",49("."),"::")')
            call print_ptb_results(env%unit, results, set%verbose, allocated(self%solvation))
            write (env%unit, outfmt) "add. restraining  ", efix, "Eh   "
            write (env%unit, outfmt) "total charge      ", sum(chk%wfn%q), "e    "
            if (set%verbose) then
               write (env%unit, '(9x,"::",49("."),"::")')
               write (env%unit, outfmt) "atomisation energy", results%e_atom, "Eh   "
            end if
         end if
         write (env%unit, '(9x,53(":"))')
         write (env%unit, '(a)')
      end if

   end subroutine singlepoint

   subroutine print_ptb_results(iunit, res, verbose, lsolv)
      integer, intent(in) :: iunit ! file handle (usually output_unit=6)
      type(scc_results), intent(in) :: res
      logical, intent(in) :: verbose, lsolv
      write (iunit, outfmt) "SCC energy        ", res%e_elec, "Eh   "
      write (iunit, outfmt) "-> isotropic ES   ", res%e_es, "Eh   "
      write (iunit, outfmt) "-> anisotropic ES ", res%e_aes, "Eh   "
      write (iunit, outfmt) "-> anisotropic XC ", res%e_axc, "Eh   "
      write (iunit, outfmt) "-> dispersion     ", res%e_disp, "Eh   "
      if (lsolv) then
         write (iunit, outfmt) "-> Gsolv          ", res%g_solv, "Eh   "
         write (iunit, outfmt) "   -> Gelec       ", res%g_born, "Eh   "
         write (iunit, outfmt) "   -> Gsasa       ", res%g_sasa, "Eh   "
         write (iunit, outfmt) "   -> Ghb         ", res%g_hb, "Eh   "
         write (iunit, outfmt) "   -> Gshift      ", res%g_shift, "Eh   "
      end if
      write (iunit, outfmt) "repulsion energy  ", res%e_rep, "Eh   "
   end subroutine print_ptb_results

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
         & nspin=1, kt=calc%etemp*kt)

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
