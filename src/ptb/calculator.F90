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

!> Extended tight binding calculator
module xtb_ptb_calculator
   use xtb_mctc_accuracy, only: wp
   use xtb_type_basisset, only: TBasisset
   use xtb_type_calculator, only: TCalculator
   use xtb_type_data
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule, only: TMolecule
   use xtb_type_param, only: scc_parameter, chrg_parameter, &
   & TPTBParameter
   use xtb_type_pcem
   use xtb_type_restart, only: TRestart
   use xtb_ptb_data, only: TPTBData
   use xtb_setparam
   use xtb_fixparam
   use xtb_scf, only: scf
   use xtb_peeq, only: peeq
   use xtb_embedding, only: read_pcem
   use xtb_basis, only: newBasisset
   use xtb_mctc_systools, only: rdpath
   use xtb_readparam, only: readParam
   use xtb_paramset, only: use_parameterset
   use xtb_chargemodel, only: new_charge_model_2019
   use xtb_disp_ncoord, only: ncoord_erf
   use xtb_eeq, only: eeq_chrgeq
   use xtb_iniq, only: iniqcn
   use xtb_scc_core, only: iniqshell

   use mctc_io, only: structure_type, new
   use tblite_basis_type, only: basis_type

   use xtb_ptb_param, only: initPTB, ptbGlobals
   implicit none

   private

   public :: TPTBCalculator, newPTBCalculator, newWavefunction

   !> Calculator interface for PTB method
   type, extends(TCalculator) :: TPTBCalculator

      !> Structure type
      type(structure_type), allocatable :: struc

      !> Tight binding basis set
      type(TBasisset), allocatable :: basis

      !> Parametrisation data base
      type(TPTBData), allocatable :: ptbData

      !> Electronic temperature
      real(wp) :: etemp

      !> Maximum number of cycles for SCC convergence
      integer :: maxiter

      !> External potential
      type(tb_pcem) :: pcem

   contains

      !> Perform PTB single point calculation
      procedure :: singlepoint

      !> Write informative printout
      procedure :: writeInfo

   end type TPTBCalculator

   character(len=*), private, parameter :: outfmt = &
                                           '(9x,"::",1x,a,f23.12,1x,a,1x,"::")'

contains

   subroutine newPTBCalculator(env, mol, calc, accuracy)

      character(len=*), parameter :: source = 'xtb_ptb_calculator_newPTBCalculator'

      type(TEnvironment), intent(inout) :: env

      type(TMolecule), intent(in) :: mol

      type(TPTBCalculator), intent(out) :: calc

      real(wp), intent(in), optional :: accuracy

      character(len=:), allocatable :: filename
      type(TPTBParameter) :: globpar
      integer :: ich
      logical :: exist, okbas
      logical :: exitRun

! #if WITH_TBLITE
      type(structure_type) :: struc

      call new(struc, mol%at, mol%xyz, mol%chrg, mol%uhf, mol%lattice)
      if (allocated(mol%pdb)) struc%pdb = mol%pdb
      if (allocated(mol%sdf)) struc%sdf = mol%sdf
      calc%struc = struc

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
      allocate(calc%ptbData)
      call open_file(ich, filename, 'r')
      exist = ich /= -1
      if (exist) then
         error stop "Parameter file not supported yet."
         call close_file(ich)
      else ! no parameter file, check if we have one compiled into the code
         call initPTB(calc%ptbData)
      end if

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not load parameters", source)
         return
      end if

      !> set up the basis set for the tb-Hamiltonian
      ! allocate (calc%basis)
      ! call newBasisset(calc%ptbData, mol%n, mol%at, calc%basis, okbas)
      ! if (.not. okbas) then
      !    call env%error('basis set could not be setup completely', source)
      !    return
      ! end if

      !> check for external point charge field
      ! if (allocated(set%pcem_file)) then
      !    call open_file(ich, set%pcem_file, 'r')
      !    if (ich /= -1) then
      !       call read_pcem(ich, env, calc%pcem, calc%xtbData%coulomb)
      !       call close_file(ich)
      !    end if
      ! end if
! #else
    call feature_not_implemented(env)
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

      integer :: i, ich
      integer :: mode_sp_run = 1
      real(wp) :: efix
      logical :: inmol
      logical, parameter :: ccm = .true.
      logical :: exitRun

      call mol%update

      energy = 0.0_wp
      gradient(:, :) = 0.0_wp
      sigma(:, :) = 0.0_wp
      hlgap = 0.0_wp
      efix = 0.0_wp

      ! ------------------------------------------------------------------------
      !  actual calculation
      ! select case (self%ptbData%level)
      ! case (1, 2)
      !    ! if (allocated(self%solvation)) then
      !    !    call newSolvationModel(self%solvation, env, solvation, mol%at)
      !    ! end if
      !    call scf(env, mol, chk%wfn, self%basis, self%pcem, self%xtbData, solvation, &
      !       &   hlgap, self%etemp, self%maxiter, printlevel, restart, .true., &
      !       &   self%accuracy, energy, gradient, results)

      ! case (0)
      !    if (allocated(self%solvation)) then
      !       allocate (gbsa)
      !       call newBornModel(self%solvation, env, gbsa, mol%at)
      !    end if
      !    call peeq &
      !       & (env, mol, chk%wfn, self%basis, self%xtbData, gbsa, hlgap, self%etemp, &
      !       &  printlevel, .true., ccm, self%accuracy, energy, gradient, sigma, results)
      ! end select

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Electronic structure method terminated", source)
         return
      end if

      ! ------------------------------------------------------------------------
      !  post processing of gradient and energy

      ! point charge embedding gradient file
      if (allocated(set%pcem_grad) .and. self%pcem%n > 0) then
         call open_file(ich, set%pcem_grad, 'w')
         do i = 1, self%pcem%n
            write (ich, '(3f12.8)') self%pcem%grd(1:3, i)
         end do
         call close_file(ich)
      end if

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
      if (self%pcem%n > 0) then
         results%pcem = self%pcem
      end if

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
      use xtb_type_data
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
   subroutine newWavefunction(env, mol, calc, chk)

      implicit none
      character(len=*), parameter :: source = 'xtb_ptb_calculator_newWavefunction'
      !! Name of error producer routine
      type(TEnvironment), intent(inout) :: env
      !! Calculation environment to handle I/O stream and error log
      type(TRestart), intent(inout) :: chk
      !! Restart data wrapper for wfn and nlist
      type(TPTBCalculator), intent(in) :: calc
      !! Instance of PTB Calculator
      type(TMolecule), intent(in) :: mol
      !! Molecular structure data
      real(wp), allocatable :: cn(:)
      !! Coordination number
      type(chrg_parameter) :: chrgeq
      !! guess charges(gasteiger/goedecker/sad)
      logical :: exitRun
      !! if it is recommended to terminate the run

      associate (wfn => chk%wfn)
         allocate (cn(mol%n))
         call wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)

         !> find partial charges
         if (mol%npbc > 0) then
         !! if periodic
            wfn%q = mol%chrg/real(mol%n, wp)
            !! evenly distribute charge with the equal partial charges
         else
            if (set%guess_charges .eq. p_guess_gasteiger) then
               call iniqcn(mol%n, mol%at, mol%z, mol%xyz, nint(mol%chrg), 1.0_wp, &
                  & wfn%q, cn, calc%ptbData%level, .true.)
            else if (set%guess_charges .eq. p_guess_goedecker) then
            !! default

               call new_charge_model_2019(chrgeq, mol%n, mol%at)
               !! to get parametrized values for q (en,gam,kappa,alpha)

               call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
               !! to obtain CN
               !! (49) Extended Tight-Binding Quantum Chemistry Mehods 2020

               call eeq_chrgeq(mol, env, chrgeq, cn, wfn%q)
               !! to obtain partial charges q
               !! (47) Extended Tight-Binding Quantum Chemistry Mehods 2020

               call env%check(exitRun)
               !! to check status of environment
               if (exitRun) then
                  call env%rescue("EEQ guess failed, falling back to SAD guess", source)
                  wfn%q = mol%chrg/real(mol%n, wp)
               end if
            else
               wfn%q = mol%chrg/real(mol%n, wp)
            end if
         end if

         !> find shell charges
         ! call iniqshell(calc%ptbData, mol%n, mol%at, mol%z, calc%basis%nshell, &
         !    & wfn%q, wfn%qsh, calc%ptbData%level)

      end associate

   end subroutine newWavefunction


! #if ! WITH_TBLITE
subroutine feature_not_implemented(env)
   !> Computational environment
   type(TEnvironment), intent(inout) :: env

   call env%error("Compiled without support for tblite library")
end subroutine feature_not_implemented
! #endif
end module xtb_ptb_calculator
