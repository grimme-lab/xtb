! This file is part of xtb.
!
! Copyright (C) 2022 Christoph Plett
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

module xtb_iff_iffprepare
   use xtb_mctc_accuracy, only: wp
   use xtb_splitparam
   use xtb_iff_calculator, only: TIFFCalculator
   implicit none

contains

   subroutine prepare_IFF(env, comb)
      use xtb_docking_param
      use xtb_type_environment, only: TEnvironment
      use xtb_type_molecule, only: TMolecule
      use xtb_setparam

      type(TEnvironment), intent(inout) :: env
      type(TMolecule), intent(in) :: comb!Combined structure of molA and molB (molA has to be first)

      character(len=*), parameter :: source = 'preperation_IFFCalculator'
      type(TMolecule) :: molA, molB
      integer, allocatable :: at(:)
      real(wp), allocatable :: xyz(:, :)
      real(wp) :: molA_e, molB_e

      if (natom_molA == 0) call env%error('No atoms of Molecule A given') !Not given as input => Crash
      if(natom_molA > comb%n) call env%error('More atoms of Molecule A than contained in structure')
      molA%n = natom_molA
      molB%n = comb%n - natom_molA
      call split_mol(molA, molB, comb)
      set%lmoinfo_fname = 'xtblmoinfA'
      call precomp(env, molA, molA_e, 1)
      set%lmoinfo_fname = 'xtblmoinfB'
      call precomp(env, molB, molB_e, 2)
      set%lmoinfo_fname = 'xtblmoinfo'

   end subroutine prepare_IFF

   subroutine precomp(env, mol, etot, mol_num)
      use xtb_type_restart, only: TRestart
      use xtb_type_calculator, only: TCalculator
      use xtb_xtb_calculator, only: TxTBCalculator
      use xtb_setparam
      use xtb_setmod, only: set_gfn
      use xtb_type_environment, only: TEnvironment
      use xtb_readin, only: xfind
      use xtb_solv_state, only: solutionState
      use xtb_type_data, only: scc_results
      use xtb_main_defaults, only: initDefaults
      use xtb_disp_ncoord, only: ncoord_gfn, ncoord_erf
      use xtb_scc_core, only: iniqshell
      use xtb_eeq, only: goedecker_chrgeq
      use xtb_main_setup, only: newCalculator
      use xtb_single, only: singlepoint
      use xtb_type_molecule, only: TMolecule
      use xtb_docking_param, only: chrg, uhf, gsolvstate_iff, pre_e_A, pre_e_B, optlvl

      implicit none

      !> Molecular structure data
      type(TMolecule), intent(inout) :: mol
      !> Calculation environment
      type(TEnvironment), intent(inout) :: env
      integer, intent(in) :: mol_num

      class(TCalculator), allocatable :: calc
      type(TRestart) :: chk
      character(len=:), allocatable :: fnv !parameter file
      logical :: restart = .false.
      real(wp) :: acc = 1.0_wp !SCF accuracy
      real(wp) :: egap
      integer :: gsolvstate !Should be gas
      real(wp), allocatable :: cn(:)
      real(wp), allocatable :: dcn(:, :, :), dq(:, :, :), g(:, :)
      real(wp) :: er
      integer, external :: ncore
      logical :: lgrad
      logical :: exist
      real(wp), intent(out) :: etot
      real(wp) :: sigma(3, 3)
      type(scc_results) :: res
      integer :: i
      integer :: extrun_tmp
      integer, allocatable ::  tmp_unit
      logical :: newdisp_tmp
      integer :: itemp = 48

      allocate (cn(mol%n), g(3, mol%n))

      set%pr_lmo = .true.
      set%silent = .true.

      !Setting up z
      do i = 1, mol%n
         mol%z(i) = mol%at(i) - ncore(mol%at(i))
         ! lanthanides without f are treated as La
         if (mol%at(i) .gt. 57 .and. mol%at(i) .lt. 72) mol%z(i) = 3
      end do

      !> Set GFN2 settings
      if (optlvl /= 'gfn1') then
         set%gfn_method = 2
         fnv = xfind('param_gfn2-xtb.txt')
         extrun_tmp = set%mode_extrun
         set%mode_extrun = p_ext_xtb
         newdisp_tmp = set%newdisp
         set%newdisp = .true.
      else
         fnv = xfind('param_gfn1-xtb.txt') !Other stuff already set
      end if
      if (mol_num .eq. 1) then
         if (chrg(1) /= 0.0_wp) mol%chrg = chrg(1)
         if (uhf(1) /= 0.0_wp) mol%uhf = uhf(1)
      elseif (mol_num .eq. 2) then
         if (chrg(2) /= 0.0_wp) mol%chrg = chrg(2)
         if (uhf(2) /= 0.0_wp) mol%uhf = uhf(2)
      end if

      call open_file(itemp, 'tmp_lmo', 'w')
      tmp_unit = env%unit
      env%unit = itemp

      !> New calculator
      call newCalculator(env, mol, calc, fnv, restart, acc)
      call env%checkpoint("Could not setup parameterisation")

!   gsolvstate = solutionState%gsolv
      call initDefaults(env, calc, mol, gsolvstate_iff)
      !> initial guess, setup wavefunction
      select type (calc)
      type is (TxTBCalculator)
         calc%etemp = set%etemp
         calc%maxiter = set%maxscciter
         call chk%wfn%allocate(mol%n, calc%basis%nshell, calc%basis%nao)
         !> EN charges and CN
         call ncoord_gfn(mol%n, mol%at, mol%xyz, cn)
         if (mol%npbc > 0) then
            chk%wfn%q = real(set%ichrg, wp)/real(mol%n, wp)
         else
            call ncoord_erf(mol%n, mol%at, mol%xyz, cn)
            call goedecker_chrgeq(mol%n,mol%at,mol%xyz,real(set%ichrg,wp),cn,dcn,chk%wfn%q,dq,er,g,&
                                  .false., .false., .false.)
            chk%wfn%q = real(set%ichrg, wp)/real(mol%n, wp)
         end if
         !> initialize shell charges from gasteiger charges
      call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell,chk%wfn%q,chk%wfn%qsh,set%gfn_method)
      end select

      call delete_file('.sccnotconverged')
      call delete_file('charges')

      call env%checkpoint("Setup for calculation failed")

      !> the SP
      call singlepoint &
         &       (env, mol, chk, calc, egap, set%etemp, set%maxscciter, 2,&
         &        exist, lgrad, acc, etot, g, sigma, res)

      set%pr_lmo = .false.

      if (mol_num .eq. 1) then
         pre_e_A = etot
      elseif (mol_num .eq. 2) then
         pre_e_B = etot
      end if

      if (optlvl /= 'gfn1') then
         set%mode_extrun = extrun_tmp
         set%newdisp = newdisp_tmp
      end if

      set%silent = .false.
      env%unit = tmp_unit
      call remove_file(itemp)
      deallocate (cn, g)

   end subroutine precomp

end module
