! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> TODO
module xtb_main_setup
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools, only : rdpath
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : init
   use xtb_type_calculator, only : TCalculator
   use xtb_type_dummycalc, only : TDummyCalculator
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param, only : TxTBParameter, chrg_parameter
   use xtb_type_restart, only : TRestart
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_readparam, only : readParam
   use xtb_paramset, only : use_parameterset
   use xtb_basis, only : newBasisset
   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_eeq, only : eeq_chrgeq
   use xtb_iniq, only : iniqcn
   use xtb_scc_core, only : iniqshell
   use xtb_embedding, only : read_pcem
   use xtb_setparam
   use xtb_disp_ncoord
   use xtb_chargemodel
   implicit none
   private

   public :: newCalculator, newWavefunction, addSolvationModel
   public :: newXTBCalculator, newGFFCalculator


contains


subroutine newCalculator(env, mol, calc, fname, restart, accuracy)

   character(len=*), parameter :: source = 'main_setup_newCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   class(TCalculator), allocatable, intent(out) :: calc

   character(len=*), intent(in) :: fname

   logical, intent(in) :: restart

   real(wp), intent(in) :: accuracy

   type(TxTBCalculator), allocatable :: xtb
   type(TGFFCalculator), allocatable :: gfnff
   
   logical :: exitRun

   select case(mode_extrun)
   case default
      call env%error("Unknown calculator type", source)
   case(p_ext_eht, p_ext_xtb)
      allocate(xtb)

      call newXTBCalculator(env, mol, xtb, fname, gfn_method, accuracy)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(xtb, calc)
   case(p_ext_gfnff)
      allocate(gfnff)

      call newGFFCalculator(env, mol, gfnff, fname, restart)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(gfnff, calc)
   case(p_ext_qmdff, p_ext_orca, p_ext_turbomole, p_ext_mopac)
      allocate(TDummyCalculator :: calc)
      calc%accuracy = accuracy
   end select

end subroutine newCalculator


subroutine newXTBCalculator(env, mol, calc, fname, method, accuracy)

   character(len=*), parameter :: source = 'main_setup_newXTBCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TxTBCalculator), intent(out) :: calc

   integer, intent(in), optional :: method

   real(wp), intent(in), optional :: accuracy

   character(len=*), intent(in), optional :: fname

   character(len=:), allocatable :: filename
   type(TxTBParameter) :: globpar
   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   if (present(fname)) then
      filename = fname
   else
      if (present(method)) then
         select case(method)
         case(0)
            call rdpath(env%xtbpath, 'param_gfn0-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn0-xtb.txt'
         case(1)
            call rdpath(env%xtbpath, 'param_gfn1-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn1-xtb.txt'
         case(2)
            call rdpath(env%xtbpath, 'param_gfn2-xtb.txt', filename, exist)
            if (.not.exist) filename = 'param_gfn2-xtb.txt'
         end select
      end if
   end if
   if (.not.allocated(filename)) then
      call env%error("No parameter file or parametrisation info provided", source)
      return
   end if

   if (present(accuracy)) then
      calc%accuracy = accuracy
   else
      calc%accuracy = 1.0_wp
   end if

   calc%etemp = 300.0_wp
   calc%maxiter = 250

   !> Obtain the parameter file
   allocate(calc%xtbData)
   call open_file(ich, filename, 'r')
   exist = ich /= -1
   if (exist) then
      call readParam(env, ich, globpar, calc%xtbData, .true.)
      call close_file(ich)
   else ! no parameter file, check if we have one compiled into the code
      call use_parameterset(filename, globpar, calc%xtbData, exist)
      if (.not.exist) then
         call env%error('Parameter file '//filename//' not found!', source)
         return
      end if
   endif

   if (present(method)) then
      if (method /= calc%xtbData%level) then
         call env%error("Requested method does not match loaded method", source)
         return
      end if
   end if

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not load parameters", source)
      return
   end if

   !> set up the basis set for the tb-Hamiltonian
   allocate(calc%basis)
   call newBasisset(calc%xtbData, mol%n, mol%at, calc%basis, okbas)
   if (.not.okbas) then
      call env%error('basis set could not be setup completely', source)
      return
   end if

   !> check for external point charge field
   if (allocated(pcem_file)) then
      call open_file(ich, pcem_file, 'r')
      if (ich /= -1) then
         call read_pcem(ich, env, calc%pcem, calc%xtbData%coulomb)
         call close_file(ich)
      end if
   end if

end subroutine newXTBCalculator

subroutine newGFFCalculator(env, mol, calc, fname, restart, version)
   use xtb_gfnff_param
   use xtb_gfnff_setup, only : gfnff_setup
   use xtb_disp_dftd4, only : newD3Model

   character(len=*), parameter :: source = 'main_setup_newGFFCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TGFFCalculator), intent(out) :: calc

   character(len=*), intent(in) :: fname

   logical, intent(in) :: restart

   integer, intent(in), optional :: version

   type(TxTBParameter) :: globpar
   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   if (present(version)) then
      calc%version = version
   else
      calc%version = gffVersion%angewChem2020_2
   end if

   call calc%topo%zero
   calc%update = .true.
   ! global accuracy factor similar to acc in xtb used in SCF
   calc%accuracy = 0.1_wp
   if (mol%n > 10000) then
      calc%accuracy = 2.0_wp
   end if

   !> Obtain the parameter file
   call open_file(ich, fname, 'r')
   exist = ich /= -1
   if (exist) then
      call gfnff_read_param(ich, calc%param)
      call close_file(ich)
   else ! no parameter file, try to load internal version
      call gfnff_load_param(calc%version, calc%param, exist)
      if (.not.exist) then
         call env%error('Parameter file '//fname//' not found!', source)
         return
      end if
   endif

   call newD3Model(calc%topo%dispm, mol%n, mol%at)

   call gfnff_setup(env, verbose, restart, mol, &
      & calc%gen, calc%param, calc%topo, calc%accuracy, calc%version)

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not create force field calculator", source)
      return
   end if

end subroutine newGFFCalculator

subroutine newWavefunction(env, mol, calc, chk)
   type(TEnvironment), intent(inout) :: env
   type(TRestart), intent(inout) :: chk
   type(TxTBCalculator), intent(in) :: calc
   type(TMolecule), intent(in) :: mol
   call newWavefunction_(env, mol, calc, chk%wfn)
end subroutine newWavefunction

subroutine newWavefunction_(env, mol, calc, wfn)
   character(len=*), parameter :: source = 'main_setup_newWavefunction'
   type(TEnvironment), intent(inout) :: env
   type(TWavefunction), intent(inout) :: wfn
   type(TxTBCalculator), intent(in) :: calc
   type(TMolecule), intent(in) :: mol
   real(wp), allocatable :: cn(:)
   type(chrg_parameter) :: chrgeq
   logical :: exitRun

   allocate(cn(mol%n))
   call wfn%allocate(mol%n,calc%basis%nshell,calc%basis%nao)

   wfn%nel = idint(sum(mol%z)) - nint(mol%chrg)
   wfn%nopen = mol%uhf
   if(wfn%nopen == 0 .and. mod(wfn%nel,2) /= 0) wfn%nopen=1

   if (mol%npbc > 0) then
      wfn%q = mol%chrg/real(mol%n,wp)
   else
      if (guess_charges.eq.p_guess_gasteiger) then
         call iniqcn(mol%n,wfn%nel,mol%at,mol%z,mol%xyz,nint(mol%chrg),1.0_wp, &
            & wfn%q,cn,calc%xtbData%level,.true.)
      else if (guess_charges.eq.p_guess_goedecker) then
         call new_charge_model_2019(chrgeq,mol%n,mol%at)
         call ncoord_erf(mol%n,mol%at,mol%xyz,cn)
         call eeq_chrgeq(mol,env,chrgeq,cn,wfn%q)

         call env%check(exitRun)
         if (exitRun) then
            call env%rescue("EEQ guess failed, falling back to SAD guess", source)
            wfn%q = mol%chrg/real(mol%n,wp)
         end if
      else
         wfn%q = mol%chrg/real(mol%n,wp)
      end if
   end if

   call iniqshell(calc%xtbData,mol%n,mol%at,mol%z,calc%basis%nshell, &
      & wfn%q,wfn%qsh,calc%xtbData%level)
end subroutine newWavefunction_


subroutine addSolvationModel(env, calc, input)
   type(TEnvironment), intent(inout) :: env
   class(TCalculator), intent(inout) :: calc
   type(TSolvInput), intent(in) :: input
   integer :: level

   level = 0
   select type(calc)
   type is(TxTBCalculator)
      level = calc%xtbData%level
   end select

   if (allocated(input%solvent)) then
      calc%lSolv = input%solvent /= 'none' .and. input%solvent /= 'gas' &
         & .and. input%solvent /= 'vac'
   else
      calc%lSolv = .false.
   end if

   if (calc%lSolv) then
      allocate(calc%solvation)
      call init(calc%solvation, env, input, level)
   endif

end subroutine addSolvationModel

end module xtb_main_setup
