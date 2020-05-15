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
   use xtb_type_calculator, only : TCalculator
   use xtb_type_dummycalc, only : TDummyCalculator
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_param, only : TxTBParameter
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_readparam, only : readParam
   use xtb_paramset, only : use_parameterset
   use xtb_basis, only : newBasisset
   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_gfnff_calculator, only : TGFFCalculator
   use xtb_eeq, only : goedecker_chrgeq
   use xtb_iniq, only : iniqcn
   use xtb_scc_core, only : iniqshell
   use xtb_setparam
   use xtb_disp_ncoord
   implicit none
   private

   public :: newCalculator
   public :: newXTBCalculator, newGFFCalculator


contains


subroutine newCalculator(env, mol, calc, fname)

   character(len=*), parameter :: source = 'main_setup_newCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   class(TCalculator), allocatable, intent(out) :: calc

   character(len=*), intent(in) :: fname

   type(TxTBCalculator), allocatable :: xtb
   type(TGFFCalculator), allocatable :: gfnff
   
   logical :: exitRun

   select case(mode_extrun)
   case default
      call env%error("Unknown calculator type", source)
   case(p_ext_eht, p_ext_xtb)
      allocate(xtb)

      call newXTBCalculator(env, mol, xtb, fname)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(xtb, calc)
   case(p_ext_gfnff)
      allocate(gfnff)

      call newGFFCalculator(env, mol, gfnff, fname)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(gfnff, calc)
   case(p_ext_qmdff, p_ext_orca, p_ext_turbomole, p_ext_mopac)
      allocate(TDummyCalculator :: calc)
   end select

end subroutine newCalculator


subroutine newXTBCalculator(env, mol, calc, fname)

   character(len=*), parameter :: source = 'main_setup_newXTBCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TxTBCalculator), intent(out) :: calc

   character(len=*), intent(in) :: fname

   type(TxTBParameter) :: globpar
   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   !> Obtain the parameter file
   allocate(calc%xtbData)
   call open_file(ich, fname, 'r')
   exist = ich /= -1
   if (exist) then
      call readParam(env, ich, globpar, calc%xtbData, .true.)
      call close_file(ich)
   else ! no parameter file, check if we have one compiled into the code
      call use_parameterset(fname, globpar, calc%xtbData, exist)
      if (.not.exist) then
         call env%error('Parameter file '//fname//' not found!', source)
         return
      end if
   endif

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

end subroutine newXTBCalculator

subroutine newGFFCalculator(env, mol, calc, fname)
   use xtb_gfnff_param

   character(len=*), parameter :: source = 'main_setup_newGFFCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   type(TGFFCalculator), intent(out) :: calc

   character(len=*), intent(in) :: fname

   type(TxTBParameter) :: globpar
   integer :: ich
   logical :: exist, okbas
   logical :: exitRun

   !> Obtain the parameter file
   call open_file(ich, fname, 'r')
   exist = ich /= -1
   if (exist) then
      call gfnff_read_param(ich, calc%param)
      call close_file(ich)
   else ! no parameter file
      call env%error('Parameter file '//fname//' not found!', source)
      return
   endif

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Could not load parameters", source)
      return
   end if

end subroutine newGFFCalculator

end module xtb_main_setup
