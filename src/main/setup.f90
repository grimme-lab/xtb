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
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : init
   use xtb_extern_orca, only : TOrcaCalculator, newOrcaCalculator
   use xtb_extern_mopac, only : TMopacCalculator, newMopacCalculator
   use xtb_extern_turbomole, only : TTMCalculator, newTMCalculator
   use xtb_extern_driver, only : TDriverCalculator, newDriverCalculator
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, newTBLiteCalculator
   use xtb_type_calculator, only : TCalculator
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   use xtb_type_restart, only : TRestart
   use xtb_type_wavefunction, only : TWavefunction
   use xtb_xtb_calculator, only : TxTBCalculator, newXTBcalculator, newWavefunction
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_iff_calculator, only : TIFFCalculator, newIFFCalculator
   use xtb_iff_data, only : TIFFData
   use xtb_oniom, only : TOniomCalculator, newOniomCalculator, oniom_input
   use xtb_setparam
   implicit none
   private

   public :: newCalculator, newWavefunction, addSolvationModel
   public :: newXTBCalculator


contains


subroutine newCalculator(env, mol, calc, fname, restart, accuracy, input, iff_data, tblite_input)

   character(len=*), parameter :: source = 'main_setup_newCalculator'

   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(in) :: mol

   class(TCalculator), allocatable, intent(out) :: calc

   character(len=*), intent(in) :: fname

   logical, intent(in) :: restart

   real(wp), intent(in) :: accuracy

   type(oniom_input), intent(in), optional :: input

   type(TIFFData), intent(in), optional, allocatable :: iff_data

   !> Input for TBLite calculator
   type(TTBLiteInput), intent(in), optional :: tblite_input

   type(TxTBCalculator), allocatable :: xtb
   type(TTBLiteCalculator), allocatable :: tblite
   type(TGFFCalculator), allocatable :: gfnff
   type(TIFFCalculator), allocatable :: iff
   type(TOrcaCalculator), allocatable :: orca
   type(TMopacCalculator), allocatable :: mopac
   type(TTMCalculator), allocatable :: turbo
   type(TOniomCalculator), allocatable :: oniom
   type(TDriverCalculator), allocatable :: driver
   
   logical :: exitRun

   select case(set%mode_extrun)
   case default
      call env%error("Unknown calculator type", source)
   case(p_ext_eht, p_ext_xtb)
      allocate(xtb)

      call newXTBCalculator(env, mol, xtb, fname, set%gfn_method, accuracy)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(xtb, calc)
   case(p_ext_tblite)
      if (.not.present(tblite_input)) then
         call env%error("TBLite calculator requires input", source)
         return
      end if
      allocate(tblite)

      call newTBLiteCalculator(env, mol, tblite, tblite_input)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(tblite, calc)
   case(p_ext_gfnff)
      allocate(gfnff)

      call newGFFCalculator(env, mol, gfnff, fname, restart)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(gfnff, calc)
   case(p_ext_iff)
      if (.not.present(iff_data)) then
         call env%error("IFF calculator requires input", source)
         return
      end if
      allocate(iff)

      if (.not. allocated(iff_data)) then
         call env%error("IFF Data not present for Calculator", source)
      end if
      
      call newIFFCalculator(env, mol, iff_data, iff)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not construct new calculator", source)
         return
      end if

      call move_alloc(iff, calc)

   case(p_ext_orca)
      allocate(orca)
      call newOrcaCalculator(orca, env, set%ext_orca)
      call move_alloc(orca, calc)
   case(p_ext_mopac)
      allocate(mopac)
      call newMopacCalculator(mopac, env, set%ext_mopac)
      call move_alloc(mopac, calc)
   case(p_ext_turbomole)
      allocate(turbo)
      call newTMCalculator(turbo, set%extcode, set%extmode)
      call move_alloc(turbo, calc)
   case(p_ext_driver)
      allocate(driver)
      call newDriverCalculator(driver, env, set%ext_driver)
      call move_alloc(driver, calc)
   case(p_ext_oniom)
      if (.not.present(input)) then
         call env%error("ONIOM calculator requires input", source)
         return
      end if
      allocate(oniom)
      call newOniomCalculator(oniom, env, mol, input)
      call move_alloc(oniom, calc)
   end select

end subroutine newCalculator



subroutine addSolvationModel(env, calc, input)
   type(TEnvironment), intent(inout) :: env
   class(TCalculator), intent(inout) :: calc
   type(TSolvInput), intent(in) :: input
   integer :: level

   level = 0
   select type(calc)
   type is(TxTBCalculator)
      level = calc%xtbData%level
   type is(TOniomCalculator)
      select type(xtb => calc%real_low)
      type is(TxTBCalculator)
         level = xtb%xtbData%level
      end select
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
