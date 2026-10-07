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

module test_oniom
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed, &
      & skip_test
   use xtb_features, only : get_xtb_feature
   implicit none
   private

   public :: collect_oniom

contains

!> Collect all exported unit tests
subroutine collect_oniom(testsuite)
   type(unittest_type),allocatable,intent(out) :: testsuite(:)
      !! array of objects(unit tests)

   testsuite = [ &
      new_unittest("calculateCharge", test_oniom_calculateCharge), &
         !! function that registers a new unit test -> result unittest_type object
      new_unittest("cutbond", test_oniom_cutbond), &
      new_unittest("singlepoint", test_oniom_singlepoint), &
      new_unittest("calculateCharge-tblite", test_oniom_calculateCharge_tblite), &
      new_unittest("singlepoint-tblite", test_oniom_singlepoint_tblite), &
      new_unittest("singlepoint-tblite-alpb", test_oniom_singlepoint_tblite_alpb), &
      new_unittest("singlepoint-tblite-gbsa", test_oniom_singlepoint_tblite_gbsa) &
      ]

end subroutine collect_oniom 

!---------------------------------------------
! Unit test for automatic charge determination
!---------------------------------------------
subroutine test_oniom_calculateCharge(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_basisset
   use xtb_type_environment
   use xtb_setparam
   
   use xtb_oniom, only : TOniomCalculator, calculateCharge, newOniomCalculator, oniom_input
   use xtb_xtb_calculator, only : TxTBCalculator,newWavefunction
   use xtb_io_writer 
   use xtb_mctc_filetypes, only : filetype
   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message
   real(wp), parameter :: thr = 1.0e-7_wp
   !> Molecular structure data
   integer, parameter :: nat = 17
      !! 1,3 pentadiene-1-cation + h2o
   integer, parameter :: at(nat) = [6,6,6,6,6,1,1,1,1,1,1,1,1,1, 8,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape(&
      [ 7.11179039031184_wp, -0.65533236978132_wp,  0.64998029924754_wp, &
      & 9.54574522537207_wp,  0.30133566481849_wp,  0.45943898636196_wp, &
      & 5.13954887436106_wp,  0.33928833942508_wp, -0.70036877110532_wp, &
      &11.69126318297630_wp, -0.71972100226778_wp,  1.75001053205592_wp, &
      & 2.55237710171780_wp, -0.58098675992530_wp, -0.59290041363421_wp, &
      & 6.76811222294413_wp, -2.23038959587740_wp,  1.90836157103276_wp, &
      & 9.86172384803294_wp,  1.88061446715386_wp, -0.80779818687663_wp, &
      & 5.50044858469288_wp,  1.92224580974752_wp, -1.95079140534105_wp, &
      &11.22134447610320_wp, -1.94366600251515_wp,  3.32992090742500_wp, &
      &12.66603076756040_wp, -1.91198818602751_wp,  0.30482336978588_wp, &
      &13.10010470005330_wp,  0.69862396283326_wp,  2.23170047213733_wp, &
      & 1.27494163620898_wp,  0.97527170197973_wp, -0.13473721665480_wp, &
      & 1.96940212226880_wp, -1.18150778329094_wp, -2.48288874184834_wp, &
      & 2.29849092954603_wp, -2.12388855284911_wp,  0.73684528757996_wp, &
      &14.20884549186960_wp, -3.85859246461132_wp, -2.12437520158188_wp, &
      &13.93153897438870_wp, -5.65462186833272_wp, -2.15399113818471_wp, &
      &15.98310331300670_wp, -3.59190426260830_wp, -2.41348990450454_wp  &
      & ], shape(xyz))
   
   type(TMolecule)  :: mol
   type(TRestart)   :: chk1,chk2
   type(TEnvironment)       :: env
   type(TOniomCalculator)   :: calc1,calc2
   type(oniom_input)        :: oniom_input1, oniom_input2

   
   integer ::  innchrg, innchrg2
   
   call init(env)
      !! construct calculation environment

   call init(mol, at, xyz)
      !! construct molecular structure type
      !! interface to initMoleculeNumbers
   mol%chrg = 1.0_wp
   !> allocate calculation  settings
   oniom_input1%first_arg = 'turbomole:gfn2'
   oniom_input2%first_arg = 'orca:gfn2'
   oniom_input1%second_arg = '1-14'
   oniom_input2%second_arg = '15-17'
   
   call newOniomCalculator(calc1,env,mol,oniom_input1)
   select type(xtb => calc1%real_low)
   type is(TxTBCalculator)
      call chk1%wfn%allocate(mol%n,xtb%basis%nshell,xtb%basis%nao)
      call newWavefunction(env,mol,xtb,chk1)
   end select
   innchrg=calculateCharge(calc1,env,mol,chk1)
   call check_(error,calc1%method_low,2)
   call check_(error,calc1%method_high,5)
   call check_(error,chk1%wfn%q(6),0.159091489555113_wp, thr=thr)
   call check_(error,chk1%wfn%q(17),0.340373819431617_wp, thr=thr)
   call check_(error,innchrg,1)

   call newOniomCalculator(calc2,env,mol,oniom_input2)
   select type(xtb => calc2%real_low)
   type is(TxTBCalculator)
      call chk2%wfn%allocate(mol%n,xtb%basis%nshell,xtb%basis%nao)
      call newWavefunction(env,mol,xtb,chk2)
   end select

   innchrg2=calculateCharge(calc2,env,mol,chk2)
   call check_(error,calc2%method_low,2)
   call check_(error,calc2%method_high,4)
   call check_(error,chk2%wfn%q(6),0.159091489555113_wp, thr=thr)
   call check_(error,chk2%wfn%q(17),0.340373819431617_wp, thr=thr)
   call check_(error,innchrg2,0)

endsubroutine test_oniom_calculateCharge

!---------------------------------------------
! Unit test for cutbond
!---------------------------------------------
subroutine test_oniom_cutbond(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : aatoau

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_basisset
   use xtb_type_environment
   use xtb_setparam
   
   use xtb_oniom, only : TOniomCalculator, calculateCharge, newOniomCalculator, oniom_input
   use xtb_xtb_calculator, only : TxTBCalculator,newWavefunction
   use xtb_io_writer 
   use xtb_mctc_filetypes, only : filetype
   use xtb_type_data, only : scc_results

   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message
   real(wp), parameter :: thr = 1.0e-7_wp
      !! molecular structure data
   integer, parameter :: nat = 8
      !! ethane 
   integer, parameter :: at(nat) = [6,6,1,1,1,1,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape(&
      [-6.42056194812604_wp,  4.42017706707601_wp,  0.00000514267503_wp, &
      &-3.54497465898797_wp,  4.42017905349380_wp, -0.00000041470841_wp, &
      &-7.14339394244295_wp,  4.96696898911297_wp,  1.84484217131586_wp, &
      &-7.14339651464075_wp,  2.54910578073934_wp, -0.44887528991923_wp, &
      &-7.14340046326272_wp,  5.74445490934900_wp, -1.39594753088679_wp, &
      &-2.82213644603571_wp,  3.09590294303751_wp,  1.39595405551176_wp, &
      &-2.82214050550548_wp,  6.29125104633740_wp,  0.44887759570039_wp, &
      &-2.82214295551184_wp,  3.87338472886870_wp, -1.84483683242736_wp &
      & ], shape(xyz))
   
   type(TMolecule)  :: mol, inner_mol
   type(TRestart)   :: chk
   type(TEnvironment)       :: env
   type(TOniomCalculator)   :: calc
   type(oniom_input)        :: oniom_input3
   type(scc_results)        :: results
   real(wp), allocatable :: jacobian(:,:)
   integer, allocatable  :: idx2(:)
   real(wp) :: energy, hlgap, sigma(3,3)
   real(wp), allocatable :: gradient(:,:)


   call init(env)
      !! To initialize environment
   call init(mol,at,xyz)
      !! construct molecular structure type
      !! interface to initMoleculeNumbers
   mol%chrg = 0.0_wp
   allocate(gradient(3,mol%n))

   oniom_input3%first_arg = "gfn2:gfn2" 
   oniom_input3%second_arg = '1,3-5'

   call newOniomCalculator(calc,env,mol,oniom_input3)
   
   select type(xtb => calc%real_low)
   type is(TxTBCalculator)
      call chk%wfn%allocate(mol%n,xtb%basis%nshell,xtb%basis%nao)
      call newWavefunction(env,mol,xtb,chk)
   end select
   
   !> Check newOniomCalculator
   call check_(error,calc%method_low,2)
   call check_(error,calc%method_high,2)
   call check_(error,calc%idx(4),5)
   
   !> check newWavefunction
   call check_(error,chk%wfn%q(1),-0.258519131601850_wp, thr=thr)

   !> SP for whole molecule
   call calc%real_low%singlepoint(env,mol,chk,0,.false.,energy,gradient,sigma,hlgap,results)
   
   call check_(error,energy, -7.336370550119_wp, thr=thr)
   call check_(error,hlgap, 15.501740190914_wp, thr=thr)
   call check_(error,norm2(gradient), 0.000327194200_wp, thr=thr)
   
   !> cutbond
   call calc%cutbond(env,mol,chk,calc%topo,inner_mol,jacobian,idx2)
   
   call check_(error,idx2(5),2)
   call check_(error,inner_mol%xyz(1,5),-4.38055107024201_wp,thr=thr)
   call check_(error,inner_mol%xyz(2,5),4.42017847630110_wp,thr=thr)
   call check_(error,inner_mol%xyz(3,5),0.00000120013337_wp,thr=thr)
   
end subroutine test_oniom_cutbond

!---------------------------------------------
! Unit test ONIOM sp 
!---------------------------------------------
subroutine test_oniom_singlepoint(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : autoaa

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_basisset
   use xtb_type_environment
   use xtb_setparam
   
   use xtb_oniom, only : TOniomCalculator, calculateCharge, newOniomCalculator, oniom_input
   use xtb_xtb_calculator, only : TxTBCalculator,newWavefunction
   use xtb_io_writer 
   use xtb_mctc_filetypes, only : filetype
   use xtb_type_data, only : scc_results
   use xtb_gfnff_calculator, only :TGFFCalculator


   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message
   real(wp), parameter :: thr = 1.0e-7_wp
   integer, parameter :: nat = 12
      !! water cluster
   integer, parameter :: at(nat) = [8,1,1, 8,1,1, 8,1,1, 8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape([&
      &-2.75237178376284_wp, 2.43247309226225_wp,-0.01392519847964_wp, &
      &-0.93157260886974_wp, 2.79621404458590_wp,-0.01863384029005_wp, &
      &-3.43820531288547_wp, 3.30583608421060_wp, 1.42134539425148_wp, &
      &-2.43247309226225_wp,-2.75237178376284_wp, 0.01392519847964_wp, &
      &-2.79621404458590_wp,-0.93157260886974_wp, 0.01863384029005_wp, &
      &-3.30583608421060_wp,-3.43820531288547_wp,-1.42134539425148_wp, &
      & 2.75237178376284_wp,-2.43247309226225_wp,-0.01392519847964_wp, &
      & 0.93157260886974_wp,-2.79621404458590_wp,-0.01863384029005_wp, &
      & 3.43820531288547_wp,-3.30583608421060_wp, 1.42134539425148_wp, &
      & 2.43247309226225_wp, 2.75237178376284_wp, 0.01392519847964_wp, &
      & 2.79621404458590_wp, 0.93157260886974_wp, 0.01863384029005_wp, &
      & 3.30583608421060_wp, 3.43820531288547_wp,-1.42134539425148_wp], shape(xyz))
   
   type(TMolecule)  :: mol, inner_mol
   type(TRestart)   :: chk
   type(TEnvironment)       :: env
   type(TOniomCalculator)   :: calc
   type(oniom_input)        :: oniom_input4
   type(scc_results)        :: results
   real(wp),allocatable :: jacobian(:,:)
   integer,allocatable  :: idx_orig(:)
   real(wp) :: energy, hlgap, sigma(3,3)
   real(wp), allocatable :: gradient(:,:)
   integer :: io

   call init(env)
      !! To initialize environment
   call init(mol,at,xyz)
      !! construct molecular structure type
      !! interface to initMoleculeNumbers
   mol%chrg = 0.0_wp
   allocate(gradient(3,mol%n))
   call open_file(io,"w.coord","w")
   call writeMolecule(mol,io,filetype%tmol)
   call close_file(io)
   oniom_input4%second_arg = '1-6'

   call newOniomCalculator(calc,env,mol,oniom_input4)
   
   select type(xtb => calc%real_low)
   type is(TGFFCalculator)
      call check_(error,xtb%topo%qa(1),-0.570632781827558_wp, thr=thr)
      call check_(error,xtb%topo%qa(2),0.285316390913779_wp, thr=thr)
   end select
   
   !> Check newOniomCalculator
   call check_(error,calc%method_low,3)
   call check_(error,calc%method_high,2)
   call check_(error,mol%at(calc%idx(4)),8)
   

   !> SP for whole molecule
   call calc%singlepoint(env,mol,chk,0,.false.,energy,gradient,sigma,hlgap,results)
   
   call check_(error,energy, -10.831044042094_wp, thr=thr)
   call check_(error,norm2(gradient), 0.034482361663_wp, thr=thr)
   !call check_(error,hlgap, 15.501740190914_wp, thr=thr)
   
   !> cutbond
   !call calc%cutbond(env,mol,chk,calc%topo,inner_mol,jacobian,idx_orig)
   
   !call check_(error,idx_orig(5),2)
   !call check_(error,inner_mol%xyz(1,5),-4.38055107024201_wp,thr=thr)
   !call check_(error,inner_mol%xyz(2,5),4.42017847630110_wp,thr=thr)
   !call check_(error,inner_mol%xyz(3,5),0.00000120013337_wp,thr=thr)

end subroutine test_oniom_singlepoint

!---------------------------------------------
! Unit test for automatic charge determination
! with tblite, compared to the native backend
!---------------------------------------------
subroutine test_oniom_calculateCharge_tblite(error)
   use xtb_mctc_accuracy, only : wp

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   
   use xtb_oniom, only : TOniomCalculator, calculateCharge, newOniomCalculator, oniom_input
   use xtb_xtb_calculator, only : TxTBCalculator, newWavefunction
   use xtb_tblite_calculator, only : TTBLiteCalculator, TTBLiteInput, newTBLiteWavefunction
   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message
   real(wp), parameter :: thr = 2.0e-3_wp
      !! EEQ guess charges of native and tblite backend differ slightly
   !> Molecular structure data
   integer, parameter :: nat = 17
      !! 1,3 pentadiene-1-cation + h2o
   integer, parameter :: at(nat) = [6,6,6,6,6,1,1,1,1,1,1,1,1,1, 8,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape(&
      [ 7.11179039031184_wp, -0.65533236978132_wp,  0.64998029924754_wp, &
      & 9.54574522537207_wp,  0.30133566481849_wp,  0.45943898636196_wp, &
      & 5.13954887436106_wp,  0.33928833942508_wp, -0.70036877110532_wp, &
      &11.69126318297630_wp, -0.71972100226778_wp,  1.75001053205592_wp, &
      & 2.55237710171780_wp, -0.58098675992530_wp, -0.59290041363421_wp, &
      & 6.76811222294413_wp, -2.23038959587740_wp,  1.90836157103276_wp, &
      & 9.86172384803294_wp,  1.88061446715386_wp, -0.80779818687663_wp, &
      & 5.50044858469288_wp,  1.92224580974752_wp, -1.95079140534105_wp, &
      &11.22134447610320_wp, -1.94366600251515_wp,  3.32992090742500_wp, &
      &12.66603076756040_wp, -1.91198818602751_wp,  0.30482336978588_wp, &
      &13.10010470005330_wp,  0.69862396283326_wp,  2.23170047213733_wp, &
      & 1.27494163620898_wp,  0.97527170197973_wp, -0.13473721665480_wp, &
      & 1.96940212226880_wp, -1.18150778329094_wp, -2.48288874184834_wp, &
      & 2.29849092954603_wp, -2.12388855284911_wp,  0.73684528757996_wp, &
      &14.20884549186960_wp, -3.85859246461132_wp, -2.12437520158188_wp, &
      &13.93153897438870_wp, -5.65462186833272_wp, -2.15399113818471_wp, &
      &15.98310331300670_wp, -3.59190426260830_wp, -2.41348990450454_wp  &
      & ], shape(xyz))
   
   type(TMolecule) :: mol
   type(TRestart) :: chk_xtb, chk_tblite
   type(TEnvironment) :: env
   type(TOniomCalculator) :: calc_xtb, calc_tblite
   type(oniom_input) :: input
   type(TTBLiteInput) :: tblite
   integer :: iat, innchrg_xtb, innchrg_tblite

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   call init(env)
   call init(mol, at, xyz)
   mol%chrg = 1.0_wp

   input%first_arg = 'gfn2:gfn2'
   input%second_arg = '1-14'

   ! native backend !
   call newOniomCalculator(calc_xtb, env, mol, input)
   select type(xtb => calc_xtb%real_low)
   type is(TxTBCalculator)
      call newWavefunction(env, mol, xtb, chk_xtb)
   class default
      call test_failed(error, "Native ONIOM calculator does not use native xTB")
      return
   end select
   innchrg_xtb = calculateCharge(calc_xtb, env, mol, chk_xtb)

   ! tblite backend !
   input%tblite = .true.
   call newOniomCalculator(calc_tblite, env, mol, input, tblite)
   select type(xtb => calc_tblite%real_low)
   type is(TTBLiteCalculator)
      call newTBLiteWavefunction(env, mol, xtb, chk_tblite)
   class default
      call test_failed(error, "tblite ONIOM calculator does not use tblite")
      return
   end select
   innchrg_tblite = calculateCharge(calc_tblite, env, mol, chk_tblite)

   call check_(error, innchrg_xtb, 1)
   if (allocated(error)) return
   call check_(error, innchrg_tblite, innchrg_xtb)
   if (allocated(error)) return

   ! both backends start from EEQ guess charges !
   do iat = 1, mol%n
      call check_(error, chk_tblite%wfn%q(iat), chk_xtb%wfn%q(iat), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_oniom_calculateCharge_tblite

!---------------------------------------------
! Unit test for ONIOM sp with tblite in gas phase,
! compared to the native backend
!---------------------------------------------
subroutine test_oniom_singlepoint_tblite(error)
   use xtb_mctc_accuracy, only : wp

   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message

   call oniom_singlepoint_tblite(error, -7.370829949440_wp)

end subroutine test_oniom_singlepoint_tblite

!---------------------------------------------
! Unit test for ONIOM sp with tblite in ALPB(water),
! compared to the native backend
!---------------------------------------------
subroutine test_oniom_singlepoint_tblite_alpb(error)
   use xtb_mctc_accuracy, only : wp

   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message

   call oniom_singlepoint_tblite(error, -7.368810571421_wp, "alpb")

end subroutine test_oniom_singlepoint_tblite_alpb

!---------------------------------------------
! Unit test for ONIOM sp with tblite in GBSA(water),
! compared to the native backend
!---------------------------------------------
subroutine test_oniom_singlepoint_tblite_gbsa(error)
   use xtb_mctc_accuracy, only : wp

   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message

   call oniom_singlepoint_tblite(error, -7.369155170153_wp, "gbsa")

end subroutine test_oniom_singlepoint_tblite_gbsa

!---------------------------------------------
! ONIOM sp with native and tblite backend,
! optionally with implicit solvation
!---------------------------------------------
subroutine oniom_singlepoint_tblite(error, reference, solvation_model)
   use xtb_mctc_accuracy, only : wp

   use xtb_type_environment
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data, only : scc_results
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel
   
   use xtb_main_setup, only : addSolvationModel
   use xtb_oniom, only : TOniomCalculator, newOniomCalculator, oniom_input
   use xtb_tblite_calculator, only : TTBLiteInput, TTBLiteSolvationInput

   type(error_type), allocatable, intent(out) :: error
      !! error message type, with stat amd message
   real(wp), intent(in) :: reference
      !! ONIOM energy of the native backend
   character(len=*), intent(in), optional :: solvation_model
      !! solvation model ('alpb' or 'gbsa') with water, gas phase if absent
   real(wp), parameter :: thr = 1.0e-7_wp
   real(wp), parameter :: thr2 = 1.0e-6_wp
   integer, parameter :: nat = 8
      !! ethane 
   integer, parameter :: at(nat) = [6,6,1,1,1,1,1,1]
   real(wp), parameter :: xyz(3,nat) = reshape(&
      [-6.42056194812604_wp,  4.42017706707601_wp,  0.00000514267503_wp, &
      &-3.54497465898797_wp,  4.42017905349380_wp, -0.00000041470841_wp, &
      &-7.14339394244295_wp,  4.96696898911297_wp,  1.84484217131586_wp, &
      &-7.14339651464075_wp,  2.54910578073934_wp, -0.44887528991923_wp, &
      &-7.14340046326272_wp,  5.74445490934900_wp, -1.39594753088679_wp, &
      &-2.82213644603571_wp,  3.09590294303751_wp,  1.39595405551176_wp, &
      &-2.82214050550548_wp,  6.29125104633740_wp,  0.44887759570039_wp, &
      &-2.82214295551184_wp,  3.87338472886870_wp, -1.84483683242736_wp &
      & ], shape(xyz))
   
   type(TMolecule) :: mol, inner_xtb, inner_tblite
   type(TRestart) :: chk_xtb, chk_tblite
   type(TEnvironment) :: env
   type(TOniomCalculator) :: calc_xtb, calc_tblite
   type(oniom_input) :: input
   type(TTBLiteInput) :: tblite
   type(scc_results) :: results
   real(wp), allocatable :: jacobian_xtb(:, :), jacobian_tblite(:, :)
   integer, allocatable :: idx_xtb(:), idx_tblite(:)
   real(wp) :: energy_xtb, energy_tblite, hlgap, sigma(3, 3)
   real(wp), allocatable :: gradient_xtb(:, :), gradient_tblite(:, :)
   integer :: iat, ic, jc

   if (.not.get_xtb_feature('tblite')) then
      call skip_test(error, "xtb not compiled with tblite support")
      return
   end if

   call init(env)
   call init(mol, at, xyz)
   mol%chrg = 0.0_wp
   allocate(gradient_xtb(3, mol%n), gradient_tblite(3, mol%n))

   ! C-C bond is cut, all three calculations use GFN-xTB !
   input%first_arg = "gfn2:gfn1"
   input%second_arg = '1,3-5'

   ! native backend !
   call newOniomCalculator(calc_xtb, env, mol, input)
   if (present(solvation_model)) then
      ! same settings as --alpb and --gbsa in the main program !
      call addSolvationModel(env, calc_xtb, TSolvInput(solvent="water", &
         & alpb=solvation_model == "alpb", &
         & kernel=merge(gbKernel%p16, gbKernel%still, solvation_model == "alpb")))
      tblite%solvation = TTBLiteSolvationInput(solvation_model=solvation_model, &
         & solvent="water")
   end if
   call calc_xtb%singlepoint(env, mol, chk_xtb, 0, .false., energy_xtb, &
      & gradient_xtb, sigma, hlgap, results)

   ! solvation has to be applied, the reference differs from the gas phase !
   call check_(error, energy_xtb, reference, thr=thr)
   if (allocated(error)) return

   ! tblite backend !
   input%tblite = .true.
   call newOniomCalculator(calc_tblite, env, mol, input, tblite)
   call calc_tblite%singlepoint(env, mol, chk_tblite, 0, .false., energy_tblite, &
      & gradient_tblite, sigma, hlgap, results)

   call check_(error, energy_tblite, energy_xtb, thr=thr)
   if (allocated(error)) return
   do iat = 1, mol%n
      do ic = 1, 3
         call check_(error, gradient_tblite(ic, iat), gradient_xtb(ic, iat), thr=thr2)
         if (allocated(error)) return
      end do
   end do

   ! same linked atoms for both backends !
   call calc_xtb%cutbond(env, mol, chk_xtb, calc_xtb%topo, inner_xtb, &
      & jacobian_xtb, idx_xtb)
   call calc_tblite%cutbond(env, mol, chk_tblite, calc_tblite%topo, inner_tblite, &
      & jacobian_tblite, idx_tblite)

   call check_(error, size(idx_tblite), size(idx_xtb))
   if (allocated(error)) return
   call check_(error, size(idx_tblite), 5)
   if (allocated(error)) return
   do iat = 1, size(idx_xtb)
      call check_(error, idx_tblite(iat), idx_xtb(iat))
      if (allocated(error)) return
      do ic = 1, 3
         call check_(error, inner_tblite%xyz(ic, iat), inner_xtb%xyz(ic, iat), thr=thr)
         if (allocated(error)) return
      end do
   end do
   do jc = 1, size(jacobian_xtb, 2)
      do ic = 1, size(jacobian_xtb, 1)
         call check_(error, jacobian_tblite(ic, jc), jacobian_xtb(ic, jc), thr=thr)
         if (allocated(error)) return
      end do
   end do

end subroutine oniom_singlepoint_tblite


end module test_oniom
