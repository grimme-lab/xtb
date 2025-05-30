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

module test_gfnff
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_gfnff

contains

!> Collect all exported unit tests
subroutine collect_gfnff(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("sp", test_gfnff_sp), &
      new_unittest("hb", test_gfnff_hb), &
      new_unittest("gbsa", test_gfnff_gbsa), &
      new_unittest("mindless", test_gfnff_mindless_basic), &
      new_unittest("mindless-solvation", test_gfnff_mindless_solvation), &
      new_unittest("scaleup", test_gfnff_scaleup), &
      new_unittest("pdb", test_gfnff_pdb), &
      new_unittest("sdf", test_gfnff_sdf), &
      new_unittest("pbc", test_gfnff_pbc), &
      new_unittest("Ln_An", test_gfnff_LnAn_H) &
      ]

end subroutine collect_gfnff


subroutine test_gfnff_sp(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_solv_gbsa
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_gfnff_neighbourlist
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 8
   integer, parameter :: at(nat) = [7,15,9,9,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[1.50040286526241_wp,-2.88219140061585_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 1.90142164792207_wp, 0.00000000000000_wp, &
      &-0.02649585010919_wp,-5.05651406856634_wp, 0.00000000000000_wp, &
      &-1.39762821979929_wp, 4.65236211997835_wp, 0.00000000000000_wp, &
      & 2.62205170116282_wp,-3.14316635901963_wp, 1.53958066744940_wp, &
      &-1.46489869067775_wp, 0.78885483581631_wp, 1.94964934855945_wp, &
      & 2.62205170116282_wp,-3.14316635901963_wp,-1.53958066744940_wp, &
      &-1.46489869067775_wp, 0.78885483581631_wp,-1.94964934855945_wp],&
      & shape(xyz))
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc
   type(TGFFNeighbourList) :: nlist
   type(TBorn), allocatable :: solvation

   real(wp) :: etot
   real(wp), allocatable :: g(:,:)
   real(wp) :: sigma(3, 3)
   character(len=:),allocatable :: fnv
   integer  :: ipar
   real(wp) :: efield(3) = [0.0_wp, 0.0_wp, 0.0_wp]

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,6)
   call check_(error, calc%topo%nangl,6)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call gfnff_eg(env,mol,gff_print,mol%n,nint(mol%chrg),mol%at,mol%xyz,sigma, &
      & g,etot,res_gff,calc%param,calc%topo,calc%neigh,nlist,efield,solvation,.true.,calc%version, &
      & calc%accuracy)

   call check_(error, res_gff%e_total,-0.76480130317838_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.06237477492373_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.74131049663951_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.00633910404059_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.00004724445432_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.05070333390156_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.00224146422313_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.03086605590295_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.00003142616658_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.00776698664545_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.00000000000000_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_sp

subroutine test_gfnff_hb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,1,8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[-5.78520874132429_wp,-1.92328475821000_wp,-0.02944611115854_wp, &
      & -5.57801768832583_wp, 0.17912532844037_wp, 0.72444143178660_wp, &
      & -4.27822256938673_wp,-2.74845397256109_wp,-1.13038073598642_wp, &
      & -7.47879539136783_wp,-2.97570121473211_wp, 0.39488815557786_wp, &
      & -0.83005338399036_wp, 2.43458470560665_wp,-0.78566331969245_wp, &
      & -0.74201439536855_wp, 4.04199055249898_wp, 0.09144422329636_wp, &
      & -2.44679415487233_wp, 1.69392751177087_wp,-0.27417668699116_wp],&
      & shape(xyz))
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TRestart)      :: chk
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc

   real(wp) :: etot, sigma(3,3), gap
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '---', .false.)

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,5)
   call check_(error, calc%topo%nangl,4)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call calc%singlepoint(env, mol, chk, 1, .false., etot, g, sigma, gap, res_gff)

   call check_(error, res_gff%e_total,-0.949706677118_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.001152720923_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.856707643513_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.000579711773_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.000000008811_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.152313816530_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.001251669186_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.066881023899_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.006894292337_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.000000000000_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_hb

subroutine test_gfnff_gbsa(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_solv_input
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_gfnff_param
   use xtb_gfnff_setup
   use xtb_gfnff_eg
   use xtb_gfnff_ini
   use xtb_setparam
   use xtb_setmod
   use xtb_disp_dftd3param
   use xtb_disp_dftd4
   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [6,8,1,1,8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[-5.78520874132429_wp,-1.92328475821000_wp,-0.02944611115854_wp, &
      & -5.57801768832583_wp, 0.17912532844037_wp, 0.72444143178660_wp, &
      & -4.27822256938673_wp,-2.74845397256109_wp,-1.13038073598642_wp, &
      & -7.47879539136783_wp,-2.97570121473211_wp, 0.39488815557786_wp, &
      & -0.83005338399036_wp, 2.43458470560665_wp,-0.78566331969245_wp, &
      & -0.74201439536855_wp, 4.04199055249898_wp, 0.09144422329636_wp, &
      & -2.44679415487233_wp, 1.69392751177087_wp,-0.27417668699116_wp],&
      & shape(xyz))
   type(scc_options),parameter :: opt = scc_options( &
      &  prlevel = 2, maxiter = 30, acc = 1.0_wp, etemp = 300.0_wp, grad = .true.,&
      &  solvent = "h2o")
   logical, parameter :: restart = .false.

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TRestart)      :: chk
   type(scc_results)   :: res_gff
   type(TGFFCalculator) :: calc

   real(wp) :: etot, sigma(3,3), gap
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   call delete_file('charges')
   call newGFFCalculator(env, mol, calc, '---', .false.)
   call addSolvationModel(env, calc, TSolvInput(solvent=opt%solvent, alpb=.false., kernel=gbKernel%still))

   call env%checkpoint("GFN-FF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call check_(error, calc%neigh%nbond,5)
   call check_(error, calc%topo%nangl,4)
   call check_(error, calc%topo%ntors,1)

   g = 0.0_wp
   gff_print=.true.

   call calc%singlepoint(env, mol, chk, 1, .false., etot, g, sigma, gap, res_gff)

   call check_(error, res_gff%e_total,-0.964158677062_wp, thr=thr)
   call check_(error, res_gff%gnorm,   0.013624276205_wp, thr=thr)
   call check_(error, res_gff%e_bond, -0.856707643513_wp, thr=thr)
   call check_(error, res_gff%e_angl,  0.000579711773_wp, thr=thr)
   call check_(error, res_gff%e_tors,  0.000000008811_wp, thr=thr)
   call check_(error, res_gff%e_es,   -0.150043166563_wp, thr=thr*10)
   call check_(error, res_gff%e_disp, -0.001251669186_wp, thr=thr)
   call check_(error, res_gff%e_rep,   0.066881023899_wp, thr=thr)
   call check_(error, res_gff%e_hb,   -0.006894292337_wp, thr=thr)
   call check_(error, res_gff%e_xb,   -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%e_batm, -0.000000000000_wp, thr=thr)
   call check_(error, res_gff%g_solv, -0.016722649876_wp, thr=thr)
   call check_(error, res_gff%g_sasa,  0.000126368690_wp, thr=thr)
   call check_(error, res_gff%g_hb,   -0.009238122476_wp, thr=thr)
   call check_(error, res_gff%g_born, -0.009468339217_wp, thr=thr)
   call check_(error, res_gff%g_shift, 0.001857443126_wp, thr=thr)

   call mol%deallocate

end subroutine test_gfnff_gbsa


subroutine test_gfnff_mindless_basic(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   real(wp), parameter :: ref_energies(10) = &
      &[-1.6235670601725_wp, -1.2974060907872_wp, -1.5825072926857_wp, &
      & -1.6086171378897_wp, -1.5093596077875_wp, -1.6076220027918_wp, &
      & -1.7328634195448_wp, -1.8875339867396_wp, -1.3924035489143_wp, &
      & -1.9583702712389_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.11714711890640_wp, 0.08933168067298_wp, 0.15687133018871_wp, &
      & 0.09820451651462_wp, 0.08460010429134_wp, 0.08787739425161_wp, &
      & 0.12463658172704_wp, 0.10062734775717_wp, 0.06347506656236_wp, &
      & 0.09445561400996_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_mindless_basic


subroutine test_gfnff_mindless_solvation(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   character(len=*), parameter :: solvents(10) = [character(len=20) ::&
      & "h2o", "chcl3", "thf", "acetonitrile", "toluene", &
      & "ch2cl2", "ether", "methanol", "cs2", "dmso"]
   real(wp), parameter :: ref_energies(10) = &
      &[-1.6540916426629_wp, -1.3143056510098_wp, -1.6141899469716_wp, &
      & -1.6204572392623_wp, -1.5288738676799_wp, -1.6293089943213_wp, &
      & -1.7562513711414_wp, -1.9034123090085_wp, -1.4068010735061_wp, &
      & -1.9722601751413_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.11594122130461_wp, 0.08871796411903_wp, 0.15190830479600_wp, &
      & 0.09873901233954_wp, 0.08427306878219_wp, 0.08670875610818_wp, &
      & 0.11948602245812_wp, 0.10022284156198_wp, 0.06288313258789_wp, &
      & 0.09564569352042_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0, kernel=gbKernel%still))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_mindless_solvation


subroutine test_gfnff_scaleup(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: molecules(5) = [character(len=20) ::&
      & "caffeine", "rivaroxaban", "grubbs", "remdesivir", "taxol"]
   character(len=*), parameter :: solvents(5) = [character(len=20) ::&
      & "h2o", "acetonitrile", "toluene", "ether", "dmso"]
   real(wp), parameter :: ref_energies(5) = &
      &[-4.6919926039901_wp, -8.8807760138817_wp, -13.40715596616603_wp, &
      & -13.822994859730_wp, -20.577952329212_wp]
   real(wp), parameter :: ref_gnorms(5) = &
      &[0.05947676640487_wp, 0.09522104624089_wp, 0.1753769495334539_wp, &
      & 0.12496592222660_wp, 0.19366599743810_wp]

   call init(env)
   do iMol = 1, 5

      call getMolecule(mol, trim(molecules(iMol)))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      call addSolvationModel(env, calc, TSolvInput(solvent=trim(solvents(iMol)), &
         & alpb=mod(iMol, 2)==0, kernel=gbKernel%still))

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_scaleup


subroutine test_gfnff_pdb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   real(wp), parameter :: ref_energies(3) = &
      &[-12.358885277947_wp, -12.442589851697_wp, -12.442131341448_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.15855661051730_wp, 0.15148957207346_wp, 0.15151167424491_wp]

   call init(env)
   do iMol = 1, 3

      call getMolecule(mol, 'pdb-4qxx')

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='h2o', &
            & alpb=iMol==3, kernel=gbKernel%still))
      end if

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_pdb


subroutine test_gfnff_sdf(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
   use xtb_main_setup, only : addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   use xtb_setparam, only : set

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)
   real(wp), parameter :: ref_energies(3) = &
      &[-0.98330642628373_wp, -1.0826001974262_wp, -1.0826319137928_wp]
   real(wp), parameter :: ref_gnorms(3) = &
      &[0.11515550863614e-2_wp, 0.59760221346330e-2_wp, 0.59645861418074e-2_wp]

   call init(env)
   do iMol = 1, 3

      call getMolecule(mol, 'bug332')
      set%ichrg = nint(mol%chrg)

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call delete_file('charges')
      call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
      if (iMol > 1) then
         call addSolvationModel(env, calc, TSolvInput(solvent='h2o', &
            & alpb=iMol==3, kernel=gbKernel%p16))
      end if

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
         & hl_gap, res)

      call env%check(exitRun)
      call check_(error, .not.exitRun)
      if (exitRun) exit

      call check_(error, energy, ref_energies(iMol), thr=thr)
      call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)

   end do

end subroutine test_gfnff_sdf


subroutine test_gfnff_pbc(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp
   real(wp), parameter :: thr2 = 3.0e-3_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

    integer :: iMol, i
    logical :: exitRun
    real(wp) :: energy, energy2, hl_gap, sigma(3, 3)
    real(wp), allocatable :: gradient(:, :)
    real(wp), allocatable :: xyztmp(:,:), lattmp(:,:)
    integer, allocatable :: attmp(:)
    character(len=4), allocatable :: symtmp(:)
    logical, parameter :: pbc(3) = [.true., .true., .true. ]
    ! structures from X23, mcVOL22, and GFN-FF for Ln/An paper
    character(len=*), parameter :: pbc_strucs(3) = [&
       & "x06_b", "mcv15", "Th_15"]
    !> references for original GFN-FF
    real(wp), parameter :: ref_energies(3) = &
        & [-9.522300429916_wp, -11.059826732607_wp, -46.265672914_wp ]
    real(wp), parameter :: ref_gnorms(3) = &
        & [ 0.083513313043_wp, 0.083870455567_wp, 1.6127503408_wp ]
    real(wp), parameter :: ref_snorm(2) = &
        & [ 1.410826672_wp, 0.419545156_wp ]
    !> references for mcGFN-FF
    real(wp), parameter :: mcref_energies(2) = &
        & [-9.419561502589_wp, -10.850342398282_wp ]
    real(wp), parameter :: mcref_gnorms(2) = &
        & [ 0.071130586474_wp, 0.104285414859_wp ]
    real(wp), parameter :: mcref_snorm(2) = &
        & [ 0.839554942_wp, 0.729174866_wp ]

    call init(env)
    do iMol = 1, 2

       !  Calucaltions with original parameterization  !
       ! load molecule from molstock
       call getMolecule(mol, pbc_strucs(iMol))
       ! reset energy, gradient, sigma
       energy=0.0_wp
       if (allocated(gradient)) deallocate(gradient)
       allocate(gradient(3, mol%n))
       sigma = 0.0_wp
       ! setup new calculator
       call delete_file('charges')
       ! original angewChem2020_2 version is default
       call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! run single point calculation
       call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
          & hl_gap, res)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! check energy, gradient and sigma versus reference calculation
       call check_(error, energy, ref_energies(iMol), thr=thr)
       call check_(error, norm2(gradient), ref_gnorms(iMol), thr=thr)
       call check_(error, norm2(sigma), ref_snorm(iMol), thr=thr)

       !  Calucaltions with mcGFN-FF parameterization  !
       ! reset energy, gradient, sigma
       energy=0.0_wp
       if (allocated(gradient)) deallocate(gradient)
       allocate(gradient(3, mol%n))
       sigma = 0.0_wp
       ! setup new calculator
       call delete_file('charges')
       ! mcGFN-FF is version=4
       call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! run single point calculation
       call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
          & hl_gap, res)
       call env%check(exitRun)
       call check_(error, .not.exitRun)
       if (exitRun) exit
       ! check energy, gradient and sigma versus reference calculation
       call check_(error, energy, mcref_energies(iMol), thr=thr)
       call check_(error, norm2(gradient), mcref_gnorms(iMol), thr=thr)
       call check_(error, norm2(sigma), mcref_snorm(iMol), thr=thr)
    end do

    ! check GFN-FF calculation on periodic system with actinide (Th)
    ! load molecule "Th_15" from molstock
    call getMolecule(mol, pbc_strucs(3))
    ! reset energy, gradient, sigma
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! original angewChem2020_2 version is default
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! check energy and gradient versus reference calculation
    call check_(error, energy, ref_energies(3), thr=thr)
    call check_(error, norm2(gradient), ref_gnorms(3), thr=thr)


    !  check super cell scaling  !
    ! load molecule from molstock
    call getMolecule(mol, "mcv15")
    ! reset energy, gradient, sigma
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! mcGFN-FF is version=4
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! build 2x1x1 supercell
    ! coordinates
    if (allocated(xyztmp)) deallocate(xyztmp)
    allocate(xyztmp(3,mol%n))
    xyztmp = mol%xyz
    deallocate(mol%xyz)
    allocate(mol%xyz(3,2*mol%n), source=0.0_wp)
    mol%xyz(:,1:mol%n)=xyztmp
    do i=1, mol%n
       mol%xyz(:,i+mol%n) = xyztmp(:,i)+mol%lattice(:,1)
    enddo
    ! atom types
    allocate(attmp(mol%n), source=0)
    attmp=mol%at
    deallocate(mol%at)
    allocate(mol%at(2*mol%n))
    mol%at(1:mol%n)=attmp
    mol%at(mol%n+1:2*mol%n)=attmp
    ! symbols
    allocate(symtmp(2*mol%n))
    symtmp = mol%sym
    deallocate(mol%sym)
    allocate(character(len=len(symtmp)) :: mol%sym(2*mol%n))
    mol%sym(1:mol%n)= symtmp
    mol%sym(mol%n+1:2*mol%n) = symtmp
    ! lattice
    mol%lattice(:,1) = 2.0_wp*mol%lattice(:,1)
    ! number of atoms
    mol%n = 2*mol%n
    ! init mol
    deallocate(symtmp, xyztmp)
    allocate(symtmp(mol%n))
    allocate(xyztmp(3,mol%n), source=0.0_wp)
    allocate(lattmp(3,3), source=0.0_wp)
    symtmp=mol%sym
    xyztmp=mol%xyz
    lattmp=mol%lattice
    call init(mol, symtmp, xyztmp, 0.0_wp, 0, lattmp, pbc)

    energy2=0.0_wp ! energy of supercell
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    sigma = 0.0_wp
    ! setup new calculator
    call delete_file('charges')
    ! mcGFN-FF is version=4
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false., 4)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy2, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)

    ! scale down energy of supercell and compare
    energy2 = energy2/2.0_wp
    call check_(error, energy2, energy, thr=thr2)

end subroutine test_gfnff_pbc

subroutine test_gfnff_LnAn_H(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp ! threshold for energy
   real(wp), parameter :: qthr = 1.0e-4_wp ! threshold for charge

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TGFFCalculator) :: calc
   type(scc_results) :: res

    integer :: i
    integer, parameter :: nat = 33
    logical :: exitRun
    real(wp) :: energy, charges(nat), hl_gap, sigma(3,3)
    real(wp), allocatable :: gradient(:, :)
    real(wp), allocatable :: xyztmp(:,:), lattmp(:,:)
    integer, allocatable :: attmp(:)
    character(len=4), allocatable :: symtmp(:)
    logical, parameter :: pbc(3) = [.true., .true., .true. ]
    ! structures from X23, mcVOL22, and GFN-FF for Ln/An paper
    character(len=*), parameter :: struc(1) = ["Ce_0a7745"]
    !> references for original GFN-FF
    real(wp), parameter :: ref_charges(nat) = [&
      &  0.93714798_wp, -0.39685427_wp, -0.00104516_wp, -0.01815703_wp, -0.26122075_wp, &
      & -0.01774628_wp, -0.00287629_wp, -0.36261674_wp,  0.19471344_wp,  0.20497869_wp, &
      &  0.09972845_wp,  0.08148308_wp,  0.10674954_wp,  0.09486206_wp,  0.17232582_wp, &
      &  0.11414769_wp,  0.09722964_wp,  0.10512679_wp,  0.08798849_wp,  0.22945410_wp, &
      &  0.15310280_wp,  0.88629154_wp, -0.48033324_wp, -0.48478265_wp, -0.49522170_wp, &
      & -0.14264511_wp, -0.12548232_wp, -0.01061975_wp, -0.40598965_wp,  0.13119126_wp, &
      &  0.11822500_wp,  0.13303034_wp,  0.25781424_wp &
      & ]
    real(wp), parameter :: ref_energy = -3.550331926678 ! reference energy
    real(wp), parameter :: ref_gnorm = 0.213662542141 ! reference gradient norm
    ! number of neighbors and sum of neighbor indices for atom 1 (Ce)
    integer, parameter :: ref_numnb_1 = 9
    integer, parameter :: ref_sumnb_1 = 177 ! 2+5+8+22+24+25+26+27+29 +9
    ! number of neighbors and sum of neighbor indices for atom 21 (H)
    integer, parameter :: ref_numnb_21 = 1
    integer, parameter :: ref_sumnb_21 = 9 ! 8 +1
    ! number of neighbors and sum of neighbor indices for atom 33 (H)
    integer, parameter :: ref_numnb_33 = 1
    integer, parameter :: ref_sumnb_33 = 30 ! 29 +1



    call init(env)

    ! check GFN-FF calculation on periodic system with actinide (Th)
    ! load molecule from molstock
    call getMolecule(mol, struc(1))
    ! reset energy and gradient
    energy=0.0_wp
    if (allocated(gradient)) deallocate(gradient)
    allocate(gradient(3, mol%n))
    ! setup new calculator
    call delete_file('charges')
    ! original angewChem2020_2 version is default
    call newGFFCalculator(env, mol, calc, '.param_gfnff.xtb', .false.)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! run single point calculation
    call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
       & hl_gap, res)
    call env%check(exitRun)
    call check_(error, .not.exitRun)
    ! check energy and gradient versus reference calculation
    call check_(error, energy, ref_energy, thr=thr)
    call check_(error, norm2(gradient), ref_gnorm, thr=thr)
    ! get charges from SP calculation
    charges = chk%nlist%q ! these charges are in the gfnff_charges file
    ! check charges of all atoms
    do i=1, nat
      ! compare to reference
      call check_(error, charges(i), ref_charges(i), thr=qthr)
    enddo
    ! The Ln/An-H bonds in GFN-FF are charge dependent
    ! But it is compared in the gfnff_ini where only topo%qa is available
    ! Therefore we check the neighbor list here additionally
    call check_(error, calc%neigh%nb(42,1,1), ref_numnb_1) ! check numnb of atom 1
    call check_(error, sum(calc%neigh%nb(:,1,1)), ref_sumnb_1) ! check sumnb of atom 1
    call check_(error, calc%neigh%nb(42,21,1), ref_numnb_21) ! check numnb of atom 21
    call check_(error, sum(calc%neigh%nb(:,21,1)), ref_sumnb_21) ! check sumnb of atom 21
    call check_(error, calc%neigh%nb(42,33,1), ref_numnb_33) ! check numnb of atom 33
    call check_(error, sum(calc%neigh%nb(:,33,1)), ref_sumnb_33) ! check sumnb of atom 33

end subroutine test_gfnff_LnAn_H

end module test_gfnff
