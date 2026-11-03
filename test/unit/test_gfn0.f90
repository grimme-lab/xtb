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

module test_gfn0
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, &
      & test_failed, to_string
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   !> Tolerance for comparing parameters
   real(wp), parameter :: thr_param = 1.0e-12_wp

   !> Compare optional per-element parameter arrays
   interface check_param
      module procedure :: check_param_real1
      module procedure :: check_param_real2
      module procedure :: check_param_int1
      module procedure :: check_param_int2
   end interface check_param

   public :: collect_gfn0

contains

!> Collect all exported unit tests
subroutine collect_gfn0(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("builtin-param", test_gfn0_builtin_param), &
      new_unittest("sp", test_gfn0_sp), &
      new_unittest("api", test_gfn0_api), &
      new_unittest("api-srb", test_gfn0_api_srb), &
      new_unittest("mindless-basic", test_gfn0_mindless_basic), &
      new_unittest("mindless-solvation", test_gfn0_mindless_solvation) &
      ]

end subroutine collect_gfn0


subroutine test_gfn0_sp(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_mctc_systools
   use xtb_type_options

   use xtb_type_molecule
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment
   use xtb_solv_gbsa

   use xtb_setparam
   use xtb_basis
   use xtb_peeq
   use xtb_readparam
   use xtb_paramset

   use xtb_xtb_data
   use xtb_xtb_gfn0

   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-7_wp
   real(wp),parameter :: thr2 = 1.0e-5_wp
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
   real(wp),parameter :: et = 300.0_wp
   integer, parameter :: maxiter = 50
   integer, parameter :: prlevel = 2
   logical, parameter :: lgrad = .true.
   logical, parameter :: restart = .false.
   real(wp),parameter :: acc = 1.0_wp

   type(TEnvironment)  :: env
   type(TMolecule)     :: mol
   type(scc_results)     :: res
   type(TBasisset)     :: basis
   type(TWavefunction) :: wfn
   type(TxTBData) :: xtbData
   type(TBorn), allocatable :: gbsa

   real(wp) :: etot,egap,sigma(3,3)
   real(wp), allocatable :: g(:,:)
   character(len=:),allocatable :: fnv
   integer  :: ipar

   type(TxTBParameter) :: globpar
   logical  :: okpar,okbas,exist,diff

   call init(env)

   set%gfn_method = 0

   call init(mol, at, xyz)

   wfn%nel = idint(sum(mol%z))
   wfn%nopen = 0

   allocate( g(3,mol%n), source = 0.0_wp )
 
   ! use the built-in parameters, no parameter file required
   call use_parameterset('param_gfn0-xtb.txt',globpar,xtbData,okpar)
   call check_(error, okpar)
   if (allocated(error)) return

   call newBasisset(xtbData,mol%n,mol%at,basis,okbas)
   call check_(error, okbas)

   call check_(error, basis%nshell,17)
   call check_(error, basis%nbf,   30)
   call check_(error, basis%nao,   29)

   call wfn%allocate(mol%n,basis%nshell,basis%nao)
   wfn%q = mol%chrg/real(mol%n,wp)

   g = 0.0_wp

   call peeq(env,mol,wfn,basis,xtbData,gbsa, &
      &   egap,et,prlevel,lgrad,.false.,acc,etot,g,sigma,res)

   call check_(error, res%converged)

   call check_(error, res%e_total,-15.271927490420_wp, thr=thr)
   call check_(error, res%gnorm,0.55600868499186E-01_wp, thr=thr)
   ! value in electron volt
   call check_(error, res%hl_gap, 4.4910991783546_wp, thr=1.0e-4_wp)
   call check_(error, res%e_elec,-15.238855850629_wp, thr=thr)
   ! es and aes are not welldefined at this accuracy level
   call check_(error, res%e_es,   -0.091088841288_wp, thr=thr*10)
   call check_(error, res%e_disp, -0.001331408911_wp, thr=thr)
   call check_(error, res%e_rep,   0.071830283799_wp, thr=thr)
   call check_(error, res%e_xb,   -0.012481673391_wp, thr=thr)

   call mol%deallocate
   call wfn%deallocate
   call basis%deallocate

end subroutine test_gfn0_sp

subroutine test_gfn0_api(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_data
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-7_wp
   integer, parameter :: nat = 7
   integer, parameter :: at(nat) = [51,1,1,1,1,16,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      [-0.33839323285714_wp,-2.52035427985714_wp, 0.00000000000000_wp, &
      & 0.87126950614286_wp,-4.36564525285714_wp, 2.28275498100000_wp, &
      & 0.87126950614286_wp,-4.36564525285714_wp,-2.28275498100000_wp, &
      &-3.07540243785714_wp,-4.12867550385714_wp, 0.00000000000000_wp, &
      & 2.16321847414286_wp, 6.15526246914286_wp, 0.00000000000000_wp, &
      &-0.33839323285714_wp, 5.87053934414286_wp, 0.00000000000000_wp, &
      &-0.15356858285714_wp, 3.35451847614286_wp, 0.00000000000000_wp], shape(xyz))
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .false., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TRestart) :: chk
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=0)
   call env%checkpoint("failed setup")
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 5.5384030361036_wp, thr=thr)
   call check_(error, energy,-8.6908532561691_wp, thr=thr)
   call check_(error, norm2(gradient),0.24709483673929E-01_wp, thr=thr)

   call check_(error, gradient(1,3), 0.12460501676405E-02_wp, thr=thr)
   call check_(error, gradient(2,1),-0.20307939298595E-01_wp, thr=thr)
   call check_(error, gradient(1,6),-0.29075269403059E-02_wp, thr=thr)
   call check_(error, gradient(3,5), 0.00000000000000E+00_wp, thr=thr)

end subroutine test_gfn0_api

subroutine test_gfn0_api_srb(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout

   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_data
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_environment

   use xtb_pbc_tools

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   real(wp),parameter :: thr = 1.0e-7_wp
   integer, parameter :: nat = 24
   integer, parameter :: at(nat) = [6,7,6,7,6,6,6,8,7,6,8,7,6,6, &
      &                             1,1,1,1,1,1,1,1,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
      &[ 2.02799738646442_wp,  0.09231312124713_wp, -0.14310895950963_wp, &
      &  4.75011007621000_wp,  0.02373496014051_wp, -0.14324124033844_wp, &
      &  6.33434307654413_wp,  2.07098865582721_wp, -0.14235306905930_wp, &
      &  8.72860718071825_wp,  1.38002919517619_wp, -0.14265542523943_wp, &
      &  8.65318821103610_wp, -1.19324866489847_wp, -0.14231527453678_wp, &
      &  6.23857175648671_wp, -2.08353643730276_wp, -0.14218299370797_wp, &
      &  5.63266886875962_wp, -4.69950321056008_wp, -0.13940509630299_wp, &
      &  3.44931709749015_wp, -5.48092386085491_wp, -0.14318454855466_wp, &
      &  7.77508917214346_wp, -6.24427872938674_wp, -0.13107140408805_wp, &
      & 10.30229550927022_wp, -5.39739796609292_wp, -0.13672168520430_wp, &
      & 12.07410272485492_wp, -6.91573621641911_wp, -0.13666499342053_wp, &
      & 10.70038521493902_wp, -2.79078533715849_wp, -0.14148379504141_wp, &
      & 13.24597858727017_wp, -1.76969072232377_wp, -0.14218299370797_wp, &
      &  7.40891694074004_wp, -8.95905928176407_wp, -0.11636933482904_wp, &
      &  1.38702118184179_wp,  2.05575746325296_wp, -0.14178615122154_wp, &
      &  1.34622199478497_wp, -0.86356704498496_wp,  1.55590600570783_wp, &
      &  1.34624089204623_wp, -0.86133716815647_wp, -1.84340893849267_wp, &
      &  5.65596919189118_wp,  4.00172183859480_wp, -0.14131371969009_wp, &
      & 14.67430918222276_wp, -3.26230980007732_wp, -0.14344911021228_wp, &
      & 13.50897177220290_wp, -0.60815166181684_wp,  1.54898960808727_wp, &
      & 13.50780014200488_wp, -0.60614855212345_wp, -1.83214617078268_wp, &
      &  5.41408424778406_wp, -9.49239668625902_wp, -0.11022772492007_wp, &
      &  8.31919801555568_wp, -9.74947502841788_wp,  1.56539243085954_wp, &
      &  8.31511620712388_wp, -9.76854236502758_wp, -1.79108242206824_wp],&
      &  shape(xyz))
   type(peeq_options),parameter :: opt = peeq_options( &
      &  prlevel = 2, ccm = .false., acc = 1.0_wp, etemp = 300.0_wp, grad = .true. )

   type(TMolecule)    :: mol
   type(TEnvironment) :: env
   type(TRestart) :: chk
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! setup the environment variables
   call init(env)

   call init(mol, at, xyz)

   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=0)
   call env%checkpoint("failed setup")
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   call check_(error, hl_gap, 3.1192454818777_wp, thr=thr)
   call check_(error, energy,-40.908850483831_wp, thr=thr)
   call check_(error, norm2(gradient),0.83573853247451E-01_wp, thr=thr)

   call check_(error, gradient(1, 3), 0.29870640452823E-01_wp, thr=thr)
   call check_(error, gradient(2, 1), 0.67498165159819E-03_wp, thr=thr)
   call check_(error, gradient(1, 6),-0.54145817368542E-02_wp, thr=thr)
   call check_(error, gradient(3, 5),-0.39879097860049E-04_wp, thr=thr)
   call check_(error, gradient(1, 8),-0.63393710255449E-02_wp, thr=thr)
   call check_(error, gradient(2, 6), 0.19478157220852E-01_wp, thr=thr)
   call check_(error, gradient(1,20), 0.10070932785359E-03_wp, thr=thr)
   call check_(error, gradient(3,10),-0.54834076982200E-04_wp, thr=thr)

end subroutine test_gfn0_api_srb


subroutine test_gfn0_mindless_basic(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator, newXTBCalculator, newWavefunction

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
   type(scc_results) :: res

   integer :: iMol
   logical :: exitRun
   real(wp) :: energy, hl_gap, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :)

   character(len=*), parameter :: mindless(10) = [&
      & "mindless01", "mindless02", "mindless03", "mindless04", "mindless05", &
      & "mindless06", "mindless07", "mindless08", "mindless09", "mindless10"]
   real(wp), parameter :: ref_energies(10) = &
      &[-28.497588487382_wp, -23.254292745575_wp, -23.060461949336_wp, &
      & -21.822909541168_wp, -27.551996317693_wp, -18.814024099835_wp, &
      & -31.195804737400_wp, -27.317394865162_wp, -20.258074680634_wp, &
      & -25.051348432577_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.05890968017315_wp, 0.04174928560148_wp, 0.05768263390495_wp, &
      & 0.06591502692791_wp, 0.08382457849620_wp, 0.09995151813704_wp, &
      & 0.08343577915160_wp, 0.08653792253034_wp, 0.06328938755566_wp, &
      & 0.06269022672416_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[4.9355819131481_wp, 0.5817816809422_wp, 2.4534064694863_wp, &
      & 0.8490082384623_wp, 1.4811843296726_wp, 0.9873933294559_wp, &
      & 4.4138881278461_wp, 1.1871661530314_wp, 1.9816408428452_wp, &
      & 0.3432085685632_wp]

   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, method=0)
      call newWavefunction(env, mol, calc, chk)

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
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn0_mindless_basic


subroutine test_gfn0_mindless_solvation(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_test_molstock, only : getMolecule

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_pcem
   use xtb_type_data, only : scc_results
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_restart, only : TRestart

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction, addSolvationModel
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel

   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: thr = 1.0e-8_wp

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TxTBCalculator) :: calc
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
      &[-28.510632057897_wp, -23.266701495294_wp, -23.086125389488_wp, &
      & -21.822303675300_wp, -27.566928099653_wp, -18.832541288587_wp, &
      & -31.204815411606_wp, -27.323931755616_wp, -20.270904021114_wp, &
      & -25.053546425391_wp]
   real(wp), parameter :: ref_gnorms(10) = &
      &[0.05788277333095_wp, 0.04117813956764_wp, 0.05838285007764_wp, &
      & 0.06675003613700_wp, 0.08237342790744_wp, 0.10215598499722_wp, &
      & 0.08394783095013_wp, 0.086507389465479_wp, 0.06283764829270_wp, &
      & 0.06263572703682_wp]
   real(wp), parameter :: ref_hlgaps(10) = &
      &[4.9494566269479_wp, 0.5809055820838_wp, 2.4549057899374_wp, &
      & 0.8522898232117_wp, 1.4827704526244_wp, 0.9710672600740_wp, &
      & 4.3955271628719_wp, 1.1868022513933_wp, 1.9777106058585_wp, &
      & 0.3654557534774_wp]


   call init(env)
   do iMol = 1, 10

      call getMolecule(mol, mindless(iMol))

      if (allocated(gradient)) deallocate(gradient)
      allocate(gradient(3, len(mol)))

      call newXTBCalculator(env, mol, calc, method=0)
      call newWavefunction(env, mol, calc, chk)
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
      call check_(error, hl_gap, ref_hlgaps(iMol), thr=thr)

   end do

end subroutine test_gfn0_mindless_solvation

!> The built-in GFN0-xTB parameters must reproduce the parameter file
subroutine test_gfn0_builtin_param(error)
   use xtb_mctc_systools, only : rdpath
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_param, only : TxTBParameter
   use xtb_xtb_data, only : TxTBData
   use xtb_readparam, only : readParam
   use xtb_paramset, only : use_parameterset

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   character(len=*), parameter :: pfile = 'param_gfn0-xtb.txt'
   type(TEnvironment) :: env
   type(TxTBParameter) :: globpar
   type(TxTBData) :: builtin, reference
   character(len=:), allocatable :: fname
   logical :: exist, exitRun
   integer :: ipar

   call init(env)

   call use_parameterset(pfile, globpar, builtin, exist)
   call check_(error, exist, "Built-in GFN0-xTB parameters are not available")
   if (allocated(error)) return

   call rdpath(env%xtbpath, pfile, fname, exist)
   if (.not.exist) then
      call test_failed(error, "Parameter file "//pfile//" not found in XTBPATH")
      return
   end if
   call open_file(ipar, fname, 'r')
   call readParam(env, ipar, globpar, reference, .true.)
   call close_file(ipar)
   call env%check(exitRun)
   call check_(error, .not.exitRun, "Reading the parameter file failed")
   if (allocated(error)) return

   call check_xtb_data(error, builtin, reference)

end subroutine test_gfn0_builtin_param


!> Compare all parameters of two xTB data sets for the elements of the
!> built-in data
subroutine check_xtb_data(error, actual, expected)
   use xtb_xtb_data, only : TxTBData

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Parametrisation to check
   type(TxTBData), intent(in) :: actual

   !> Reference parametrisation
   type(TxTBData), intent(in) :: expected

   integer :: nElem

   nElem = size(actual%nShell)

   call check_(error, actual%name, expected%name, more="name")
   if (allocated(error)) return
   call check_(error, actual%level, expected%level, more="level")
   if (allocated(error)) return
   call check_(error, actual%ipeashift, expected%ipeashift, more="ipeashift", &
      & thr=thr_param)
   if (allocated(error)) return
   call check_param(error, "nShell", actual%nShell, expected%nShell, nElem)
   if (allocated(error)) return

   associate(act => actual%repulsion, ref => expected%repulsion)
      call check_(error, act%kExp, ref%kExp, more="repulsion%kExp", thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%kExpLight, ref%kExpLight, &
         & more="repulsion%kExpLight", thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%rExp, ref%rExp, more="repulsion%rExp", thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%enScale, ref%enScale, more="repulsion%enScale", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%cutoff, ref%cutoff, more="repulsion%cutoff", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_param(error, "repulsion%alpha", act%alpha, ref%alpha, nElem)
      if (allocated(error)) return
      call check_param(error, "repulsion%zeff", act%zeff, ref%zeff, nElem)
      if (allocated(error)) return
      call check_param(error, "repulsion%electronegativity", &
         & act%electronegativity, ref%electronegativity, nElem)
      if (allocated(error)) return
   end associate

   associate(act => actual%coulomb, ref => expected%coulomb)
      call check_(error, act%gExp, ref%gExp, more="coulomb%gExp", thr=thr_param)
      if (allocated(error)) return
      call check_param(error, "coulomb%chemicalHardness", &
         & act%chemicalHardness, ref%chemicalHardness, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%shellHardness", act%shellHardness, &
         & ref%shellHardness, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%thirdOrderAtom", act%thirdOrderAtom, &
         & ref%thirdOrderAtom, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%thirdOrderShell", act%thirdOrderShell, &
         & ref%thirdOrderShell, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%chargeWidth", act%chargeWidth, &
         & ref%chargeWidth, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%electronegativity", &
         & act%electronegativity, ref%electronegativity, nElem)
      if (allocated(error)) return
      call check_param(error, "coulomb%kCN", act%kCN, ref%kCN, nElem)
      if (allocated(error)) return
   end associate

   associate(act => actual%dispersion, ref => expected%dispersion)
      call check_(error, act%dpar%s6, ref%dpar%s6, more="dispersion%s6", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%s8, ref%dpar%s8, more="dispersion%s8", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%s10, ref%dpar%s10, more="dispersion%s10", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%s9, ref%dpar%s9, more="dispersion%s9", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%a1, ref%dpar%a1, more="dispersion%a1", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%a2, ref%dpar%a2, more="dispersion%a2", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%dpar%alp, ref%dpar%alp, more="dispersion%alp")
      if (allocated(error)) return
      call check_(error, act%g_a, ref%g_a, more="dispersion%g_a", thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%g_c, ref%g_c, more="dispersion%g_c", thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%wf, ref%wf, more="dispersion%wf", thr=thr_param)
      if (allocated(error)) return
   end associate

   call check_(error, allocated(actual%srb), allocated(expected%srb), &
      & more="allocation of srb")
   if (allocated(error)) return
   if (allocated(actual%srb)) then
      call check_(error, actual%srb%shift, expected%srb%shift, &
         & more="srb%shift", thr=thr_param)
      if (allocated(error)) return
      call check_(error, actual%srb%prefactor, expected%srb%prefactor, &
         & more="srb%prefactor", thr=thr_param)
      if (allocated(error)) return
      call check_(error, actual%srb%steepness, expected%srb%steepness, &
         & more="srb%steepness", thr=thr_param)
      if (allocated(error)) return
      call check_(error, actual%srb%enScale, expected%srb%enScale, &
         & more="srb%enScale", thr=thr_param)
      if (allocated(error)) return
   end if
   call check_(error, allocated(actual%halogen), allocated(expected%halogen), &
      & more="allocation of halogen")
   if (allocated(error)) return
   call check_(error, allocated(actual%multipole), &
      & allocated(expected%multipole), more="allocation of multipole")
   if (allocated(error)) return

   associate(act => actual%hamiltonian, ref => expected%hamiltonian)
      call check_real_matrix(error, "hamiltonian%kScale", act%kScale, &
         & ref%kScale, size(ref%kScale, 2))
      if (allocated(error)) return
      call check_(error, act%kDiff, ref%kDiff, more="hamiltonian%kDiff", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_real_matrix(error, "hamiltonian%enScale", act%enScale, &
         & ref%enScale, size(ref%enScale, 2))
      if (allocated(error)) return
      call check_(error, act%enScale4, ref%enScale4, more="hamiltonian%enScale4", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_(error, act%wExp, ref%wExp, more="hamiltonian%wExp", &
         & thr=thr_param)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%principalQuantumNumber", &
         & act%principalQuantumNumber, ref%principalQuantumNumber, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%angShell", act%angShell, &
         & ref%angShell, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%valenceShell", act%valenceShell, &
         & ref%valenceShell, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%numberOfPrimitives", &
         & act%numberOfPrimitives, ref%numberOfPrimitives, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%slaterExponent", act%slaterExponent, &
         & ref%slaterExponent, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%selfEnergy", act%selfEnergy, &
         & ref%selfEnergy, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%referenceOcc", act%referenceOcc, &
         & ref%referenceOcc, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%kCN", act%kCN, ref%kCN, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%electronegativity", &
         & act%electronegativity, ref%electronegativity, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%atomicRad", act%atomicRad, &
         & ref%atomicRad, nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%shellPoly", act%shellPoly, &
         & ref%shellPoly, nElem)
      if (allocated(error)) return
      ! pair parameters are indexed by elements in both dimensions
      call check_param(error, "hamiltonian%pairParam", act%pairParam, &
         & ref%pairParam, nElem, nrow=nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%kQShell", act%kQShell, ref%kQShell, &
         & nElem)
      if (allocated(error)) return
      call check_param(error, "hamiltonian%kQAtom", act%kQAtom, ref%kQAtom, &
         & nElem)
      if (allocated(error)) return
   end associate

end subroutine check_xtb_data


!> Compare the first n entries of two optional per-element real arrays
subroutine check_param_real1(error, name, actual, expected, n)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   real(wp), allocatable, intent(in) :: actual(:)

   !> Reference values
   real(wp), allocatable, intent(in) :: expected(:)

   !> Number of entries to compare
   integer, intent(in) :: n

   call check_(error, allocated(actual), allocated(expected), &
      & more="allocation of "//name)
   if (allocated(error) .or. .not.allocated(actual)) return

   call check_real_matrix(error, name, reshape(actual, [1, size(actual)]), &
      & reshape(expected, [1, size(expected)]), n)

end subroutine check_param_real1


!> Compare the first n columns of two optional per-element real arrays
subroutine check_param_real2(error, name, actual, expected, n, nrow)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   real(wp), allocatable, intent(in) :: actual(:, :)

   !> Reference values
   real(wp), allocatable, intent(in) :: expected(:, :)

   !> Number of columns to compare
   integer, intent(in) :: n

   !> Number of rows to compare, by default both arrays must have the same rows
   integer, intent(in), optional :: nrow

   call check_(error, allocated(actual), allocated(expected), &
      & more="allocation of "//name)
   if (allocated(error) .or. .not.allocated(actual)) return

   if (present(nrow)) then
      call check_(error, min(size(actual, 1), size(expected, 1)) >= nrow, &
         & more="too few rows in "//name)
      if (allocated(error)) return
      call check_real_matrix(error, name, actual(:nrow, :), expected(:nrow, :), n)
   else
      call check_real_matrix(error, name, actual, expected, n)
   end if

end subroutine check_param_real2


!> Compare the first n entries of two optional per-element integer arrays
subroutine check_param_int1(error, name, actual, expected, n)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   integer, allocatable, intent(in) :: actual(:)

   !> Reference values
   integer, allocatable, intent(in) :: expected(:)

   !> Number of entries to compare
   integer, intent(in) :: n

   call check_(error, allocated(actual), allocated(expected), &
      & more="allocation of "//name)
   if (allocated(error) .or. .not.allocated(actual)) return

   call check_int_matrix(error, name, reshape(actual, [1, size(actual)]), &
      & reshape(expected, [1, size(expected)]), n)

end subroutine check_param_int1


!> Compare the first n columns of two optional per-element integer arrays
subroutine check_param_int2(error, name, actual, expected, n)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   integer, allocatable, intent(in) :: actual(:, :)

   !> Reference values
   integer, allocatable, intent(in) :: expected(:, :)

   !> Number of columns to compare
   integer, intent(in) :: n

   call check_(error, allocated(actual), allocated(expected), &
      & more="allocation of "//name)
   if (allocated(error) .or. .not.allocated(actual)) return

   call check_int_matrix(error, name, actual, expected, n)

end subroutine check_param_int2


!> Compare the first n columns of two real arrays with the same number of rows
subroutine check_real_matrix(error, name, actual, expected, n)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   real(wp), intent(in) :: actual(:, :)

   !> Reference values
   real(wp), intent(in) :: expected(:, :)

   !> Number of columns to compare
   integer, intent(in) :: n

   integer :: ii, jj

   call check_(error, size(actual, 1), size(expected, 1), &
      & more="number of rows in "//name)
   if (allocated(error)) return
   call check_(error, min(size(actual, 2), size(expected, 2)) >= n, &
      & more="too few entries in "//name)
   if (allocated(error)) return

   do jj = 1, n
      do ii = 1, size(actual, 1)
         call check_(error, actual(ii, jj), expected(ii, jj), thr=thr_param, &
            & more=name//"("//to_string(ii)//", "//to_string(jj)//")")
         if (allocated(error)) return
      end do
   end do

end subroutine check_real_matrix


!> Compare the first n columns of two integer arrays with the same number of rows
subroutine check_int_matrix(error, name, actual, expected, n)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Name of the compared parameter
   character(len=*), intent(in) :: name

   !> Values to check
   integer, intent(in) :: actual(:, :)

   !> Reference values
   integer, intent(in) :: expected(:, :)

   !> Number of columns to compare
   integer, intent(in) :: n

   integer :: ii, jj

   call check_(error, size(actual, 1), size(expected, 1), &
      & more="number of rows in "//name)
   if (allocated(error)) return
   call check_(error, min(size(actual, 2), size(expected, 2)) >= n, &
      & more="too few entries in "//name)
   if (allocated(error)) return

   do jj = 1, n
      do ii = 1, size(actual, 1)
         call check_(error, actual(ii, jj), expected(ii, jj), &
            & more=name//"("//to_string(ii)//", "//to_string(jj)//")")
         if (allocated(error)) return
      end do
   end do

end subroutine check_int_matrix


end module test_gfn0
