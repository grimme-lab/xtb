! This file is part of xtb.
!
! Copyright (C) 2021 Sebastian Ehlert
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

module test_hessian
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : autoaa
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction
   use xtb_modelhessian, only : mh_swart
   use xtb_type_setvar, only : modhess_setvar
   use xtb_o1numhess, only : gen_displdir, get_neighbor_list, adj_list, gen_local_hessian, lr_loop, swart
   implicit none
   private

   public :: collect_hessian


contains

function get_fixture_path(filename) result(fullpath)
   character(len=*), intent(in) :: filename
   character(len=:), allocatable :: fullpath
   character(len=1024) :: base_dir
   integer :: length, stat
   
   ! Get variable; status 0 means it exists
   call get_environment_variable("TEST_FIXTURES_DIR", base_dir, length, stat)
   
   if (stat == 0) then
      ! We effectively perform: "path/to/dir" // "/" // "file.txt"
      ! On Windows this creates: "C:/path/to/dir/file.txt" (Valid)
      ! On Linux this creates: "/path/to/dir/file.txt" (Valid)
      fullpath = trim(base_dir) // '/' // filename
   else
      ! Fallback for manual running (assumes running from build dir or similar)
      fullpath = "../test/unit/fixtures/" // filename
   end if
end function get_fixture_path

!> Collect all exported unit tests
subroutine collect_hessian(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1_hessian", test_gfn1_hessian), &
      new_unittest("gfn2_hessian", test_gfn2_hessian), &
      new_unittest("nblist_o1numhess", test_nblist_o1numhess), &
      new_unittest("modified_swart", test_modified_swart), &
      new_unittest("o1numhess_gen_displdir", test_gendispldir_o1numhess), &
      new_unittest("gfn1_o1numhess", test_o1numhess_gfn1), &
      new_unittest("gfn2_o1numhess", test_o1numhess_gfn2), &
      new_unittest("linear_h2o_gfn1_o1numhess", test_o1numhess_linear_h2o_gfn1), &
      new_unittest("linear_h2o_gfn2_o1numhess", test_o1numhess_linear_h2o_gfn2) &
      ]

end subroutine collect_hessian

subroutine test_gfn1_hessian(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: dipgrad_ref(3, 3*nat) = reshape([ &
      & -1.013452580143500E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  3.527708416572683E-11_wp, -0.580105714422270E+00_wp, -4.440892098500626E-10_wp, &
      &  1.640601522247707E-11_wp, -1.577916091825909E-09_wp, -0.719923066050399E+00_wp, &
      &  0.506726290141942E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  5.168504173660994E-11_wp,  0.290827888051893E+00_wp,  7.276368407804767E-02_wp, &
      & -3.527728552279606E-11_wp,  0.159475763089199E+00_wp,  0.359943068684032E+00_wp, &
      &  0.506726290062044E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  0.000000000000000E+00_wp,  0.290827888619521E+00_wp, -7.276368707564984E-02_wp, &
      &  0.000000000000000E+00_wp, -0.159475761295180E+00_wp,  0.359943068239943E+00_wp],&
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      &  9.642724315151719E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -4.821349241627352E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -4.821375073522908E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -5.183734810089242E-12_wp,  0.654328674106580E+00_wp, -8.411850317864566E-10_wp, &
      &  4.912366132678069E-12_wp, -0.327164337503372E+00_wp,  0.241508066094075E+00_wp, &
      &  2.713686774111726E-13_wp, -0.327164336592530E+00_wp, -0.241508065247398E+00_wp, &
      &  1.269389491007523E-11_wp, -1.003642487147619E-09_wp,  0.396956797942315E+00_wp, &
      & -1.357110357831859E-11_wp,  0.193314810701155E+00_wp, -0.198478399246260E+00_wp, &
      &  8.772086682433594E-13_wp, -0.193314809690279E+00_wp, -0.198478398680934E+00_wp, &
      & -4.821349763250659E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  4.418957886905253E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  4.023918763456385E-06_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  7.510234345234422E-12_wp, -0.327069271103985E+00_wp,  0.193419608180277E+00_wp, &
      & -8.658575901216439E-12_wp,  0.348666988822030E+00_wp, -0.217428743767287E+00_wp, &
      &  1.148341555982017E-12_wp, -2.159771772581407E-02_wp,  2.400913557787203E-02_wp, &
      &  5.183777743453684E-12_wp,  0.241435744915371E+00_wp, -0.198481745303595E+00_wp, &
      & -4.912396810694304E-12_wp, -0.217376911834996E+00_wp,  0.188386308225154E+00_wp, &
      & -2.713809327593799E-13_wp, -2.405883307999757E-02_wp,  1.009543707641887E-02_wp, &
      & -4.821379744720639E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  4.023966094020255E-06_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  4.418983135319131E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  0.000000000000000E+00_wp, -0.327069270111905E+00_wp, -0.193419606975667E+00_wp, &
      &  0.000000000000000E+00_wp, -2.159771757181671E-02_wp, -2.400913565342059E-02_wp, &
      &  0.000000000000000E+00_wp,  0.348666987689540E+00_wp,  0.217428742639103E+00_wp, &
      &  0.000000000000000E+00_wp, -0.241435744067920E+00_wp, -0.198481744759780E+00_wp, &
      &  0.000000000000000E+00_wp,  2.405883312826490E-02_wp,  1.009543692497955E-02_wp, &
      &  0.000000000000000E+00_wp,  0.217376910941702E+00_wp,  0.188386307847977E+00_wp],&
      &  shape(hessian_ref))
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i,j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:), dipgrad(:,:), hessian(:,:)
   integer, allocatable :: list(:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3,mol%n), dipgrad(3,3*mol%n), hessian(3*mol%n,3*mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=1)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   dipgrad = 0.0_wp
   hessian = 0.0_wp
   list = [(i, i = 1, mol%n)]
   call calc%hessian(env, mol, chk, list, step, hessian, dipgrad)

   do i = 1, size(dipgrad_ref, 2)
      do j = 1, size(dipgrad_ref, 1)
         call check(error, dipgrad(j, i), dipgrad_ref(j, i), thr=thr)
      end do
   end do

   do i = 1, size(hessian_ref, 2)
      do j = 1, size(hessian_ref, 1)
         call check(error, hessian(j, i), hessian_ref(j, i), thr=thr)
      end do
   end do

end subroutine test_gfn1_hessian

subroutine test_gfn2_hessian(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,   -0.00000000077760_wp,    0.18829790750029_wp, &
      & 0.00000000000000_wp,    1.45987612440076_wp,   -0.88615189669760_wp, &
      &-0.00000000000000_wp,   -1.45987612362316_wp,   -0.88615189608629_wp],&
      & shape(xyz))
   real(wp), parameter :: dipgrad_ref(3, 3*nat) = reshape([ &
      & -0.811313733750829E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -3.379448427312665E-11_wp, -0.412116147045635E+00_wp, -2.220446049250313E-09_wp, &
      & -1.688661825390276E-11_wp,  1.035790186933606E-09_wp, -0.483878709767183E+00_wp, &     
      &  0.405648009573818E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  7.146633102892883E-11_wp,  0.205575204598719E+00_wp,  9.503822501200077E-02_wp, &
      & -1.829769461273858E-14_wp,  0.148069825228005E+00_wp,  0.243425242296702E+00_wp, &     
      &  0.405648009946286E+00_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -4.559890805138046E-14_wp,  0.205575207063050E+00_wp, -9.503822462342271E-02_wp, &
      &  3.356458019834560E-11_wp, -0.148069829195679E+00_wp,  0.243425244239592E+00_wp],&
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      & -1.939596096290860E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  9.697686009794110E-06_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  9.698274953126006E-06_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -1.978173339857878E-11_wp,  0.612198809002427E+00_wp,  2.697691516787698E-09_wp, &
      &  2.481358892927353E-11_wp, -0.306099403183025E+00_wp,  0.225852148702132E+00_wp, &
      & -5.031855530694753E-12_wp, -0.306099405820178E+00_wp, -0.225852151393925E+00_wp, &
      & -1.024802643153815E-11_wp,  2.163320452807998E-09_wp,  0.399810948644134E+00_wp, &
      &  1.255439924419533E-11_wp,  0.184154284413794E+00_wp, -0.199905473616475E+00_wp, &
      & -2.306372812657182E-12_wp, -0.184154286575537E+00_wp, -0.199905475030756E+00_wp, &
      & -6.559552217481455E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  5.182604735260440E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  1.376947482221124E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  3.321776245993920E-11_wp, -0.305714839568990E+00_wp,  0.183399699337614E+00_wp, &
      & -3.597681672480956E-11_wp,  0.339222824678649E+00_wp, -0.204443891042213E+00_wp, &
      &  2.759054264870358E-12_wp, -3.350798510798255E-02_wp,  2.104419170521866E-02_wp, &
      & -2.747668766007676E-13_wp,  0.225612945682292E+00_wp, -0.199989952327133E+00_wp, &
      &  2.050024431283972E-13_wp, -0.204892631012658E+00_wp,  0.183497773686177E+00_wp, &
      &  6.976443347237028E-14_wp, -2.072031467084191E-02_wp,  1.649217862695682E-02_wp, &
      & -6.559491291679194E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  1.376944591474809E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      &  5.182546700204432E-05_wp,  0.000000000000000E+00_wp,  0.000000000000000E+00_wp, &
      & -2.703865614441994E-13_wp, -0.305714842513751E+00_wp, -0.183399702402830E+00_wp, &
      &  2.719647608323040E-13_wp, -3.350798515432542E-02_wp, -2.104419150402123E-02_wp, &
      & -1.578199388104547E-15_wp,  0.339222827661363E+00_wp,  0.204443893905696E+00_wp, &
      &  2.029210950346453E-11_wp, -0.225612948054816E+00_wp, -0.199989953443651E+00_wp, &
      & -2.499966571670633E-11_wp,  2.072031476318205E-02_wp,  1.649217857736473E-02_wp, &
      &  4.707556213241799E-12_wp,  0.204892633286792E+00_wp,  0.183497774863884E+00_wp],&
      &  shape(hessian_ref))

   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i,j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:), dipgrad(:,:), hessian(:,:)
   integer, allocatable :: list(:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3,mol%n), dipgrad(3,3*mol%n), hessian(3*mol%n,3*mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2)
   call newWavefunction(env, mol, calc, chk)

   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   dipgrad = 0.0_wp
   hessian = 0.0_wp
   list = [(i, i = 1, mol%n)]
   call calc%hessian(env, mol, chk, list, step, hessian, dipgrad)

   do i = 1, size(dipgrad_ref, 2)
      do j = 1, size(dipgrad_ref, 1)
         call check(error, dipgrad(j, i), dipgrad_ref(j, i), thr=thr)
      end do
   end do

   do i = 1, size(hessian_ref, 2)
      do j = 1, size(hessian_ref, 1)
         call check(error, hessian(j, i), hessian_ref(j, i), thr=thr)
      end do
   end do

end subroutine test_gfn2_hessian

subroutine test_nblist_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 54
   real(wp), parameter :: thr = 1.0e-4_wp
   character(len=*), parameter :: sym(nat) = ['O', 'O', 'O', 'O', 'O', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'H']
   integer, parameter :: at(nat) = [8, 8, 8, 8, 8, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      &-3.06876197244764_wp,   0.19448556443832_wp, -1.94335516402362_wp, &
      & 0.48672317087736_wp,  -2.37414977416329_wp,  0.74640403964754_wp, &
      &-5.01392286848523_wp,   0.06500518477614_wp,  0.99438680911247_wp, &
      &-5.92926707825204_wp,  -2.14964638902524_wp, -0.09679054786651_wp, &
      & 6.46124291527335_wp,  -0.36657109114947_wp, -0.97517530474911_wp, &
      &-0.99293297479044_wp,   1.21448854000979_wp, -0.25678568019737_wp, &
      &-1.88187467277990_wp,   0.02222909651218_wp,  0.13250075340634_wp, &
      & 0.40066850753014_wp,   1.13354627833840_wp,  0.35239280297367_wp, &
      & 1.06360061679260_wp,  -0.18376426171793_wp, -0.10101865874150_wp, &
      &-3.21584077209344_wp,   0.45371531410135_wp, -0.55983996021412_wp, &
      & 2.50598996635889_wp,  -0.36661297069103_wp,  0.44518293352787_wp, &
      &-1.85548086158727_wp,   2.43127007566346_wp,  0.07859124220201_wp, &
      &-1.24286384401417_wp,  -1.24793995135794_wp, -0.42270418964592_wp, &
      &-3.28684067393579_wp,   1.96830601326335_wp, -0.25245092803743_wp, &
      & 1.26790003151346_wp,   2.30531243156587_wp, -0.10078672435379_wp, &
      & 0.15407031860948_wp,  -1.37873895338713_wp,  0.15531287229906_wp, &
      &-2.07585892781203_wp,  -0.14181953375097_wp,  1.63925566896310_wp, &
      & 3.31294973541281_wp,   0.88258965214860_wp,  0.16968950515907_wp, &
      & 2.66729325199697_wp,   2.19367639453615_wp,  0.49736188521583_wp, &
      & 3.14556493096524_wp,  -1.55982985828490_wp, -0.28809822699537_wp, &
      & 2.51422839817214_wp,  -0.60483849487020_wp,  1.96446208259101_wp, &
      &-4.41921521111641_wp,  -0.28742042420872_wp,  0.00085168234654_wp, &
      & 4.65019175446183_wp,  -1.65135422072961_wp, -0.07463306598685_wp, &
      & 4.57013308085424_wp,   0.84147357516259_wp, -0.28491650710015_wp, &
      & 5.34019312592334_wp,  -0.38176061789806_wp, -0.51245960465392_wp, &
      &-4.85468325039185_wp,  -1.53235168144890_wp, -0.73618634233715_wp, &
      &-0.88536000505000_wp,   1.17775196980949_wp, -1.34698146609481_wp, &
      & 0.32451699483681_wp,   1.15980925138336_wp,  1.44214229435536_wp, &
      & 1.14530302264227_wp,  -0.12136956510912_wp, -1.19525159153613_wp, &
      &-1.57456838689741_wp,   3.30265656519013_wp, -0.50994722436245_wp, &
      &-1.76218625627092_wp,   2.69093642014312_wp,  1.13270257354206_wp, &
      &-1.18927493326531_wp,  -1.18835746298682_wp, -1.51046637912891_wp, &
      &-1.80784940235908_wp,  -2.13948164603414_wp, -0.14598752625224_wp, &
      &-3.96611761537147_wp,   2.13411373620025_wp,  0.58290206904660_wp, &
      &-3.68591892958296_wp,   2.49809302066445_wp, -1.11882955134477_wp, &
      & 0.81387859146875_wp,   3.24645064346949_wp,  0.21243822059035_wp, &
      & 1.33729877444188_wp,   2.31000542053058_wp, -1.19098039484345_wp, &
      &-1.12118685681340_wp,  -0.32162061906065_wp,  2.12639045628194_wp, &
      &-2.54463670293864_wp,   0.72913296780967_wp,  2.08567527228266_wp, &
      &-2.70911977624777_wp,  -1.00248434251682_wp,  1.84262485728548_wp, &
      & 2.59391711447919_wp,   2.27953745210943_wp,  1.58574891610354_wp, &
      & 3.29343891115649_wp,   3.01461931811614_wp,  0.14480636755127_wp, &
      & 2.95983104730468_wp,  -1.45350131071300_wp, -1.36092075616257_wp, &
      & 2.67353385789435_wp,  -2.48394831537700_wp,  0.04336940312963_wp, &
      & 3.53564571272084_wp,  -0.60633642649981_wp,  2.33669404875740_wp, &
      & 1.97081004589295_wp,   0.17618341572015_wp,  2.48911941746077_wp, &
      & 2.05619026864024_wp,  -1.56170278103109_wp,  2.19420595639869_wp, &
      &-3.73057965515078_wp,   0.69737538781391_wp, -2.43068227224802_wp, &
      & 5.06838675345177_wp,  -2.49659445557115_wp, -0.62352909086624_wp, &
      & 4.88048000733017_wp,  -1.80335421693271_wp,  0.98284255140168_wp, &
      & 5.10835341462563_wp,   1.76010665629942_wp, -0.47938272166668_wp, &
      &-4.00084962479015_wp,  -2.21956028206792_wp, -0.80500973828872_wp, &
      &-5.13568412213777_wp,  -1.26351703151787_wp, -1.76363779186583_wp, &
      &-6.13475894704690_wp,  -1.61304366767513_wp,  0.68825272793312_wp],&
      & shape(xyz))
   real(wp), parameter :: vdw_radii(1:103) = [ &
      2.886_wp, 2.362_wp, 2.451_wp, 2.745_wp, 4.083_wp, 3.851_wp, 3.66_wp, 3.5_wp, 3.364_wp, &
      3.243_wp, 2.983_wp, 3.021_wp, 4.499_wp, 4.295_wp, 4.147_wp, 4.035_wp, 3.947_wp, 3.868_wp, &
      3.812_wp, 3.399_wp, 3.295_wp, 3.175_wp, 3.144_wp, 3.023_wp, 2.961_wp, 2.912_wp, 2.872_wp, &
      2.834_wp, 3.495_wp, 2.763_wp, 4.383_wp, 4.28_wp, 4.23_wp, 4.205_wp, 4.189_wp, 4.141_wp, &
      4.114_wp, 3.641_wp, 3.345_wp, 3.124_wp, 3.165_wp, 3.052_wp, 2.998_wp, 2.963_wp, 2.929_wp, &
      2.899_wp, 3.148_wp, 2.848_wp, 4.463_wp, 4.392_wp, 4.42_wp, 4.47_wp, 4.5_wp, 4.404_wp, &
      4.517_wp, 3.703_wp, 3.522_wp, 3.556_wp, 3.606_wp, 3.575_wp, 3.547_wp, 3.52_wp, 3.493_wp, &
      3.368_wp, 3.451_wp, 3.428_wp, 3.409_wp, 3.391_wp, 3.374_wp, 3.355_wp, 3.64_wp, 3.141_wp, &
      3.17_wp, 3.069_wp, 2.954_wp, 3.12_wp, 2.84_wp, 2.754_wp, 3.293_wp, 2.705_wp, 4.347_wp, &
      4.297_wp, 4.37_wp, 4.709_wp, 4.75_wp, 4.765_wp, 4.9_wp, 3.677_wp, 3.478_wp, 3.396_wp, &
      3.424_wp, 3.395_wp, 3.424_wp, 3.424_wp, 3.381_wp, 3.326_wp, 3.339_wp, 3.313_wp, 3.299_wp, &
      3.286_wp, 3.274_wp, 3.248_wp, 3.236_wp &
   ] / 2.0_wp / autoaa
   real(wp), allocatable :: distmat_ref(:,:)
   real(wp), parameter :: dmax = 1.0_wp, eps = 1.0e-8_wp, eps2 = 1.0e-15_wp

   type(adj_list), allocatable :: nblist_ref(:), nblist(:)
   real(wp), allocatable :: distmat(:, :)
   real(wp) :: dist
   integer :: N, i, j, iunit, nnb
   character(len=2048) :: line
   character(len=:), allocatable :: fixture_path
   logical :: is_same

   N = 3 * nat

   ! Read neighbor list from fixture file
   fixture_path = get_fixture_path("nblist_ref")
   allocate(nblist_ref(N))
   open(newunit=iunit, file=fixture_path, status='old')
   do i = 1, N
      read(iunit, '(a)') line
      ! Count commas to determine number of neighbors
      nnb = 1
      do j = 1, len_trim(line)
         if (line(j:j) == ',') nnb = nnb + 1
      end do
      allocate(nblist_ref(i)%neighbors(nnb))
      read(line, *) nblist_ref(i)%neighbors
   end do
   close(iunit)

   ! Read distance matrix from fixture file
   fixture_path = get_fixture_path("distmat_ref")
   allocate(distmat_ref(N, N))
   open(newunit=iunit, file=fixture_path, status='old')
   read(iunit, *) distmat_ref
   close(iunit)

   allocate(distmat(N, N))
   distmat = 0.0_wp
   do i = 1, nat
      do j = i, nat
         ! effective distmat
         dist = norm2(xyz(:, i) - xyz(:, j)) - vdw_radii(at(i)) - vdw_radii(at(j))
         distmat(3 * i - 2:3 * i, 3 * j - 2:3 * j) = dist
         distmat(3 * j - 2:3 * j, 3 * i - 2:3 * i) = dist
      end do
   end do

   call get_neighbor_list(distmat, dmax, nblist)

   if (any(abs(distmat - distmat_ref) > thr)) then
      call test_failed(error, "Distance Matrix not matching")
      print *, "--- Distance Matrix ---"
      do i = 1, N
         print '(*(F21.14))', distmat(i, :) 
      end do

      print *, "--- Ref. Distance Matrix ---"
      do i = 1, N
         print '(*(F21.14))', distmat_ref(i, :) 
      end do
   end if

   is_same = .true.
   do i = 1, N
      ! 1. Check if allocation status matches (both allocated or both not)
      if (allocated(nblist_ref(i)%neighbors) .neqv. allocated(nblist(i)%neighbors)) then
         is_same = .false.
      end if

      ! 2. If allocated, check size and values
      if (allocated(nblist_ref(i)%neighbors)) then

         ! Check Size
         if (size(nblist_ref(i)%neighbors) /= size(nblist(i)%neighbors)) then
            is_same = .false.
         end if

         ! Check Values using ALL() intrinsic
         if (.not. all(nblist_ref(i)%neighbors == nblist(i)%neighbors)) then
            is_same = .false.
         end if
      end if
   end do

   if (.not. is_same) then
      call test_failed(error, "Adjacency list not matching")
      print *, "--- Adjacency List ---"
      do i = 1, N
         if (allocated(nblist(i)%neighbors)) then
            ! Print the node index, then all neighbors on one line
            print '("Node ", I2, ": ", *(I4))', i, nblist(i)%neighbors
         else
            print '("Node ", I2, ": <empty>")', i
         end if
      end do

      print *, "--- Ref. Adjacency List ---"
      do i = 1, N
         if (allocated(nblist_ref(i)%neighbors)) then
            ! Print the node index, then all neighbors on one line
            print '("Node ", I2, ": ", *(I4))', i, nblist_ref(i)%neighbors
         else
            print '("Node ", I2, ": <empty>")', i
         end if
      end do
   end if

end subroutine test_nblist_o1numhess

subroutine test_modified_swart(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   integer, parameter :: at(nat) = [8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   ! swart model hessian from o1numhess utils
   real(wp), parameter :: h0(9, 9) = reshape([&
      & 1.1092305614e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp,  4.4184751144e-01_wp, -2.7552743027e-10_wp, &
      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.6310770462e-01_wp, &
      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.6310770435e-01_wp, &
      & 0.0000000000e+00_wp, -2.7552743027e-10_wp,  2.9417034746e-01_wp, &
      & 0.0000000000e+00_wp,  1.4319582461e-01_wp, -1.4708517387e-01_wp, &
      & 0.0000000000e+00_wp, -1.4319582434e-01_wp, -1.4708517359e-01_wp, &
      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 1.6724709350e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.4319582461e-01_wp, &
      & 0.0000000000e+00_wp,  2.3323207511e-01_wp, -1.5315176462e-01_wp, &
      & 0.0000000000e+00_wp, -1.2308319237e-02_wp,  9.9559400089e-03_wp, &
      & 0.0000000000e+00_wp,  1.6310770462e-01_wp, -1.4708517387e-01_wp, &
      & 0.0000000000e+00_wp, -1.5315176462e-01_wp,  1.3375368579e-01_wp, &
      & 0.0000000000e+00_wp, -9.9559400048e-03_wp,  1.3331488088e-02_wp, &
      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 1.6724709315e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.4319582434e-01_wp, &
      & 0.0000000000e+00_wp, -1.2308319237e-02_wp, -9.9559400048e-03_wp, &
      & 0.0000000000e+00_wp,  2.3323207480e-01_wp,  1.5315176434e-01_wp, &
      & 0.0000000000e+00_wp, -1.6310770435e-01_wp, -1.4708517359e-01_wp, &
      & 0.0000000000e+00_wp,  9.9559400089e-03_wp,  1.3331488088e-02_wp, &
      & 0.0000000000e+00_wp,  1.5315176434e-01_wp,  1.3375368550e-01_wp],&
      & shape(h0))

   integer :: i
   real(wp), allocatable :: hess_out(:, :)

   allocate(hess_out(3*nat, 3*nat))
   hess_out = 0.0_wp
   call swart(xyz, at, hess_out)

   if (any(abs(hess_out - h0) > thr)) then
      call test_failed(error, "Hessians do not match")

      print *, "--- hessian ---"
      do i = 1, 3*nat
         print '(*(F21.14))', hess_out(i, :) 
      end do

      print *, "--- Ref. hessian ---"
      do i = 1, 3*nat
         print '(*(F21.14))', h0(i, :) 
      end do
   end if
end subroutine test_modified_swart

subroutine test_gendispldir_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   ! swart model hessian from o1numhess utils
   real(wp), parameter :: h0(9, 9) = reshape([&
      & 1.1092305614e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp,  4.4184751144e-01_wp, -2.7552743027e-10_wp, &
      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.6310770462e-01_wp, &
      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.6310770435e-01_wp, &
      & 0.0000000000e+00_wp, -2.7552743027e-10_wp,  2.9417034746e-01_wp, &
      & 0.0000000000e+00_wp,  1.4319582461e-01_wp, -1.4708517387e-01_wp, &
      & 0.0000000000e+00_wp, -1.4319582434e-01_wp, -1.4708517359e-01_wp, &
      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 1.6724709350e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.4319582461e-01_wp, &
      & 0.0000000000e+00_wp,  2.3323207511e-01_wp, -1.5315176462e-01_wp, &
      & 0.0000000000e+00_wp, -1.2308319237e-02_wp,  9.9559400089e-03_wp, &
      & 0.0000000000e+00_wp,  1.6310770462e-01_wp, -1.4708517387e-01_wp, &
      & 0.0000000000e+00_wp, -1.5315176462e-01_wp,  1.3375368579e-01_wp, &
      & 0.0000000000e+00_wp, -9.9559400048e-03_wp,  1.3331488088e-02_wp, &
      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 1.6724709315e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.4319582434e-01_wp, &
      & 0.0000000000e+00_wp, -1.2308319237e-02_wp, -9.9559400048e-03_wp, &
      & 0.0000000000e+00_wp,  2.3323207480e-01_wp,  1.5315176434e-01_wp, &
      & 0.0000000000e+00_wp, -1.6310770435e-01_wp, -1.4708517359e-01_wp, &
      & 0.0000000000e+00_wp,  9.9559400089e-03_wp,  1.3331488088e-02_wp, &
      & 0.0000000000e+00_wp,  1.5315176434e-01_wp,  1.3375368550e-01_wp],&
      & shape(h0))
   ! Ref from o1numhess python implementation using the above h0
   real(wp), parameter :: displdir_ref(9, 3) = reshape([&
      & 0.00000000000000000e+00_wp, 5.85920171106577004e-01_wp, 1.51794937958024002e-01_wp, &
      & 0.00000000000000000e+00_wp, 4.38891753068939995e-01_wp, -2.48878406404109010e-01_wp, &
      & 0.00000000000000000e+00_wp, 5.10352159263066002e-01_wp, 3.44427193104884989e-01_wp, &
      & 2.83390078813432008e-01_wp, -1.26365402678647992e-01_wp, 5.36139605082272009e-01_wp, &
      & -2.09750793632672975e-17_wp, -2.71810289210553018e-01_wp, -4.85072467210535029e-01_wp, &
      & 4.19501587265346009e-17_wp, -3.21710376864937989e-01_wp, 4.51223308073935011e-01_wp, &
      & 5.04407285834402014e-01_wp, 4.29600509606577008e-02_wp, -1.42605077919505999e-01_wp, &
      & 5.98879533835269020e-01_wp, 2.40065553767858998e-02_wp, 4.97688397830703991e-02_wp, &
      & 5.98879533835269020e-01_wp, 2.84749410320055008e-02_wp, -4.70536852490466000e-02_wp],&
      & shape(displdir_ref))   
   real(wp), parameter :: dmax = 1.0_wp, eps = 1.0e-8_wp, eps2 = 1.0e-15_wp

   type(adj_list), allocatable :: neighborlist(:)
   real(wp), allocatable :: distmat(:, :), displdir(:, :)
   integer, allocatable :: nbcounts(:)
   integer :: i, j, N, max_nb, ndispl0, ndispl_final

   N = 3 * nat

   ! setup distmat
   allocate(distmat(N, N))
   do i = 1, N
      do j = 1, N
         distmat(i, j) = abs(i - j)
      end do
   end do

   ! compute neighbor list
   call get_neighbor_list(distmat, dmax, neighborlist)
   allocate(nbcounts(N))
   max_nb = 0
   do i = 1, N
      nbcounts(i) = size(neighborlist(i)%neighbors)
      if (nbcounts(i) > max_nb) max_nb = nbcounts(i)
   end do
   
   ! populate displdir
   ndispl0 = 0
   allocate(displdir(N, N))
   call gen_displdir(N, ndispl0, h0, max_nb, neighborlist, nbcounts, eps, eps2, displdir, ndispl_final)

   if (ndispl_final /= 3) then
      call test_failed(error, "Number of displacements not matching")
      print *, "ndispl_final: ", ndispl_final
   end if

   if (any(abs(displdir_ref - displdir(:, :ndispl_final)) > thr)) then
      call test_failed(error, "Displacements not matching")
      print *, "--- displdir ---"
      do i = 1, N
         print '(*(F21.14))', displdir(i, :ndispl_final) 
      end do

      print *, "--- Ref. displdir ---"
      do i = 1, N
         print '(*(F21.14))', displdir_ref(i, :) 
      end do
   end if

end subroutine test_gendispldir_o1numhess

subroutine test_o1numhess_gfn1(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      & 9.6427214452304158E-005_wp,  7.9724999781439290E-011_wp, -1.3840003654029018E-007_wp, &
      &-4.8213460516503025E-005_wp, -6.7504265780829648E-008_wp,  6.9205909101946085E-008_wp, & 
      &-4.8213753935792019E-005_wp,  6.7417854867650525E-008_wp,  6.9189718349509533E-008_wp, & 
      & 7.9724999781439290E-011_wp,  0.65431356973123023_wp,     -7.4720523570619281E-007_wp, & 
      & 4.4325798795663372E-011_wp, -0.32715720199319237_wp,      0.24150286357164832_wp,     & 
      &-1.2405079857710117E-010_wp, -0.32715636773669671_wp,     -0.24150211636497051_wp,     & 
      &-1.3840003654029018E-007_wp, -7.4720523570619281E-007_wp,  0.39702315113905740_wp,     & 
      & 3.6004504315957250E-008_wp,  0.19343120065992384_wp,     -0.19851179291786120_wp,     & 
      & 1.0239553222433293E-007_wp,- 0.19343045344911730_wp,     -0.19851135822773272_wp,     & 
      &-4.8213460516503025E-005_wp,  4.4325798795663372E-011_wp,  3.6004504315957250E-008_wp, & 
      & 4.4189538089979658E-005_wp,  1.9027949582154744E-009_wp, -1.8014006716533711E-008_wp, & 
      & 4.0239224265213695E-006_wp, -1.9457835743300046E-009_wp, -1.7987696743810155E-008_wp, & 
      &-6.7504265780829648E-008_wp, -0.32715720199319237_wp,      0.19343120065992384_wp,     & 
      & 1.9027949582154744E-009_wp,  0.34874368456677368_wp,     -0.21746695535952912_wp,     & 
      & 6.5601470822614196E-008_wp, -2.1586482575562844E-002_wp,  2.4035754703519811E-002_wp, & 
      & 6.9205909101946085E-008_wp,  0.24150286357164832_wp,     -0.19851179291786120_wp,     & 
      &-1.8014006716533711E-008_wp, -0.21746695535952912_wp,      0.18842607949434392_wp,     & 
      &-5.1191902385412367E-008_wp, -2.4035908212783213E-002_wp,  1.0085713431492066E-002_wp, & 
      &-4.8213753935792019E-005_wp, -1.2405079857710117E-010_wp,  1.0239553222433293E-007_wp, & 
      & 4.0239224265213695E-006_wp,  6.5601470822614196E-008_wp, -5.1191902385412367E-008_wp, & 
      & 4.4189831509263550E-005_wp, -6.5472071293320507E-008_wp, -5.1202021605699378E-008_wp, & 
      & 6.7417854867650525E-008_wp, -0.32715636773669671_wp,     -0.19343045344911730_wp,     & 
      &-1.9457835743300046E-009_wp, -2.1586482575562844E-002_wp, -2.4035908212783213E-002_wp, & 
      &-6.5472071293320507E-008_wp,  0.34874285031290009_wp,      0.21746636165654390_wp,     & 
      & 6.9189718349509533E-008_wp, -0.24150211636497051_wp,     -0.19851135822773272_wp,     & 
      &-1.7987696743810155E-008_wp,  2.4035754703519811E-002_wp,  1.0085713431492066E-002_wp, & 
      &-5.1202021605699378E-008_wp,  0.21746636165654390_wp,      0.18842564479480242_wp],    & 
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 2.01474475e-6_wp
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :)
   real(wp) :: final_err
   integer :: N, i

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N))
   call calc%odlrhessian(env, mol, chk, step, hessian, final_err)

   if (any(abs(hessian - hessian_ref) > thr1)) then
      call test_failed(error, "Final Hessians do not match")

      print *, "--- hessian ---"
      do i = 1, N
         print '(*(F21.14))', hessian(i, :) 
      end do

      print *, "--- Ref. hessian ---"
      do i = 1, N
         print '(*(F21.14))', hessian_ref(i, :) 
      end do
   end if

   if (abs(final_err - final_err_ref) > thr2) then
      call test_failed(error, "Final error does not match")
      print *, "--- Final error ---"
      print '(*(F21.14))', final_err
      print *, "--- Ref. final error ---"
      print '(*(F21.14))', final_err_ref
   end if
end subroutine test_o1numhess_gfn1

subroutine test_o1numhess_gfn2(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      &-7.8530493188148887E-004_wp, -2.0253444638258027E-011_wp,  -1.4017991789878334E-007_wp, &
      & 3.9265259955460915E-004_wp, -6.4571737815115928E-008_wp,   7.0089204184450145E-008_wp, &
      & 3.9265233232687208E-004_wp,  6.4590834777436623E-008_wp,   7.0091291955493636E-008_wp, &
      &-2.0253444638258027E-011_wp,  0.61315118431143367_wp,      -7.7689695353273194E-007_wp, &
      & 5.4589296151605569E-011_wp, -0.30657602838911235_wp,       0.22663338413584239_wp,     &
      &-3.4335851513348699E-011_wp, -0.30657515591650192_wp,      -0.22663260724061299_wp,     &
      &-1.4017991789878334E-007_wp, -7.7689695353273194E-007_wp,   0.40096663617601619_wp,     &
      & 3.8327586404134452E-008_wp,  0.18510543401408652_wp,      -0.20048354755593945_wp,     &
      & 1.0185233149464886E-007_wp, -0.18510465711519453_wp,      -0.20048308862513586_wp,     &
      & 3.9265259955460915E-004_wp,  5.4589296151605569E-011_wp,   3.8327586404134452E-008_wp, &
      &-5.0343577500715762E-004_wp,  1.7015555656545986E-011_wp,  -1.9169941125483937E-008_wp, &
      & 1.1078317545255174E-004_wp, -7.1677131952690313E-011_wp,  -1.9152730228802068E-008_wp, &
      &-6.4571737815115928E-008_wp, -0.30657602838911235_wp,       0.18510543401408652_wp,     &
      & 1.7015555656545986E-011_wp,  0.34071042144520453_wp,      -0.20586933199381116_wp,     &
      & 6.4554722259459369E-008_wp, -3.4134393059175372E-002_wp,   2.0763897984405841E-002_wp, &
      & 7.0089204184450145E-008_wp,  0.22663338413584239_wp,      -0.20048354755593945_wp,     &
      &-1.9169941125483937E-008_wp, -0.20586933199381116_wp,       0.18359573693614070_wp,     &
      &-5.0919263058966194E-008_wp, -2.0764052139551202E-002_wp,   1.6887810619957113E-002_wp, &
      & 3.9265233232687208E-004_wp, -3.4335851513348699E-011_wp,   1.0185233149464886E-007_wp, &
      & 1.1078317545255174E-004_wp,  6.4554722259459369E-008_wp,  -5.0919263058966194E-008_wp, &
      &-5.0343550777941942E-004_wp, -6.4519157645483930E-008_wp,  -5.0938561726691561E-008_wp, &
      & 6.4590834777436623E-008_wp, -0.30657515591650192_wp,      -0.18510465711519453_wp,     &
      &-7.1677131952690313E-011_wp, -3.4134393059175372E-002_wp,  -2.0764052139551202E-002_wp, &
      &-6.4519157645483930E-008_wp,  0.34070954897294109_wp,       0.20586870925178866_wp,     &
      & 7.0091291955493636E-008_wp, -0.22663260724061299_wp,      -0.20048308862513586_wp,     &
      &-1.9152730228802068E-008_wp,  2.0763897984405841E-002_wp,   1.6887810619957113E-002_wp, &
      &-5.0938561726691561E-008_wp,  0.20586870925178866_wp,       0.18359527801007947_wp],    &
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 2.10070777e-6_wp
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :)
   real(wp) :: final_err
   integer :: N, i

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N))
   call calc%odlrhessian(env, mol, chk, step, hessian, final_err)

   if (any(abs(hessian - hessian_ref) > thr1)) then
      call test_failed(error, "Final Hessians do not match")

      print *, "--- hessian ---"
      do i = 1, N
         print '(*(F21.14))', hessian(i, :) 
      end do

      print *, "--- Ref. hessian ---"
      do i = 1, N
         print '(*(F21.14))', hessian_ref(i, :) 
      end do
   end if

   if (abs(final_err - final_err_ref) > thr2) then
      call test_failed(error, "Final error does not match")
      print *, "--- Final error ---"
      print '(*(F21.14))', final_err
      print *, "--- Ref. final error ---"
      print '(*(F21.14))', final_err_ref
   end if
end subroutine test_o1numhess_gfn2

subroutine test_o1numhess_linear_h2o_gfn1(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,  0.00000000000000_wp,  0.00000000000000_wp, &
      & 0.00000000000000_wp,  0.00000000000000_wp, -1.81075448577205_wp, &
      & 0.00000000000000_wp,  0.00000000000000_wp,  1.81075448676713_wp],&
      & shape(xyz))
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :), aux(:), freq(:)
   real(wp) :: final_err
   integer :: N, i, lwork, info

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N))
   call calc%odlrhessian(env, mol, chk, step, hessian, final_err)
   allocate(freq(N))
   lwork  = 1 + 6*N + 2*N**2
   allocate(aux(lwork))
   call dsyev ('V', 'U', N, hessian, N, freq, aux, lwork, info)
   if (count(abs(freq) < 1.0e-4_wp) /= 3) then
      call test_failed(error, "Linear H2O should have exactly three ~0 freqs")
   end if

   if (count(freq < -1.0e-10_wp) /= 2) then
      call test_failed(error, "Linear H2O should have exactly two negative freqs")
   end if

   ! freqs around -0.3
   if (freq(1) > -0.3_wp .or. freq(2) > -0.3_wp) then
      call test_failed(error, "First two freqs should be negative")
   end if
end subroutine test_o1numhess_linear_h2o_gfn1

subroutine test_o1numhess_linear_h2o_gfn2(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,  0.00000000000000_wp,  0.00000000000000_wp, &
      & 0.00000000000000_wp,  0.00000000000000_wp, -1.81075448577205_wp, &
      & 0.00000000000000_wp,  0.00000000000000_wp,  1.81075448676713_wp],&
      & shape(xyz))
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :), aux(:), freq(:)
   real(wp) :: final_err
   integer :: N, i, lwork, info

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N))
   call calc%odlrhessian(env, mol, chk, step, hessian, final_err)
   allocate(freq(N))
   lwork  = 1 + 6*N + 2*N**2
   allocate(aux(lwork))
   call dsyev ('V', 'U', N, hessian, N, freq, aux, lwork, info)

   if (count(abs(freq) < 1.0e-4_wp) /= 3) then
      call test_failed(error, "Linear H2O should have exactly three ~0 freqs")
   end if

   if (count(freq < -1.0e-10_wp) /= 2) then
      call test_failed(error, "Linear H2O should have exactly two negative freqs")
   end if

   ! freqs around -0.3
   if (freq(1) > -0.3_wp .or. freq(2) > -0.3_wp) then
      call test_failed(error, "First two freqs should be negative")
   end if
end subroutine test_o1numhess_linear_h2o_gfn2

! NOTE: the following is not used anymore because get_gradient_derivs is now private
! but might be useful for future debugging

!!> Utility subroutine to setup neighbor list and displdir for o1numhess tests
!subroutine setup_o1numhess_test(N, displdir, ndispl_final, h0, eps, eps2)
!   integer, intent(in) :: N
!   real(wp), allocatable, intent(out) :: displdir(:, :)
!   integer, intent(out) :: ndispl_final
!   real(wp), intent(in) :: h0(:, :)
!   real(wp), intent(in) :: eps, eps2
!
!   type(adj_list), allocatable :: neighborlist(:)
!   real(wp), allocatable :: distmat(:, :)
!   integer, allocatable :: nbcounts(:)
!   integer :: i, j, max_nb, ndispl0
!
!   ! setup distmat
!   allocate(distmat(N, N))
!   do i = 1, N
!      do j = 1, N
!         distmat(i, j) = abs(i - j)
!      end do
!   end do
!
!   ! compute neighbor list
!   call get_neighbor_list(distmat, 1.0_wp, neighborlist)
!   allocate(nbcounts(N))
!   max_nb = 0
!   do i = 1, N
!      nbcounts(i) = size(neighborlist(i)%neighbors)
!      if (nbcounts(i) > max_nb) max_nb = nbcounts(i)
!   end do
!
!   ! populate displdir
!   ndispl0 = 0
!   allocate(displdir(N, N))
!   call gen_displdir(N, ndispl0, h0, max_nb, neighborlist, nbcounts, eps, eps2, displdir, ndispl_final)
!
!end subroutine setup_o1numhess_test
!
!!> Utility subroutine to calculate local and final hessian for o1numhess tests
!subroutine calculate_o1numhess_hessian(mol, env, calc, chk, displdir, ndispl_final, hessian, hessian_local, step)
!   type(TMolecule), intent(inout) :: mol
!   type(TEnvironment), intent(inout) :: env
!   type(TxTBCalculator), intent(inout) :: calc
!   type(TRestart), intent(inout) :: chk
!   real(wp), intent(in) :: displdir(:, :)
!   integer, intent(in) :: ndispl_final
!   real(wp), allocatable, intent(out) :: hessian(:, :), hessian_local(:, :)
!   real(wp), intent(in) :: step
!
!   real(wp), allocatable :: tmp_grad(:, :), g0(:), g(:, :), distmat(:, :)
!   real(wp) :: energy, sigma(3, 3), egap
!   type(scc_results) :: res
!   real(wp) :: final_err
!   integer :: i, j, N, ndispl0
!
!   N = 3 * mol%n
!   ndispl0 = 0
!
!   ! calculate gradient derivs
!   allocate(tmp_grad(3, mol%n))
!   call calc%singlepoint(env, mol, chk, -1, .false., energy, tmp_grad, sigma, egap, res)
!   g0 = reshape(tmp_grad,[N])
!   allocate(g(N, ndispl_final))
!   g = 0.0_wp
!
!   call calc%get_gradient_derivs(env, step, ndispl0, ndispl_final, displdir, mol, chk, g0, .false., g)
!
!   ! setup distmat (needed for gen_local_hessian)
!   allocate(distmat(N, N))
!   do i = 1, N
!      do j = 1, N
!         distmat(i, j) = abs(i - j)
!      end do
!   end do
!
!   ! calculate local hessian
!   allocate(hessian(N, N))
!   hessian = 0.0_wp
!   call gen_local_hessian(env, ndispl_final, distmat, displdir, g, 1.0_wp, hessian)
!   hessian_local = hessian
!   call lr_loop(env, ndispl_final, g, hessian, displdir, final_err)
!
!   ! print *, "final err", final_err
!
!end subroutine calculate_o1numhess_hessian
!
!subroutine test_o1numhess_gfn1(error)
!   type(error_type), allocatable, intent(out) :: error
!   integer, parameter :: nat = 3
!   real(wp),parameter :: thr = 1.0e-6_wp
!   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
!   real(wp), parameter :: xyz(3, nat) = reshape([&
!      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
!      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
!      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
!      & shape(xyz))
!   ! swart model hessian from o1numhess utils
!   real(wp), parameter :: h0(9, 9) = reshape([&
!      & 1.1092305614e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp,  4.4184751144e-01_wp, -2.7552743027e-10_wp, &
!      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.6310770462e-01_wp, &
!      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.6310770435e-01_wp, &
!      & 0.0000000000e+00_wp, -2.7552743027e-10_wp,  2.9417034746e-01_wp, &
!      & 0.0000000000e+00_wp,  1.4319582461e-01_wp, -1.4708517387e-01_wp, &
!      & 0.0000000000e+00_wp, -1.4319582434e-01_wp, -1.4708517359e-01_wp, &
!      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 1.6724709350e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.4319582461e-01_wp, &
!      & 0.0000000000e+00_wp,  2.3323207511e-01_wp, -1.5315176462e-01_wp, &
!      & 0.0000000000e+00_wp, -1.2308319237e-02_wp,  9.9559400089e-03_wp, &
!      & 0.0000000000e+00_wp,  1.6310770462e-01_wp, -1.4708517387e-01_wp, &
!      & 0.0000000000e+00_wp, -1.5315176462e-01_wp,  1.3375368579e-01_wp, &
!      & 0.0000000000e+00_wp, -9.9559400048e-03_wp,  1.3331488088e-02_wp, &
!      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 1.6724709315e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.4319582434e-01_wp, &
!      & 0.0000000000e+00_wp, -1.2308319237e-02_wp, -9.9559400048e-03_wp, &
!      & 0.0000000000e+00_wp,  2.3323207480e-01_wp,  1.5315176434e-01_wp, &
!      & 0.0000000000e+00_wp, -1.6310770435e-01_wp, -1.4708517359e-01_wp, &
!      & 0.0000000000e+00_wp,  9.9559400089e-03_wp,  1.3331488088e-02_wp, &
!      & 0.0000000000e+00_wp,  1.5315176434e-01_wp,  1.3375368550e-01_wp],&
!      & shape(h0))
!   real(wp), parameter :: hessian_local_ref(9, 9) = reshape([&
!      & 3.51981175325187e-04_wp, -3.76511734388845e-03_wp, 3.81178115691305e-04_wp,&
!      & -2.37598469686780e-04_wp, 2.52570715380943e-03_wp, -3.32541413989657e-04_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & -3.76511734388845e-03_wp, -6.72401397966386e-02_wp, -2.11729004550184e-01_wp,&
!      & 1.32115390027701e-02_wp, 8.06252321926636e-03_wp, 7.39146674965733e-03_wp,&
!      & 1.17788921786219e-04_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & 3.81178115691305e-04_wp, -2.11729004550184e-01_wp, 4.63279981243663e-01_wp,&
!      & 2.00861682675189e-02_wp, 1.75152612759093e-01_wp, -5.17419081415414e-03_wp,&
!      & 9.72230292022307e-04_wp, 8.49401740708410e-04_wp, 0.00000000000000e+00_wp,&
!      & -2.37598469686780e-04_wp, 1.32115390027701e-02_wp, 2.00861682675189e-02_wp,&
!      & 2.49132380283098e-03_wp, -1.15993237268750e-02_wp, 2.38421110728877e-02_wp,&
!      & -1.70620005390799e-14_wp, 5.64969948975101e-04_wp, -1.06083536905192e-04_wp,&
!      & 2.52570715380943e-03_wp, 8.06252321926636e-03_wp, 1.75152612759093e-01_wp,&
!      & -1.15993237268750e-02_wp, -5.36618040845405e-02_wp, -1.22992832629143e-01_wp,&
!      & -1.28817650172161e-02_wp, -1.41730858229258e-04_wp, 4.94454130056628e-03_wp,&
!      & -3.32541413989657e-04_wp, 7.39146674965733e-03_wp, -5.17419081415414e-03_wp,&
!      & 2.38421110728877e-02_wp, -1.22992832629143e-01_wp, 2.84196776383033e-01_wp,&
!      & 1.79323845243361e-02_wp, 1.49966320866984e-01_wp, 3.13392912791258e-03_wp,&
!      & 0.00000000000000e+00_wp, 1.17788921786219e-04_wp, 9.72230292022307e-04_wp,&
!      & -1.70620005390799e-14_wp, -1.28817650172161e-02_wp, 1.79323845243361e-02_wp,&
!      & 4.76557582555649e-05_wp, 8.95077909660287e-03_wp, 1.54823774005412e-02_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 8.49401740708410e-04_wp,&
!      & 5.64969948975101e-04_wp, -1.41730858229258e-04_wp, 1.49966320866984e-01_wp,&
!      & 8.95077909660287e-03_wp, 7.32108733232320e-02_wp, 8.32366313319442e-02_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & -1.06083536905192e-04_wp, 4.94454130056628e-03_wp, 3.13392912791258e-03_wp,&
!      & 1.54823774005412e-02_wp, 8.32366313319442e-02_wp, -9.34713030397501e-02_wp],&
!      & shape(hessian_local_ref))
!   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
!      & 3.58438579043206e-04_wp, -3.81039356197628e-03_wp, 4.03399313822165e-04_wp,&
!      & -2.82476978351356e-04_wp, 3.84811786546152e-03_wp, -7.63006730622276e-04_wp,&
!      & 2.31299048839153e-04_wp, 4.84032295451731e-04_wp, 1.32189016912669e-04_wp,&
!      & -3.81039356197628e-03_wp, -6.72927252587342e-02_wp, -2.11638325546983e-01_wp,&
!      & 1.33755408386769e-02_wp, 9.29003639664687e-03_wp, 1.13194001517243e-02_wp,&
!      & 1.68329797547482e-04_wp, 3.50224644185952e-04_wp, -1.00445342525072e-03_wp,&
!      & 4.03399313822165e-04_wp, -2.11638325546983e-01_wp, 4.63339799742256e-01_wp,&
!      & 2.00825620440430e-02_wp, 1.78712614527116e-01_wp, -6.07451375400960e-03_wp,&
!      & 1.49362842611287e-03_wp, 2.00250342612253e-03_wp, 5.14406818716446e-04_wp,&
!      & -2.82476978351356e-04_wp, 1.33755408386769e-02_wp, 2.00825620440430e-02_wp,&
!      & 2.47276791227794e-03_wp, -1.16796282245497e-02_wp, 2.43427547421269e-02_wp,&
!      & -1.85558905700834e-05_wp, 7.76482115758251e-04_wp, -2.82293395959265e-04_wp,&
!      & 3.84811786546152e-03_wp, 9.29003639664687e-03_wp, 1.78712614527116e-01_wp,&
!      & -1.16796282245497e-02_wp, -5.37238742371844e-02_wp, -1.23025151971481e-01_wp,&
!      & -1.32197048151783e-02_wp, -2.37180362396234e-04_wp, 7.63607499039816e-03_wp,&
!      & -7.63006730622276e-04_wp, 1.13194001517243e-02_wp, -6.07451375400960e-03_wp,&
!      & 2.43427547421269e-02_wp, -1.23025151971481e-01_wp, 2.84273693850240e-01_wp,&
!      & 1.79561859722289e-02_wp, 1.52928583030049e-01_wp, 3.54908290642412e-03_wp,&
!      & 2.31299048839153e-04_wp, 1.68329797547482e-04_wp, 1.49362842611287e-03_wp,&
!      & -1.85558905700834e-05_wp, -1.32197048151783e-02_wp, 1.79561859722289e-02_wp,&
!      & 2.90998677025632e-05_wp, 8.85720749101923e-03_wp, 1.57516020166363e-02_wp,&
!      & 4.84032295451731e-04_wp, 3.50224644185952e-04_wp, 2.00250342612253e-03_wp,&
!      & 7.76482115758251e-04_wp, -2.37180362396234e-04_wp, 1.52928583030049e-01_wp,&
!      & 8.85720749101923e-03_wp, 7.31255602924301e-02_wp, 8.32609561079180e-02_wp,&
!      & 1.32189016912669e-04_wp, -1.00445342525072e-03_wp, 5.14406818716446e-04_wp,&
!      & -2.82293395959265e-04_wp, 7.63607499039816e-03_wp, 3.54908290642412e-03_wp,&
!      & 1.57516020166363e-02_wp, 8.32609561079180e-02_wp, -9.33774159826147e-02_wp],&
!      & shape(hessian_ref))
!   real(wp), parameter :: eps = 1.0e-8_wp, eps2 = 1.0e-15_wp, step = 1.0e-6_wp
!
!   type(TMolecule) :: mol
!   type(TRestart) :: chk
!   type(TEnvironment) :: env
!   type(TxTBCalculator) :: calc
!   real(wp), allocatable :: displdir(:, :), hessian(:, :), hessian_local(:, :)
!   integer :: i, N, ndispl_final
!
!   call init(env)
!   call init(mol, sym, xyz)
!
!   call setup_o1numhess_test(3*mol%n, displdir, ndispl_final, h0, eps, eps2)
!
!   ! use GFN1
!   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-8_wp)
!   call newWavefunction(env, mol, calc, chk)
!
!   call calculate_o1numhess_hessian(mol, env, calc, chk, displdir, ndispl_final, hessian, hessian_local, step)
!
!   ! compare
!   N = 3 * nat
!   if (any(abs(hessian_local - hessian_local_ref) > thr)) then
!      call test_failed(error, "Local Hessians do not match")
!
!      print *, "--- hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_local(i, :) 
!      end do
!
!      print *, "--- Ref. hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_local_ref(i, :) 
!      end do
!   end if
!
!   if (any(abs(hessian - hessian_ref) > thr)) then
!      call test_failed(error, "Final Hessians do not match")
!
!      print *, "--- hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian(i, :) 
!      end do
!
!      print *, "--- Ref. hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_ref(i, :) 
!      end do
!   end if
!end subroutine test_o1numhess_gfn1
!
!subroutine test_o1numhess_gfn2(error)
!   type(error_type), allocatable, intent(out) :: error
!   integer, parameter :: nat = 3
!   real(wp),parameter :: thr = 1.0e-6_wp
!   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
!   real(wp), parameter :: xyz(3, nat) = reshape([&
!      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
!      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
!      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
!      & shape(xyz))
!   ! swart model hessian from o1numhess utils
!   real(wp), parameter :: h0(9, 9) = reshape([&
!      & 1.1092305614e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp,  4.4184751144e-01_wp, -2.7552743027e-10_wp, &
!      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.6310770462e-01_wp, &
!      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.6310770435e-01_wp, &
!      & 0.0000000000e+00_wp, -2.7552743027e-10_wp,  2.9417034746e-01_wp, &
!      & 0.0000000000e+00_wp,  1.4319582461e-01_wp, -1.4708517387e-01_wp, &
!      & 0.0000000000e+00_wp, -1.4319582434e-01_wp, -1.4708517359e-01_wp, &
!      &-5.5461528244e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 1.6724709350e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp, -2.2092375587e-01_wp,  1.4319582461e-01_wp, &
!      & 0.0000000000e+00_wp,  2.3323207511e-01_wp, -1.5315176462e-01_wp, &
!      & 0.0000000000e+00_wp, -1.2308319237e-02_wp,  9.9559400089e-03_wp, &
!      & 0.0000000000e+00_wp,  1.6310770462e-01_wp, -1.4708517387e-01_wp, &
!      & 0.0000000000e+00_wp, -1.5315176462e-01_wp,  1.3375368579e-01_wp, &
!      & 0.0000000000e+00_wp, -9.9559400048e-03_wp,  1.3331488088e-02_wp, &
!      &-5.5461527896e-04_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      &-1.1178556526e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 1.6724709315e-03_wp,  0.0000000000e+00_wp,  0.0000000000e+00_wp, &
!      & 0.0000000000e+00_wp, -2.2092375557e-01_wp, -1.4319582434e-01_wp, &
!      & 0.0000000000e+00_wp, -1.2308319237e-02_wp, -9.9559400048e-03_wp, &
!      & 0.0000000000e+00_wp,  2.3323207480e-01_wp,  1.5315176434e-01_wp, &
!      & 0.0000000000e+00_wp, -1.6310770435e-01_wp, -1.4708517359e-01_wp, &
!      & 0.0000000000e+00_wp,  9.9559400089e-03_wp,  1.3331488088e-02_wp, &
!      & 0.0000000000e+00_wp,  1.5315176434e-01_wp,  1.3375368550e-01_wp],&
!      & shape(h0))
!   real(wp), parameter :: hessian_local_ref(9, 9) = reshape([&
!      & 3.19677112844234e-04_wp, -3.32342131309818e-03_wp, -2.67790616248500e-05_wp,&
!      & -1.93324060431286e-04_wp, 2.34887629241175e-03_wp, -2.97382602124213e-04_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & -3.32342131309818e-03_wp, -6.30870143425176e-02_wp, -1.99496626493542e-01_wp,&
!      & 1.20902490299375e-02_wp, 7.47363090613361e-03_wp, 6.85160128840844e-03_wp,&
!      & 8.90280217865149e-05_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & -2.67790616248500e-05_wp, -1.99496626493542e-01_wp, 4.63216871841495e-01_wp,&
!      & 1.92444959775377e-02_wp, 1.62320366204895e-01_wp, -4.80311579353320e-03_wp,&
!      & 9.25923730998759e-04_wp, 6.23997695176231e-04_wp, 0.00000000000000e+00_wp,&
!      & -1.93324060431286e-04_wp, 1.20902490299375e-02_wp, 1.92444959775377e-02_wp,&
!      & 2.33820095994906e-03_wp, -1.06506356763412e-02_wp, 2.25289791347748e-02_wp,&
!      & -1.55391471185750e-14_wp, 5.73518635910525e-04_wp, -6.78104372076033e-05_wp,&
!      & 2.34887629241175e-03_wp, 7.47363090613361e-03_wp, 1.62320366204895e-01_wp,&
!      & -1.06506356763412e-02_wp, -5.41413149509843e-02_wp, -1.20812305434173e-01_wp,&
!      & -1.30893376495030e-02_wp, -1.41183939391319e-04_wp, 4.92546103224633e-03_wp,&
!      & -2.97382602124213e-04_wp, 6.85160128840844e-03_wp, -4.80311579353320e-03_wp,&
!      & 2.25289791347748e-02_wp, -1.20812305434173e-01_wp, 2.79913227822072e-01_wp,&
!      & 1.87438497579491e-02_wp, 1.49387622567730e-01_wp, 3.12183574957329e-03_wp,&
!      & 0.00000000000000e+00_wp, 8.90280217865149e-05_wp, 9.25923730998759e-04_wp,&
!      & -1.55391471185750e-14_wp, -1.30893376495030e-02_wp, 1.87438497579491e-02_wp,&
!      & -4.79319030380975e-06_wp, 8.88421644367883e-03_wp, 1.65261717538513e-02_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 6.23997695176231e-04_wp,&
!      & 5.73518635910525e-04_wp, -1.41183939391319e-04_wp, 1.49387622567730e-01_wp,&
!      & 8.88421644367883e-03_wp, 7.07044219384617e-02_wp, 8.37522763859865e-02_wp,&
!      & 0.00000000000000e+00_wp, 0.00000000000000e+00_wp, 0.00000000000000e+00_wp,&
!      & -6.78104372076033e-05_wp, 4.92546103224633e-03_wp, 3.12183574957329e-03_wp,&
!      & 1.65261717538513e-02_wp, 8.37522763859865e-02_wp, -1.01307468842741e-01_wp],&
!      & shape(hessian_local_ref))
!   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
!      & 3.21202548348861e-04_wp, -3.36143522890183e-03_wp, -1.61434670765408e-05_wp,&
!      & -2.33476898312511e-04_wp, 3.58311191231894e-03_wp, -6.80597145923440e-04_wp,&
!      & 2.25424943548716e-04_wp, 3.90896402245910e-04_wp, 8.85821968433675e-05_wp,&
!      & -3.36143522890183e-03_wp, -6.31369977795672e-02_wp, -1.99406199899265e-01_wp,&
!      & 1.22397867307581e-02_wp, 8.60841380834760e-03_wp, 1.04867482257833e-02_wp,&
!      & 1.10716609923727e-04_wp, 2.44005464184140e-04_wp, -7.27804274818003e-04_wp,&
!      & -1.61434670765408e-05_wp, -1.99406199899265e-01_wp, 4.63276927090022e-01_wp,&
!      & 1.92372026481329e-02_wp, 1.65623784950145e-01_wp, -5.64376804607281e-03_wp,&
!      & 1.41862921621685e-03_wp, 1.48865713775440e-03_wp, 4.03179144287053e-04_wp,&
!      & -2.33476898312511e-04_wp, 1.22397867307581e-02_wp, 1.92372026481329e-02_wp,&
!      & 2.32016668873703e-03_wp, -1.07232346626970e-02_wp, 2.30053559824699e-02_wp,&
!      & -1.80342712296196e-05_wp, 7.98663420130482e-04_wp, -1.95601212615882e-04_wp,&
!      & 3.58311191231894e-03_wp, 8.60841380834760e-03_wp, 1.65623784950145e-01_wp,&
!      & -1.07232346626970e-02_wp, -5.42028871409179e-02_wp, -1.20845949315161e-01_wp,&
!      & -1.34237233887938e-02_wp, -2.35976097335490e-04_wp, 7.60816561498387e-03_wp,&
!      & -6.80597145923440e-04_wp, 1.04867482257833e-02_wp, -5.64376804607281e-03_wp,&
!      & 2.30053559824699e-02_wp, -1.20845949315161e-01_wp, 2.79989829800501e-01_wp,&
!      & 1.87696470230527e-02_wp, 1.52336761147164e-01_wp, 3.53512397590933e-03_wp,&
!      & 2.25424943548716e-04_wp, 1.10716609923727e-04_wp, 1.41862921621685e-03_wp,&
!      & -1.80342712296196e-05_wp, -1.34237233887938e-02_wp, 1.87696470230527e-02_wp,&
!      & -2.28274615199358e-05_wp, 8.79966116458134e-03_wp, 1.68157017730541e-02_wp,&
!      & 3.90896402245910e-04_wp, 2.44005464184140e-04_wp, 1.48865713775440e-03_wp,&
!      & 7.98663420130482e-04_wp, -2.35976097335490e-04_wp, 1.52336761147164e-01_wp,&
!      & 8.79966116458134e-03_wp, 7.06197606988344e-02_wp, 8.37783234293568e-02_wp,&
!      & 8.85821968433675e-05_wp, -7.27804274818003e-04_wp, 4.03179144287053e-04_wp,&
!      & -1.95601212615882e-04_wp, 7.60816561498387e-03_wp, 3.53512397590933e-03_wp,&
!      & 1.68157017730541e-02_wp, 8.37783234293568e-02_wp, -1.01213366096153e-01_wp],&
!      & shape(hessian_ref))
!   real(wp), parameter :: eps = 1.0e-8_wp, eps2 = 1.0e-15_wp, step = 1.0e-6_wp
!
!   type(TMolecule) :: mol
!   type(TRestart) :: chk
!   type(TEnvironment) :: env
!   type(TxTBCalculator) :: calc
!   real(wp), allocatable :: displdir(:, :), hessian(:, :), hessian_local(:, :)
!   real(wp) :: final_err
!   integer :: i, N, ndispl_final
!
!   N = 3 * nat
!   call init(env)
!   call init(mol, sym, xyz)
!
!   call setup_o1numhess_test(3*mol%n, displdir, ndispl_final, h0, eps, eps2)
!
!   ! use GFN2
!   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-8_wp)
!   call newWavefunction(env, mol, calc, chk)
!
!   call calculate_o1numhess_hessian(mol, env, calc, chk, displdir, ndispl_final, hessian, hessian_local, step)
!
!   ! compare
!   if (any(abs(hessian_local - hessian_local_ref) > thr)) then
!      call test_failed(error, "Local Hessians do not match")
!
!      print *, "--- hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_local(i, :) 
!      end do
!
!      print *, "--- Ref. hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_local_ref(i, :) 
!      end do
!   end if
!
!   if (any(abs(hessian - hessian_ref) > thr)) then
!      call test_failed(error, "Final Hessians do not match")
!
!      print *, "--- hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian(i, :) 
!      end do
!
!      print *, "--- Ref. hessian ---"
!      do i = 1, N
!         print '(*(F21.14))', hessian_ref(i, :) 
!      end do
!   end if
!end subroutine test_o1numhess_gfn2

end module test_hessian
