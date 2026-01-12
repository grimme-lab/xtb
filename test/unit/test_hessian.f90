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
   use xtb_o1numhess, only : gen_displdir, get_neighbor_list, adj_list
   implicit none
   private

   public :: collect_hessian


contains

!> Collect all exported unit tests
subroutine collect_hessian(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1", test_gfn1_hessian), &
      new_unittest("gfn2", test_gfn2_hessian), &
      new_unittest("gfn1_o1numhess", test_gfn1_o1numhess), &
      new_unittest("gfn2_o1numhess", test_gfn2_o1numhess), &
      new_unittest("o1numhess_nblist", test_nblist_o1numhess) &
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

! TODO: linear h2o test
! NOTE: need to get o1numhess from python implementation for reference
! TODO: o1numhess for gfn1
subroutine test_gfn1_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,   -0.00000000077760_wp,    0.18829790750029_wp, &
      & 0.00000000000000_wp,    1.45987612440076_wp,   -0.88615189669760_wp, &
      &-0.00000000000000_wp,   -1.45987612362316_wp,   -0.88615189608629_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      &  9.92358267447000E-03_wp, -1.59939748593000E-02_wp, -1.87048285136700E-02_wp, &
      & -2.59359330698000E-03_wp,  3.64688203467000E-03_wp, -5.39801379054000E-03_wp, &
      &  6.25238023315000E-03_wp, -3.27295611700000E-03_wp,  1.42168933961900E-02_wp, &
      & -1.59939748593000E-02_wp,  3.42218900459050E-01_wp, -3.08980536676550E-01_wp, &
      & -1.04138037574800E-01_wp, -1.25477379828200E-02_wp,  5.59504259788000E-03_wp, &
      & -7.47376889756000E-03_wp, -9.06559557738000E-03_wp, -1.59622384225000E-03_wp, &
      & -1.87048285136700E-02_wp, -3.08980536676550E-01_wp,  4.02905416993150E-01_wp, &
      & -3.03167484390390E-01_wp,  1.74873176732200E-01_wp, -1.63205233915600E-02_wp, &
      & -5.61303887266000E-03_wp, -1.87754058917000E-02_wp,  2.50371617394000E-03_wp, &
      & -2.59359330698000E-03_wp, -1.04138037574800E-01_wp, -3.03167484390390E-01_wp, &
      &  2.38883417311780E-01_wp,  3.38812242967830E-01_wp, -2.43831143274860E-01_wp, &
      &  5.19134281304400E-02_wp, -2.37579310095000E-03_wp,  1.85778177995000E-02_wp, &
      &  3.64688203467000E-03_wp, -1.25477379828200E-02_wp,  1.74873176732200E-01_wp, &
      &  3.38812242967830E-01_wp,  3.17724391774870E-01_wp, -6.05448741990820E-01_wp, &
      & -1.23828892769770E-01_wp, -7.22740379575100E-02_wp, -7.24042304863000E-03_wp, &
      & -5.39801379054000E-03_wp,  5.59504259788000E-03_wp, -1.63205233915600E-02_wp, &
      & -2.43831143274860E-01_wp, -6.05448741990820E-01_wp,  8.93127897294960E-01_wp, &
      &  2.36343834116960E-01_wp, -3.47959068635070E-01_wp,  6.11533252314200E-02_wp, &
      &  6.25238023315000E-03_wp, -7.47376889756000E-03_wp, -5.61303887266000E-03_wp, &
      &  5.19134281304400E-02_wp, -1.23828892769770E-01_wp,  2.36343834116960E-01_wp, &
      &  1.76329986096260E-01_wp, -2.50214065151940E-01_wp, -5.79648528100500E-02_wp, &
      & -3.27295611700000E-03_wp, -9.06559557738000E-03_wp, -1.87754058917000E-02_wp, &
      & -2.37579310095000E-03_wp, -7.22740379575100E-02_wp, -3.47959068635070E-01_wp, &
      & -2.50214065151940E-01_wp,  5.82460760279690E-01_wp,  1.54646889977970E-01_wp, &
      &  1.42168933961900E-02_wp, -1.59622384225000E-03_wp,  2.50371617394000E-03_wp, &
      &  1.85778177995000E-02_wp, -7.24042304863000E-03_wp,  6.11533252314200E-02_wp, &
      & -5.79648528100500E-02_wp,  1.54646889977970E-01_wp, -6.48551074732600E-02_wp],&
      &  shape(hessian_ref))

   ! real(wp), parameter :: step = 1.0e-6_wp
   !
   ! type(TMolecule) :: mol
   ! type(TRestart) :: chk
   ! type(TEnvironment) :: env
   ! type(scc_results) :: res
   ! type(TxTBCalculator) :: calc
   !
   ! integer :: i,j
   ! real(wp) :: energy, sigma(3, 3)
   ! real(wp) :: hl_gap
   ! real(wp),allocatable :: gradient(:,:), hessian(:,:), g(:, :), displdir(:, :), dipgrad(:, :)
   ! integer, allocatable :: list(:)
   !
   ! call init(env)
   ! call init(mol, sym, xyz)
   !
   ! allocate(gradient(3,mol%n), dipgrad(3,3*mol%n), hessian(3*mol%n,3*mol%n))
   ! energy = 0.0_wp
   ! gradient = 0.0_wp
   !
   ! call newXTBCalculator(env, mol, calc, method=1)
   ! call newWavefunction(env, mol, calc, chk)
   !
   ! call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
   !    & hl_gap, res)
   !
   ! dipgrad = 0.0_wp
   ! hessian = 0.0_wp
   ! list = [(i, i = 1, mol%n)]
   ! allocate(displdir(3*mol%n, 3*mol%n), g(3*mol%n, 3*mol%n))
   ! call calc%odlrhessian(env, mol, chk, list, step, displdir, g, hessian)
   !
   ! do i = 1, size(hessian_ref, 2)
   !    do j = 1, size(hessian_ref, 1)
   !       call check(error, hessian(j, i), hessian_ref(j, i), thr=thr)
   !    end do
   ! end do


end subroutine test_gfn1_o1numhess

! TODO: o1numhess for gfn2
subroutine test_gfn2_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,   -0.00000000077760_wp,    0.18829790750029_wp, &
      & 0.00000000000000_wp,    1.45987612440076_wp,   -0.88615189669760_wp, &
      &-0.00000000000000_wp,   -1.45987612362316_wp,   -0.88615189608629_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      &  1.62933952901000E-02_wp, -1.80625576050600E-02_wp, -2.47987822331700E-02_wp, &
      & -2.17470476835000E-03_wp,  3.08355509349000E-03_wp, -5.77375836794000E-03_wp, &
      &  6.25362285030000E-03_wp, -3.11443648298000E-03_wp,  1.35631483207000E-02_wp, &
      & -1.80625576050600E-02_wp,  3.19925401982700E-01_wp, -2.91364996984180E-01_wp, &
      & -9.36910289643000E-02_wp, -1.03103620156000E-02_wp,  5.09713565000000E-03_wp, &
      & -5.81344232255000E-03_wp, -8.63433650894000E-03_wp, -1.05664939542000E-03_wp, &
      & -2.47987822331700E-02_wp, -2.91364996984180E-01_wp,  4.14873085252300E-01_wp, &
      & -2.90696693193510E-01_wp,  1.45179613487700E-01_wp, -2.01830124826900E-02_wp, &
      & -5.31400396812000E-03_wp, -1.83529708442400E-02_wp,  1.38466608361000E-03_wp, &
      & -2.17470476835000E-03_wp, -9.36910289643000E-02_wp, -2.90696693193510E-01_wp, &
      &  2.28144675991550E-01_wp,  3.29092219214610E-01_wp, -2.44282490185450E-01_wp, &
      &  5.08526187694800E-02_wp, -2.31969703816000E-03_wp,  1.79600150151300E-02_wp, &
      &  3.08355509349000E-03_wp, -1.03103620156000E-02_wp,  1.45179613487700E-01_wp, &
      &  3.29092219214610E-01_wp,  3.34219270726620E-01_wp, -5.98300437953140E-01_wp, &
      & -1.13867596488980E-01_wp, -7.01285294822000E-02_wp, -7.14316983316000E-03_wp, &
      & -5.77375836794000E-03_wp,  5.09713565000000E-03_wp, -2.01830124826900E-02_wp, &
      & -2.44282490185450E-01_wp, -5.98300437953140E-01_wp,  8.82616691941500E-01_wp, &
      &  2.42819312635370E-01_wp, -3.42472382781960E-01_wp,  6.15935432766700E-02_wp, &
      &  6.25362285030000E-03_wp, -5.81344232255000E-03_wp, -5.31400396812000E-03_wp, &
      &  5.08526187694800E-02_wp, -1.13867596488980E-01_wp,  2.42819312635370E-01_wp, &
      &  1.80449827323740E-01_wp, -2.51959773923450E-01_wp, -7.12477083861000E-02_wp, &
      & -3.11443648298000E-03_wp, -8.63433650894000E-03_wp, -1.83529708442400E-02_wp, &
      & -2.31969703816000E-03_wp, -7.01285294822000E-02_wp, -3.42472382781960E-01_wp, &
      & -2.51959773923450E-01_wp,  5.76128840973710E-01_wp,  1.50731990717000E-01_wp, &
      &  1.35631483207000E-02_wp, -1.05664939542000E-03_wp,  1.38466608361000E-03_wp, &
      &  1.79600150151300E-02_wp, -7.14316983316000E-03_wp,  6.15935432766700E-02_wp, &
      & -7.12477083861000E-02_wp,  1.50731990717000E-01_wp, -5.42020080753600E-02_wp],&
      &  shape(hessian_ref))
end subroutine test_gfn2_o1numhess


subroutine test_nblist_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp), parameter :: thr = 1.0e-4_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   integer, parameter :: at(nat) = [8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,   -0.00000000077760_wp,    0.18829790750029_wp, &
      & 0.00000000000000_wp,    1.45987612440076_wp,   -0.88615189669760_wp, &
      &-0.00000000000000_wp,   -1.45987612362316_wp,   -0.88615189608629_wp],&
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
   real(wp), parameter :: distmat_ref(9, 9) = reshape([ &
      & -6.61404096000_wp,   -6.61404096000_wp,   -6.61404096000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -6.61404096000_wp,   -6.61404096000_wp,   -6.61404096000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -6.61404096000_wp,   -6.61404096000_wp,   -6.61404096000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp, &
      & -4.22125042000_wp,   -4.22125042000_wp,   -4.22125042000_wp, &
      & -2.53399696000_wp,   -2.53399696000_wp,   -2.53399696000_wp, &
      & -5.45374920000_wp,   -5.45374920000_wp,   -5.45374920000_wp],&
      & shape(distmat_ref))
   real(wp), parameter :: dmax = 1.0_wp, eps = 1.0e-8_wp, eps2 = 1.0e-15_wp

   type(adj_list), allocatable :: nblist_ref(:), nblist(:)
   real(wp), allocatable :: distmat(:, :)
   real(wp) :: dist
   integer :: N, i, j
   logical :: is_same

   N = 3 * nat

   allocate(nblist_ref(N))
   do i = 1, N
      allocate(nblist_ref(i)%neighbors(9))
      nblist_ref(i)%neighbors = [(j, j=1, 9)]
   end do

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

subroutine test_gendispldir_o1numhess(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,   -0.00000000077760_wp,    0.18829790750029_wp, &
      & 0.00000000000000_wp,    1.45987612440076_wp,   -0.88615189669760_wp, &
      &-0.00000000000000_wp,   -1.45987612362316_wp,   -0.88615189608629_wp],&
      & shape(xyz))

   ! real(wp), allocatable :: displdir0(:, :), displdir(:, :), h0(:, :)
   ! integer :: N, i
   !
   ! N = 3 * nat
   !
   ! allocate(displdir0(N, N), displdir(N, N), h0(N, N))
   ! h0 = 0.0_wp
   ! do i = 1, N
   !    h0(i, i) = 1.0_wp
   ! end do

   ! call gen_displdir(n, 0, h0, displdir0, max_nb, nblist, nbcounts, eps, eps2, displdir, ndispl_final)
end subroutine test_gendispldir_o1numhess

end module test_hessian
