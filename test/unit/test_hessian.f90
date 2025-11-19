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
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction
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
      new_unittest("gfn2", test_gfn2_hessian) &
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
! TODO: linear h2o test for numhess
! TODO: o1numhess for gfn1
! TODO: o1numhess for gfn2
! NOTES: 
! - can probably reuse the arrays from above for testing (need to know the error threshold for new hessian)
end module test_hessian
