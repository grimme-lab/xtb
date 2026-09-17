! This file is part of xtb.
!
! Copyright (C) 2021 Sebastian Ehlert
! Copyright (C) 2026 Leopold M. Seidler
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
   implicit none
   private

   public :: collect_hessian


contains

!> Collect all exported unit tests
subroutine collect_hessian(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("gfn1_hessian", test_gfn1_hessian), &
      new_unittest("gfn2_hessian", test_gfn2_hessian), &
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

subroutine test_o1numhess_gfn1(error)
   type(error_type), allocatable, intent(out) :: error
   integer, parameter :: nat = 3
   real(wp),parameter :: thr1 = 1.4901161193847656e-08_wp, thr2 = 1.4901161193847656e-08_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      & 9.642722735617289E-05_wp,  2.150578923697825E-12_wp,  1.071019141435998E-12_wp, &
      &-4.821347480851584E-05_wp, -4.062528712426575E-12_wp,  2.583698911064602E-13_wp, &
      &-4.821375254766774E-05_wp,  1.911949788728749E-12_wp, -1.329389032542458E-12_wp, &
      & 2.150578923697825E-12_wp,  6.543135697355346E-01_wp, -7.229475687954884E-07_wp, &
      &-1.642628093963891E-12_wp, -3.271571836987667E-01_wp,  2.415028515059434E-01_wp, &
      &-5.079508297339344E-13_wp, -3.271563860451477E-01_wp, -2.415021285496890E-01_wp, &
      & 1.071019141435998E-12_wp, -7.229475687954884E-07_wp,  3.970231512770018E-01_wp, &
      &-6.195795301597118E-13_wp,  1.934311885548390E-01_wp, -1.985118423833666E-01_wp, &
      &-4.514396112762865E-13_wp, -1.934304656060189E-01_wp, -1.985113088943410E-01_wp, &
      &-4.821347480851584E-05_wp, -1.642628093963891E-12_wp, -6.195795301597118E-13_wp, &
      & 4.418955528562922E-05_wp,  3.506229450266824E-12_wp, -2.965819713964842E-13_wp, &
      & 4.023919522892954E-06_wp, -1.863601356302934E-12_wp,  9.161615015561960E-13_wp, &
      &-4.062528712426575E-12_wp, -3.271571836987667E-01_wp,  1.934311885548390E-01_wp, &
      & 3.506229450266824E-12_wp,  3.487436663442405E-01_wp, -2.174669866056496E-01_wp, &
      & 5.562992621597500E-13_wp, -2.158648263991621E-02_wp,  2.403579804767703E-02_wp, &
      & 2.583698911064602E-13_wp,  2.415028515059434E-01_wp, -1.985118423833666E-01_wp, &
      &-2.965819713964842E-13_wp, -2.174669866056496E-01_wp,  1.884261289483365E-01_wp, &
      & 3.821208029002418E-14_wp, -2.403586490401287E-02_wp,  1.008571343858931E-02_wp, &
      &-4.821375254766774E-05_wp, -5.079508297339344E-13_wp, -4.514396112762865E-13_wp, &
      & 4.023919522892954E-06_wp,  5.562992621597500E-13_wp,  3.821208029002418E-14_wp, &
      & 4.418983302477912E-05_wp, -4.834843242581555E-14_wp,  4.132275309862622E-13_wp, &
      & 1.911949788728749E-12_wp, -3.271563860451477E-01_wp, -1.934304656060189E-01_wp, &
      &-1.863601356302934E-12_wp, -2.158648263991621E-02_wp, -2.403586490401287E-02_wp, &
      &-4.834843242581555E-14_wp,  3.487428686878863E-01_wp,  2.174663305044797E-01_wp, &
      &-1.329389032542458E-12_wp, -2.415021285496890E-01_wp, -1.985113088943410E-01_wp, &
      & 9.161615015561960E-13_wp,  2.403579804767703E-02_wp,  1.008571343858931E-02_wp, &
      & 4.132275309862622E-13_wp,  2.174663305044797E-01_wp,  1.884255954528983E-01_wp],&
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 1.53703515069416350E-06_wp
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :), dipgrad_dummy(:, :)
   real(wp) :: final_err
   integer :: N, i

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N))
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)

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
   real(wp),parameter :: thr1 = 1.4901161193847656e-08_wp, thr2 = 1.4901161193847656e-08_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp,    0.00000000034546_wp,    0.18900383618455_wp, &
      & 0.00000000000000_wp,    1.45674735348811_wp,   -0.88650486059828_wp, &
      &-0.00000000000000_wp,   -1.45674735383357_wp,   -0.88650486086986_wp],&
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      &-7.853049280457453E-04_wp,  1.955977244785751E-12_wp,  1.247080264997499E-12_wp, &
      & 3.926525943422639E-04_wp, -1.210600324495050E-12_wp,  9.850364921710785E-14_wp, &
      & 3.926523337034877E-04_wp, -7.453769202907016E-13_wp, -1.345583914214606E-12_wp, &
      & 1.955977244785751E-12_wp,  6.131511844074383E-01_wp, -7.530713906089167E-07_wp, &
      &-4.508671918596937E-13_wp, -3.065760101726123E-01_wp,  2.266333722977527E-01_wp, &
      &-1.505110052926058E-12_wp, -3.065751742284646E-01_wp, -2.266326192307660E-01_wp, &
      & 1.247080264997499E-12_wp, -7.530713906089167E-07_wp,  4.009666360184508E-01_wp, &
      &-8.063227290784491E-13_wp,  1.851054221071441E-01_wp, -2.004835959201852E-01_wp, &
      &-4.407575359190495E-13_wp, -1.851046690338150E-01_wp, -2.004830401052821E-01_wp, &
      & 3.926525943422639E-04_wp, -4.508671918596937E-13_wp, -8.063227290784491E-13_wp, &
      &-5.034357727917302E-04_wp,  9.315986645004943E-13_wp,  2.367249284567081E-13_wp, &
      & 1.107831784494622E-04_wp, -4.807314726408010E-13_wp,  5.695978006217410E-13_wp, &
      &-1.210600324495050E-12_wp, -3.065760101726123E-01_wp,  1.851054221071441E-01_wp, &
      & 9.315986645004943E-13_wp,  3.407104032218639E-01_wp, -2.058693632980307E-01_wp, &
      & 2.790016599945554E-13_wp, -3.413439305260558E-02_wp,  2.076394119763031E-02_wp, &
      & 9.850364921710785E-14_wp,  2.266333722977527E-01_wp, -2.004835959201852E-01_wp, &
      & 2.367249284567081E-13_wp, -2.058693632980307E-01_wp,  1.835957853733738E-01_wp, &
      &-3.352285776738158E-13_wp, -2.076400899834286E-02_wp,  1.688781054869399E-02_wp, &
      & 3.926523337034877E-04_wp, -1.505110052926058E-12_wp, -4.407575359190495E-13_wp, &
      & 1.107831784494622E-04_wp,  2.790016599945554E-13_wp, -3.352285776738158E-13_wp, &
      &-5.034355121529521E-04_wp,  1.226108392931502E-12_wp,  7.759861135928650E-13_wp, &
      &-7.453769202907016E-13_wp, -3.065751742284646E-01_wp, -1.851046690338150E-01_wp, &
      &-4.807314726408010E-13_wp, -3.413439305260558E-02_wp, -2.076400899834286E-02_wp, &
      & 1.226108392931502E-12_wp,  3.407095672780631E-01_wp,  2.058686780298182E-01_wp, &
      &-1.345583914214606E-12_wp, -2.266326192307660E-01_wp, -2.004830401052821E-01_wp, &
      & 5.695978006217410E-13_wp,  2.076394119763031E-02_wp,  1.688781054869399E-02_wp, &
      & 7.759861135928650E-13_wp,  2.058686780298182E-01_wp,  1.835952295617220E-01_wp],&
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 1.27700129644838372E-06_wp
   real(wp), parameter :: step = 1.0e-6_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(TxTBCalculator) :: calc
   real(wp), allocatable :: hessian(:, :), dipgrad_dummy(:, :)
   real(wp) :: final_err
   integer :: N, i

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N))
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)

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
   real(wp), allocatable :: hessian(:, :), aux(:), freq(:), dipgrad_dummy(:, :)
   real(wp) :: final_err
   integer :: N, i, lwork, info

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N))
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)
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

   ! freqs around -0.16
   if (freq(1) > -0.16_wp .or. freq(2) > -0.16_wp) then
      call test_failed(error, "First two freqs should be negative")
      print *, freq
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
   real(wp), allocatable :: hessian(:, :), aux(:), freq(:), dipgrad_dummy(:, :)
   real(wp) :: final_err
   integer :: N, i, lwork, info

   N = 3 * nat
   call init(env)
   call init(mol, sym, xyz)

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-8_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N))
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)
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

end module test_hessian
