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

!> Hessian, compliance and redundant-internal behaviour locks
!>
!> Numerical and O1NumHess (ODLR) Hessians of water are pinned against reference
!> tables and against their eigenvalue spectra, the GFN1 compliance constants of
!> water are pinned as a reference matrix, and the analytic B-matrix of a
!> non-collinear redundant coordinate set is validated against central finite
!> differences.
module test_hessian
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_mctc_convert, only : autoaa
   use xtb_mctc_constants, only : pi
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_param
   use xtb_type_data
   use xtb_type_environment

   use xtb_xtb_calculator, only : TxTBCalculator
   use xtb_main_setup, only : newXTBCalculator, newWavefunction
   use xtb_compliance, only : compliance_driver, compute_compliance
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_mctc_param, only : covalent_radius_d3
   use xtb_internals_graph, only : graph_type, init
   use xtb_internals_redundant, only : redundant_type, init
   use xtb_internals_type, only : coord_bond, coord_angle, coord_dihedral, coord_linbend
   use xtb_bmatrix, only : get_bmatrix
   use xtb_o1numhess, only : lr_loop
   use xtb_freq_project, only : projectHessian
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
      new_unittest("linear_h2o_gfn2_o1numhess", test_o1numhess_linear_h2o_gfn2), &
      new_unittest("rectangular_lr_loop", test_rectangular_lr_loop), &
      new_unittest("compliance", test_compliance), &
      new_unittest("compliance_water", test_compliance_water), &
      new_unittest("covalent_neighbour_list", test_covalent_neighbour_list), &
      new_unittest("redundant_bmatrix_fd", test_redundant_bmatrix_fd), &
      new_unittest("linear_chain_compliance", test_linear_chain_compliance), &
      new_unittest("near_linear_bmatrix_fd", test_near_linear_bmatrix_fd), &
      new_unittest("compliance_driver_redundant", test_compliance_driver_redundant), &
      new_unittest("linear_projection", test_linear_projection), &
      new_unittest("fixed_projection", test_fixed_projection), &
      new_unittest("fixed_linear_projection", test_fixed_linear_projection) &
      ]

end subroutine collect_hessian

!> A rectangular displacement set must produce a symmetric correction without
!> overrunning the GEMM workspace; final_err must describe the returned matrix.
subroutine test_rectangular_lr_loop(error)
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: scale = 1.0e-4_wp
   real(wp), parameter :: displdir(2, 1) = reshape([1.0_wp, 0.0_wp], [2, 1])
   type(TEnvironment) :: env
   real(wp) :: gradient(2, 1), hessian(2, 2), final_err, residual

   call init(env)
   gradient(:, 1) = [scale, 2.0_wp*scale]
   hessian = 0.0_wp
   call lr_loop(env, 1, gradient, hessian, displdir, final_err)

   call check(error, hessian(1, 1)/scale, 1.0_wp, thr=1.0e-4_wp)
   if (allocated(error)) return
   call check(error, hessian(2, 1)/scale, 2.0_wp, thr=1.0e-4_wp)
   if (allocated(error)) return
   call check(error, hessian(1, 2), hessian(2, 1), thr=scale*1.0e-12_wp)
   if (allocated(error)) return
   call check(error, hessian(2, 2), 0.0_wp, thr=scale*1.0e-12_wp)
   if (allocated(error)) return

   residual = norm2(gradient - matmul(hessian, displdir))/norm2(gradient)
   call check(error, final_err, residual, thr=1.0e-14_wp)
end subroutine test_rectangular_lr_loop

!> Pins the GFN1 numerical Hessian and dipole gradient of water against the
!> reference tables below.
subroutine test_gfn1_hessian(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000034546_wp, 0.18900383618455_wp, &
      & 0.00000000000000_wp, 1.45674735348811_wp, -0.88650486059828_wp, &
      &-0.00000000000000_wp, -1.45674735383357_wp, -0.88650486086986_wp], &
      & shape(xyz))
   real(wp), parameter :: dipgrad_ref(3, 3*nat) = reshape([ &
      & -1.013449328123924E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  3.527708416572683E-11_wp, -0.580988540908312E+00_wp, -4.440892098500626E-10_wp, &
      &  1.640601522247707E-11_wp, -1.577916091825909E-09_wp, -0.712523062581760E+00_wp, &
      &  0.506725097155187E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  5.168504173660994E-11_wp, 0.290494797892338E+00_wp, 7.508966683134055E-02_wp, &
      & -3.527728552279606E-11_wp, 0.159642686552089E+00_wp, 0.356260853040524E+00_wp, &
      &  0.506725096971555E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  0.000000000000000E+00_wp, 0.290494797707754E+00_wp, -7.508966686984309E-02_wp, &
      &  0.000000000000000E+00_wp, -0.159642686648137E+00_wp, 0.356260854150747E+00_wp], &
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      &  10.071650499390740E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -5.035811363445364E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -5.035839135945918E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -5.183734810089242E-12_wp, 0.654325256090513E+00_wp, -8.411850317864566E-10_wp, &
      &  4.912366132678069E-12_wp, -0.327162628574232E+00_wp, 0.241514196911523E+00_wp, &
      &  2.713686774111726E-13_wp, -0.327162627516280E+00_wp, -0.241514195869865E+00_wp, &
      &  1.269389491007523E-11_wp, -1.003642487147619E-09_wp, 0.397015891111508E+00_wp, &
      & -1.357110357831859E-11_wp, 0.193433637903868E+00_wp, -0.198507945853971E+00_wp, &
      &  8.772086682433594E-13_wp, -0.193433636862216E+00_wp, -0.198507945257537E+00_wp, &
      & -5.009881575075736E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  4.607945338911357E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  4.023918763456385E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  7.510234345234422E-12_wp, -0.327162715268037E+00_wp, 0.193442787219949E+00_wp, &
      & -8.658575901216439E-12_wp, 0.348749304550916E+00_wp, -0.217478532934599E+00_wp, &
      &  1.148341555982017E-12_wp, -2.158658928287718E-02_wp, 2.403578845504295E-02_wp, &
      &  5.183777743453684E-12_wp, 0.241504825916428E+00_wp, -0.198508418786384E+00_wp, &
      & -4.912396810694304E-12_wp, -0.217469144437701E+00_wp, 0.188422842621929E+00_wp, &
      & -2.713809327593799E-13_wp, -2.403568147872955E-02_wp, 1.008567370418267E-02_wp, &
      & -5.009909348870555E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  4.023966094020255E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  4.607973112642175E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  0.000000000000000E+00_wp, -0.327162714210147E+00_wp, -0.193442786178324E+00_wp, &
      &  0.000000000000000E+00_wp, -2.158658928289722E-02_wp, -2.403578852147744E-02_wp, &
      &  0.000000000000000E+00_wp, 0.348749303493043E+00_wp, 0.217478531886452E+00_wp, &
      &  0.000000000000000E+00_wp, -0.241504824874788E+00_wp, -0.198508418190000E+00_wp, &
      &  0.000000000000000E+00_wp, 2.403568148528039E-02_wp, 1.008567395243448E-02_wp, &
      &  0.000000000000000E+00_wp, 0.217469143389508E+00_wp, 0.188422842025550E+00_wp], &
      &  shape(hessian_ref))
   real(wp), parameter :: step = 5.0e-3_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i, j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp), allocatable :: gradient(:, :), dipgrad(:, :), hessian(:, :)
   integer, allocatable :: list(:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3, mol%n), dipgrad(3, 3*mol%n), hessian(3*mol%n, 3*mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-4_wp)
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
         if (allocated(error)) return
      end do
   end do

   do i = 1, size(hessian_ref, 2)
      do j = 1, size(hessian_ref, 1)
         call check(error, hessian(j, i), hessian_ref(j, i), thr=thr)
         if (allocated(error)) return
      end do
   end do

end subroutine test_gfn1_hessian

!> Pins the GFN2 numerical Hessian and dipole gradient of water against the
!> reference tables below.
subroutine test_gfn2_hessian(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr = 1.0e-7_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, -0.00000000077760_wp, 0.18829790750029_wp, &
      & 0.00000000000000_wp, 1.45987612440076_wp, -0.88615189669760_wp, &
      &-0.00000000000000_wp, -1.45987612362316_wp, -0.88615189608629_wp], &
      & shape(xyz))
   real(wp), parameter :: dipgrad_ref(3, 3*nat) = reshape([ &
      & -0.811304874337466E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -3.379448427312665E-11_wp, -0.407381159139852E+00_wp, -2.220446049250313E-09_wp, &
      & -1.688661825390276E-11_wp, 1.035790186933606E-09_wp, -0.470822851681985E+00_wp, &
      &  0.405653028129658E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  7.146633102892883E-11_wp, 0.203691696290942E+00_wp, 10.003961418175722E-02_wp, &
      & -1.829769461273858E-14_wp, 0.148641975785757E+00_wp, 0.235410623840382E+00_wp, &
      &  0.405653028514606E+00_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -4.559890805138046E-14_wp, 0.203691696462815E+00_wp, -10.003961419398077E-02_wp, &
      &  3.356458019834560E-11_wp, -0.148641975709096E+00_wp, 0.235410623986088E+00_wp], &
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      & -10.350003623222118E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  51.749727501138594E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  51.750308731077851E-06_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -1.978173339857878E-11_wp, 0.609790159587238E+00_wp, 2.697691516787698E-09_wp, &
      &  2.481358892927353E-11_wp, -0.304895078542458E+00_wp, 0.224446027814586E+00_wp, &
      & -5.031855530694753E-12_wp, -0.304895081044780E+00_wp, -0.224446030224059E+00_wp, &
      & -1.024802643153815E-11_wp, 2.163320452807998E-09_wp, 0.397834712878787E+00_wp, &
      &  1.255439924419533E-11_wp, 0.182997022810964E+00_wp, -0.198917355722765E+00_wp, &
      & -2.306372812657182E-12_wp, -0.182997025220330E+00_wp, -0.198917357156021E+00_wp, &
      & 5.210681441559375E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  0.206210280389948E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  -5.416516874286603E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  3.321776245993920E-11_wp, -0.304895066033135E+00_wp, 0.183006624733882E+00_wp, &
      & -3.597681672480956E-11_wp, 0.338691636342253E+00_wp, -0.203726326929642E+00_wp, &
      &  2.759054264870358E-12_wp, -3.379657030911816E-02_wp, 2.071975479991245E-02_wp, &
      & -2.747668766007676E-13_wp, 0.224436311742210E+00_wp, -0.198917773260492E+00_wp, &
      &  2.050024431283972E-13_wp, -0.203716602155418E+00_wp, 0.182074825008421E+00_wp, &
      &  6.976443347237028E-14_wp, -2.071975470580370E-02_wp, 1.684294825206964E-02_wp, &
      & 5.210739564922810E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  -5.416514608773755E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      &  0.206152156648532E-05_wp, 0.000000000000000E+00_wp, 0.000000000000000E+00_wp, &
      & -2.703865614441994E-13_wp, -0.304895068535454E+00_wp, -0.183006627143470E+00_wp, &
      &  2.719647608323040E-13_wp, -3.379657030920542E-02_wp, -2.071975532206421E-02_wp, &
      & -1.578199388104547E-15_wp, 0.338691638844660E+00_wp, 0.203726329337441E+00_wp, &
      &  2.029210950346453E-11_wp, -0.224436314151628E+00_wp, -0.198917774693790E+00_wp, &
      & -2.499966571670633E-11_wp, 2.071975470580370E-02_wp, 1.684294825207260E-02_wp, &
      &  4.707556213241799E-12_wp, 0.203716604563118E+00_wp, 0.182074826441720E+00_wp], &
      &  shape(hessian_ref))

   real(wp), parameter :: step = 5.0e-3_wp

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   integer :: i, j
   real(wp) :: energy, sigma(3, 3)
   real(wp) :: hl_gap
   real(wp), allocatable :: gradient(:, :), dipgrad(:, :), hessian(:, :)
   integer, allocatable :: list(:)

   call init(env)
   call init(mol, sym, xyz)

   allocate(gradient(3, mol%n), dipgrad(3, 3*mol%n), hessian(3*mol%n, 3*mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-4_wp)
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
         if (allocated(error)) return
      end do
   end do

   do i = 1, size(hessian_ref, 2)
      do j = 1, size(hessian_ref, 1)
         call check(error, hessian(j, i), hessian_ref(j, i), thr=thr)
         if (allocated(error)) return
      end do
   end do

end subroutine test_gfn2_hessian

!> Pins the GFN1 O1NumHess (ODLR) Hessian of water and its final error
!> estimate against the reference values below.
subroutine test_o1numhess_gfn1(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr1 = 1.4901161193847656e-08_wp, thr2 = 1.4901161193847656e-08_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000034546_wp, 0.18900383618455_wp, &
      & 0.00000000000000_wp, 1.45674735348811_wp, -0.88650486059828_wp, &
      &-0.00000000000000_wp, -1.45674735383357_wp, -0.88650486086986_wp], &
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      & 9.642722910000000E-05_wp, -2.920000000000000E-12_wp, -5.610000000000000E-12_wp, &
      &-4.821347568000000E-05_wp,  2.240000000000000E-12_wp,  1.730000000000000E-12_wp, &
      &-4.821375342000000E-05_wp,  6.900000000000000E-13_wp,  3.890000000000000E-12_wp, &
      &-2.920000000000000E-12_wp,  6.543061464970300E-01_wp, -7.350380400000000E-07_wp, &
      & 2.770000000000000E-12_wp, -3.271534756103600E-01_wp,  2.415001173009400E-01_wp, &
      & 1.500000000000000E-13_wp, -3.271526708881400E-01_wp, -2.414993822564200E-01_wp, &
      &-5.610000000000000E-12_wp, -7.350380400000000E-07_wp,  3.970231511250400E-01_wp, &
      & 4.930000000000000E-12_wp,  1.934311946184400E-01_wp, -1.985118468089900E-01_wp, &
      & 6.900000000000000E-13_wp, -1.934304595822700E-01_wp, -1.985113043193300E-01_wp, &
      &-4.821347568000000E-05_wp,  2.770000000000000E-12_wp,  4.930000000000000E-12_wp, &
      & 4.418955591000000E-05_wp, -1.980000000000000E-12_wp, -1.440000000000000E-12_wp, &
      & 4.023919770000000E-06_wp, -8.000000000000000E-13_wp, -3.490000000000000E-12_wp, &
      & 2.240000000000000E-12_wp, -3.271534756103600E-01_wp,  1.934311946184400E-01_wp, &
      &-1.980000000000000E-12_wp,  3.487418139707800E-01_wp, -2.174656208231900E-01_wp, &
      &-2.600000000000000E-13_wp, -2.158833835586000E-02_wp,  2.403442620289000E-02_wp, &
      & 1.730000000000000E-12_wp,  2.415001173009400E-01_wp, -1.985118468089900E-01_wp, &
      &-1.440000000000000E-12_wp, -2.174656208231900E-01_wp,  1.884251218631200E-01_wp, &
      &-2.900000000000000E-13_wp, -2.403449647736000E-02_wp,  1.008672494990000E-02_wp, &
      &-4.821375342000000E-05_wp,  1.500000000000000E-13_wp,  6.900000000000000E-13_wp, &
      & 4.023919770000000E-06_wp, -2.600000000000000E-13_wp, -2.900000000000000E-13_wp, &
      & 4.418983365000000E-05_wp,  1.100000000000000E-13_wp, -4.000000000000000E-13_wp, &
      & 6.900000000000000E-13_wp, -3.271526708881400E-01_wp, -1.934304595822700E-01_wp, &
      &-8.000000000000000E-13_wp, -2.158833835586000E-02_wp, -2.403449647736000E-02_wp, &
      & 1.100000000000000E-13_wp,  3.487410092408900E-01_wp,  2.174649560550100E-01_wp, &
      & 3.890000000000000E-12_wp, -2.414993822564200E-01_wp, -1.985113043193300E-01_wp, &
      &-3.490000000000000E-12_wp,  2.403442620289000E-02_wp,  1.008672494990000E-02_wp, &
      &-4.000000000000000E-13_wp,  2.174649560550100E-01_wp,  1.884245793686800E-01_wp], &
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 1.56104730000000000E-06_wp
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

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-4_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N), source=0.0_wp)
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)

   if (any(abs(hessian - hessian_ref) > thr1)) then
      call test_failed(error, "Final Hessians do not match")

      print *, "--- hessian ---"
      do i = 1, N
         print "(*(F21.14))", hessian(i, :)
      end do

      print *, "--- Ref. hessian ---"
      do i = 1, N
         print "(*(F21.14))", hessian_ref(i, :)
      end do
      return
   end if

   if (abs(final_err - final_err_ref) > thr2) then
      call test_failed(error, "Final error does not match")
      print *, "--- Final error ---"
      print "(*(F21.14))", final_err
      print *, "--- Ref. final error ---"
      print "(*(F21.14))", final_err_ref
   end if
end subroutine test_o1numhess_gfn1

!> Pins the GFN2 O1NumHess (ODLR) Hessian of water and its final error
!> estimate against the reference values below.
subroutine test_o1numhess_gfn2(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr1 = 1.4901161193847656e-08_wp, thr2 = 1.4901161193847656e-08_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000034546_wp, 0.18900383618455_wp, &
      & 0.00000000000000_wp, 1.45674735348811_wp, -0.88650486059828_wp, &
      &-0.00000000000000_wp, -1.45674735383357_wp, -0.88650486086986_wp], &
      & shape(xyz))
   real(wp), parameter :: hessian_ref(9, 9) = reshape([&
      &-7.853052671800001E-04_wp,  8.700000000000000E-13_wp, -7.300000000000000E-13_wp, &
      & 3.926527639100000E-04_wp, -8.000000000000000E-14_wp,  6.800000000000000E-13_wp, &
      & 3.926525032700000E-04_wp, -7.900000000000000E-13_wp,  4.000000000000000E-14_wp, &
      & 8.700000000000000E-13_wp,  6.131573464383200E-01_wp, -7.350901300000000E-07_wp, &
      &-6.100000000000000E-13_wp, -3.065790827913500E-01_wp,  2.266356332227100E-01_wp, &
      &-2.600000000000000E-13_wp, -3.065782636468700E-01_wp, -2.266348981325700E-01_wp, &
      &-7.300000000000000E-13_wp, -7.350901300000000E-07_wp,  4.009666311434500E-01_wp, &
      & 6.800000000000000E-13_wp,  1.851054123830400E-01_wp, -2.004835868517400E-01_wp, &
      & 5.000000000000000E-14_wp, -1.851046772928600E-01_wp, -2.004830442917700E-01_wp, &
      & 3.926527639100000E-04_wp, -6.100000000000000E-13_wp,  6.800000000000000E-13_wp, &
      &-5.034358995600001E-04_wp, -3.000000000000000E-14_wp, -5.600000000000000E-13_wp, &
      & 1.107831356500000E-04_wp,  6.400000000000000E-13_wp, -1.100000000000000E-13_wp, &
      &-8.000000000000000E-14_wp, -3.065790827913500E-01_wp,  1.851054123830400E-01_wp, &
      &-3.000000000000000E-14_wp,  3.407119356534000E-01_wp, -2.058704902859900E-01_wp, &
      & 1.100000000000000E-13_wp, -3.413285286200000E-02_wp,  2.076507790304000E-02_wp, &
      & 6.800000000000000E-13_wp,  2.266356332227100E-01_wp, -2.004835868517400E-01_wp, &
      &-5.600000000000000E-13_wp, -2.058704902859900E-01_wp,  1.835966136036200E-01_wp, &
      &-1.200000000000000E-13_wp, -2.076514293670000E-02_wp,  1.688697324816000E-02_wp, &
      & 3.926525032700000E-04_wp, -2.600000000000000E-13_wp,  5.000000000000000E-14_wp, &
      & 1.107831356500000E-04_wp,  1.100000000000000E-13_wp, -1.200000000000000E-13_wp, &
      &-5.034356389200000E-04_wp,  1.500000000000000E-13_wp,  7.000000000000000E-14_wp, &
      &-7.900000000000000E-13_wp, -3.065782636468700E-01_wp, -1.851046772928600E-01_wp, &
      & 6.400000000000000E-13_wp, -3.413285286200000E-02_wp, -2.076514293670000E-02_wp, &
      & 1.500000000000000E-13_wp,  3.407111165087200E-01_wp,  2.058698202294700E-01_wp, &
      & 4.000000000000000E-14_wp, -2.266348981325700E-01_wp, -2.004830442917700E-01_wp, &
      &-1.100000000000000E-13_wp,  2.076507790304000E-02_wp,  1.688697324816000E-02_wp, &
      & 7.000000000000000E-14_wp,  2.058698202294700E-01_wp,  1.835960710436500E-01_wp], &
      & shape(hessian_ref))
   real(wp), parameter :: final_err_ref = 4.02917678000000000E-06_wp
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

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-4_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N), source=0.0_wp)
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)

   if (any(abs(hessian - hessian_ref) > thr1)) then
      call test_failed(error, "Final Hessians do not match")

      print *, "--- hessian ---"
      do i = 1, N
         print "(*(F21.14))", hessian(i, :)
      end do

      print *, "--- Ref. hessian ---"
      do i = 1, N
         print "(*(F21.14))", hessian_ref(i, :)
      end do
      return
   end if

   if (abs(final_err - final_err_ref) > thr2) then
      call test_failed(error, "Final error does not match")
      print *, "--- Final error ---"
      print "(*(F21.14))", final_err
      print *, "--- Ref. final error ---"
      print "(*(F21.14))", final_err_ref
   end if
end subroutine test_o1numhess_gfn2

!> Pins the GFN1 O1NumHess Hessian of linear water through its eigenvalue
!> spectrum: three near-zero and two negative frequencies around -0.16.
subroutine test_o1numhess_linear_h2o_gfn1(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, -1.81075448577205_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 1.81075448676713_wp], &
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

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-4_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N), source=0.0_wp)
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)
   allocate(freq(N))
   lwork = 1 + 6 * N + 2 * N**2
   allocate(aux(lwork))
   call dsyev("V", "U", N, hessian, N, freq, aux, lwork, info)
   call check(error, info, 0)
   if (allocated(error)) return
   if (count(abs(freq) < 1.0e-4_wp) /= 3) then
      call test_failed(error, "Linear H2O should have exactly three ~0 freqs")
      return
   end if

   if (count(freq < -1.0e-10_wp) /= 2) then
      call test_failed(error, "Linear H2O should have exactly two negative freqs")
      return
   end if

   ! freqs around -0.16
   if (freq(1) > -0.16_wp .or. freq(2) > -0.16_wp) then
      call test_failed(error, "First two freqs should be negative")
      print *, freq
      return
   end if
end subroutine test_o1numhess_linear_h2o_gfn1

!> Pins the GFN2 O1NumHess Hessian of linear water through its eigenvalue
!> spectrum: three near-zero and two negative frequencies around -0.3.
subroutine test_o1numhess_linear_h2o_gfn2(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   real(wp), parameter :: thr1 = 1.0e-9_wp, thr2 = 1.0e-5_wp
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000000000_wp, 0.00000000000000_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, -1.81075448577205_wp, &
      & 0.00000000000000_wp, 0.00000000000000_wp, 1.81075448676713_wp], &
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

   call newXTBCalculator(env, mol, calc, method=2, accuracy=1.0e-4_wp)
   call newWavefunction(env, mol, calc, chk)

   allocate(hessian(N, N), dipgrad_dummy(3, N), source=0.0_wp)
   call calc%hessian(env, mol, chk, [(i, i=1, mol%n)], step, hessian, dipgrad_dummy, odlr=.true., final_err=final_err)
   allocate(freq(N))
   lwork = 1 + 6 * N + 2 * N**2
   allocate(aux(lwork))
   call dsyev("V", "U", N, hessian, N, freq, aux, lwork, info)
   call check(error, info, 0)
   if (allocated(error)) return

   if (count(abs(freq) < 1.0e-4_wp) /= 3) then
      call test_failed(error, "Linear H2O should have exactly three ~0 freqs")
      return
   end if

   if (count(freq < -1.0e-10_wp) /= 2) then
      call test_failed(error, "Linear H2O should have exactly two negative freqs")
      return
   end if

   ! freqs around -0.3
   if (freq(1) > -0.3_wp .or. freq(2) > -0.3_wp) then
      call test_failed(error, "First two freqs should be negative")
      return
   end if
end subroutine test_o1numhess_linear_h2o_gfn2

!> Pins the compliance constants of a two-atom harmonic model against their
!> analytic values, including the redundant B-matrix with duplicated rows.
subroutine test_compliance(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 2
   integer, parameter :: ndim = 3 * nat
   real(wp), parameter :: force_constant = 2.0_wp
   real(wp), parameter :: thr = 1.0e-10_wp
   integer, parameter :: pairs(2, 1) = reshape([1, 2], [2, 1])
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp], shape(xyz))

   type(TEnvironment) :: env
   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals
   integer :: i, j
   logical :: failed
   real(wp), allocatable :: bmat(:, :), compliance(:, :), hessian(:, :)
   real(wp) :: redundant_bmat(2, ndim), redundant_compliance(2, 2)

   call init(env)
   env%unit = 0
   neigh_list = build_neigh_list(nat, pairs)
   call init(graph, neigh_list)
   call init(internals, graph, xyz)
   call check(error, internals%ncoords, 1)
   if (allocated(error)) return

   allocate(bmat(internals%ncoords, ndim), compliance(internals%ncoords, internals%ncoords), &
      & hessian(ndim, ndim))
   call get_bmatrix(internals, xyz, bmat)

   do i = 1, ndim
      do j = 1, ndim
         hessian(i, j) = force_constant * bmat(1, i) * bmat(1, j)
      end do
   end do

   call compute_compliance(env, hessian, bmat, xyz, nat, internals%ncoords, compliance)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Compliance calculation failed")
      return
   end if
   call check(error, bmat(1, 1), -1.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, bmat(1, 4), 1.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, compliance(1, 1), 1.0_wp/force_constant, thr=thr)
   if (allocated(error)) return

   redundant_bmat(1, :) = bmat(1, :)
   redundant_bmat(2, :) = bmat(1, :)
   call compute_compliance(env, hessian, redundant_bmat, xyz, nat, 2, redundant_compliance)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Compliance calculation failed")
      return
   end if
   call check(error, redundant_compliance(1, 1), 1.0_wp/force_constant, thr=thr)
   if (allocated(error)) return
   call check(error, redundant_compliance(1, 2), 1.0_wp/force_constant, thr=thr)
   if (allocated(error)) return
   call check(error, redundant_compliance(2, 1), 1.0_wp/force_constant, thr=thr)
   if (allocated(error)) return
   call check(error, redundant_compliance(2, 2), 1.0_wp/force_constant, thr=thr)

end subroutine test_compliance

!> Pins the GFN1 compliance constants of water at the default Hessian step
!> against the reference matrix below.
subroutine test_compliance_water(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   integer, parameter :: ndim = 3 * nat
   real(wp), parameter :: step = 5.0e-3_wp
   integer, parameter :: at(nat) = [8, 1, 1]
   character(len=*), parameter :: sym(nat) = ["O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.00000000000000_wp, 0.00000000034546_wp, 0.18900383618455_wp, &
      & 0.00000000000000_wp, 1.45674735348811_wp, -0.88650486059828_wp, &
      &-0.00000000000000_wp, -1.45674735383357_wp, -0.88650486086986_wp], &
      & shape(xyz))
   real(wp), parameter :: thr = 1.0e-10_wp
   !> GFN1 compliance constants of water, in (E_h/Bohr^2)^{-1}, in the
   !> redundant ordering bond(1,2), bond(1,3), angle(2,1,3):
   !> the two OH bonds are symmetry equivalent, so C(1,1) == C(2,2)
   real(wp), parameter :: ref(3, 3) = reshape([&
      &  2.0178418732549255_wp, 3.9503605639448727e-02_wp, -3.8479908315707656e-01_wp, &
      &  3.9503605639448727e-02_wp, 2.0178418810941556_wp, -3.8479908613342922e-01_wp, &
      & -3.8479908315707656e-01_wp, -3.8479908613342922e-01_wp, 8.3576726219491668_wp], &
      & shape(ref))

   type(TMolecule) :: mol
   type(TRestart) :: chk
   type(TEnvironment) :: env
   type(scc_results) :: res
   type(TxTBCalculator) :: calc

   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals
   integer :: i, j
   logical :: failed
   integer, allocatable :: list(:)
   real(wp) :: energy, sigma(3, 3), hl_gap
   real(wp), allocatable :: gradient(:, :), dipgrad(:, :), hessian(:, :)
   real(wp), allocatable :: bmat(:, :), compliance(:, :)

   call init(env)
   env%unit = 0
   call init(mol, sym, xyz)

   allocate(gradient(3, nat), dipgrad(3, ndim), hessian(ndim, ndim))
   energy = 0.0_wp
   gradient = 0.0_wp

   call newXTBCalculator(env, mol, calc, method=1, accuracy=1.0e-4_wp)
   call newWavefunction(env, mol, calc, chk)
   call calc%singlepoint(env, mol, chk, 2, .false., energy, gradient, sigma, &
      & hl_gap, res)

   dipgrad = 0.0_wp
   hessian = 0.0_wp
   list = [(i, i = 1, nat)]
   call calc%hessian(env, mol, chk, list, step, hessian, dipgrad)

   call init(neigh_list, nat)
   call neigh_list%generate_covalent(at, xyz)
   call init(graph, neigh_list)
   call init(internals, graph, xyz)
   call check(error, internals%ncoords, 3)
   if (allocated(error)) return

   allocate(bmat(internals%ncoords, ndim), compliance(internals%ncoords, internals%ncoords))
   call get_bmatrix(internals, xyz, bmat)
   call compute_compliance(env, hessian, bmat, xyz, nat, internals%ncoords, compliance)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Compliance calculation failed")
      return
   end if

   do i = 1, internals%ncoords
      do j = 1, internals%ncoords
         call check(error, compliance(i, j), ref(i, j), thr=thr)
         if (allocated(error)) return
      end do
   end do

end subroutine test_compliance_water

!> Pins a redundant-coordinate B-matrix against central finite differences.
subroutine test_redundant_bmatrix_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   integer, parameter :: ncoord = 6
   integer, parameter :: pairs(2, 3) = reshape([1, 2, 2, 3, 3, 4], [2, 3])
   integer, parameter :: expected_kind(ncoord) = [ &
      & coord_bond, coord_bond, coord_bond, coord_angle, coord_angle, coord_dihedral]
   integer, parameter :: expected_atoms(4, ncoord) = reshape([ &
      & 1, 2, 0, 0, &
      & 2, 3, 0, 0, &
      & 3, 4, 0, 0, &
      & 1, 2, 3, 0, &
      & 2, 3, 4, 0, &
      & 1, 2, 3, 4], [4, ncoord])
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 0.0_wp, 0.0_wp, &
      & 2.5_wp, 1.0_wp, 0.0_wp, &
      & 3.0_wp, 1.0_wp, 1.0_wp], shape(xyz))
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), parameter :: thr = 1.0e-6_wp

   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals, right, left
   real(wp) :: coord(3, nat), dq
   real(wp), allocatable :: bmat(:, :)
   integer :: i, ic, ii, k

   neigh_list = build_neigh_list(nat, pairs)
   call init(graph, neigh_list)
   call init(internals, graph, xyz)
   call check(error, internals%ncoords, ncoord)
   if (allocated(error)) return

   do k = 1, ncoord
      call check(error, internals%kind(k), expected_kind(k))
      if (allocated(error)) return
      do i = 1, 4
         call check(error, internals%atoms(i, k), expected_atoms(i, k))
         if (allocated(error)) return
      end do
   end do

   allocate(bmat(ncoord, 3*nat))
   call get_bmatrix(internals, xyz, bmat)

   do i = 1, nat
      do ic = 1, 3
         ii = 3 * (i - 1) + ic
         coord = xyz
         coord(ic, i) = coord(ic, i) + step
         call init(right, graph, coord)
         coord = xyz
         coord(ic, i) = coord(ic, i) - step
         call init(left, graph, coord)
         do k = 1, ncoord
            dq = right%q(k) - left%q(k)
            if (internals%kind(k) == coord_dihedral) then
               dq = modulo(dq + pi, 2.0_wp*pi) - pi
            end if
            call check(error, bmat(k, ii), 0.5_wp*dq/step, thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_redundant_bmatrix_fd

!> A straight four-atom chain has three stretches and both transverse bend
!> directions at each centre, but no defined torsion.
subroutine test_linear_chain_compliance(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   integer, parameter :: ndim = 3 * nat
   integer, parameter :: ncoord = 7
   integer, parameter :: pairs(2, 3) = reshape([1, 2, 2, 3, 3, 4], [2, 3])
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 2.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 0.0_wp, 0.0_wp], shape(xyz))
   real(wp), parameter :: thr = 1.0e-9_wp

   type(TEnvironment) :: env
   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals
   integer :: i, j
   logical :: failed
   real(wp) :: expected_bmat(ncoord, ndim), hessian(ndim, ndim)
   real(wp), allocatable :: bmat(:, :), compliance(:, :)

   call init(env)
   env%unit = 0
   neigh_list = build_neigh_list(nat, pairs)
   call init(graph, neigh_list)
   call init(internals, graph, xyz)
   call check(error, internals%ncoords, ncoord)
   if (allocated(error)) return
   call check(error, internals%ndihedral, 0)
   if (allocated(error)) return
   do i = 1, 3
      call check(error, internals%kind(i), coord_bond)
      if (allocated(error)) return
   end do
   do i = 4, ncoord
      call check(error, internals%kind(i), coord_linbend)
      if (allocated(error)) return
   end do

   expected_bmat = 0.0_wp
   do i = 1, 3
      expected_bmat(i, 3*(i - 1) + 1) = -1.0_wp
      expected_bmat(i, 3*i + 1) = 1.0_wp
   end do
   expected_bmat(4, [3, 6, 9]) = [1.0_wp, -2.0_wp, 1.0_wp]
   expected_bmat(5, [2, 5, 8]) = [1.0_wp, -2.0_wp, 1.0_wp]
   expected_bmat(6, [6, 9, 12]) = [1.0_wp, -2.0_wp, 1.0_wp]
   expected_bmat(7, [5, 8, 11]) = [1.0_wp, -2.0_wp, 1.0_wp]
   hessian = matmul(transpose(expected_bmat), expected_bmat)

   allocate(bmat(ncoord, ndim), compliance(ncoord, ncoord))
   call get_bmatrix(internals, xyz, bmat)
   call compute_compliance(env, hessian, bmat, xyz, nat, ncoord, compliance)
   call env%check(failed)
   if (failed) then
      call test_failed(error, "Compliance calculation failed")
      return
   end if
   do i = 1, ncoord
      do j = 1, ncoord
         if (i == j) then
            call check(error, compliance(i, j), 1.0_wp, thr=thr)
         else
            call check(error, compliance(i, j), 0.0_wp, thr=thr)
         end if
         if (allocated(error)) return
      end do
   end do

end subroutine test_linear_chain_compliance


!> Near-linear coordinates use the reference frame even when their B rows are
!> evaluated at a displaced geometry.
subroutine test_near_linear_bmatrix_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   integer, parameter :: ncoord = 4
   integer, parameter :: pairs(2, 2) = reshape([1, 2, 2, 3], [2, 2])
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), parameter :: thr = 1.0e-7_wp
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -1.0_wp, 0.0_wp, 0.0_wp, &
      &  0.0_wp, 0.0_wp, 0.0_wp, &
      &  cos(pi/180.0_wp), sin(pi/180.0_wp), 0.0_wp], shape(xyz))

   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals
   real(wp) :: coord(3, nat), right(3, nat), left(3, nat)
   real(wp) :: v1(3), v2(3), qref, qright, qleft
   real(wp), allocatable :: bmat(:, :)
   integer :: i, ic, ii, k

   neigh_list = build_neigh_list(nat, pairs)
   call init(graph, neigh_list)
   call init(internals, graph, xyz)
   call check(error, internals%ncoords, ncoord)
   if (allocated(error)) return
   do k = 3, 4
      call check(error, internals%kind(k), coord_linbend)
      if (allocated(error)) return
      v1 = xyz(:, 1) - xyz(:, 2)
      v2 = xyz(:, 3) - xyz(:, 2)
      qref = dot_product(internals%frame(:, k), v1/norm2(v1) + v2/norm2(v2))
      call check(error, internals%q(k), qref, thr=thr)
      if (allocated(error)) return
   end do

   coord = xyz
   coord(2, 1) = coord(2, 1) + 0.02_wp
   allocate(bmat(ncoord, 3*nat))
   call get_bmatrix(internals, coord, bmat)
   do i = 1, nat
      do ic = 1, 3
         ii = 3 * (i - 1) + ic
         right = coord
         right(ic, i) = right(ic, i) + step
         left = coord
         left(ic, i) = left(ic, i) - step
         do k = 3, 4
            v1 = right(:, 1) - right(:, 2)
            v2 = right(:, 3) - right(:, 2)
            qright = dot_product(internals%frame(:, k), v1/norm2(v1) + v2/norm2(v2))
            v1 = left(:, 1) - left(:, 2)
            v2 = left(:, 3) - left(:, 2)
            qleft = dot_product(internals%frame(:, k), v1/norm2(v1) + v2/norm2(v2))
            call check(error, bmat(k, ii), 0.5_wp*(qright - qleft)/step, thr=thr)
            if (allocated(error)) return
         end do
      end do
   end do

end subroutine test_near_linear_bmatrix_fd


!> The driver writes two bonds and one angle for water, none for separated atoms.
subroutine test_compliance_driver_redundant(error)
   type(error_type), allocatable, intent(out) :: error

   type(TEnvironment) :: env
   type(TMolecule) :: mol
   real(wp) :: xyz(3, 3), hess(9, 9)
   integer :: unit, dat_unit, ios, i, icase, counts(3), expected(3)
   logical :: failed
   character(256) :: line
   character(4) :: kind

   xyz(:, 1) = 0.0_wp
   xyz(:, 2) = [0.0_wp, 1.8_wp, 0.0_wp]
   xyz(:, 3) = [1.7_wp, -0.5_wp, 0.0_wp]
   call init(mol, [8, 1, 1], xyz)
   hess = 0.0_wp
   do i = 1, 9
      hess(i, i) = 1.0_wp
   end do
   call init(env)
   open(newunit=unit, status="scratch")
   env%unit = unit

   expected = [3, 2, 1] ! Total, bond and angle records.
   do icase = 1, 2
      if (icase == 2) then
         mol%xyz = 10.0_wp * mol%xyz
         expected = 0
      end if
      call compliance_driver(env, mol, hess)
      call env%check(failed)
      call check(error, .not. failed)
      if (allocated(error)) exit

      counts = 0
      open(newunit=dat_unit, file="compliance.dat", status="old")
      do
         read(dat_unit, "(a)", iostat=ios) line
         if (ios /= 0) exit
         if (index(line, "# coord ") /= 1) cycle
         read(line(8:), *) kind
         counts(1) = counts(1) + 1
         if (kind == "bond") counts(2) = counts(2) + 1
         if (kind == "ang") counts(3) = counts(3) + 1
      end do
      close(dat_unit, status="delete")
      do i = 1, 3
         call check(error, counts(i), expected(i))
         if (allocated(error)) exit
      end do
      if (allocated(error)) exit
   end do
   close(unit)

end subroutine test_compliance_driver_redundant


!> Locks covalent neighbour-list generation: canonical half-list storage,
!> replacement semantics, and the supported-cutoff tolerance margin for
!> getNeighs.
subroutine test_covalent_neighbour_list(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   integer, parameter :: at(nat) = [8, 1, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.8_wp, 0.0_wp, &
      & 1.7_wp, -0.5_wp, 0.0_wp, &
      & 10.0_wp, 0.0_wp, 0.0_wp], shape(xyz))
   integer, parameter :: deg_ref(nat) = [2, 1, 1, 0]
   real(wp), parameter :: dist2_ref(2) = [3.14_wp, 3.24_wp]
   real(wp), parameter :: thr = 1.0e-10_wp

   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   integer :: neighs(nat), i

   ! Empty list must have zero supported cutoff.
   call init(neigh_list, 0)
   call neigh_list%generate_covalent(at(:0), xyz(:, :0))
   call check(error, neigh_list%cutoff, 0.0_wp)
   if (allocated(error)) return

   ! O/H/H plus one isolated H: two O-H bonds, no H-H bond
   call init(neigh_list, nat)
   call neigh_list%generate_covalent(at, xyz)
   call check(error, neigh_list%neighs(1), 2)
   if (allocated(error)) return
   call check(error, all(neigh_list%neighs(2:4) == 0), .true.)
   if (allocated(error)) return
   call check(error, all(neigh_list%image == [(i, i = 1, nat)]), .true.)
   if (allocated(error)) return
   call check(error, all(neigh_list%iNeigh(0, :) == [(i, i = 1, nat)]), .true.)
   if (allocated(error)) return
   call check(error, all(abs(neigh_list%dist2(1:2, 1) - dist2_ref) < thr), .true.)
   if (allocated(error)) return

   call init(graph, neigh_list)
   do i = 1, nat
      call check(error, graph%neighs(i), deg_ref(i))
      if (allocated(error)) return
   end do
   call check(error, graph%parent(4), 0)
   if (allocated(error)) return

   ! Replacement semantics: regenerating the same geometry must replace
   ! (not append to) the previous list
   call init(neigh_list, nat)
   call neigh_list%generate_covalent(at, xyz)
   call neigh_list%generate_covalent(at, xyz)
   call check(error, neigh_list%neighs(1), 2)
   if (allocated(error)) return
   call check(error, all(neigh_list%neighs(2:4) == 0), .true.)
   if (allocated(error)) return

   ! Default getNeighs must return every stored neighbour, also when the
   ! bond sits infinitesimally below the largest pair threshold
   call init(neigh_list, 3)
   call neigh_list%generate_covalent([1, 1, 1], &
      & reshape([0.0_wp, 0.0_wp, 0.0_wp, &
      &          0.5_wp, 0.0_wp, 0.0_wp, &
      &          2.0_wp*covalent_radius_d3(1) - 1.0e-12_wp, 0.0_wp, 0.0_wp], &
      &         [3, 3]))
   call neigh_list%getNeighs(neighs(:3))
   call check(error, neighs(1), 2)
   if (allocated(error)) return
   call check(error, neighs(2), 1)
   if (allocated(error)) return
   call check(error, neighs(3), 0)

end subroutine test_covalent_neighbour_list


!> Projecting a linear saddle point must keep both imaginary bending modes
!> and remove exactly three translations and two rotations.
subroutine test_linear_projection(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   character(len=*), parameter :: sym(nat) = ["H", "O", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 1.8_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, -1.8_wp], shape(xyz))
   real(wp), parameter :: thr = 1.0e-10_wp

   type(TMolecule) :: mol
   real(wp), allocatable :: hessian(:, :), eig(:)
   integer :: nModes

   call init(mol, sym, xyz)
   hessian = linear_model_hessian(mol)

   call check(error, mol%linear)
   if (allocated(error)) return

   call projectHessian(hessian, mol, .true., .true., nModes=nModes)
   call check(error, nModes, 5)
   if (allocated(error)) return

   call get_eigenvalues(error, hessian, eig)
   if (allocated(error)) return
   call check(error, count(abs(eig) < thr), 5)
   if (allocated(error)) return
   call check(error, count(eig < -thr), 2)

end subroutine test_linear_projection


!> Fixing atoms of a non-linear molecule removes their coordinates and leaves
!> only rotations about axes through all fixed atoms. The vibrations must
!> match the mass-weighted Hessian with the fixed atoms deleted.
subroutine test_fixed_projection(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 4
   character(len=*), parameter :: sym(nat) = ["O", "O", "H", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.00_wp, 1.38_wp, -0.12_wp, &
      & 0.00_wp, -1.38_wp, -0.12_wp, &
      & 1.62_wp, 1.78_wp, 0.95_wp, &
      & -1.62_wp, -1.78_wp, 0.95_wp], shape(xyz))

   type(TMolecule) :: mol
   real(wp), allocatable :: hessian(:, :)

   call init(mol, sym, xyz)
   hessian = spring_model_hessian(mol)

   call check(error, .not.mol%linear)
   if (allocated(error)) return

   ! one fixed atom: three rotations about the fixed atom remain
   call check_fixed_projection(error, mol, hessian, [1], 6)
   if (allocated(error)) return

   ! two fixed atoms: only the rotation about their connecting axis remains
   call check_fixed_projection(error, mol, hessian, [1, 2], 7)
   if (allocated(error)) return

   ! three non-collinear fixed atoms: no rotation remains
   call check_fixed_projection(error, mol, hessian, [1, 3, 4], 9)

end subroutine test_fixed_projection


!> Fixing atoms of a linear molecule leaves only the rotations about the fixed
!> atoms perpendicular to the molecular axis. The center of mass is off the
!> fixed atom, so rotations about it would wrongly be projected as well.
subroutine test_fixed_linear_projection(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: nat = 3
   character(len=*), parameter :: sym(nat) = ["H", "O", "H"]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & 0.0_wp, 0.0_wp, 1.8_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, -2.1_wp], shape(xyz))

   type(TMolecule) :: mol
   real(wp), allocatable :: hessian(:, :)

   call init(mol, sym, xyz)
   hessian = linear_model_hessian(mol)

   call check(error, mol%linear)
   if (allocated(error)) return

   ! fixed central atom: two rotations about the central atom remain
   call check_fixed_projection(error, mol, hessian, [2], 5)
   if (allocated(error)) return

   ! two fixed atoms: no rotation remains
   call check_fixed_projection(error, mol, hessian, [1, 3], 6)

end subroutine test_fixed_linear_projection


!> Project a mass-weighted Hessian with fixed atoms and compare its spectrum
!> against the mass-weighted Hessian with the fixed atoms deleted.
subroutine check_fixed_projection(error, mol, hessian, fixed, nref)

   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Mass-weighted model Hessian invariant under rigid-body motions
   real(wp), intent(in) :: hessian(:, :)

   !> Fixed atoms
   integer, intent(in) :: fixed(:)

   !> Expected number of projected modes
   integer, intent(in) :: nref

   real(wp), parameter :: thr = 1.0e-10_wp
   real(wp), allocatable :: projected(:, :), eig(:), eigref(:)
   integer, allocatable :: free(:)
   integer :: nModes, iat, ic

   projected = hessian
   call projectHessian(projected, mol, .false., .true., fixedAtoms=fixed, &
      & nModes=nModes)
   call check(error, nModes, nref)
   if (allocated(error)) return

   ! coordinates of all atoms that are not fixed
   free = [((3*(iat-1)+ic, ic = 1, 3), iat = 1, mol%n)]
   free = pack(free, [(all(fixed /= iat), all(fixed /= iat), &
      & all(fixed /= iat), iat = 1, mol%n)])

   call get_eigenvalues(error, projected, eig)
   if (allocated(error)) return
   call get_eigenvalues(error, hessian(free, free), eigref)
   if (allocated(error)) return

   ! deleted Cartesian coordinates are not part of the reduced Hessian
   call check(error, count(abs(eig) < thr), nref)
   if (allocated(error)) return
   call check(error, count(abs(eigref) < thr), nref - 3*size(fixed))
   if (allocated(error)) return

   eig = pack(eig, abs(eig) >= thr)
   eigref = pack(eigref, abs(eigref) >= thr)
   call check(error, size(eig), size(eigref))
   if (allocated(error)) return
   call check(error, maxval(abs(eig - eigref)), 0.0_wp, thr=thr)

end subroutine check_fixed_projection


!> Mass-weighted Hessian of harmonic springs between all atom pairs at their
!> rest lengths, which is exactly invariant under rigid-body motions.
function spring_model_hessian(mol) result(hessian)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Mass-weighted Hessian
   real(wp), allocatable :: hessian(:, :)

   real(wp) :: vec(3), block(3, 3), force
   integer :: iat, jat, ii, jj

   allocate(hessian(3*mol%n, 3*mol%n), source=0.0_wp)
   do iat = 1, mol%n
      do jat = 1, iat - 1
         ! distinct force constants avoid degenerate eigenvalues
         force = 0.1_wp * real(iat + 2*jat, wp)
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         vec(:) = vec / norm2(vec)
         block(:, :) = force * spread(vec, 1, 3) * spread(vec, 2, 3)
         ii = 3*(iat-1)
         jj = 3*(jat-1)
         hessian(ii+1:ii+3, ii+1:ii+3) = hessian(ii+1:ii+3, ii+1:ii+3) + block
         hessian(jj+1:jj+3, jj+1:jj+3) = hessian(jj+1:jj+3, jj+1:jj+3) + block
         hessian(ii+1:ii+3, jj+1:jj+3) = hessian(ii+1:ii+3, jj+1:jj+3) - block
         hessian(jj+1:jj+3, ii+1:ii+3) = hessian(jj+1:jj+3, ii+1:ii+3) - block
      end do
   end do

   call mass_weight(mol, hessian)

end function spring_model_hessian


!> Mass-weighted Hessian of a linear triatomic along z with stretching springs
!> and a negative bending force constant, i.e. a second-order saddle point
!> like linear water. It is exactly invariant under rigid-body motions.
function linear_model_hessian(mol) result(hessian)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Mass-weighted Hessian
   real(wp), allocatable :: hessian(:, :)

   real(wp), parameter :: kbend = -0.05_wp
   real(wp), parameter :: kstretch(3, 3) = reshape([ &
      & 0.00_wp, 0.50_wp, 0.01_wp, &
      & 0.50_wp, 0.00_wp, 0.45_wp, &
      & 0.01_wp, 0.45_wp, 0.00_wp], [3, 3])
   real(wp) :: r1, r2, bend(3)
   integer :: iat, jat, ic

   r1 = abs(mol%xyz(3, 1) - mol%xyz(3, 2))
   r2 = abs(mol%xyz(3, 3) - mol%xyz(3, 2))
   ! derivative of the bending angle with respect to the x or y coordinates
   bend(:) = [1.0_wp/r1, -(1.0_wp/r1 + 1.0_wp/r2), 1.0_wp/r2]

   allocate(hessian(9, 9), source=0.0_wp)
   do iat = 1, 3
      do jat = 1, 3
         do ic = 1, 2
            hessian(3*(iat-1)+ic, 3*(jat-1)+ic) = kbend * bend(iat) * bend(jat)
         end do
         if (iat == jat) then
            hessian(3*iat, 3*jat) = sum(kstretch(:, iat))
         else
            hessian(3*iat, 3*jat) = -kstretch(jat, iat)
         end if
      end do
   end do

   call mass_weight(mol, hessian)

end function linear_model_hessian


!> Transform a Cartesian Hessian to mass-weighted coordinates
subroutine mass_weight(mol, hessian)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Hessian, mass-weighted on exit
   real(wp), intent(inout) :: hessian(:, :)

   real(wp), allocatable :: invsqrtmass(:)
   integer :: ii, jj

   invsqrtmass = [(1.0_wp/sqrt(mol%atmass((ii-1)/3+1)), ii = 1, 3*mol%n)]
   do ii = 1, size(hessian, 2)
      do jj = 1, size(hessian, 1)
         hessian(jj, ii) = hessian(jj, ii) * invsqrtmass(jj) * invsqrtmass(ii)
      end do
   end do

end subroutine mass_weight


!> Eigenvalues of a symmetric matrix in ascending order
subroutine get_eigenvalues(error, matrix, eig)

   !> Failure report, allocated when the diagonalization fails.
   type(error_type), allocatable, intent(out) :: error

   !> Symmetric matrix
   real(wp), intent(in) :: matrix(:, :)

   !> Eigenvalues in ascending order
   real(wp), allocatable, intent(out) :: eig(:)

   real(wp), allocatable :: work(:), scratch(:, :)
   integer :: ndim, info

   ndim = size(matrix, 1)
   scratch = matrix
   allocate(eig(ndim), work(3*ndim))
   call dsyev("N", "U", ndim, scratch, ndim, eig, work, size(work), info)
   call check(error, info, 0)

end subroutine get_eigenvalues


!> Build a symmetric neighbour list from undirected atom pairs.
function build_neigh_list(nat, pairs) result(neigh_list)
   !> Number of atoms.
   integer, intent(in) :: nat
   !> Undirected atom pairs, dimension (2, number of pairs).
   integer, intent(in) :: pairs(:, :)
   !> Neighbour list with both adjacency columns filled.
   type(TNeighbourList) :: neigh_list

   integer :: i, j, k

   call init(neigh_list, nat)
   neigh_list%neighs = 0
   neigh_list%iNeigh = 0
   neigh_list%image = [(i, i = 1, nat)]
   do k = 1, size(pairs, 2)
      i = pairs(1, k)
      j = pairs(2, k)
      neigh_list%neighs(i) = neigh_list%neighs(i) + 1
      neigh_list%iNeigh(neigh_list%neighs(i), i) = j
      neigh_list%neighs(j) = neigh_list%neighs(j) + 1
      neigh_list%iNeigh(neigh_list%neighs(j), j) = i
   end do

end function build_neigh_list

end module test_hessian
