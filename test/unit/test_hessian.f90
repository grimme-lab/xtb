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
      & -1.013452580143E+00_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -2.398345095286E-15_wp, -5.801057146788E-01_wp,  4.440892098501E-10_wp, &
      &  1.886892359658E-11_wp,  1.204196255680E-09_wp, -7.199230668276E-01_wp, &
      &  5.067262901419E-01_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -5.168463665021E-11_wp,  2.908278878443E-01_wp,  7.276368485520E-02_wp, &
      &  0.000000000000E+00_wp,  1.594757638403E-01_wp,  3.599430671297E-01_wp, &
      &  5.067262899828E-01_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  1.887034647398E-11_wp,  2.908278878482E-01_wp, -7.276368718667E-02_wp, &
      &  3.527704322735E-11_wp, -1.594757620026E-01_wp,  3.599430659085E-01_wp],&
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      &  9.642724315145E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -4.821349241637E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -4.821375073506E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -1.193695511930E-16_wp,  6.543286738602E-01_wp, -1.383245460448E-09_wp, &
      &  2.351384649776E-16_wp, -3.271643375731E-01_wp,  2.415080663264E-01_wp, &
      & -1.157689137746E-16_wp, -3.271643362758E-01_wp, -2.415080649501E-01_wp, &
      & -7.979405595065E-11_wp, -8.671336189204E-10_wp,  3.969567976825E-01_wp, &
      &  9.617770183336E-11_wp,  1.933148104688E-01_wp, -1.984783990326E-01_wp, &
      & -1.638364588271E-11_wp, -1.933148095978E-01_wp, -1.984783986206E-01_wp, &
      & -1.362023380178E-03_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  4.418957886947E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  1.317833801308E-03_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -4.323271238214E-11_wp, -3.270692709850E-01_wp,  1.934196079808E-01_wp, &
      &  4.723871392038E-11_wp,  3.486669886682E-01_wp, -2.174287435611E-01_wp, &
      & -4.006001538234E-12_wp, -2.159771769511E-02_wp,  2.400913555821E-02_wp, &
      &  0.000000000000E+00_wp,  2.414357451224E-01_wp, -1.984817451162E-01_wp, &
      &  0.000000000000E+00_wp, -2.173769120375E-01_wp,  1.883863080853E-01_wp, &
      &  0.000000000000E+00_wp, -2.405883308764E-02_wp,  1.009543702395E-02_wp, &
      & -1.362023712098E-03_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  1.317833873830E-03_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  4.418983826786E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -7.979413971062E-11_wp, -3.270692701994E-01_wp, -1.934196070410E-01_wp, &
      &  9.617754838868E-11_wp, -2.159771765280E-02_wp, -2.400913556735E-02_wp, &
      & -1.638340867806E-11_wp,  3.486669878673E-01_wp,  2.174287426208E-01_wp, &
      & -1.828072800336E-11_wp, -2.414357440668E-01_wp, -1.984817441066E-01_wp, &
      &  2.446948471031E-11_wp,  2.405883338662E-02_wp,  1.009543671832E-02_wp, &
      & -6.188756706950E-12_wp,  2.173769106779E-01_wp,  1.883863073862E-01_wp],&
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
      & -8.113137336986E-01_wp,  5.356654180818E-10_wp,  0.000000000000E+00_wp, &
      & -4.581749290942E-11_wp, -4.121161465351E-01_wp, -4.996003610813E-10_wp, &
      &  5.570500822590E-16_wp, -7.377395210573E-10_wp, -4.838787099892E-01_wp, &
      &  4.056480094735E-01_wp, -2.200999955020E-09_wp,  1.665334536938E-09_wp, &
      & -3.988074658262E-11_wp,  2.055752062158E-01_wp,  9.503822412382E-02_wp, &
      & -3.969738621316E-11_wp,  1.480698268022E-01_wp,  2.434252440731E-01_wp, &
      &  4.056480099166E-01_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -6.052638931112E-11_wp,  2.055752051281E-01_wp, -9.503822262502E-02_wp, &
      &  1.674180515478E-11_wp, -1.480698285344E-01_wp,  2.434252415195E-01_wp],&
      &  shape(dipgrad_ref))
   real(wp), parameter :: hessian_ref(3*nat, 3*nat) = reshape([ &
      & -1.939578553890E-05_wp,  6.507635536734E-12_wp,  2.039452049081E-10_wp, &
      &  9.697594542869E-06_wp,  5.901447950110E-11_wp, -1.346714623499E-10_wp, &
      &  9.698190996037E-06_wp, -6.773553072603E-11_wp, -7.145569943037E-11_wp, &
      & -1.648255831783E-11_wp,  6.121988089582E-01_wp,  2.181719702902E-09_wp, &
      &  1.784659444500E-11_wp, -3.060994033735E-01_wp,  2.258521490017E-01_wp, &
      & -1.364036127169E-12_wp, -3.060994055863E-01_wp, -2.258521511869E-01_wp, &
      &  1.330615415687E-14_wp,  2.485255512685E-09_wp,  3.998109489554E-01_wp, &
      & -1.478661148333E-14_wp,  1.841542843085E-01_wp, -1.999054737643E-01_wp, &
      &  1.480457326447E-15_wp, -1.841542867939E-01_wp, -1.999054751878E-01_wp, &
      & -1.071178290455E-02_wp,  6.630362607819E-11_wp, -1.553932763715E-10_wp, &
      &  5.182610736054E-05_wp, -3.686287386451E-11_wp,  1.200855549981E-10_wp, &
      &  1.065995679719E-02_wp, -2.949707535518E-11_wp,  3.970551643549E-11_wp, &
      & -2.623112661213E-11_wp, -3.057148399989E-01_wp,  1.833996998222E-01_wp, &
      &  2.414732886753E-11_wp,  3.392228251096E-01_wp, -2.044438914287E-01_wp, &
      &  2.083797744602E-12_wp, -3.350798510988E-02_wp,  2.104419160463E-02_wp, &
      & -2.832253701727E-11_wp,  2.256129456682E-01_wp, -1.999899521472E-01_wp, &
      &  2.493748844816E-11_wp, -2.048926309476E-01_wp,  1.834977735563E-01_wp, &
      &  3.385048569112E-12_wp, -2.072031472126E-02_wp,  1.649217859268E-02_wp, &
      & -1.071178231982E-02_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  1.065995680557E-02_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      &  5.182551425220E-05_wp,  0.000000000000E+00_wp,  0.000000000000E+00_wp, &
      & -4.798341604929E-11_wp, -3.057148423748E-01_wp, -1.833997022848E-01_wp, &
      &  5.114468347305E-11_wp, -3.350798504577E-02_wp, -2.104419155262E-02_wp, &
      & -3.161267423758E-12_wp,  3.392228274191E-01_wp,  2.044438938336E-01_wp, &
      &  2.174439509959E-11_wp, -2.256129480723E-01_wp, -1.999899531846E-01_wp, &
      & -1.789465301866E-11_wp,  2.072031483689E-02_wp,  1.649217845321E-02_wp, &
      & -3.849742080923E-12_wp,  2.048926332339E-01_wp,  1.834977747311E-01_wp],&
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

end module test_hessian
