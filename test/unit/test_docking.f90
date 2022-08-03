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

module test_docking
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_docking

contains

!> Collect all exported unit tests
subroutine collect_docking(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("eth_wat", test_dock_eth_wat)]!, &
!      new_unittest("ellips", test_iff_ellips), &
!      ]

end subroutine collect_docking


subroutine test_dock_eth_wat(error)
   use xtb_type_environment, only: TEnvironment, init
   use xtb_mctc_accuracy, only: wp
   use xtb_type_environment, only: TEnvironment
   use xtb_type_molecule
   use xtb_setparam, only: initrand
   use xtb_mctc_systools
   use xtb_setmod
   use xtb_docking_set_module
   use xtb_docking_param
   use xtb_iff_iffini, only: init_iff
   use xtb_iff_iffprepare, only: precomp
   use xtb_iff_iffenergy
   use xtb_docking_search_nci, only: docking_search
   use xtb_mctc_systools
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-6_wp
   real(wp) :: icoord0(6),icoord(6) 
   integer :: n1 
   integer :: at1(9)
   real(wp) :: xyz1(3,9) 
   integer :: n2 
   integer :: at2(3)
   real(wp) :: xyz2(3,3)
   integer,parameter  :: n = 12
   integer :: at(n)
   real(wp) :: xyz(3,n)
   real(wp) :: molA_e, molB_e

   integer :: nlmo1, nlmo2
   real(wp), allocatable :: rlmo1(:, :)
   real(wp), allocatable :: q1(:)
   real(wp), allocatable :: qdr1(:)
   real(wp), allocatable ::xyzdr1(:, :)
   real(wp), allocatable :: cn1(:)
   real(wp), allocatable :: z1(:)
   real(wp), allocatable :: alp1(:)
   real(wp), allocatable :: qct1(:, :)
   integer, allocatable :: lmo1(:)
   real(wp), allocatable :: rlmo2(:, :)
   real(wp), allocatable :: q2(:)
   real(wp), allocatable :: qdr2(:)
   real(wp), allocatable ::xyzdr2(:, :)
   real(wp), allocatable :: cn2(:)
   real(wp), allocatable :: z2(:)
   real(wp), allocatable :: alp2(:)
   real(wp), allocatable :: qct2(:, :)
   integer, allocatable :: lmo2(:)
   real(wp), allocatable :: c6ab(:, :)
   real(wp), allocatable :: alpab(:, :)
   real(wp), allocatable :: cprob(:)
   real(wp), allocatable :: q(:)
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: alp0(:)
   real(wp), allocatable :: gab1(:, :)
   real(wp), allocatable :: gab2(:, :)
   real(wp), allocatable :: den1(:, :, :)
   real(wp), allocatable :: den2(:, :, :)
   real(wp), allocatable :: qcm1(:)
   real(wp), allocatable :: qcm2(:)
   integer, allocatable :: neigh(:, :)


   logical, parameter :: restart = .false.

   type(TMolecule)     :: molA,molB,comb
   real(wp) :: r(3),e
   integer :: i, j, k
   type(TEnvironment)  :: env
   character(len=:), allocatable :: fnam

   logical  :: exist

   n1 = 9
   at1 = [6,6,1,1,1,1,8,1,1]
   xyz1 = reshape(&
     &[-10.23731657240310_wp, 4.49526140160057_wp,-0.01444418438295_wp, &
     &  -8.92346521409264_wp, 1.93070071439448_wp,-0.06395255870166_wp, &
     & -11.68223026601517_wp, 4.59275319594179_wp,-1.47773402464476_wp, &
     &  -8.87186386250594_wp, 5.99323351202269_wp,-0.34218865473603_wp, &
     & -11.13173905821587_wp, 4.79659470137971_wp, 1.81087156124027_wp, &
     &  -8.00237319881732_wp, 1.64827398274591_wp,-1.90766384850460_wp, &
     & -10.57950417459641_wp,-0.08678859323640_wp, 0.45375344239209_wp, &
     &  -7.47634978432660_wp, 1.83563324786661_wp, 1.40507877109963_wp, &
     & -11.96440476139715_wp,-0.02296384600667_wp,-0.72783557254661_wp],&
     & shape(xyz1))

   n2 = 3
   at2 = [8,1,1]
   xyz2 = reshape(&
     &[-14.55824225787638_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -12.72520790897730_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -15.16924740842229_wp,-0.86379381534203_wp,-0.15285994688912_wp],&
     & shape(xyz2))

   call init(env)
   call init(molA,at1,xyz1)
   call init(molB,at2,xyz2)

   call initrand

   call set_iff_param
   fnam = 'xtblmoinfo'
   set%pr_local = .false.

   call precomp(env, molA, molA_e, 1)
   call check_(error, molA_e,-11.3943358674_wp, thr=thr)

   call rd0(1, trim(fnam), n1, nlmo1)
   allocate (lmo1(10*n1), source=0)
   allocate (rlmo1(4, 10*n1), q1(n1),&
   &cn1(n1), alp1(n1), qct1(n1, 2), qdr1(n1), xyzdr1(3, n1),&
   &z1(n1), den1(2, 4, n1), gab1(n1, n1), qcm1(n1), cprob(n1), source=0.0_wp)
   call rd(trim(fnam), 1, n1, xyz1, at1, nlmo1, lmo1, rlmo1, q1, qct1)
   call delete_file(trim(fnam))

   call precomp(env, molB, molB_e, 2)

   call rd0(2, trim(fnam), n2, nlmo2)
   allocate (lmo2(10*n2), source=0)
   allocate (rlmo2(4, 10*n2), q2(n2),&
   & cn2(n2), alp2(n2), qct2(n2, 2), qdr2(n2), xyzdr2(3, n2),&
   & z2(n2), den2(2, 4, n2), gab2(n2, n2), qcm2(n2), source=0.0_wp)
   call rd(trim(fnam), 2, n2, xyz2, at2, nlmo2, lmo2, rlmo2, q2, qct2)
   call delete_file(trim(fnam))

   call env%checkpoint("LMO computation failed")

   call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
   do i = 1, 3
      molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
   end do
   call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
   do i = 1, 3
      molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
   end do

   allocate (neigh(0:n, n), source=0)
   allocate (q(n), c6ab(n, n), alp0(n), cn(n),&
   &        alpab(n2, n1), source=0.0_wp)

   call init_iff(env, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab, z1, z2,&
      &          cprob, nlmo1, nlmo2, lmo1, lmo2,&
      &          qdr1, qdr2, rlmo1, rlmo2, cn1, cn2, alp1, alp2, alpab,&
      &          den1, den2, gab1, gab2, qcm1, qcm2, n, at, xyz, q, icoord, icoord0,&
      &          .false.)

   call check_(error, icoord(1) ,-4.5265_wp, thr=1.0e-3_wp)


   call env%checkpoint("Initializing xtb-IFF failed.")

   call iff_e(env, n, n1, n2, at1, at2, neigh, xyz1, xyz2, q1, q2, c6ab, z1, z2,&
                    & nlmo1, nlmo2, lmo1, lmo2, rlmo1, rlmo2,&
                    & qdr1, qdr2, cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
                    & den1, den2, gab1, gab2,&
                    & .true., 0, e, icoord)
   call check_(error, e,0.07745330953243718_wp, thr=thr)

   call env%checkpoint("xtb-IFF sp computation failed")

   call set_optlvl(env) !Sets the optimization level for search in global parameters
!   call docking_search(env, molA, molB, n, n1, n2, at1, at2, neigh, xyz1,&
!                 & xyz2, q1, q2, c6ab, z1, z2,&
!                 &nlmo1, nlmo2, lmo1, lmo2, rlmo1, rlmo2,&
!                 &qdr1, qdr2, cn1, cn2, alp1, alp2, alpab, qct1, qct2,&
!                 &den1, den2, gab1, gab2, molA_e, molB_e,&
!                 &cprob, e, icoord, comb)

!   call env%checkpoint("Docking algorithm failed")

!   call check_(error, calc%topo%nbond,6)
!
!   call check_(error, res_gff%e_total,-0.76480130317838_wp, thr=thr)

   call molA%deallocate
   call molB%deallocate

end subroutine test_dock_eth_wat

end module test_docking
