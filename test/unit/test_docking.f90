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
      new_unittest("dock_gfn2_eth_wat", test_dock_eth_wat_gfn2), &
      new_unittest("dock_gfn2_wat_wat_wall", test_dock_wat_wat_gfn2_wall), &
      new_unittest("dock_gfn2_wat_wat_attpot", test_dock_wat_wat_gfn2_attpot), &
      new_unittest("dock_gfnff_wat_wat", test_dock_wat_wat_gfnff) &
      ]

end subroutine collect_docking


subroutine test_dock_eth_wat_gfn2(error)
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
   use xtb_iff_data, only: TIFFData
   type(error_type), allocatable, intent(out) :: error
   !> All important variables stored here
   type(TIFFData) :: iff_data
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
   call init(molA, at1, xyz1)
   call init(molB, at2, xyz2)
   call iff_data%allocateIFFData(molA%n, molB%n)

   call initrand

   call set_iff_param
   fnam = 'xtblmoinfo'
   set%pr_local = .false.

   call precomp(env, iff_data, molA, molA_e, 1)
   call check_(error, molA_e,-11.39433586739173_wp, thr=thr)
   call precomp(env, iff_data, molB, molB_e, 2)
   call check_(error, molB_e,-5.070201841808753_wp, thr=thr)

   call env%checkpoint("LMO computation failed")

   call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
   do i = 1, 3
      molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
   end do
   call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
   do i = 1, 3
      molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
   end do

   call init_iff(env, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
      &          iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, &
      &          iff_data%q2, iff_data%c6ab, iff_data%z1, iff_data%z2,&
      &          iff_data%cprob, iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2,&
      &          iff_data%qdr1, iff_data%qdr2, iff_data%rlmo1, iff_data%rlmo2,&
      &          iff_data%cn1, iff_data%cn2, iff_data%alp1, iff_data%alp2, iff_data%alpab,&
      &          iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, iff_data%qcm1,&
      &          iff_data%qcm2, iff_data%n, iff_data%at, iff_data%xyz, iff_data%q, icoord, icoord0,&
      &          .false.)

   call check_(error, icoord(1) ,-4.5265_wp, thr=1.0e-3_wp)

   call env%checkpoint("Initializing xtb-IFF failed.")

   call iff_e(env, iff_data%n, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
             &  iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, iff_data%q2,&
             & iff_data%c6ab, iff_data%z1, iff_data%z2,&
             & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2, &
             & iff_data%rlmo1, iff_data%rlmo2,&
             & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1, &
             & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2, &
             & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, &
             & set%verbose, 0, e, icoord)

   call check_(error, e,.7745330953199690E-01_wp, thr=thr)

   call env%checkpoint("xtb-IFF sp computation failed")

   optlvl="gfn2"
   call set_optlvl(env) !Sets the optimization level for search in global parameters

   !Settings for fast
   maxparent = 30
   maxgen = 3    
   mxcma = 250   
   stepr = 4.0   
   stepa = 60    
   n_opt = 2

   call docking_search(env, molA, molB, iff_data%n, iff_data%n1, iff_data%n2,&
                 & iff_data%at1, iff_data%at2, iff_data%neigh, iff_data%xyz1,&
                 & iff_data%xyz2, iff_data%q1, iff_data%q2, iff_data%c6ab,&
                 & iff_data%z1, iff_data%z2,&
                 & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1,&
                 & iff_data%lmo2, iff_data%rlmo1, iff_data%rlmo2,&
                 & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1,&
                 & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2,&
                 & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, molA_e, molB_e,&
                 & iff_data%cprob, e, icoord, comb)

   call env%checkpoint("Docking algorithm failed")

   call molA%deallocate
   call molB%deallocate

end subroutine test_dock_eth_wat_gfn2


subroutine test_dock_wat_wat_gfn2_wall(error)
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
   use xtb_iff_data, only: TIFFData
   use xtb_sphereparam, only: init_walls, maxwalls, set_sphere_radius,&
           & spherepot_type, p_type_polynomial, clear_walls

   type(error_type), allocatable, intent(out) :: error
   !> All important variables stored here
   type(TIFFData) :: iff_data
   real(wp),parameter :: thr = 1.0e-6_wp
   real(wp) :: icoord0(6),icoord(6) 
   integer :: n1 
   integer :: at1(3)
   real(wp) :: xyz1(3,3) 
   integer :: n2 
   integer :: at2(3)
   real(wp) :: xyz2(3,3)
   integer,parameter  :: n = 6
   integer :: at(n)
   real(wp) :: xyz(3,n)
   real(wp) :: molA_e, molB_e

   logical, parameter :: restart = .false.

   type(TMolecule)     :: molA,molB,comb
   real(wp) :: r(3),e
   integer :: i, j, k
   type(TEnvironment)  :: env
   character(len=:), allocatable :: fnam

   logical  :: exist

   n2 = 3
   at2 = [8,1,1]
   xyz2 = reshape(&
     &[-14.55824225787638_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -12.72520790897730_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -15.16924740842229_wp,-0.86379381534203_wp,-0.15285994688912_wp],&
     & shape(xyz2))

   n1=n2
   at1=at2
   xyz1=xyz2

   qcg=.true.
   call set_gbsa(env, 'solvent', "h2o")
   call init(env)

   !Molecular stuff:
   call init(molA, at1, xyz1)
   call init(molB, at2, xyz2)
   call iff_data%allocateIFFData(molA%n, molB%n)

   call initrand

   call set_iff_param
   fnam = 'xtblmoinfo'
   set%pr_local = .false.

   call precomp(env, iff_data, molA, molA_e, 1)
!   call check_(error, molA_e,-5.084810260314862_wp, thr=thr)
   call precomp(env, iff_data, molB, molB_e, 2)
   call check_(error, molB_e,-5.08481026031486_wp, thr=thr)

   call env%checkpoint("LMO computation failed")

   call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
   do i = 1, 3
      molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
   end do
   call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
   do i = 1, 3
      molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
   end do

   call init_iff(env, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
      &          iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, &
      &          iff_data%q2, iff_data%c6ab, iff_data%z1, iff_data%z2,&
      &          iff_data%cprob, iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2,&
      &          iff_data%qdr1, iff_data%qdr2, iff_data%rlmo1, iff_data%rlmo2,&
      &          iff_data%cn1, iff_data%cn2, iff_data%alp1, iff_data%alp2, iff_data%alpab,&
      &          iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, iff_data%qcm1,&
      &          iff_data%qcm2, iff_data%n, iff_data%at, iff_data%xyz, iff_data%q, icoord, icoord0,&
      &          .false.)

   call check_(error, icoord(1) ,0.0000_wp, thr=1.0e-3_wp)

   !Wall potential stuff:
   maxwalls=3
   call init_walls
   spherepot_type = p_type_polynomial
   call set_sphere_radius([6.7,4.7,3.8],[0.0,0.0,0.0],3,[1,2,3]) !"solute"
   call set_sphere_radius([6.9,5.7,5.2],[0.0,0.0,0.0]) !"all"

   call env%checkpoint("Initializing xtb-IFF failed.")

   call iff_e(env, iff_data%n, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
             &  iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, iff_data%q2,&
             & iff_data%c6ab, iff_data%z1, iff_data%z2,&
             & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2, &
             & iff_data%rlmo1, iff_data%rlmo2,&
             & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1, &
             & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2, &
             & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, &
             & set%verbose, 0, e, icoord)

   call env%checkpoint("xtb-IFF sp computation failed")

   optlvl="gfn2"
   call set_optlvl(env) !Sets the optimization level for search in global parameters

   !Settings for fast
   maxparent = 30
   maxgen = 3    
   mxcma = 250   
   stepr = 4.0   
   stepa = 60    
   n_opt = 2

   call docking_search(env, molA, molB, iff_data%n, iff_data%n1, iff_data%n2,&
                 & iff_data%at1, iff_data%at2, iff_data%neigh, iff_data%xyz1,&
                 & iff_data%xyz2, iff_data%q1, iff_data%q2, iff_data%c6ab,&
                 & iff_data%z1, iff_data%z2,&
                 & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1,&
                 & iff_data%lmo2, iff_data%rlmo1, iff_data%rlmo2,&
                 & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1,&
                 & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2,&
                 & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, molA_e, molB_e,&
                 & iff_data%cprob, e, icoord, comb)
   call env%checkpoint("Docking algorithm failed")

   call molA%deallocate
   call molB%deallocate
   call clear_walls
   maxwalls=0

end subroutine test_dock_wat_wat_gfn2_wall

subroutine test_dock_wat_wat_gfn2_attpot(error)
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
   use xtb_iff_data, only: TIFFData
   use xtb_sphereparam, only: init_walls, maxwalls, set_sphere_radius,&
           & spherepot_type, p_type_polynomial, clear_walls

   type(error_type), allocatable, intent(out) :: error
   !> All important variables stored here
   type(TIFFData) :: iff_data
   real(wp),parameter :: thr = 1.0e-6_wp
   real(wp) :: icoord0(6),icoord(6) 
   integer :: n1 
   integer :: at1(3)
   real(wp) :: xyz1(3,3) 
   integer :: n2 
   integer :: at2(3)
   real(wp) :: xyz2(3,3)
   integer,parameter  :: n = 6
   integer :: at(n)
   real(wp) :: xyz(3,n)
   real(wp) :: molA_e, molB_e

   logical, parameter :: restart = .false.

   type(TMolecule)     :: molA,molB,comb
   real(wp) :: r(3),e
   integer :: i, j, k
   type(TEnvironment)  :: env
   character(len=:), allocatable :: fnam

   logical  :: exist

   n2 = 3
   at2 = [8,1,1]
   xyz2 = reshape(&
     &[-14.55824225787638_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -12.72520790897730_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -15.16924740842229_wp,-0.86379381534203_wp,-0.15285994688912_wp],&
     & shape(xyz2))

   n1=n2
   at1=at2
   xyz1=xyz2

   qcg=.true.
   call set_gbsa(env, 'solvent', "h2o")
   call init(env)

   !Molecular stuff:
   call init(molA, at1, xyz1)
   call init(molB, at2, xyz2)
   call iff_data%allocateIFFData(molA%n, molB%n)

   call initrand

   call set_iff_param
   fnam = 'xtblmoinfo'
   set%pr_local = .false.

   call precomp(env, iff_data, molA, molA_e, 1)
!   call check_(error, molA_e,-5.084810260314862_wp, thr=thr)
   call precomp(env, iff_data, molB, molB_e, 2)
   call check_(error, molB_e,-5.08481026031486_wp, thr=thr)

   call env%checkpoint("LMO computation failed")

   call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
   do i = 1, 3
      molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
   end do
   call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
   do i = 1, 3
      molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
   end do

   call init_iff(env, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
      &          iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, &
      &          iff_data%q2, iff_data%c6ab, iff_data%z1, iff_data%z2,&
      &          iff_data%cprob, iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2,&
      &          iff_data%qdr1, iff_data%qdr2, iff_data%rlmo1, iff_data%rlmo2,&
      &          iff_data%cn1, iff_data%cn2, iff_data%alp1, iff_data%alp2, iff_data%alpab,&
      &          iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, iff_data%qcm1,&
      &          iff_data%qcm2, iff_data%n, iff_data%at, iff_data%xyz, iff_data%q, icoord, icoord0,&
      &          .false.)

   call check_(error, icoord(1) ,0.0000_wp, thr=1.0e-3_wp)

   directed_type = p_atom_qcg
   allocate(directedset%expo(2), source=0.0_wp)
   allocate(directedset%val(1), source=0.0_wp)
   directedset%val(1)=0.011928
   directedset%expo(1)=0.01
   directedset%expo(2)=2.000
   directedset%n=3
   directedset%atoms=[1,2,3]

   call env%checkpoint("Initializing xtb-IFF failed.")

   call iff_e(env, iff_data%n, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
             &  iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, iff_data%q2,&
             & iff_data%c6ab, iff_data%z1, iff_data%z2,&
             & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2, &
             & iff_data%rlmo1, iff_data%rlmo2,&
             & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1, &
             & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2, &
             & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, &
             & set%verbose, 0, e, icoord)

   call env%checkpoint("xtb-IFF sp computation failed")

   optlvl="gfn2"
   call set_optlvl(env) !Sets the optimization level for search in global parameters

   !Settings for fast
   maxparent = 30
   maxgen = 3    
   mxcma = 250   
   stepr = 4.0   
   stepa = 60    
   n_opt = 2

   call docking_search(env, molA, molB, iff_data%n, iff_data%n1, iff_data%n2,&
                 & iff_data%at1, iff_data%at2, iff_data%neigh, iff_data%xyz1,&
                 & iff_data%xyz2, iff_data%q1, iff_data%q2, iff_data%c6ab,&
                 & iff_data%z1, iff_data%z2,&
                 & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1,&
                 & iff_data%lmo2, iff_data%rlmo1, iff_data%rlmo2,&
                 & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1,&
                 & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2,&
                 & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, molA_e, molB_e,&
                 & iff_data%cprob, e, icoord, comb)
   call env%checkpoint("Docking algorithm failed")

   call molA%deallocate
   call molB%deallocate
   call directedset%deallocate

end subroutine test_dock_wat_wat_gfn2_attpot


subroutine test_dock_wat_wat_gfnff(error)
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
   use xtb_iff_data, only: TIFFData
   use xtb_mctc_global, only : persistentEnv

   type(error_type), allocatable, intent(out) :: error
   !> All important variables stored here
   type(TIFFData) :: iff_data
   real(wp),parameter :: thr = 1.0e-6_wp
   real(wp) :: icoord0(6),icoord(6) 
   integer :: n1 
   integer :: at1(3)
   real(wp) :: xyz1(3,3) 
   integer :: n2 
   integer :: at2(3)
   real(wp) :: xyz2(3,3)
   integer,parameter  :: n = 6
   integer :: at(n)
   real(wp) :: xyz(3,n)
   real(wp) :: molA_e, molB_e

   logical, parameter :: restart = .false.

   type(TMolecule)     :: molA,molB,comb
   real(wp) :: r(3),e
   integer :: i, j, k
   type(TEnvironment)  :: env
   character(len=:), allocatable :: fnam

   logical  :: exist

   n2 = 3
   at2 = [8,1,1]
   xyz2 = reshape(&
     &[-14.55824225787638_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -12.72520790897730_wp, 0.85763330814882_wp, 0.00000000000000_wp, &
     & -15.16924740842229_wp,-0.86379381534203_wp,-0.15285994688912_wp],&
     & shape(xyz2))

   n1=n2
   at1=at2
   xyz1=xyz2

   call set_gbsa(env, 'solvent', "h2o")
   !GFN-FF 
   call set_exttyp('ff')

!   deallocate(persistentEnv%whoami)
!   deallocate(persistentEnv%home)
!   deallocate(persistentEnv%path)
!   deallocate(persistentEnv%xtbpath)
!   deallocate(persistentEnv%xtbhome)
!   deallocate(persistentEnv%io%log)
!   persistentEnv%io%count=0


   call init(env)
   call init(molA, at1, xyz1)
   call init(molB, at2, xyz2)
   call iff_data%allocateIFFData(molA%n, molB%n)

   call initrand

   call set_iff_param
   fnam = 'xtblmoinfo'
   set%pr_local = .false.

   call precomp(env, iff_data, molA, molA_e, 1)
!   call check_(error, molA_e,-5.084810260314862_wp, thr=thr)
   call precomp(env, iff_data, molB, molB_e, 2)
   call check_(error, molB_e,-5.08481026031486_wp, thr=thr)

   call env%checkpoint("LMO computation failed")

   call cmadock(molA%n, molA%n, molA%at, molA%xyz, r)
   do i = 1, 3
      molA%xyz(i, 1:molA%n) = molA%xyz(i, 1:molA%n) - r(i)
   end do
   call cmadock(molB%n, molB%n, molB%at, molB%xyz, r)
   do i = 1, 3
      molB%xyz(i, 1:molB%n) = molB%xyz(i, 1:molB%n) - r(i)
   end do

   call init_iff(env, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
      &          iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, &
      &          iff_data%q2, iff_data%c6ab, iff_data%z1, iff_data%z2,&
      &          iff_data%cprob, iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2,&
      &          iff_data%qdr1, iff_data%qdr2, iff_data%rlmo1, iff_data%rlmo2,&
      &          iff_data%cn1, iff_data%cn2, iff_data%alp1, iff_data%alp2, iff_data%alpab,&
      &          iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, iff_data%qcm1,&
      &          iff_data%qcm2, iff_data%n, iff_data%at, iff_data%xyz, iff_data%q, icoord, icoord0,&
      &          .false.)

   call check_(error, icoord(1) ,0.0000_wp, thr=1.0e-3_wp)

   call env%checkpoint("Initializing xtb-IFF failed.")

   call iff_e(env, iff_data%n, iff_data%n1, iff_data%n2, iff_data%at1, iff_data%at2,&
             &  iff_data%neigh, iff_data%xyz1, iff_data%xyz2, iff_data%q1, iff_data%q2,&
             & iff_data%c6ab, iff_data%z1, iff_data%z2,&
             & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1, iff_data%lmo2, &
             & iff_data%rlmo1, iff_data%rlmo2,&
             & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1, &
             & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2, &
             & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, &
             & set%verbose, 0, e, icoord)

   call env%checkpoint("xtb-IFF sp computation failed")

   optlvl = "gfnff"
   call set_optlvl(env) !Sets the optimization level for search in global parameters

   !Settings for fast
   maxparent = 30
   maxgen = 3    
   mxcma = 250   
   stepr = 4.0   
   stepa = 60    
   n_opt = 2

   call docking_search(env, molA, molB, iff_data%n, iff_data%n1, iff_data%n2,&
                 & iff_data%at1, iff_data%at2, iff_data%neigh, iff_data%xyz1,&
                 & iff_data%xyz2, iff_data%q1, iff_data%q2, iff_data%c6ab,&
                 & iff_data%z1, iff_data%z2,&
                 & iff_data%nlmo1, iff_data%nlmo2, iff_data%lmo1,&
                 & iff_data%lmo2, iff_data%rlmo1, iff_data%rlmo2,&
                 & iff_data%qdr1, iff_data%qdr2, iff_data%cn1, iff_data%cn2, iff_data%alp1,&
                 & iff_data%alp2, iff_data%alpab, iff_data%qct1, iff_data%qct2,&
                 & iff_data%den1, iff_data%den2, iff_data%gab1, iff_data%gab2, molA_e, molB_e,&
                 & iff_data%cprob, e, icoord, comb)
   call env%checkpoint("Docking algorithm failed")

   call molA%deallocate
   call molB%deallocate

end subroutine test_dock_wat_wat_gfnff

end module test_docking
