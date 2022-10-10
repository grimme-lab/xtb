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

module test_iff
   use testdrive, only : new_unittest, unittest_type, error_type, check_ => check, test_failed
   implicit none
   private

   public :: collect_iff

contains

!> Collect all exported unit tests
subroutine collect_iff(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("iff_sp", test_iff_sp) &
      ]

end subroutine collect_iff


subroutine test_iff_sp(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_systools
   use xtb_type_environment
   use xtb_type_options
   use xtb_type_molecule
   use xtb_type_data
   use xtb_setparam
   use xtb_setmod
   use xtb_docking_param
   use xtb_type_restart, only: TRestart
   use xtb_iff_calculator, only : TIFFCalculator, newIFFCalculator
   use xtb_iff_data, only : TIFFData
   use xtb_iff_iffprepare, only : prepare_IFF
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [8,1,1,8,1,1]
   real(wp),parameter :: xyz(3,nat) = reshape(&
   &[-0.04394313126682_wp,-0.15162589575685_wp,-0.12104386666899_wp, &
   &  1.57021196791551_wp, 0.39852223676357_wp,-0.74238980424452_wp, &
   & -0.52096572705612_wp,-1.61453119575215_wp,-1.08142949421996_wp, &
   & -3.35527313086215_wp, 3.63311912367937_wp,-1.97688072772302_wp, &
   & -4.33218397586988_wp, 4.42399249914350_wp,-0.67195282333972_wp, &
   & -2.30284007732108_wp, 2.38679302539156_wp,-1.14422985849389_wp],&
   & shape(xyz))

   type(TMolecule)     :: mol
   type(TEnvironment)  :: env
   type(TIFFCalculator) :: calc
   type(TIFFData) :: iff_data

   real(wp) :: etot,egap
   real(wp), allocatable :: g(:,:)
   type(TRestart) :: chk
   real(wp) :: sigma(3,3)
   type(scc_results) :: res


   logical  :: exist

   call init(env)
   call init(mol,at,xyz)

   natom_arg = '1-3'
   call prepare_IFF(env, mol, iff_data)
   call env%checkpoint("Could not generate electronic properties for IFF")

   call newIFFCalculator(env, mol, iff_data, calc)

   call env%checkpoint("xtb-IFF parameter setup failed")

   allocate( g(3,mol%n), source = 0.0_wp )
 
   call calc%singlepoint(env, mol, chk, 2, exist, etot, g, sigma, egap, res)
   write(*,*) 'WW', etot

   call check_(error, etot,-0.002628051951328639_wp, thr=thr)

   call mol%deallocate

end subroutine test_iff_sp

end module test_iff
