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

module test_pbc_tools
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   implicit none
   private

   public :: collect_pbc_tools

contains

!> Collect all exported unit tests
subroutine collect_pbc_tools(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("cutoff", test_pbc_tools_cutoff), &
      new_unittest("convert", test_pbc_tools_convert) &
      ]

end subroutine collect_pbc_tools


subroutine test_pbc_tools_convert(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_pbc_tools
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: cellpar(6) = &
      [9.09903131_wp, 9.09903131_wp, 30.46049560_wp, 90.0_wp, 90.0_wp, 120.0_wp]
   real(wp),parameter :: lattice(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice))
   real(wp) :: dlat(3,3),rlat(3,3),volume,cpar(6),center(3)

   call cell_to_dlat(cellpar,dlat)
   call check(error, cellpar(1),dlat(1,1), thr=thr)
   call check(error, dlat(1,3),-13.648544412579_wp, thr=thr)
   call check(error, dlat(3,3), 26.878966813429_wp, thr=thr)

   call cell_to_rlat(cellpar,rlat)
   call check(error, rlat(1,1), 0.69053343077018_wp, thr=thr)
   call check(error, rlat(2,1),-0.96832302603374_wp, thr=thr)
   call check(error, rlat(3,2), 0.19327597283376_wp, thr=thr)

   volume = cell_to_dvol(cellpar)
   call check(error, volume,1292.0766773144_wp, thr=thr)

   call dlat_to_cell(lattice,cpar)
   call check(error, cpar(1),cpar(2), thr=thr)
   call check(error, cpar(3),7.2993338808807_wp, thr=thr)
   call check(error, cpar(4),1.0471975511966_wp, thr=thr)

   call dlat_to_rlat(lattice,rlat)
   call check(error, rlat(1,1), 1.0542467445047_wp, thr=thr)
   call check(error, rlat(1,3),-.35141558150155_wp, thr=thr)
   call check(error, rlat(2,2), .99395336277746_wp, thr=thr)

   volume = dlat_to_dvol(lattice)
   call check(error, volume, 275.00126402469_wp, thr=thr)

   center = get_center_dlat(lattice)
   call check(error, center(1),2.9799405783945_wp, thr=thr)
   call check(error, center(2),4.2142723810314_wp, thr=thr)

end subroutine test_pbc_tools_convert

subroutine test_pbc_tools_cutoff(error)
   use xtb_mctc_accuracy, only : wp
   use xtb_pbc_tools
   use xtb_pbc
   type(error_type), allocatable, intent(out) :: error
   real(wp),parameter :: thr = 1.0e-10_wp
   real(wp),parameter :: lattice_1(3,3) = reshape(&
      &[9.0990313100000_wp,      0.0000000000000_wp,      0.0000000000000_wp,    &
      & 7.4082581428274_wp,      5.2829993440840_wp,      0.0000000000000_wp,    &
      &-13.648544412579_wp,     -4.3680854682659_wp,      26.878966813429_wp],   &
      & shape(lattice_1))
   real(wp),parameter :: lattice_2(3,3) = reshape(&
      &[5.9598811567890_wp,      2.1071361905157_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      6.3214085715472_wp,      3.6496669404404_wp,    &
      & 0.0000000000000_wp,      0.0000000000000_wp,      7.2993338808807_wp],   &
      & shape(lattice_2))
   real(wp),parameter :: rthr_1 =  40.0_wp ** 2
   real(wp),parameter :: rthr_2 =  70.0_wp ** 2
   real(wp),parameter :: rthr_3 = 100.0_wp ** 2
   integer :: rep(3)

   call get_realspace_cutoff(lattice_1,rthr_1,rep)
   call check(error, rep(1), 8)
   call check(error, rep(3), 2)

   call get_realspace_cutoff(lattice_1,rthr_2,rep)
   call check(error, rep(2),14)
   call check(error, rep(3), 3)

   call get_realspace_cutoff(lattice_1,rthr_3,rep)
   call check(error, rep(1),20)
   call check(error, rep(2),20)
   call check(error, rep(3), 4)

   call get_realspace_cutoff(lattice_2,rthr_1,rep)
   call check(error, rep(1),rep(2))
   call check(error, rep(3), 7)

   call get_realspace_cutoff(lattice_2,rthr_2,rep)
   call check(error, rep(1),rep(3))
   call check(error, rep(2),12)

   call get_realspace_cutoff(lattice_2,rthr_3,rep)
   call check(error, rep(3),rep(2))
   call check(error, rep(1),17)

end subroutine test_pbc_tools_cutoff

end module test_pbc_tools
