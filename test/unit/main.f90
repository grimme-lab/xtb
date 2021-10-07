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

!> Driver for unit testing
program tester
   use, intrinsic :: iso_fortran_env, only : error_unit
   use testdrive, only : run_testsuite, new_testsuite, testsuite_type, &
      & select_suite, run_selected, get_argument
   use test_coulomb, only : collect_coulomb
   use test_dftd3, only : collect_dftd3
   use test_dftd4, only : collect_dftd4
   use test_eeq, only : collect_eeq
   use test_geometry_reader, only : collect_geometry_reader
   use test_hessian, only : collect_hessian
   use test_molecule, only : collect_molecule
   use test_thermo, only : collect_thermo
   use test_wsc, only : collect_wsc
   implicit none
   integer :: stat, is
   character(len=:), allocatable :: suite_name, test_name
   type(testsuite_type), allocatable :: testsuites(:)
   character(len=*), parameter :: fmt = '("#", *(1x, a))'

   call mctc_init('test',10,.true.)

   stat = 0

   testsuites = [ &
      new_testsuite("coulomb", collect_coulomb), &
      new_testsuite("dftd3", collect_dftd3), &
      new_testsuite("dftd4", collect_dftd4), &
      new_testsuite("eeq", collect_eeq), &
      new_testsuite("geometry-reader", collect_geometry_reader), &
      new_testsuite("hessian", collect_hessian), &
      new_testsuite("molecule", collect_molecule), &
      new_testsuite("thermo", collect_thermo), &
      new_testsuite("wsc", collect_wsc) &
      ]

   call get_argument(1, suite_name)
   call get_argument(2, test_name)

   if (allocated(suite_name)) then
      is = select_suite(testsuites, suite_name)
      if (is > 0 .and. is <= size(testsuites)) then
         if (allocated(test_name)) then
            write(error_unit, fmt) "Suite:", testsuites(is)%name
            call run_selected(testsuites(is)%collect, test_name, error_unit, stat)
            if (stat < 0) then
               error stop 1
            end if
         else
            write(error_unit, fmt) "Testing:", testsuites(is)%name
            call run_testsuite(testsuites(is)%collect, error_unit, stat)
         end if
      else
         write(error_unit, fmt) "Available testsuites"
         do is = 1, size(testsuites)
            write(error_unit, fmt) "-", testsuites(is)%name
         end do
         error stop 1
      end if
   else
      do is = 1, size(testsuites)
         write(error_unit, fmt) "Testing:", testsuites(is)%name
         call run_testsuite(testsuites(is)%collect, error_unit, stat)
      end do
   end if

   if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop 1
   end if


end program tester

