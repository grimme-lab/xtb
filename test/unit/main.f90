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
   use testdrive, only : new_testsuite, testsuite_type, select_suite, run_selected, &
      & get_argument, unittest_type, collect_interface, error_type
   use test_atomlist, only : collect_atomlist
   use test_coordinationnumber, only : collect_coordinationnumber
   use test_coulomb, only : collect_coulomb
   use test_dftd3, only : collect_dftd3
   use test_dftd4, only : collect_dftd4
   use test_docking, only : collect_docking
   use test_eeq, only : collect_eeq
   use test_gfn0, only : collect_gfn0
   use test_gfn1, only : collect_gfn1
   use test_gfn2, only : collect_gfn2
   use test_gfnff, only : collect_gfnff
   use test_hessian, only : collect_hessian
   use test_iff, only : collect_iff
   use test_latticepoint, only : collect_latticepoint
   use test_molecule, only : collect_molecule
   use test_pbc_tools, only : collect_pbc_tools
   use test_peeq, only : collect_peeq
   use test_repulsion, only : collect_repulsion
   use test_symmetry, only : collect_symmetry
   use test_tblite, only : collect_tblite
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
      new_testsuite("atomlist", collect_atomlist), &
      new_testsuite("coordinationnumber", collect_coordinationnumber), &
      new_testsuite("coulomb", collect_coulomb), &
      new_testsuite("dftd3", collect_dftd3), &
      new_testsuite("dftd4", collect_dftd4), &
      new_testsuite("docking", collect_docking), &
      new_testsuite("eeq", collect_eeq), &
      new_testsuite("gfn0", collect_gfn0), &
      new_testsuite("gfn1", collect_gfn1), &
      new_testsuite("gfn2", collect_gfn2), &
      new_testsuite("gfnff", collect_gfnff), &
      new_testsuite("hessian", collect_hessian), &
      new_testsuite("iff", collect_iff), &
      new_testsuite("latticepoint", collect_latticepoint), &
      new_testsuite("molecule", collect_molecule), &
      new_testsuite("pbc-tools", collect_pbc_tools), &
      new_testsuite("peeq", collect_peeq), &
      new_testsuite("repulsion", collect_repulsion), &
      new_testsuite("symmetry", collect_symmetry), &
      new_testsuite("tblite", collect_tblite), &
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


contains

   !> Driver for testsuite
   subroutine run_testsuite(collect, unit, stat)

      !> Collect tests
      procedure(collect_interface) :: collect

      !> Unit for IO
      integer, intent(in) :: unit

      !> Number of failed tests
      integer, intent(inout) :: stat

      type(unittest_type), allocatable :: testsuite(:)
      integer :: it

      call collect(testsuite)

      do it = 1, size(testsuite)
         !$omp critical(testdrive_testsuite)
         write(unit, '(1x, 3(1x, a), 1x, "(", i0, "/", i0, ")")') &
            & "Starting", testsuite(it)%name, "...", it, size(testsuite)
         !$omp end critical(testdrive_testsuite)
         call run_unittest(testsuite(it), unit, stat)
      end do

   end subroutine run_testsuite

   !> Run a selected unit test
   subroutine run_unittest(test, unit, stat)

      !> Unit test
      type(unittest_type), intent(in) :: test

      !> Unit for IO
      integer, intent(in) :: unit

      !> Number of failed tests
      integer, intent(inout) :: stat

      type(error_type), allocatable :: error
      character(len=:), allocatable :: message

      call test%test(error)
      if (.not.test_skipped(error) .and. allocated(error) .neqv. test%should_fail) then
         stat = stat + 1
      end if
      call make_output(message, test, error)
      !$omp critical(testdrive_testsuite)
      write(unit, '(a)') message
      !$omp end critical(testdrive_testsuite)
      if (allocated(error)) then
         call clear_error(error)
      end if

   end subroutine run_unittest

   !> Create output message for test (this procedure is pure and therefore cannot launch tests)
   pure subroutine make_output(output, test, error)

      !> Output message for display
      character(len=:), allocatable, intent(out) :: output

      !> Unit test
      type(unittest_type), intent(in) :: test

      !> Error handling
      type(error_type), intent(in), optional :: error

      character(len=:), allocatable :: label
      character(len=*), parameter :: indent = repeat(" ", 7) // repeat(".", 3) // " "

      if (test_skipped(error)) then
         output = indent // test%name // " [SKIPPED]" &
            & // new_line("a") // "  Message: " // error%message
         return
      end if

      if (present(error) .neqv. test%should_fail) then
         if (test%should_fail) then
            label = " [UNEXPECTED PASS]"
         else
            label = " [FAILED]"
         end if
      else
         if (test%should_fail) then
            label = " [EXPECTED FAIL]"
         else
            label = " [PASSED]"
         end if
      end if
      output = indent // test%name // label
      if (present(error)) then
         output = output // new_line("a") // "  Message: " // error%message
      end if
   end subroutine make_output

   pure function test_skipped(error) result(is_skipped)

      !> Error handling
      type(error_type), intent(in), optional :: error

      !> Test was skipped
      logical :: is_skipped

      is_skipped = .false.
      if (present(error)) then
         is_skipped = error%stat == 77
      end if

   end function test_skipped

   !> Clear error type after it has been handled.
   subroutine clear_error(error)

      !> Error handling
      type(error_type), intent(inout) :: error

      if (error%stat /= 0) then
         error%stat = 0
      end if

      if (allocated(error%message)) then
         deallocate(error%message)
      end if

   end subroutine clear_error

end program tester

