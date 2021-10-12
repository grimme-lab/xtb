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

module test_atomlist
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stderr
   use xtb_type_atomlist
   implicit none
   private

   public :: collect_atomlist

   integer, parameter :: atoms(9) = [3,1,1,5,8,1,1,2,5]
   logical, parameter :: lpar(9) = [.true., .false., .true., .true., .true., &
      &                             .false., .false., .true., .false.]
   integer, parameter :: ipar(5) = [1, 3, 4, 5, 8]
   character(len=*), parameter :: cpar = '1,3-5,8'

contains

!> Collect all exported unit tests
subroutine collect_atomlist(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("defaults", test_atomlist_defaults), &
      new_unittest("constr1", test_atomlist_constr1), &
      new_unittest("constr2", test_atomlist_constr2), &
      new_unittest("constr3", test_atomlist_constr3), &
      new_unittest("constr4", test_atomlist_constr4), &
      new_unittest("manip1", test_atomlist_manip1), &
      new_unittest("manip2", test_atomlist_manip2) &
      ]

end subroutine collect_atomlist


subroutine test_atomlist_defaults(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   call check(error, atl%get_truth() .eqv. .true.)
   call atl%switch_truth
   call check(error, atl%get_truth() .eqv. .false.)
   call check(error, size(atl), 0)

end subroutine test_atomlist_defaults

subroutine test_atomlist_constr1(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   atl = TAtomList(list=lpar, truth=.true.)
   call check(error, size(atl), 9)
   call check(error, len(atl), 5)

end subroutine test_atomlist_constr1

subroutine test_atomlist_constr2(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   atl = TAtomList(list=lpar, truth=.false.)
   call check(error, size(atl), 9)
   call check(error, len(atl), 4)

end subroutine test_atomlist_constr2

subroutine test_atomlist_constr3(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   atl = TAtomList(list=ipar, truth=.true.)
   call check(error, size(atl), 8)
   call check(error, len(atl), 5)
   call atl%to_list(list)
   call check(error, all(list == ipar))

end subroutine test_atomlist_constr3

subroutine test_atomlist_constr4(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   atl = TAtomList(list=cpar, truth=.false.)
   call check(error, size(atl), 8)
   call check(error, len(atl), 5)
   call atl%to_string(string)
   call check(error, string, cpar)
   call atl%new

end subroutine test_atomlist_constr4

subroutine test_atomlist_manip1(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   call atl%new(lpar)
   call atl%resize(9)
   call atl%switch_truth
   call atl%to_string(string)
   call check(error, string, '2,6-7,9')
   call atl%gather(atoms, list)
   call check(error, all(list == [1,1,1,5]))
   call check(error, size(list), len(atl))
   call atl%new

end subroutine test_atomlist_manip1

subroutine test_atomlist_manip2(error)
   type(error_type), allocatable, intent(out) :: error
   type(TAtomList) :: atl
   character(len=:), allocatable :: string
   integer, allocatable :: list(:)

   atl = TAtomList(list=lpar, truth=.false., delimiter=' ', skip=':')
   call atl%to_string(string)
   call check(error, string, '2 6:7 9')
   call atl%switch_truth
   call atl%to_list(list)
   call atl%switch_truth
   call atl%add(list)
   call check(error, len(atl), size(atl))
   call atl%to_string(string)
   call check(error, string, '1:9')

end subroutine test_atomlist_manip2

end module test_atomlist
