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
      new_unittest("manip2", test_atomlist_manip2), &
      new_unittest("fix_free", test_fix_free), &
      new_unittest("constrain_free", test_constrain_free) &
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

!> Fix all atoms except the ones selected with the free key in $fix
subroutine test_fix_free(error)
   type(error_type), allocatable, intent(out) :: error

   call check_selection(error, 'fix', 'free: 2-3', [1, 4, 5, 6])
   if (allocated(error)) return
   ! elements can be freed by symbol and mixed with atom numbers
   call check_selection(error, 'fix', 'free: H', [1, 4])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: h, 4', [1])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: H,O', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: 1-6', [integer ::])
   if (allocated(error)) return
   ! a trailing comma is ignored, like for all other keys
   call check_selection(error, 'fix', 'free: H,', [1, 4])
   if (allocated(error)) return
   ! elements not present in the molecule free nothing, but are reported
   call check_selection(error, 'fix', 'free: N', [1, 2, 3, 4, 5, 6], &
      & warning=.true.)
   if (allocated(error)) return
   ! repeated free keywords are combined
   call check_selection(error, 'fix', 'free: H; free: O', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: 2; free: 5-6', [1, 3, 4])
   if (allocated(error)) return
   ! atoms selected by other keys stay fixed, independent of the order
   call check_selection(error, 'fix', 'free: 4-6; atoms: 5', [1, 2, 3, 5])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'atoms: 5; free: 4-6', [1, 2, 3, 5])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: H; elements: H', [1, 2, 3, 4, 5, 6])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'elements: H; free: H', [1, 2, 3, 4, 5, 6])
   if (allocated(error)) return
   ! atoms not present in the molecule are rejected
   call check_selection(error, 'fix', 'atoms: 1-8', [integer ::], warning=.true.)
   if (allocated(error)) return
   call check_selection(error, 'fix', 'atoms: 1-8; free: 2', [1, 3, 4, 5, 6], &
      & warning=.true.)
   if (allocated(error)) return
   ! every entry must be an atom, a range of atoms or an element symbol
   call check_rejected(error, 'fix', 'free: 1,7')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: 0')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: 2-')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: H,Xx')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free:')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: H,,O')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: ,H')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: all')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: H 4')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: 2 3')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: C1')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: H2O')
   if (allocated(error)) return
   call check_rejected(error, 'fix', 'free: H; free: Xx')

end subroutine test_fix_free

!> Constrain all atoms except the ones selected with the free key in $constrain
subroutine test_constrain_free(error)
   type(error_type), allocatable, intent(out) :: error

   call check_selection(error, 'constrain', 'free: H', [1, 4])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'free: H; free: O', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'free: 4-6; atoms: 5', [1, 2, 3, 5])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'atoms: 5; free: 4-6', [1, 2, 3, 5])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'free: 2-5; atoms: 3', [1, 3, 6])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'atoms: 2,7', [integer ::], &
      & warning=.true.)
   if (allocated(error)) return
   ! a single constrained atom is rejected, the potential needs pairs of atoms
   call check_selection(error, 'constrain', 'free: 2-6', [integer ::], &
      & warning=.true.)
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'atoms: 3', [integer ::], &
      & warning=.true.)
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'atoms: 1; free: 2-6', [integer ::], &
      & warning=.true.)
   if (allocated(error)) return
   call check_rejected(error, 'constrain', 'free: Xx')
   if (allocated(error)) return
   call check_rejected(error, 'constrain', 'free:')

end subroutine test_constrain_free

!> Read a $fix or $constrain group for a water dimer and compare the selected
!> atoms with the reference
subroutine check_selection(error, group, input, ref, warning)
   use xtb_type_environment, only : TEnvironment
   type(error_type), allocatable, intent(out) :: error
   !> Data group, either fix or constrain
   character(len=*), intent(in) :: group
   !> Lines of the data group separated by semicolons
   character(len=*), intent(in) :: input
   !> Expected fixed or constrained atoms
   integer, intent(in) :: ref(:)
   !> A warning is expected, default is no message at all
   logical, intent(in), optional :: warning

   type(TEnvironment) :: env
   integer, allocatable :: selected(:)
   logical :: expect_warning, failed

   call read_selection(group, input, env, selected)

   expect_warning = .false.
   if (present(warning)) expect_warning = warning
   failed = any(env%log(:env%nLog)%error)
   if (failed) then
      call test_failed(error, &
      & "Unexpected error for $"//group//" input '"//input//"'")
      return
   end if
   if ((env%nLog > 0) .neqv. expect_warning) then
      call test_failed(error, &
      & "Unexpected warnings for $"//group//" input '"//input//"'")
      return
   end if

   if (size(selected) == size(ref)) then
      if (all(selected == ref)) return
   end if
   call test_failed(error, &
   & "Wrong selection for $"//group//" input '"//input//"'")

end subroutine check_selection

!> Read a $fix or $constrain group for a water dimer and check that the
!> input is rejected with an error
subroutine check_rejected(error, group, input)
   use xtb_type_environment, only : TEnvironment
   type(error_type), allocatable, intent(out) :: error
   !> Data group, either fix or constrain
   character(len=*), intent(in) :: group
   !> Lines of the data group separated by semicolons
   character(len=*), intent(in) :: input

   type(TEnvironment) :: env
   integer, allocatable :: selected(:)
   logical :: failed

   call read_selection(group, input, env, selected)

   failed = any(env%log(:env%nLog)%error)
   if (.not.failed) then
      call test_failed(error, &
      & "Missing error for $"//group//" input '"//input//"'")
   end if

end subroutine check_rejected

!> Write the lines of a data group to an input file, read it with
!> read_userdata for a water dimer and return the fixed or constrained atoms
subroutine read_selection(group, input, env, selected)
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init
   use xtb_fixparam, only : fixset
   use xtb_scanparam, only : potset
   use xtb_constrain_param, only : read_userdata
   !> Data group, either fix or constrain
   character(len=*), intent(in) :: group
   !> Lines of the data group separated by semicolons
   character(len=*), intent(in) :: input
   !> Calculation environment holding the raised warnings and errors
   type(TEnvironment), intent(out) :: env
   !> Fixed or constrained atoms
   integer, allocatable, intent(out) :: selected(:)

   character(len=*), parameter :: fname = 'test-atomlist-free.inp'
   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [8, 1, 1, 8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
   & 0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, -0.4_wp, 1.7_wp, 0.0_wp, &
   & 5.5_wp, 0.0_wp, 0.0_wp, 6.0_wp, 1.7_wp, 0.0_wp, 6.0_wp, -0.8_wp, 1.5_wp], &
   & shape(xyz))
   type(TMolecule) :: mol
   integer :: unit, ist, ie

   call init(env)
   call init(mol, at, xyz)
   call fixset%deallocate
   call potset%pos%deallocate

   open(newunit=unit, file=fname, status='replace', action='write')
   write(unit, '(a)') '$'//group
   ist = 1
   do while (ist <= len(input))
      ie = index(input(ist:), ';')
      if (ie == 0) then
         ie = len(input) + 1
      else
         ie = ist + ie - 1
      end if
      write(unit, '(3x,a)') trim(adjustl(input(ist:ie-1)))
      ist = ie + 1
   end do
   write(unit, '(a)') '$end'
   close(unit)

   ! keep the summary of the molecule out of the test output
   open(newunit=env%unit, status='scratch')
   call read_userdata(fname, env, mol)
   close(env%unit)
   open(newunit=unit, file=fname, status='old')
   close(unit, status='delete')

   allocate(selected(0))
   if (group == 'fix') then
      if (fixset%n > 0) selected = fixset%atoms(:fixset%n)
   else
      if (potset%pos%n > 0) selected = potset%pos%atoms(:potset%pos%n)
   end if
   call fixset%deallocate
   call potset%pos%deallocate

end subroutine read_selection

end module test_atomlist
