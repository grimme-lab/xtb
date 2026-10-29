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
   ! elements not present in the molecule free nothing
   call check_selection(error, 'fix', 'free: N', [1, 2, 3, 4, 5, 6])
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
   ! invalid lists are rejected without changing the selection
   call check_selection(error, 'fix', 'free: 1,7', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: H,Xx', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free:', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: H,,O', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'fix', 'free: H; free: Xx', [1, 4])

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
   call check_selection(error, 'constrain', 'free: Xx', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'free:', [integer ::])
   if (allocated(error)) return
   call check_selection(error, 'constrain', 'atoms: 2,7', [integer ::])

end subroutine test_constrain_free

!> Process the lines of a $fix or $constrain group, separated by semicolons,
!> for a water dimer and compare the selected atoms with the reference
subroutine check_selection(error, group, input, ref)
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule, only : TMolecule, init
   use xtb_type_identitymap, only : TIdentityMap, init
   use xtb_fixparam, only : fixset
   use xtb_scanparam, only : potset
   use xtb_constrain_param, only : set_fix, set_constr, apply_free_atoms
   type(error_type), allocatable, intent(out) :: error
   !> Data group, either fix or constrain
   character(len=*), intent(in) :: group
   !> Lines of the data group separated by semicolons
   character(len=*), intent(in) :: input
   !> Expected fixed or constrained atoms
   integer, intent(in) :: ref(:)

   integer, parameter :: nat = 6
   integer, parameter :: at(nat) = [8, 1, 1, 8, 1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([&
      & 0.0_wp, 0.0_wp, 0.0_wp, 1.8_wp, 0.0_wp, 0.0_wp, -0.4_wp, 1.7_wp, 0.0_wp, &
      & 5.5_wp, 0.0_wp, 0.0_wp, 6.0_wp, 1.7_wp, 0.0_wp, 6.0_wp, -0.8_wp, 1.5_wp], &
      & shape(xyz))
   type(TEnvironment) :: env
   type(TMolecule) :: mol
   type(TIdentityMap) :: idMap
   character(len=:), allocatable :: line
   integer, allocatable :: selected(:)
   integer :: ist, ie, ic

   call init(env)
   call init(mol, at, xyz)
   call init(idMap, mol)
   call fixset%deallocate
   call potset%pos%deallocate

   ist = 1
   do while (ist <= len(input))
      ie = index(input(ist:), ';')
      if (ie == 0) then
         ie = len(input) + 1
      else
         ie = ist + ie - 1
      end if
      line = input(ist:ie-1)
      ist = ie + 1
      ic = index(line, ':')
      if (group == 'fix') then
         call set_fix(env, trim(adjustl(line(:ic-1))), trim(adjustl(line(ic+1:))), &
            & nat, at, idMap, xyz)
      else
         call set_constr(env, trim(adjustl(line(:ic-1))), trim(adjustl(line(ic+1:))), &
            & nat, at, idMap, xyz)
      end if
   end do
   call apply_free_atoms

   if (group == 'fix') then
      selected = fixset%atoms(:fixset%n)
   else
      selected = potset%pos%atoms(:potset%pos%n)
   end if
   call fixset%deallocate
   call potset%pos%deallocate

   if (size(selected) == size(ref)) then
      if (all(selected == ref)) return
   end if
   call test_failed(error, "Wrong selection for $"//group//" input '"//input//"'")

end subroutine check_selection

end module test_atomlist
