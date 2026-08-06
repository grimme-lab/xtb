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

module test_random
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_setparam, only : set, initrand, expand_random_seed
   implicit none
   private

   public :: collect_random

contains

!> Collect tests for scalar RNG seed expansion and replay.
subroutine collect_random(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("seed expansion", test_seed_expansion), &
      new_unittest("explicit seed replay", test_explicit_seed), &
      new_unittest("automatic seed replay", test_automatic_seed), &
      new_unittest("samerand alias", test_samerand) &
      ]
end subroutine collect_random

!> Ensure expansion is deterministic and independent of output-vector length.
subroutine test_seed_expansion(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: short_seed(2), long_seed(8)
   integer, parameter :: expected(8) = [ &
      2075653, 1409598201, 1842888923, 728608805, &
      1335939236, 336425193, 309152689, 245587716]

   call expand_random_seed(42, short_seed)
   call expand_random_seed(42, long_seed)

   call check(error, all(short_seed == long_seed(:size(short_seed))))
   call check(error, all(long_seed == expected))
end subroutine test_seed_expansion

!> Ensure repeated initialization with an explicit scalar replays the sequence.
subroutine test_explicit_seed(error)
   type(error_type), allocatable, intent(out) :: error
   real :: first(4), second(4)
   integer :: old_randseed, nseed
   integer, allocatable :: old_seed_vector(:)
   logical :: old_randseed_set, old_samerand

   old_randseed = set%randseed
   old_randseed_set = set%randseed_set
   old_samerand = set%samerand
   call random_seed(size=nseed)
   allocate(old_seed_vector(nseed))
   call random_seed(get=old_seed_vector)

   set%randseed = 314159
   set%randseed_set = .true.
   set%samerand = .false.
   call initrand
   call random_number(first)
   call initrand
   call random_number(second)

   call check(error, all(first == second))

   set%randseed = old_randseed
   set%randseed_set = old_randseed_set
   set%samerand = old_samerand
   call random_seed(put=old_seed_vector)
end subroutine test_explicit_seed

!> Ensure the automatically selected scalar can reconstruct the generated state.
subroutine test_automatic_seed(error)
   type(error_type), allocatable, intent(out) :: error
   real :: first(4), second(4)
   integer :: old_randseed, selected_seed, nseed
   integer, allocatable :: old_seed_vector(:)
   logical :: old_randseed_set, old_samerand

   old_randseed = set%randseed
   old_randseed_set = set%randseed_set
   old_samerand = set%samerand
   call random_seed(size=nseed)
   allocate(old_seed_vector(nseed))
   call random_seed(get=old_seed_vector)

   set%randseed_set = .false.
   set%samerand = .false.
   call initrand
   selected_seed = set%randseed
   call random_number(first)
   set%randseed = selected_seed
   set%randseed_set = .true.
   call initrand
   call random_number(second)

   call check(error, all(first == second))

   set%randseed = old_randseed
   set%randseed_set = old_randseed_set
   set%samerand = old_samerand
   call random_seed(put=old_seed_vector)
end subroutine test_automatic_seed

!> Ensure the legacy `$samerand` input is equivalent to scalar seed 41.
subroutine test_samerand(error)
   type(error_type), allocatable, intent(out) :: error
   real :: samerand_values(4), seed_values(4)
   integer :: old_randseed, nseed
   integer, allocatable :: old_seed_vector(:)
   logical :: old_randseed_set, old_samerand

   old_randseed = set%randseed
   old_randseed_set = set%randseed_set
   old_samerand = set%samerand
   call random_seed(size=nseed)
   allocate(old_seed_vector(nseed))
   call random_seed(get=old_seed_vector)

   set%randseed_set = .false.
   set%samerand = .true.
   call initrand
   call random_number(samerand_values)
   call check(error, set%randseed, 41)

   set%randseed = 41
   set%randseed_set = .true.
   set%samerand = .false.
   call initrand
   call random_number(seed_values)
   call check(error, all(samerand_values == seed_values))

   set%randseed = old_randseed
   set%randseed_set = old_randseed_set
   set%samerand = old_samerand
   call random_seed(put=old_seed_vector)
end subroutine test_samerand

end module test_random
