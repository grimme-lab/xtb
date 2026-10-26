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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

module test_fermi
   use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, &
      & skip_test
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : autoev
   use xtb_mctc_constants, only : kB
   use xtb_scc_core, only : fermismear
#if WITH_TBLITE
   use tblite_wavefunction_fermi, only : get_fermi_filling
#endif
   implicit none
   private

   public :: collect_fermi

   !> Convergence threshold of the Fermi level optimization (same as in tblite)
   real(wp), parameter :: thr = min(sqrt(epsilon(1.0_wp)), 1.0e5_wp*epsilon(1.0_wp))
   !> Tolerance for quantities which are exact up to rounding
   real(wp), parameter :: thr_exact = 10.0_wp*epsilon(1.0_wp)
   integer, parameter :: norbs = 6
   real(wp), parameter :: emo(norbs) = &
      & [-0.8_wp, -0.5_wp, -0.3_wp, -0.1_wp, 0.15_wp, 0.4_wp]

contains

!> Collect all exported unit tests
subroutine collect_fermi(testsuite)
   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("electron-count", test_electron_count), &
      new_unittest("reference", test_reference), &
      new_unittest("zero-temperature", test_zero_temperature), &
      new_unittest("empty", test_empty), &
      new_unittest("full", test_full), &
      new_unittest("low-temperature", test_low_temperature), &
      new_unittest("high-temperature", test_high_temperature), &
      new_unittest("tblite", test_tblite) &
      ]

end subroutine collect_fermi


!> Number of electrons is conserved for a range of temperatures
subroutine test_electron_count(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s
   integer :: nel, it
   real(wp), parameter :: temps(4) = [50.0_wp, 300.0_wp, 5000.0_wp, 20000.0_wp]

   do nel = 1, norbs-1
      do it = 1, size(temps)
         call fermismear(.false., norbs, nel, temps(it), emo, occ, fod, ef, s)
         call check(error, sum(occ), real(nel, wp), thr=thr)
         if (allocated(error)) return
         call check(error, all(occ >= 0.0_wp .and. occ <= 1.0_wp))
         if (allocated(error)) return
         call check(error, s <= 0.0_wp)
         if (allocated(error)) return
      end do
   end do

end subroutine test_electron_count


!> Compare against the Fermi level found by an independent bisection
subroutine test_reference(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s, lo, hi, mid, bkt
   integer :: nel, iter
   real(wp), parameter :: t = 2500.0_wp

   bkt = kB * autoev * t
   nel = 3
   lo = emo(1) - 1.0_wp
   hi = emo(norbs) + 1.0_wp
   do iter = 1, 200
      mid = 0.5_wp * (lo + hi)
      if (sum(1.0_wp / (exp((emo - mid) / bkt) + 1.0_wp)) > real(nel, wp)) then
         hi = mid
      else
         lo = mid
      end if
   end do

   call fermismear(.false., norbs, nel, t, emo, occ, fod, ef, s)
   call check(error, ef, mid, thr=thr)
   if (allocated(error)) return
   call check(error, occ(1), 1.0_wp / (exp((emo(1) - mid) / bkt) + 1.0_wp), thr=thr)

end subroutine test_reference


!> Non-positive temperature gives the Aufbau occupation
subroutine test_zero_temperature(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s
   integer :: nel

   nel = 4
   occ = -1.0_wp
   call fermismear(.false., norbs, nel, 0.0_wp, emo, occ, fod, ef, s)
   call check(error, all(occ(:nel) == 1.0_wp))
   if (allocated(error)) return
   call check(error, all(occ(nel+1:) == 0.0_wp))
   if (allocated(error)) return
   call check(error, fod, 0.0_wp, thr=thr_exact)
   if (allocated(error)) return
   call check(error, s, 0.0_wp, thr=thr_exact)
   if (allocated(error)) return
   call check(error, ef, 0.5_wp * (emo(nel) + emo(nel+1)), thr=thr_exact)

end subroutine test_zero_temperature


!> Without electrons nothing is occupied
subroutine test_empty(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s

   occ = -1.0_wp
   call fermismear(.false., norbs, 0, 300.0_wp, emo, occ, fod, ef, s)
   call check(error, all(occ == 0.0_wp))
   if (allocated(error)) return
   call check(error, fod, 0.0_wp, thr=thr_exact)

end subroutine test_empty


!> Completely filled spectrum, e.g. a minimal basis set without LUMO
subroutine test_full(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s

   occ = -1.0_wp
   call fermismear(.false., norbs, norbs, 300.0_wp, emo, occ, fod, ef, s)
   call check(error, all(occ == 1.0_wp))
   if (allocated(error)) return
   call check(error, fod, 0.0_wp, thr=thr_exact)

end subroutine test_full


!> Tiny electronic temperature must not produce overflow or NaN
subroutine test_low_temperature(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s
   integer :: nel

   nel = 3
   call fermismear(.false., norbs, nel, 1.0e-3_wp, emo, occ, fod, ef, s)
   call check(error, all(occ == occ))
   if (allocated(error)) return
   call check(error, ef == ef)
   if (allocated(error)) return
   call check(error, sum(occ), real(nel, wp), thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(occ(:nel) - 1.0_wp) < thr))
   if (allocated(error)) return
   call check(error, all(occ(nel+1:) < thr))

end subroutine test_low_temperature


!> Degenerate levels at very high temperature
subroutine test_high_temperature(error)
   type(error_type), allocatable, intent(out) :: error
   real(wp) :: occ(norbs), fod, ef, s
   real(wp), parameter :: deg(norbs) = [-0.2_wp, -0.2_wp, -0.2_wp, 0.1_wp, 0.1_wp, 0.1_wp]

   call fermismear(.false., norbs, 3, 1.0e6_wp, deg, occ, fod, ef, s)
   call check(error, sum(occ), 3.0_wp, thr=thr)
   if (allocated(error)) return
   call check(error, all(abs(occ - 0.5_wp) < 1.0e-2_wp))

end subroutine test_high_temperature


!> Occupations and Fermi level agree with the tblite implementation
subroutine test_tblite(error)
   type(error_type), allocatable, intent(out) :: error

#if WITH_TBLITE
   real(wp) :: occ(norbs), occ_ref(norbs), fod, ef, ef_ref, s
   integer :: nel, it, homo
   real(wp), parameter :: temps(4) = [50.0_wp, 300.0_wp, 5000.0_wp, 20000.0_wp]

   do nel = 1, norbs-1
      do it = 1, size(temps)
         call fermismear(.false., norbs, nel, temps(it), emo, occ, fod, ef, s)
         call get_fermi_filling(real(nel, wp), kB * autoev * temps(it), emo, homo, &
            & occ_ref, ef_ref)
         call check(error, maxval(abs(occ - occ_ref)), 0.0_wp, thr=thr_exact)
         if (allocated(error)) return
         call check(error, ef, ef_ref, thr=thr_exact)
         if (allocated(error)) return
      end do
   end do

   ! Non-positive temperature and empty spectrum use the Aufbau filling
   call fermismear(.false., norbs, 3, 0.0_wp, emo, occ, fod, ef, s)
   call get_fermi_filling(3.0_wp, 0.0_wp, emo, homo, occ_ref, ef_ref)
   call check(error, maxval(abs(occ - occ_ref)), 0.0_wp, thr=thr_exact)
   if (allocated(error)) return
   call fermismear(.false., norbs, 0, 300.0_wp, emo, occ, fod, ef, s)
   call get_fermi_filling(0.0_wp, kB * autoev * 300.0_wp, emo, homo, occ_ref, ef_ref)
   call check(error, maxval(abs(occ - occ_ref)), 0.0_wp, thr=thr_exact)
#else
   call skip_test(error, "Compiled without tblite support")
#endif

end subroutine test_tblite

end module test_fermi
