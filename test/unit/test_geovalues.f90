! This file is part of xtb.
!
! Copyright (C) 2026 Leopold M. Seidler
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

!> Regression tests for internal-coordinate value and constraint-potential callers.
module test_geovalues
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_intmodes, only : xyzgeo
   use xtb_type_setvar, only : fix_setvar
   use xtb_constrainpot, only : constrain_angle, constrain_dihedral
   use xtb_scanparam, only : nconstr, fcconstr, atconstr, valconstr, clear_constr
   implicit none
   private

   public :: collect_geovalues

   real(wp), parameter :: xyz4(3, 4) = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 1.4_wp,  0.0_wp,  0.0_wp, &
      & 2.1_wp,  1.2_wp,  0.3_wp, &
      & 2.9_wp,  1.5_wp, -0.8_wp], [3, 4])

   real(wp), parameter :: ref_constrain_e = 3.6499794648900553e-2_wp
   real(wp), parameter :: ref_constrain_g(3, 4) = reshape([ &
      &  1.1113262322961571e-17_wp, -2.9120284633054924e-2_wp,  3.2317570544484905e-2_wp, &
      & -9.8267126626542817e-3_wp,   3.7751940745654286e-2_wp, -4.3915198782022297e-2_wp, &
      & -1.6442336449696515e-2_wp,   1.0185180485609371e-2_wp, -2.3752702264789828e-3_wp, &
      &  2.6269049112350784e-2_wp,  -1.8816836598208722e-2_wp,  1.3972898464016375e-2_wp], &
      & [3, 4])

   interface
      subroutine calc_angles(n, at, xyz, bond, maxbend, nbend, bend, ib)
         import :: wp
         integer, intent(in) :: n
         integer, intent(in) :: at(n)
         real(wp), intent(in) :: xyz(3, n)
         integer, intent(in) :: bond(n, n)
         integer(8), intent(in) :: maxbend
         integer, intent(out) :: nbend
         real(wp), intent(out) :: bend(maxbend)
         integer, intent(out) :: ib(3, maxbend)
      end subroutine calc_angles

      subroutine calc_torsions(n, at, xyz, bond, maxtrsn, ntrsn, trsn, it)
         import :: wp
         integer, intent(in) :: n
         integer, intent(in) :: at(n)
         real(wp), intent(in) :: xyz(3, n)
         integer, intent(in) :: bond(n, n)
         integer(8), intent(in) :: maxtrsn
         integer, intent(out) :: ntrsn
         real(wp), intent(out) :: trsn(maxtrsn)
         integer, intent(out) :: it(4, maxtrsn)
      end subroutine calc_torsions

      subroutine constrpot(nat, at, xyz, g, e)
         import :: wp
         integer :: nat, at(nat)
         real(wp) :: e
         real(wp) :: xyz(3, nat)
         real(wp) :: g(3, nat)
      end subroutine constrpot
   end interface


contains


!> Collect all exported unit tests
subroutine collect_geovalues(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("xyzgeo-values", test_xyzgeo_values), &
      & new_unittest("xyzgeo-trans-continuity", test_xyzgeo_trans_continuity), &
      & new_unittest("geometry-calcs", test_geometry_calcs), &
      & new_unittest("constrain-angle-dihedral", test_constrain_angle_dihedral), &
      & new_unittest("constrpot-angle-dihedral", test_constrpot_angle_dihedral) &
      & ]

end subroutine collect_geovalues


subroutine test_xyzgeo_values(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: na(4) = [0, 1, 2, 3]
   integer, parameter :: nb(4) = [0, 0, 1, 2]
   integer, parameter :: nc(4) = [0, 0, 0, 1]
   real(wp), parameter :: thr = 1.0e-12_wp
   real(wp) :: geo(3, 4)

   geo = 0.0_wp
   call xyzgeo(xyz4, 4, na, nb, nc, 1.0_wp, geo)

   call check(error, geo(2, 3), 2.0857772899407054_wp, thr=thr)
   if (allocated(error)) return

   call check(error, geo(2, 4), 1.8734356567038593_wp, thr=thr)
   if (allocated(error)) return

   call check(error, geo(3, 4), -2.1114389127119253_wp, thr=thr)

end subroutine test_xyzgeo_values


subroutine test_xyzgeo_trans_continuity(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: na(4) = [0, 1, 2, 3]
   integer, parameter :: nb(4) = [0, 0, 1, 2]
   integer, parameter :: nc(4) = [0, 0, 0, 1]
   real(wp), parameter :: phi_p = 179.9_wp * pi / 180.0_wp
   real(wp), parameter :: phi_m = -179.9_wp * pi / 180.0_wp
   real(wp) :: xyz_p(3, 4), xyz_m(3, 4), geo_p(3, 4), geo_m(3, 4)

   xyz_p(:, 1) = [0.0_wp, 1.0_wp, 0.0_wp]
   xyz_p(:, 2) = [0.0_wp, 0.0_wp, 0.0_wp]
   xyz_p(:, 3) = [1.0_wp, 0.0_wp, 0.0_wp]
   xyz_p(:, 4) = [1.0_wp, cos(phi_p), sin(phi_p)]

   xyz_m = xyz_p
   xyz_m(:, 4) = [1.0_wp, cos(phi_m), sin(phi_m)]

   geo_p = 0.0_wp
   geo_m = 0.0_wp
   call xyzgeo(xyz_p, 4, na, nb, nc, 180.0_wp / pi, geo_p)
   call xyzgeo(xyz_m, 4, na, nb, nc, 180.0_wp / pi, geo_m)

   call check(error, abs(geo_p(3, 4) - geo_m(3, 4)), 0.2_wp, thr=1.0e-10_wp)

end subroutine test_xyzgeo_trans_continuity


subroutine test_geometry_calcs(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz5(3, 5) = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 1.4_wp,  0.0_wp,  0.0_wp, &
      & 2.1_wp,  1.2_wp,  0.3_wp, &
      & 2.9_wp,  1.5_wp, -0.8_wp, &
      & 4.1_wp,  0.9_wp, -0.4_wp], [3, 5])
   integer, parameter :: at5(5) = [6, 6, 6, 6, 6]
   integer, parameter :: ref_ib(3, 2) = reshape([1, 2, 3, 2, 3, 4], [3, 2])
   integer, parameter :: ref_it(4, 2) = reshape([5, 4, 3, 2, 1, 2, 3, 4], [4, 2])
   real(wp), parameter :: ref_bend(2) = [2.0857772899407059_wp, 1.8734356567038595_wp]
   real(wp), parameter :: ref_trsn(2) = [-1.6438986669503857_wp, -2.1114389127119249_wp]
   real(wp), parameter :: thr = 1.0e-11_wp
   integer(8), parameter :: maxc = 10
   integer :: bond(5, 5), nbend, ntrsn, ntrsn_mir, i
   integer :: ib(3, maxc), it(4, maxc), it_mir(4, maxc)
   real(wp) :: bend(maxc), trsn(maxc), trsn_mir(maxc), xyz5_mir(3, 5)

   bond = 0
   bond(1, 1) = 1
   bond(2, 2) = 2
   bond(3, 3) = 2
   bond(4, 4) = 2
   bond(5, 5) = 1
   bond(1, 2) = 1
   bond(2, 1) = 1
   bond(2, 3) = 1
   bond(3, 2) = 1
   bond(3, 4) = 1
   bond(4, 3) = 1
   bond(4, 5) = 1
   bond(5, 4) = 1

   call calc_angles(5, at5, xyz5, bond, maxc, nbend, bend, ib)
   call calc_torsions(5, at5, xyz5, bond, maxc, ntrsn, trsn, it)

   call check(error, nbend, 2)
   if (allocated(error)) return

   call check(error, ntrsn, 2)
   if (allocated(error)) return

   call check(error, all(ib(:, 1:2) == ref_ib))
   if (allocated(error)) return

   call check(error, all(it(:, 1:2) == ref_it))
   if (allocated(error)) return

   do i = 1, 2
      call check(error, bend(i), ref_bend(i), thr=thr)
      if (allocated(error)) return

      call check(error, trsn(i), ref_trsn(i), thr=thr)
      if (allocated(error)) return
   end do

   xyz5_mir = xyz5
   xyz5_mir(3, :) = -xyz5_mir(3, :)
   call calc_torsions(5, at5, xyz5_mir, bond, maxc, ntrsn_mir, trsn_mir, it_mir)

   call check(error, ntrsn_mir, ntrsn)
   if (allocated(error)) return

   do i = 1, ntrsn
      call check(error, trsn_mir(i), -trsn(i), thr=thr)
      if (allocated(error)) return
   end do

end subroutine test_geometry_calcs


subroutine test_constrain_angle_dihedral(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: at4(4) = [6, 6, 6, 6]
   type(fix_setvar) :: fix_ang, fix_dih
   real(wp) :: e, g(3, 4)
   integer :: i, c

   fix_ang%n = 1
   fix_ang%fc = 0.05_wp
   fix_ang%atoms = [1, 2, 3]
   fix_ang%val = [1.8_wp]

   fix_dih%n = 1
   fix_dih%fc = 0.05_wp
   fix_dih%atoms = [1, 2, 3, 4]
   fix_dih%val = [-0.9_wp]

   e = 0.0_wp
   g = 0.0_wp
   call constrain_angle(fix_ang, 4, at4, xyz4, g, e)
   call constrain_dihedral(fix_dih, 4, at4, xyz4, g, e)

   call check(error, e, ref_constrain_e, thr=1.0e-10_wp)
   if (allocated(error)) return

   do i = 1, 4
      do c = 1, 3
         call check(error, g(c, i), ref_constrain_g(c, i), thr=1.0e-9_wp)
         if (allocated(error)) return
      end do
   end do

end subroutine test_constrain_angle_dihedral


subroutine test_constrpot_angle_dihedral(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: at4(4) = [6, 6, 6, 6]
   integer :: saved_nconstr, i, c
   real(wp) :: saved_fc, e, g(3, 4)

   saved_nconstr = nconstr
   saved_fc = fcconstr
   call clear_constr

   nconstr = 2
   fcconstr = 0.05_wp
   allocate(atconstr(4, 2), source=0)
   allocate(valconstr(2), source=0.0_wp)
   atconstr(:, 1) = [1, 2, 3, 0]
   valconstr(1) = 1.8_wp
   atconstr(:, 2) = [1, 2, 3, 4]
   valconstr(2) = -0.9_wp

   e = 0.0_wp
   g = 0.0_wp
   call constrpot(4, at4, xyz4, g, e)

   call clear_constr
   nconstr = saved_nconstr
   fcconstr = saved_fc

   call check(error, e, ref_constrain_e, thr=1.0e-10_wp)
   if (allocated(error)) return

   do i = 1, 4
      do c = 1, 3
         call check(error, g(c, i), ref_constrain_g(c, i), thr=1.0e-9_wp)
         if (allocated(error)) return
      end do
   end do

end subroutine test_constrpot_angle_dihedral


end module test_geovalues
