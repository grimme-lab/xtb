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

!> Unit tests for the xtb_bmatrix module (Wilson B-matrix derivatives).
!> Every analytic B row is validated against central finite differences of the
!> coordinate it differentiates; the linear-bend pair additionally gets
!> translational-invariance and orthogonality checks on its two frame rows.
module test_bmatrix
   use testdrive, only : new_unittest, unittest_type, error_type, check
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, linbend_frame, &
      & bmat_torsion, bmat_outofplane
   implicit none
   private

   real(wp), parameter :: fd_eps = 1.0e-6_wp
   real(wp), parameter :: thr = 100 * epsilon(0.0_wp)
   real(wp), parameter :: thr_loose = 1.0e-9_wp

   public :: collect_bmatrix

contains

!> Collect all exported unit tests
subroutine collect_bmatrix(testsuite)
   !> Tests exported to the harness.
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("bond-fd", test_bond_fd), &
      new_unittest("angle-fd", test_angle_fd), &
      new_unittest("angle-near-antiparallel", test_angle_near_antiparallel), &
      new_unittest("linbend-fd", test_linbend_fd), &
      new_unittest("linbend-invariance", test_linbend_invariance), &
      new_unittest("torsion-fd", test_torsion_fd), &
      new_unittest("torsion-sign", test_torsion_sign), &
      new_unittest("torsion-atom-order", test_torsion_atom_order), &
      new_unittest("outofplane-fd", test_outofplane_fd), &
      new_unittest("outofplane-degenerate", test_outofplane_degenerate) &
      ]

end subroutine collect_bmatrix


!> Pins the analytic bond row of bmat_bond against central finite differences
!> of the bond length for every Cartesian displacement.
subroutine test_bond_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz2(3, 2) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.2_wp, -0.3_wp], shape(xyz2))
   real(wp) :: bmat(6), fd(6)
   integer :: k, c

   bmat = bmat_bond(xyz2(:, 1) - xyz2(:, 2))
   fd = 0.0_wp
   do k = 1, 2
      do c = 1, 3
         fd(c+(k-1)*3) = fd_bond(xyz2, c, k)
      end do
   end do
   call check_bmat(error, bmat, fd, thr_loose)
end subroutine test_bond_fd

!> Central FD for bond B-row entry (c, k)
function fd_bond(xyz0, c, k) result(dval)
   !> Reference geometry, two atoms.
   real(wp), intent(in) :: xyz0(3, 2)
   !> Displaced Cartesian component and atom index.
   integer, intent(in) :: c, k

   real(wp) :: dval
   real(wp) :: xyz(3, 2), rp, rm

   xyz = xyz0
   xyz(c, k) = xyz(c, k) + fd_eps
   rp = bond_val(xyz)
   xyz = xyz0
   xyz(c, k) = xyz(c, k) - fd_eps
   rm = bond_val(xyz)
   dval = (rp - rm) / (2.0_wp*fd_eps)
end function fd_bond

!> Bond length as scalar function of xyz
function bond_val(xyz) result(val)
   !> Reference geometry, two atoms.
   real(wp), intent(in) :: xyz(3, 2)

   real(wp) :: val
   val = norm2(xyz(:, 2) - xyz(:, 1))
end function bond_val


!> Pins the analytic angle row of bmat_angle against central finite differences
!> of the valence angle for every Cartesian displacement.
subroutine test_angle_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! 3-atom fragment: j at origin, i and k on either side
   real(wp), parameter :: xyz3(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.5_wp, 0.8_wp, 0.1_wp], shape(xyz3))
   real(wp) :: bmat(9), fd(9)
   integer :: k, c

   ! bmat_angle takes (vec_ij, vec_kj) — i-j and k-j
   bmat = bmat_angle(xyz3(:, 1) - xyz3(:, 2), xyz3(:, 3) - xyz3(:, 2))
   fd = 0.0_wp
   do k = 1, 3
      do c = 1, 3
         fd(c+(k-1)*3) = fd_angle(xyz3, c, k)
      end do
   end do
   call check_bmat(error, bmat, fd, thr_loose)
end subroutine test_angle_fd

!> Near-antiparallel angles retain their finite limiting row and inverse-length
!> scaling even when rounding makes the normalized dot product exactly -1.
subroutine test_angle_near_antiparallel(error)
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: delta = 1.0e-9_wp
   real(wp), parameter :: scale = 8.0_wp
   real(wp), parameter :: vec1(3) = [1.0_wp, 0.0_wp, 0.0_wp]
   real(wp), parameter :: vec2(3) = [-1.0_wp, delta, 0.0_wp]
   real(wp), parameter :: expected(9) = [ &
      & 0.0_wp, -1.0_wp, 0.0_wp, &
      & 0.0_wp, 2.0_wp, 0.0_wp, &
      & 0.0_wp, -1.0_wp, 0.0_wp]
   real(wp) :: bmat(9), scaled_bmat(9)

   bmat = bmat_angle(vec1, vec2)
   scaled_bmat = bmat_angle(scale*vec1, scale*vec2)

   call check(error, all(ieee_is_finite(bmat)))
   if (allocated(error)) return
   call check_bmat(error, bmat, expected, thr_loose)
   if (allocated(error)) return
   call check_bmat(error, scaled_bmat, expected/scale, thr_loose)
end subroutine test_angle_near_antiparallel

!> Central FD for angle B-row entry (c, k)
function fd_angle(xyz0, c, k) result(dval)
   !> Reference geometry, three atoms.
   real(wp), intent(in) :: xyz0(3, 3)
   !> Displaced Cartesian component and atom index.
   integer, intent(in) :: c, k

   real(wp) :: dval
   real(wp) :: xyz(3, 3), rp, rm

   xyz = xyz0
   xyz(c, k) = xyz(c, k) + fd_eps
   rp = angle_val(xyz)
   xyz = xyz0
   xyz(c, k) = xyz(c, k) - fd_eps
   rm = angle_val(xyz)
   dval = (rp - rm) / (2.0_wp*fd_eps)
end function fd_angle

!> Valence angle (i-j-k) as scalar function of xyz
function angle_val(xyz) result(val)
   !> Reference geometry, three atoms.
   real(wp), intent(in) :: xyz(3, 3)

   real(wp) :: val
   real(wp) :: v1(3), v2(3), co
   v1 = xyz(:, 1) - xyz(:, 2)
   v2 = xyz(:, 3) - xyz(:, 2)
   co = dot_product(v1, v2) / (norm2(v1)*norm2(v2))
   co = max(-1.0_wp, min(1.0_wp, co))
   val = acos(co)
end function angle_val


!> Pins translational invariance and mutual orthogonality of the two rows a
!> linear centre contributes: the linbend_frame directions through bmat_linbend.
subroutine test_linbend_invariance(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! Near-linear 3-atom fragment: centre at the origin, atoms 1 and 3 collinear
   real(wp), parameter :: xyz_lin(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & -0.99_wp, 0.01_wp, 0.0_wp], shape(xyz_lin))
   real(wp) :: vec1(3), vec2(3), evec1(3), evec2(3), bmat(2, 9)
   real(wp) :: row_sum_x, row_sum_y, row_sum_z
   real(wp) :: dot_rows
   integer :: row, c, nonzero_count

   vec1 = xyz_lin(:, 1) - xyz_lin(:, 2)
   vec2 = xyz_lin(:, 3) - xyz_lin(:, 2)
   call linbend_frame(vec1/norm2(vec1), evec1, evec2)
   bmat(1, :) = bmat_linbend(vec1, vec2, evec1)
   bmat(2, :) = bmat_linbend(vec1, vec2, evec2)

   ! Per-component translational invariance:
   ! atom 1 contributes entries 1:3, the centre 4:6, atom 2 entries 7:9
   ! For each Cartesian component (x=1,4,7; y=2,5,8; z=3,6,9) the row sum must be 0
   do row = 1, 2
      row_sum_x = bmat(row, 1) + bmat(row, 4) + bmat(row, 7)
      row_sum_y = bmat(row, 2) + bmat(row, 5) + bmat(row, 8)
      row_sum_z = bmat(row, 3) + bmat(row, 6) + bmat(row, 9)
      call check(error, row_sum_x, 0.0_wp, thr=thr, &
         & message="linbend row x-sum nonzero")
      if (allocated(error)) return
      call check(error, row_sum_y, 0.0_wp, thr=thr, &
         & message="linbend row y-sum nonzero")
      if (allocated(error)) return
      call check(error, row_sum_z, 0.0_wp, thr=thr, &
         & message="linbend row z-sum nonzero")
      if (allocated(error)) return
   end do

   ! Both rows must be nonzero
   do row = 1, 2
      nonzero_count = 0
      do c = 1, 9
         if (abs(bmat(row, c)) > 1.0e-10_wp) nonzero_count = nonzero_count + 1
      end do
      if (nonzero_count == 0) then
         call check(error, 1.0_wp, 0.0_wp, &
            & message="linbend row is all zeros")
         return
      end if
   end do

   ! The two rows should be mutually orthogonal
   dot_rows = dot_product(bmat(1, :), bmat(2, :))
   call check(error, dot_rows, 0.0_wp, thr=thr, &
      & message="linbend rows not orthogonal")
end subroutine test_linbend_invariance


!> Pins the analytic linear-bend row of bmat_linbend: the row must be the exact
!> derivative of q = evec . v1_hat + evec . v2_hat at a near-linear centre that
!> is not aligned with a coordinate axis.
subroutine test_linbend_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! unit axis u and a perpendicular unit vector p; the centre sits at the origin
   real(wp), parameter :: u(3) = [0.6_wp, 0.8_wp, 0.0_wp]
   real(wp), parameter :: p(3) = [-0.8_wp, 0.6_wp, 0.0_wp]
   real(wp), parameter :: zax(3) = [0.0_wp, 0.0_wp, 1.0_wp]
   real(wp) :: xyz(3, 3), evec(3), bmat(9), fd(9)
   real(wp) :: dev
   integer :: k, c

   ! bond lengths 1.2 and 1.5, one degree off collinearity along p
   dev = 1.5_wp * tan(pi/180.0_wp)
   xyz = 0.0_wp
   xyz(:, 1) = -1.2_wp * u
   xyz(:, 3) = 1.5_wp * u + dev * p

   ! fixed frame: exactly perpendicular to vec1, tilted out of the u/p plane
   evec = p + 0.5_wp * zax
   evec = evec / norm2(evec)

   ! layout [atom1 (1:3), centre (4:6), atom2 (7:9)]
   bmat = bmat_linbend(xyz(:, 1) - xyz(:, 2), &
      & xyz(:, 3) - xyz(:, 2), evec)
   fd = 0.0_wp
   do k = 1, 3
      do c = 1, 3
         fd(c+(k-1)*3) = fd_linbend(xyz, c, k, evec)
      end do
   end do
   call check_bmat(error, bmat, fd, thr_loose)
end subroutine test_linbend_fd

!> Central FD for the linear-bend row entry (c, k)
function fd_linbend(xyz0, c, k, evec) result(dval)
   !> Reference geometry, three atoms.
   !> Fixed unit frame direction for the linear-bend coordinate.
   real(wp), intent(in) :: xyz0(3, 3), evec(3)
   !> Displaced Cartesian component and atom index.
   integer, intent(in) :: c, k

   real(wp) :: dval
   real(wp) :: xyz(3, 3), rp, rm

   xyz = xyz0
   xyz(c, k) = xyz(c, k) + fd_eps
   rp = linbend_val(xyz, evec)
   xyz = xyz0
   xyz(c, k) = xyz(c, k) - fd_eps
   rm = linbend_val(xyz, evec)
   dval = (rp - rm) / (2.0_wp*fd_eps)
end function fd_linbend

!> Linear-bend coordinate value as scalar function of xyz
function linbend_val(xyz, evec) result(val)
   !> Reference geometry, three atoms.
   !> Fixed unit vector defining the linear-bend coordinate.
   real(wp), intent(in) :: xyz(3, 3), evec(3)

   real(wp) :: val
   real(wp) :: v1(3), v2(3)

   v1 = xyz(:, 1) - xyz(:, 2)
   v2 = xyz(:, 3) - xyz(:, 2)
   val = dot_product(evec, v1) / norm2(v1) &
      & + dot_product(evec, v2) / norm2(v2)
end function linbend_val

!> Pins the analytic torsion row of bmat_torsion against central finite
!> differences of the dihedral angle for every Cartesian displacement.
subroutine test_torsion_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz4(3, 4) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 0.8_wp, 0.2_wp, &
      & 2.0_wp, 0.9_wp, -0.5_wp], shape(xyz4))
   real(wp) :: bt(3, 4), fd(3, 4)
   integer :: k, c

   bt = bmat_torsion(xyz4)
   fd = 0.0_wp
   do k = 1, 4
      do c = 1, 3
         fd(c, k) = fd_torsion(xyz4, c, k)
      end do
   end do
   call check_bmat_2d_t(error, bt, fd, thr_loose)
end subroutine test_torsion_fd

!> Central FD for torsion B-row entry (c, k)
function fd_torsion(xyz0, c, k) result(dval)
   !> Reference geometry, four atoms.
   real(wp), intent(in) :: xyz0(3, 4)
   !> Displaced Cartesian component and atom index.
   integer, intent(in) :: c, k

   real(wp) :: dval
   real(wp) :: xyz(3, 4), rp, rm

   xyz = xyz0
   xyz(c, k) = xyz(c, k) + fd_eps
   rp = torsion_val(xyz)
   xyz = xyz0
   xyz(c, k) = xyz(c, k) - fd_eps
   rm = torsion_val(xyz)
   dval = (rp - rm) / (2.0_wp*fd_eps)
end function fd_torsion

!> Dihedral angle as scalar function of xyz (4 atoms)
function torsion_val(xyz) result(val)
   !> Reference geometry, four atoms.
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: val
   real(wp) :: v1(3), v2(3), v3(3)
   real(wp) :: n1(3), n2(3), s, c

   v1 = xyz(:, 2) - xyz(:, 1)
   v2 = xyz(:, 3) - xyz(:, 2)
   v3 = xyz(:, 4) - xyz(:, 3)
   n1 = cross(v1, v2)
   n2 = cross(v2, v3)
   s = -norm2(v2) * dot_product(v1, n2)
   c = dot_product(n1, n2)
   val = atan2(s, c)
end function torsion_val


!> Pins the torsion sign convention: two geometries that differ only by flipping
!> atom 1 across the i-j-k plane must give dihedrals of opposite sign.
subroutine test_torsion_sign(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! Geometry A: atom 1 above the i-j-k plane
   real(wp), parameter :: xyz_a(3, 4) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.5_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 0.8_wp, 0.0_wp], shape(xyz_a))
   ! Geometry B: atom 1 below the plane (mirror in z)
   real(wp), parameter :: xyz_b(3, 4) = reshape([ &
      & 0.0_wp, 0.0_wp, -0.5_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 0.8_wp, 0.0_wp], shape(xyz_b))
   real(wp) :: tau_a, tau_b

   tau_a = torsion_val(xyz_a)
   tau_b = torsion_val(xyz_b)

   ! Dihedrals should be opposite in sign
   call check(error, tau_a, -tau_b, thr=thr, &
      & message="torsion sign mismatch: mirror geometry should flip sign")
end subroutine test_torsion_sign


!> Pins torsion invariance under atom reversal: the reversed order 4-3-2-1 gives
!> the same |tau| as 1-2-3-4.
subroutine test_torsion_atom_order(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz_fwd(3, 4) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 1.5_wp, 0.8_wp, 0.2_wp, &
      & 2.0_wp, 0.9_wp, -0.5_wp], shape(xyz_fwd))
   real(wp) :: xyz_rev(3, 4)
   real(wp) :: tau_fwd, tau_rev
   integer :: k

   ! Reverse atom order: atom 4 becomes 1, 3 becomes 2, etc.
   do k = 1, 4
      xyz_rev(:, k) = xyz_fwd(:, 5-k)
   end do

   tau_fwd = torsion_val(xyz_fwd)
   tau_rev = torsion_val(xyz_rev)

   ! |tau| should be the same for reversed ordering
   call check(error, abs(tau_fwd), abs(tau_rev), thr=thr, &
      & message="torsion |tau| differs for reversed atom order")
end subroutine test_torsion_atom_order


!> Pins the analytic out-of-plane row of bmat_outofplane against central finite
!> differences of the out-of-plane angle for every Cartesian displacement.
subroutine test_outofplane_fd(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! 4-atom fragment: atom 4 central, atom 1 out of plane, 2 and 3 in plane
   real(wp), parameter :: xyz4(3, 4) = reshape([ &
      & 0.1_wp, 0.0_wp, 0.5_wp, & ! atom 1 (out of plane, tilted)
      & 1.0_wp, 0.0_wp, 0.0_wp, & ! atom 2 (in plane)
      & -0.5_wp, 0.8_wp, 0.0_wp, & ! atom 3 (in plane)
      & 0.0_wp, 0.0_wp, 0.0_wp], & ! atom 4 (central)
      & shape(xyz4))
   real(wp) :: bt(3, 4), fd(3, 4)
   integer :: k, c

   bt = bmat_outofplane(xyz4)
   fd = 0.0_wp
   do k = 1, 4
      do c = 1, 3
         fd(c, k) = fd_outofplane(xyz4, c, k)
      end do
   end do
   call check_bmat_2d_t(error, bt, fd, thr_loose)
end subroutine test_outofplane_fd

!> Pins the degenerate-plane branch: with atoms 2, 3, 4 collinear the
!> out-of-plane angle is undefined; bmat_outofplane must return teta=0 and
!> bt=0 without crashing.
subroutine test_outofplane_degenerate(error)
   !> Failure report, allocated when the check fails.
   type(error_type), allocatable, intent(out) :: error

   ! Atoms 2, 3, 4 collinear along x-axis; atom 1 out of plane
   real(wp), parameter :: xyz_deg(3, 4) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.5_wp, & ! atom 1 (out of plane)
      & 1.0_wp, 0.0_wp, 0.0_wp, & ! atom 2 (collinear with 3, 4)
      & 2.0_wp, 0.0_wp, 0.0_wp, & ! atom 3 (collinear with 2, 4)
      & 0.0_wp, 0.0_wp, 0.0_wp], & ! atom 4 (central, collinear)
      & shape(xyz_deg))
   real(wp) :: bt(3, 4)

   bt = bmat_outofplane(xyz_deg)

   ! All B-matrix entries should be zero
   call check_bmat_2d_t(error, bt, 0.0_wp*bt, thr)
end subroutine test_outofplane_degenerate

!> Central FD for out-of-plane B-row entry (c, k)
function fd_outofplane(xyz0, c, k) result(dval)
   !> Reference geometry, four atoms with atom 4 central.
   real(wp), intent(in) :: xyz0(3, 4)
   !> Displaced Cartesian component and atom index.
   integer, intent(in) :: c, k

   real(wp) :: dval
   real(wp) :: xyz(3, 4), rp, rm

   xyz = xyz0
   xyz(c, k) = xyz(c, k) + fd_eps
   rp = outofplane_val(xyz)
   xyz = xyz0
   xyz(c, k) = xyz(c, k) - fd_eps
   rm = outofplane_val(xyz)
   dval = (rp - rm) / (2.0_wp*fd_eps)
end function fd_outofplane

!> Out-of-plane angle as scalar function of xyz (4 atoms, atom 4 central)
function outofplane_val(xyz) result(val)
   !> Reference geometry, four atoms with atom 4 central.
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: val
   real(wp) :: r42(3), r43(3), r41(3)
   real(wp) :: n(3), nn, h, r41n

   r42 = xyz(:, 2) - xyz(:, 4)
   r43 = xyz(:, 3) - xyz(:, 4)
   r41 = xyz(:, 1) - xyz(:, 4)
   n = cross(r42, r43)
   nn = norm2(n)
   if (nn < 1.0e-14_wp) then
      val = 0.0_wp
      return
   end if
   n = n / nn
   h = dot_product(r41, n)
   r41n = norm2(r41)
   if (r41n < 1.0e-14_wp) then
      val = 0.0_wp
      return
   end if
   val = -asin(max(-1.0_wp, min(1.0_wp, h/r41n)))
end function outofplane_val


!> Cross product of two 3-vectors
function cross(a, b) result(c)
   !> First vector.
   !> Second vector.
   real(wp), intent(in) :: a(3), b(3)

   real(wp) :: c(3)
   c(1) = a(2) * b(3) - a(3) * b(2)
   c(2) = a(3) * b(1) - a(1) * b(3)
   c(3) = a(1) * b(2) - a(2) * b(1)
end function cross

!> Compare two 1D B-row arrays element-wise
subroutine check_bmat(error, actual, expected, tol)
   !> Failure report, allocated when a comparison fails.
   type(error_type), allocatable, intent(out) :: error
   !> Analytic row and finite-difference reference.
   real(wp), intent(in) :: actual(:), expected(:)
   !> Absolute tolerance for the comparison.
   real(wp), intent(in) :: tol

   integer :: i

   if (size(actual) /= size(expected)) then
      call check(error, size(actual), size(expected))
      return
   end if
   do i = 1, size(actual)
      call check(error, actual(i), expected(i), thr=tol)
      if (allocated(error)) exit
   end do
end subroutine check_bmat

!> Compare two 2D B-row arrays (shape 3, N) element-wise — transposed view
subroutine check_bmat_2d_t(error, actual, expected, tol)
   !> Failure report, allocated when a comparison fails.
   type(error_type), allocatable, intent(out) :: error
   !> Analytic array and finite-difference reference, shape (3, N).
   real(wp), intent(in) :: actual(:, :), expected(:, :)
   !> Absolute tolerance for the comparison.
   real(wp), intent(in) :: tol

   integer :: i, j

   if (size(actual, 2) /= size(expected, 2)) then
      call check(error, size(actual, 2), size(expected, 2))
      return
   end if
   do j = 1, size(actual, 2)
      do i = 1, 3
         call check(error, actual(i, j), expected(i, j), thr=tol)
         if (allocated(error)) return
      end do
   end do
end subroutine check_bmat_2d_t

end module test_bmatrix
