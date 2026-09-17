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
module xtb_bmatrix
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_math, only : crossProd
   implicit none

   private
   public :: bmat_bond, bmat_angle, bmat_linbend, &
      & bmat_torsion, bmat_outofplane, oop_angle
   public :: bmat_accum_packed, bmat_accum_dense, bmat_accum_pairblock_packed

   character(len=*), parameter :: source = 'xtb_bmatrix'

   !> Tolerance for near-linear/collinear detection (radians)
   real(wp), parameter :: tol_collinear = 1.0e-13_wp
   !> Tolerance for near-zero vector norms
   real(wp), parameter :: tol_norm = 1.0e-15_wp
   !> Tolerance for cross-product magnitude (collinearity guard)
   real(wp), parameter :: tol_cross = 1.0e-10_wp
   !> Tolerance for near-linear bend detection
   real(wp), parameter :: tol_bend = 1.0e-6_wp
   !> Avoid zero division in angle derivatives (1e-24 is ~1e-12 in energy)
   real(wp), parameter :: eps_sqrt = 1.0e-24_wp

contains

!> Clamp a cosine to [-1, 1], rounding can push it outside and make acos NaN
elemental function clamp(co) result(c)
   real(wp), intent(in) :: co
   real(wp) :: c
   c = min(1.0_wp, max(-1.0_wp, co))
end function clamp

!> Wilson B matrix for bonds
pure function bmat_bond(vec) result(bmat)
   real(wp), intent(in) :: vec(3)

   real(wp) :: l, bmat(6)

   bmat = 0.0_wp
   l = norm2(vec)

   bmat(1:3) = vec(:) / l
   bmat(4:6) = -vec(:) / l
end function bmat_bond

!> Wilson B matrix for non-linear angles
pure function bmat_angle(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)

   real(wp) :: bmat(9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: dl(2, 6), dnvec(2, 3, 6), dinprod(9)
   real(wp) :: dot_n1n2
   integer :: ii

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   dl = 0.0_wp
   dl(1, 1:3) = nvec1
   dl(1, 4:6) = -nvec1
   dl(2, 1:3) = nvec2
   dl(2, 4:6) = -nvec2

   dnvec = 0.0_wp
   do ii = 1, 6
      dnvec(1, 1:3, ii) = -nvec1 * dl(1, ii) / l1
      dnvec(2, 1:3, ii) = -nvec2 * dl(2, ii) / l2
   end do
   do ii = 1, 3
      dnvec(1, ii, ii) = dnvec(1, ii, ii) + 1.0_wp / l1
      dnvec(2, ii, ii) = dnvec(2, ii, ii) + 1.0_wp / l2
      dnvec(1, ii, ii+3) = dnvec(1, ii, ii+3) - 1.0_wp / l1
      dnvec(2, ii, ii+3) = dnvec(2, ii, ii+3) - 1.0_wp / l2
   end do

   dinprod = 0.0_wp
   do ii = 1, 3
      dinprod(ii) = dot_product(dnvec(1, :, ii), nvec2)
      dinprod(ii+3) = dot_product(dnvec(1, :, ii+3), nvec2) &
         & + dot_product(dnvec(2, :, ii+3), nvec1)
      dinprod(ii+6) = dot_product(dnvec(2, :, ii), nvec1)
   end do

   dot_n1n2 = dot_product(nvec1, nvec2)
   bmat = -dinprod / sqrt(max(eps_sqrt, 1.0_wp-dot_n1n2**2))
end function bmat_angle

!> Wilson B matrix for linear angles
pure function bmat_linbend(vec1, vec2) result(bmat)
   real(wp), intent(in) :: vec1(3), vec2(3)

   real(wp) :: bmat(2, 9)
   real(wp) :: l1, l2, nvec1(3), nvec2(3)
   real(wp) :: vn(3), vn2(3), nvn
   real(wp), parameter :: xaxis(3) = [1.0_wp, 0.0_wp, 0.0_wp], yaxis(3) = [0.0_wp, 1.0_wp, 0.0_wp]

   l1 = norm2(vec1)
   l2 = norm2(vec2)
   nvec1 = vec1 / l1
   nvec2 = vec2 / l2

   vn = crossProd(vec1, vec2)
   nvn = norm2(vn)

   if (nvn < tol_norm) then
      vn = xaxis - dot_product(xaxis, vec1) / l1**2 * vec1
      nvn = norm2(vn)
      if (nvn < tol_norm) then
         vn = yaxis - dot_product(yaxis, vec1) / l1**2 * vec1
         nvn = norm2(vn)
      end if
   end if
   vn = vn / nvn

   vn2 = crossProd(vec1-vec2, vn)
   vn2 = vn2 / norm2(vn2)

   bmat = 0.0_wp
   bmat(2, 1:3) = vn / l1
   bmat(2, 7:9) = vn / l2
   bmat(2, 4:6) = -bmat(2, 1:3) - bmat(2, 7:9)
   bmat(1, 1:3) = vn2 / l1
   bmat(1, 7:9) = vn2 / l2
   bmat(1, 4:6) = -bmat(1, 1:3) - bmat(1, 7:9)
end function bmat_linbend

!> Wilson B matrix for torsions (4 atoms)
pure function bmat_torsion(xyz) result(bt)
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: bt(3, 4)
   real(wp) :: rij1, rjk1, rkl1
   real(wp) :: brij(3, 2), brjk(3, 2), brkl(3, 2)
   real(wp) :: fi2, sinfi2, cosfi2
   real(wp) :: fi3, sinfi3, cosfi3
   integer :: ix, iy, iz

   rij1 = norm2(xyz(:, 2)-xyz(:, 1))
   rjk1 = norm2(xyz(:, 3)-xyz(:, 2))
   rkl1 = norm2(xyz(:, 4)-xyz(:, 3))
   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   brkl = bmat_strtch(xyz(:, 3:4))

   fi2 = bend_angle(xyz(:, 1:3))
   sinfi2 = sin(fi2)
   cosfi2 = cos(fi2)
   fi3 = bend_angle(xyz(:, 2:4))
   sinfi3 = sin(fi3)
   cosfi3 = cos(fi3)

   do ix = 1, 3
      iy = ix + 1
      if (iy > 3) iy = iy - 3
      iz = iy + 1
      if (iz > 3) iz = iz - 3
      bt(ix, 1) = (brij(iy, 2)*brjk(iz, 2)-brij(iz, 2)*brjk(iy, 2)) &
         & / (rij1*sinfi2**2)
      bt(ix, 4) = (brkl(iy, 1)*brjk(iz, 1)-brkl(iz, 1)*brjk(iy, 1)) &
         & / (rkl1*sinfi3**2)
      bt(ix, 2) = -((rjk1-rij1*cosfi2)*bt(ix, 1) &
         & + rkl1 * cosfi3 * bt(ix, 4)) / rjk1
      bt(ix, 3) = -(bt(ix, 1)+bt(ix, 2)+bt(ix, 4))
   end do
end function bmat_torsion

!> Wilson B matrix for out-of-plane angles (4 atoms, atom 4 central)
pure function bmat_outofplane(xyz) result(bt)
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: bt(3, 4)
   real(wp) :: r2(3), r3(3)
   real(wp) :: q42, q43, e42(3), e43(3)
   real(wp) :: cosfi1, fi1
   real(wp) :: c14(3, 3), br14(3, 3)
   real(wp) :: r42(3), r43(3)
   integer :: ix, iy, iz

!  4 -> 2 (bond in plane)
   r2 = xyz(:, 2) - xyz(:, 4)
   q42 = norm2(r2)
   e42 = r2 / q42
!  4 -> 3 (bond in plane)
   r3 = xyz(:, 3) - xyz(:, 4)
   q43 = norm2(r3)
   e43 = r3 / q43

!  angle between e43 and e42
   cosfi1 = clamp(dot_product(e43, e42))
   fi1 = acos(cosfi1)

!  dirty exit when 2-3-4 collinear
   if (abs(fi1-pi) < tol_collinear) then
      bt = 0.0_wp
      return
   end if

!  first two centers
   c14(:, 1) = xyz(:, 1)
   c14(:, 2) = xyz(:, 4)

!  3rd center: cross product of r42 x r43
   r42 = xyz(:, 2) - xyz(:, 4)
   r43 = xyz(:, 3) - xyz(:, 4)
   c14(:, 3) = crossProd(r42, r43)

!  exit if 2-3-4 collinear
   if (sum(c14(:, 3)**2) < tol_cross) then
      bt = 0.0_wp
      return
   end if
   c14(:, 3) = c14(:, 3) + xyz(:, 4)

   br14 = bmat_bend(c14)

!  compute the WDC matrix
   do ix = 1, 3
      iy = mod(ix+1, 4) + (ix+1) / 4
      iz = mod(iy+1, 4) + (iy+1) / 4

      bt(ix, 1) = -br14(ix, 1)
      bt(ix, 2) = r43(iz) * br14(iy, 3) - r43(iy) * br14(iz, 3)
      bt(ix, 3) = -r42(iz) * br14(iy, 3) + r42(iy) * br14(iz, 3)
      bt(ix, 4) = -(bt(ix, 1)+bt(ix, 2)+bt(ix, 3))
   end do

   bt = -bt
end function bmat_outofplane

!> Out-of-plane angle for 4-atom fragment (atom 4 central).
!> Returns the Wilson out-of-plane angle in radians, i.e. the angle
!> between atom 1 and the plane (2,3,4) minus pi/2. Returns 0.0 when
!> atoms 2-3-4 are collinear (ill-defined).
pure function oop_angle(xyz) result(theta)
   real(wp), intent(in) :: xyz(3, 4)

   real(wp) :: theta
   real(wp) :: r2(3), r3(3), q42, q43, e42(3), e43(3)
   real(wp) :: cosfi1, fi1
   real(wp) :: r42(3), r43(3), c14(3, 3)

   r2 = xyz(:, 2) - xyz(:, 4)
   q42 = norm2(r2)
   e42 = r2 / q42
   r3 = xyz(:, 3) - xyz(:, 4)
   q43 = norm2(r3)
   e43 = r3 / q43

   cosfi1 = clamp(dot_product(e43, e42))
   fi1 = acos(cosfi1)

   if (abs(fi1-pi) < tol_collinear) then
      theta = 0.0_wp
      return
   end if

   c14(:, 1) = xyz(:, 1)
   c14(:, 2) = xyz(:, 4)
   r42 = xyz(:, 2) - xyz(:, 4)
   r43 = xyz(:, 3) - xyz(:, 4)
   c14(:, 3) = crossProd(r42, r43)

   if (sum(c14(:, 3)**2) < tol_cross) then
      theta = 0.0_wp
      return
   end if
   c14(:, 3) = c14(:, 3) + xyz(:, 4)

   theta = bend_angle(c14) - 0.5_wp*pi
end function oop_angle

!> Stretch B row for 2-atom fragment (private helper)
pure function bmat_strtch(xyz) result(b)
   real(wp), intent(in) :: xyz(3, 2)

   real(wp) :: b(3, 2)
   real(wp) :: r(3), rr

   r = xyz(:, 2) - xyz(:, 1)
   rr = norm2(r)
   b(:, 1) = -r / rr
   b(:, 2) = -b(:, 1)
end function bmat_strtch

!> Bend angle for 3-atom fragment (private helper)
pure function bend_angle(xyz) result(fir)
   real(wp), intent(in) :: xyz(3, 3)

   real(wp) :: fir
   real(wp) :: brij(3, 2), brjk(3, 2)
   real(wp) :: co, crap
   integer :: i

   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   co = 0.0_wp
   crap = 0.0_wp
   do i = 1, 3
      co = co + brij(i, 1) * brjk(i, 2)
      crap = crap + (brjk(i, 2)+brij(i, 1))**2
   end do
   co = clamp(co)
   if (sqrt(crap) < tol_bend) then
      fir = pi - asin(sqrt(crap))
   else
      fir = acos(co)
   end if
   if (abs(fir-pi) < tol_collinear) then
      fir = pi
   end if
end function bend_angle

!> Bend B row for 3-atom fragment (private helper)
pure function bmat_bend(xyz) result(bf)
   real(wp), intent(in) :: xyz(3, 3)

   real(wp) :: bf(3, 3)
   real(wp) :: brij(3, 2), brjk(3, 2)
   real(wp) :: co, crap, fir, si
   real(wp) :: rij1, rjk1
   integer :: i

   brij = bmat_strtch(xyz(:, 1:2))
   brjk = bmat_strtch(xyz(:, 2:3))
   rij1 = norm2(xyz(:, 2)-xyz(:, 1))
   rjk1 = norm2(xyz(:, 3)-xyz(:, 2))
   co = 0.0_wp
   crap = 0.0_wp
   do i = 1, 3
      co = co + brij(i, 1) * brjk(i, 2)
      crap = crap + (brjk(i, 2)+brij(i, 1))**2
   end do
   co = clamp(co)
   if (sqrt(crap) < tol_bend) then
      fir = pi - asin(sqrt(crap))
      si = sqrt(crap)
   else
      fir = acos(co)
      si = sqrt(1.0_wp-co**2)
   end if
   bf = 0.0_wp
   if (abs(fir-pi) < tol_collinear) then
      return
   end if
   do i = 1, 3
      bf(i, 1) = (co*brij(i, 1)-brjk(i, 2)) / (si*rij1)
      bf(i, 3) = (co*brjk(i, 2)-brij(i, 1)) / (si*rjk1)
      bf(i, 2) = -(bf(i, 1)+bf(i, 3))
   end do
end function bmat_bend

!> Accumulate scale * B^T B into a packed lower-triangular Cartesian Hessian.
!>
!> Layout:
!>   nat    - number of atoms in the full molecule
!>   hess   - packed array of length (3*nat)*(3*nat+1)/2
!>   atoms  - fragment atom indices (1-based global), size nfrag
!>   brow   - flattened Wilson B row, length 3*nfrag, ordered as
!>            [a1_x, a1_y, a1_z, a2_x, a2_y, a2_z, ...]
!>   scale  - scalar force constant multiplying the outer product
!>
!> Only lower-packed entries with global Cartesian j <= i are updated.
!> Caller is responsible for ensuring atoms are unique and brow has the
!> expected length; this routine does no validation.
pure subroutine bmat_accum_packed(nat, hess, atoms, brow, scale)
   integer, intent(in) :: nat
   real(wp), intent(inout) :: hess((3*nat)*(3*nat+1)/2)
   integer, intent(in) :: atoms(:)
   real(wp), intent(in) :: brow(:)
   real(wp), intent(in) :: scale

   integer :: nfrag, p, q, cp, cq, gi, gj, imax, imin
   real(wp) :: val

   nfrag = size(atoms)
   do p = 1, nfrag
      do cp = 1, 3
         gi = (atoms(p)-1)*3 + cp
         do q = 1, nfrag
            do cq = 1, 3
               gj = (atoms(q)-1)*3 + cq
               if (gj > gi) cycle
               val = scale * brow((p-1)*3+cp) * brow((q-1)*3+cq)
               imax = gi
               imin = gj
               hess(imax*(imax-1)/2 + imin) = hess(imax*(imax-1)/2 + imin) + val
            end do
         end do
      end do
   end do
end subroutine bmat_accum_packed

!> Accumulate scale * B^T B into a dense Cartesian Hessian.
!>
!> Layout:
!>   hess   - dense (3*nat) x (3*nat) array
!>   atoms  - fragment atom indices (1-based global), size nfrag
!>   brow   - flattened Wilson B row, length 3*nfrag, ordered as
!>            [a1_x, a1_y, a1_z, a2_x, a2_y, a2_z, ...]
!>   scale  - scalar force constant multiplying the outer product
!>
!> Updates both (i,j) and (j,i) entries. Caller is responsible for
!> ensuring atoms are unique and brow has the expected length; this
!> routine does no validation. Thread-safe when caller passes a
!> thread-local Hessian.
pure subroutine bmat_accum_dense(hess, atoms, brow, scale)
   real(wp), intent(inout) :: hess(:, :)
   integer, intent(in) :: atoms(:)
   real(wp), intent(in) :: brow(:)
   real(wp), intent(in) :: scale

   integer :: nfrag, p, q, cp, cq, gi, gj
   real(wp) :: val

   nfrag = size(atoms)
   do p = 1, nfrag
      do cp = 1, 3
         gi = (atoms(p)-1)*3 + cp
         do q = 1, nfrag
            do cq = 1, 3
               gj = (atoms(q)-1)*3 + cq
               val = scale * brow((p-1)*3+cp) * brow((q-1)*3+cq)
               hess(gi, gj) = hess(gi, gj) + val
            end do
         end do
      end do
   end do
end subroutine bmat_accum_dense

!> Accumulate a 2-atom Cartesian second-derivative block into a packed
!> lower-triangular Hessian.  Given a symmetric 3x3 block 'mat', stamps:
!>   (i,i) -= mat, (j,j) -= mat, (i,j) += mat
!> matching the sign pattern of a pair second derivative where moving
!> either atom in the same direction cancels and moving them apart adds.
pure subroutine bmat_accum_pairblock_packed(nat, hess, i, j, mat)
   integer, intent(in) :: nat
   real(wp), intent(inout) :: hess((3*nat)*(3*nat+1)/2)
   integer, intent(in) :: i, j
   real(wp), intent(in) :: mat(3,3)

   integer :: ic, jc, gi, gj, imax, imin

   do ic = 1, 3
      do jc = 1, ic
         gi = (i-1)*3 + ic
         gj = (i-1)*3 + jc
         imax = gi
         imin = gj
         hess(imax*(imax-1)/2 + imin) = hess(imax*(imax-1)/2 + imin) - mat(ic,jc)
         gi = (j-1)*3 + ic
         gj = (j-1)*3 + jc
         imax = gi
         imin = gj
         hess(imax*(imax-1)/2 + imin) = hess(imax*(imax-1)/2 + imin) - mat(ic,jc)
      end do
   end do
   do ic = 1, 3
      do jc = 1, 3
         gi = (i-1)*3 + ic
         gj = (j-1)*3 + jc
         imax = max(gi, gj)
         imin = min(gi, gj)
         hess(imax*(imax-1)/2 + imin) = hess(imax*(imax-1)/2 + imin) + mat(ic,jc)
      end do
   end do
end subroutine bmat_accum_pairblock_packed

end module xtb_bmatrix
