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

!> Compliance matrix driver.
!>
!> Coordinates: xtb_zmat::TZMatrix, built from the reference geometry by
!> init(zmat, n, at, xyz).  BFS connectivity from covalent radii (Alvarez
!> 2008), automatic angle/dihedral reference atoms, Decius linear-bend pairs
!> measured along the fixed perpendicular frames stored at build time.
!> zmat%bmatrix supplies the Wilson B matrix.
!>
!> C = B H^+ B^T, H projected out of translations and rotations first;
!> handles redundant coordinate sets (symmetric tops etc.).
!>
!> Equivalent to the projected-force-constant route
!>    F = G^+ B H B^T G^+ ,  G = B B^T ,  C = F^+ ,
!> but neither G nor F nor their pseudoinverses are needed: H^+ is formed
!> directly in the non-rigid subspace, which is where the rows of B live.
!>
!> Raw numerical Hessians (the ones handed over by the frequency code) carry
!> residual curvature along those directions, and its reciprocal would
!> otherwise dominate C.
!>
!> Works for any nint and any coordinate set (diatomic, linear, mixed,
!> general, redundant).  nint is passed explicitly -- no hardcoded 3N-6.
!>
!> Output: bonds -> angles -> dihedrals -> linear bends, each with C_ii,
!> 1/C_ii and the local mode frequency.  The full matrix is dumped to
!> compliance.dat: diagonal C_ii, 1/C_ii and the top-20 off-diagonal
!> couplings per coordinate, sorted by |C_ij| descending.
!>
!> Units: C in Bohr^2/Hartree, 1/C in Eh/a0^2, 1 Eh/a0^2 = 15.570 N/cm.
!>
!> Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558.
!>       J. Grunenberg, Chem. Sci. 6 (2015), 4086.
!>
!> SG (with Claude), 05/26

module xtb_compliance
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   use xtb_mctc_convert, only : autoamu
   use xtb_mctc_symbols, only : toSymbol
   use xtb_zmat_type, only : TZMatrix, init, COORD_ANGLE, COORD_DIHEDRAL, &
      & COORD_LINBEND
   implicit none
   private

   public :: compliance_driver, compute_compliance

contains

!> Print compliance constants (bonds, angles, dihedrals, linear bends) for the
!> reference geometry and dump the full matrix to compliance.dat.
subroutine compliance_driver(unit, n, at, xyz, hess, mass)
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Number of atoms.
   integer, intent(in) :: n
   !> Atomic numbers, dimension (n).
   integer, intent(in) :: at(n)
   !> Cartesian coordinates in Bohr, dimension (3, n).
   real(wp), intent(in) :: xyz(3, n)
   !> Cartesian Hessian in Hartree/Bohr^2, dimension (3*n, 3*n).
   real(wp), intent(in) :: hess(3*n, 3*n)
   !> Atomic masses in atomic mass units, dimension (n).
   real(wp), intent(in) :: mass(n)

   integer :: istat
   type(TZMatrix) :: zmat
   real(wp), allocatable :: bmat(:, :), compl(:, :)

   write(unit, *)
   write(unit, *) "             ======================================="
   write(unit, *) "             |                                     |"
   write(unit, *) "             |       compliance constants          |"
   write(unit, *) "             |                                     |"
   write(unit, *) "             ======================================="
   write(unit, *)
   write(unit, *) "Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558."

   call init(zmat, n, at, xyz)
   allocate(bmat(zmat%nint, 3*n), compl(zmat%nint, zmat%nint))
   call zmat%get_bmatrix(xyz, bmat)
   call compute_compliance(unit, hess, bmat, xyz, n, zmat%nint, compl, istat)
   if (istat /= 0) return
   call print_compl(unit, n, at, xyz, mass, zmat, compl)

end subroutine compliance_driver


!> Print diagonal elements of the compliance matrix in the order of bonds,
!> angles, dihedrals, and linear bends.
!>
!> Also calls write_compliance_dat to write the full matrix to a file.
subroutine print_compl(unit, n, at, xyz, mass, zmat, C)
   !> Z-matrix container with the internal-coordinate definitions.
   type(TZMatrix), intent(in) :: zmat
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Number of atoms.
   integer, intent(in) :: n
   !> Atomic numbers, dimension (n).
   integer, intent(in) :: at(n)
   !> Cartesian coordinates in Bohr, dimension (3, n).
   real(wp), intent(in) :: xyz(3, n)
   !> Atomic masses in atomic mass units, dimension (n).
   real(wp), intent(in) :: mass(n)
   !> Compliance matrix, dimension (nint, nint), in Bohr^2/Hartree.
   real(wp), intent(in) :: C(zmat%nint, zmat%nint)

   integer :: i, o, k
   integer :: idx_ord(zmat%nint)
   real(wp) :: q(zmat%nint), cc, mu, freq
   character(20) :: s
   ! local mode frequency: nu_a = fac * sqrt(k_a[a.u.] / mu[a.u.])
   !   k_a = 1/C_ii  in Eh/a0^2
   !   mu  = m_A*m_B/(m_A+m_B)  in amu  (NOT converted to me)
   !   fac = 1/(2*pi*c) * sqrt(Eh/(a0^2 * amu))
   !       = 219474.631 * sqrt(me/amu)
   !       = 219474.631 / sqrt(1822.888)  =  5140.487  -- WRONG if mu in me
   !   Correct: keep mu in amu, use fac below (CODATA 2018):
   !     1/(2*pi*c) * sqrt(Eh/(a0^2*amu)) = 5140.4869 cm^-1
   !   Derivation:
   !     Eh/a0^2 = 1556.89 N/m
   !     1 amu   = 1.66054e-27 kg
   !     sqrt(1556.89/1.66054e-27) / (2*pi*2.99792e10) = 5140.49 cm^-1
   real(wp), parameter :: fac = 5140.4869_wp ! cm^-1, mu must be in amu

   call zmat%get_coords(xyz, q)

   write(unit, *)
   write(unit, *) "units: Hartree, Bohr, radian"
   write(unit, *) "1 Eh/a0^2 (1/C = relaxed force constant) = 15.570 N/cm"
   write(unit, *) "local mode frequency nu_loc = 5140.487*sqrt(1/(mu[amu]*C[a.u.])) cm^-1"
   write(unit, *) "Ref.: Cremer, Kraka, Zou, J. Chem. Theory Comput. 8 (2012) 2864."

   k = 0

   ! 1) bond stretches
   ! local mode frequency (Kraka/Cremer = 1/C_ii route):
   ! nu_a = fac * sqrt(1/(mu_AB * C_ii))   [cm^-1]
   ! mu_AB = m_A*m_B/(m_A+m_B)  in atomic mass units -> convert to me
   write(unit, "(a)") &
      "     type                atoms                   coord value" // &
      "       C       1/C    nu_loc/cm-1"
   do i = 2, n
      k = k + 1;  o = zmat%qoff(i);  idx_ord(k) = o + 1
      cc = C(o+1, o+1);  s = "bond stretch"
      mu = mass(i) * mass(zmat%na(i)) / (mass(i) + mass(zmat%na(i))) * autoamu
      freq = fac * sqrt(1.0_wp/(mu*cc))
      write(unit, "(i4,1x,a14,2(a2,i3,3x),16x,4f10.2)") &
         k, s, toSymbol(at(i)), i, toSymbol(at(zmat%na(i))), zmat%na(i), &
         q(o+1), cc, 1.0_wp / cc, freq
   end do

   ! 2) valence angles
   do i = 2, n
      if (zmat%ctype(i)/=COORD_ANGLE .and. zmat%ctype(i)/=COORD_DIHEDRAL) cycle
      k = k + 1;  o = zmat%qoff(i);  idx_ord(k) = o + 2
      cc = C(o+2, o+2);  s = "angle"
      write(unit, "(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)") &
         k, s, toSymbol(at(i)), i, toSymbol(at(zmat%na(i))), zmat%na(i), &
         toSymbol(at(zmat%nb(i))), zmat%nb(i), &
         q(o+2), cc, 1.0_wp / cc
   end do

   ! 3) dihedrals
   do i = 2, n
      if (zmat%ctype(i)/=COORD_DIHEDRAL) cycle
      k = k + 1;  o = zmat%qoff(i);  idx_ord(k) = o + 3
      cc = C(o+3, o+3);  s = "dihedral"
      write(unit, "(i4,1x,a14,4(a2,i3,3x),3f10.4)") &
         k, s, toSymbol(at(i)), i, toSymbol(at(zmat%na(i))), zmat%na(i), &
         toSymbol(at(zmat%nb(i))), zmat%nb(i), toSymbol(at(zmat%nc(i))), zmat%nc(i), &
         q(o+3), cc, 1.0_wp / cc
   end do

   ! 4) linear bending coordinates
   do i = 2, n
      if (zmat%ctype(i)/=COORD_LINBEND) cycle
      o = zmat%qoff(i)
      k = k + 1;  idx_ord(k) = o + 2
      cc = C(o+2, o+2);  s = "lin.bend(e1)"
      write(unit, "(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)") &
         k, s, toSymbol(at(zmat%nb(i))), zmat%nb(i), &
         toSymbol(at(zmat%na(i))), zmat%na(i), toSymbol(at(i)), i, &
         q(o+2), cc, 1.0_wp / cc
      k = k + 1;  idx_ord(k) = o + 3
      cc = C(o+3, o+3);  s = "lin.bend(e2)"
      write(unit, "(i4,1x,a14,3(a2,i3,3x),8x,3f10.4)") &
         k, s, toSymbol(at(zmat%nb(i))), zmat%nb(i), &
         toSymbol(at(zmat%na(i))), zmat%na(i), toSymbol(at(i)), i, &
         q(o+3), cc, 1.0_wp / cc
   end do

   call write_compliance_dat(unit, n, at, zmat, C, idx_ord, k)

end subroutine print_compl


!> Write the full compliance matrix to a file named "compliance.dat".
!>
!> For each coordinate it outputs the diagonal element C_ii, its inverse 1/C_ii,
!> and the top NCOUP off-diagonal couplings |C_ij| sorted in descending order.
subroutine write_compliance_dat(unit, n, at, zmat, C, idx_ord, ncoord)
   !> Z-matrix container with the internal-coordinate definitions.
   type(TZMatrix), intent(in) :: zmat
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Number of atoms.
   integer, intent(in) :: n
   !> Number of internal coordinates.
   integer, intent(in) :: ncoord
   !> Atomic numbers, dimension (n).
   integer, intent(in) :: at(n)
   !> Compliance matrix, dimension (nint, nint), in Bohr^2/Hartree.
   real(wp), intent(in) :: C(zmat%nint, zmat%nint)
   !> Coordinate order used for the printed table, dimension (ncoord).
   integer, intent(in) :: idx_ord(ncoord)

   integer, parameter :: NCOUP = 20
   integer :: iunit, i, j, p, q, ii, jj, nc_act
   integer :: jsort(ncoord)
   real(wp) :: aval(ncoord), tmp_r
   integer :: tmp_i
   character(20) :: lbl(ncoord)

   ! Build labels in the same order as the printed coordinates.
   call build_labels(n, at, zmat, ncoord, lbl)

   open(newunit=iunit, file="compliance.dat", status="replace")

   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# compliance.dat"
   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# units: C   in Bohr^2/Hartree (a0^2/Eh)"
   write(iunit, "(a)") "#        1/C in Eh/a0^2  (relaxed force constant)"
   write(iunit, "(a)") "#        conversion: 1 Eh/a0^2 = 15.570 N/cm"
   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# Ref.: K. Brandhorst, J. Grunenberg,"
   write(iunit, "(a)") "#       Chem. Soc. Rev. 37 (2008) 1558."
   write(iunit, "(a)") "#"
   write(iunit, "(a,i6)") "# number of internal coordinates :", ncoord
   write(iunit, "(a,i4)")  "# top couplings shown per coord  :", NCOUP
   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# columns: coord_j  label_j  C_ij  [<-> coord_i label_i]"
   write(iunit, "(a)") "#"

   do i = 1, ncoord
      ii = idx_ord(i)

      write(iunit, "(a)")  ""
      write(iunit, "(a,i4,2x,a20,a,f12.6,a,f12.6)") &
         "# coord ", i, lbl(i), &
         "   C_ii=", C(ii, ii), "   1/C_ii=", 1.0_wp / C(ii, ii)

      ! diagonal entry
      write(iunit, "(2x,i4,2x,a20,2f14.6,a)") &
         i, lbl(i), C(ii, ii), 1.0_wp / C(ii, ii), "  (diagonal)"

      ! collect off-diagonal |C_ij|
      nc_act = 0
      do j = 1, ncoord
         if (j == i) cycle
         nc_act = nc_act + 1
         jsort(nc_act) = j
         aval(nc_act) = abs(C(ii, idx_ord(j)))
      end do

      ! insertion sort descending
      do p = 2, nc_act
         tmp_r = aval(p);  tmp_i = jsort(p);  q = p - 1
         do while (q >= 1 .and. aval(q) < tmp_r)
            aval(q+1) = aval(q);  jsort(q+1) = jsort(q);  q = q - 1
         end do
         aval(q+1) = tmp_r;  jsort(q+1) = tmp_i
      end do

      ! write top NCOUP couplings
      do p = 1, min(NCOUP, nc_act)
         j = jsort(p);  jj = idx_ord(j)
         if (aval(p) < 1.0e-12_wp) exit
         write(iunit, "(2x,i4,2x,a20,f14.6,a,i4,2x,a20)") &
            j, lbl(j), C(ii, jj), "  <-> ", i, lbl(i)
      end do

   end do

   write(iunit, "(a)") ""
   write(iunit, "(a)") "# end of compliance.dat"
   close(iunit)

   write(unit, *)
   write(unit, "(a,i4,a)") &
      " compliance matrix written to compliance.dat  (", ncoord, " coordinates)"

end subroutine write_compliance_dat


!> Build human-readable labels for the internal coordinates, in the order bonds,
!> angles, dihedrals, then linear bends.
subroutine build_labels(n, at, zmat, ncoord, lbl)
   !> Z-matrix container with the internal-coordinate definitions.
   type(TZMatrix), intent(in) :: zmat
   !> Number of atoms.
   integer, intent(in) :: n
   !> Number of internal coordinates.
   integer, intent(in) :: ncoord
   !> Atomic numbers, dimension (n).
   integer, intent(in) :: at(n)
   !> Coordinate labels, dimension (ncoord), truncated to 20 characters.
   character(20), intent(out) :: lbl(ncoord)

   integer :: i, k
   character(80) :: buf

   k = 0;  lbl = "??"

   do i = 2, n ! bonds
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0)") &
         "bond ", toSymbol(at(i)), i, "-", toSymbol(at(zmat%na(i))), zmat%na(i)
      lbl(k) = buf(1:20)
   end do
   do i = 2, n ! angles
      if (zmat%ctype(i)/=COORD_ANGLE .and. zmat%ctype(i)/=COORD_DIHEDRAL) cycle
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0)") &
         "ang ", toSymbol(at(i)), i, "-", toSymbol(at(zmat%na(i))), zmat%na(i), &
         "-", toSymbol(at(zmat%nb(i))), zmat%nb(i)
      lbl(k) = buf(1:20)
   end do
   do i = 2, n ! dihedrals
      if (zmat%ctype(i)/=COORD_DIHEDRAL) cycle
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0,a,a2,i0)") &
         "dih ", toSymbol(at(i)), i, "-", toSymbol(at(zmat%na(i))), zmat%na(i), "-", &
         toSymbol(at(zmat%nb(i))), zmat%nb(i), "-", toSymbol(at(zmat%nc(i))), zmat%nc(i)
      lbl(k) = buf(1:20)
   end do
   do i = 2, n ! linear bends
      if (zmat%ctype(i)/=COORD_LINBEND) cycle
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0)") &
         "lb1 ", toSymbol(at(zmat%nb(i))), zmat%nb(i), "-", &
         toSymbol(at(zmat%na(i))), zmat%na(i), "-", toSymbol(at(i)), i
      lbl(k) = buf(1:20)
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0)") &
         "lb2 ", toSymbol(at(zmat%nb(i))), zmat%nb(i), "-", &
         toSymbol(at(zmat%na(i))), zmat%na(i), "-", toSymbol(at(i)), i
      lbl(k) = buf(1:20)
   end do

end subroutine build_labels


!> Compliance constants of the reference geometry as the pseudoinverse of the
!> Hessian projected out of translations and rotations.
!>
!> C = B H^+ B^T, evaluated as Z (Z D)^T with Z = B V and D = diag(1/w_i) from
!> the eigen decomposition Hp = V W V^T of the projected, symmetrised Hessian.
subroutine compute_compliance(unit, H, B, xyz, natoms, nint, C, stat)
   !> Formatted output unit for the rank diagnostic.
   integer, intent(in) :: unit
   !> Number of atoms.
   integer, intent(in) :: natoms
   !> Number of internal coordinates, passed explicitly -- no hardcoded 3N-6.
   integer, intent(in) :: nint
   !> Cartesian Hessian in Hartree/Bohr^2, dimension (3*natoms, 3*natoms).
   real(wp), intent(in) :: H(3*natoms, 3*natoms)
   !> Wilson B matrix, dimension (nint, 3*natoms).
   real(wp), intent(in) :: B(nint, 3*natoms)
   !> Cartesian reference coordinates in Bohr, dimension (3, natoms).
   real(wp), intent(in) :: xyz(3, natoms)
   !> Compliance matrix, dimension (nint, nint), in Bohr^2/Hartree.
   real(wp), intent(out) :: C(nint, nint)
   !> Status: zero on success, the LAPACK info of the failing DSYEV otherwise.
   integer, intent(out) :: stat

   integer :: i, j, ndim, lwork, info, nrigid, nvib, rank_h
   real(wp) :: tol_h, normq, center(3), scale
   real(wp), parameter :: eps_svd = 1.0e-10_wp
   real(wp), allocatable :: Hp(:, :), Q(:, :), W(:), Z(:, :), ZD(:, :), &
      & T1(:, :), T2(:, :), work(:)

   stat = 0; ndim = 3 * natoms
   allocate(Hp(ndim, ndim), Q(ndim, 6), W(ndim), Z(nint, ndim), &
      & ZD(nint, ndim), T1(ndim, 6), T2(6, ndim))

   ! orthonormal basis of the rigid (translation + rotation) space
   center = sum(xyz, dim=2) / natoms
   scale = max(1.0_wp, maxval(abs(xyz)))
   Q = 0.0_wp
   do i = 1, 3
      Q(i::3, i) = 1.0_wp
   end do
   do j = 1, natoms
      Q(3*j-2:3*j, 4) = crossProd([1.0_wp, 0.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 5) = crossProd([0.0_wp, 1.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 6) = crossProd([0.0_wp, 0.0_wp, 1.0_wp], xyz(:, j) - center)
   end do

   ! modified Gram-Schmidt; the rotation about the molecular axis of a
   ! linear molecule produces no displacement
   nrigid = 0
   do i = 1, 6
      do j = 1, i - 1
         normq = dot_product(Q(:, i), Q(:, j))
         Q(:, i) = Q(:, i) - normq * Q(:, j)
      end do
      normq = norm2(Q(:, i))
      if (normq < 1.0e-8_wp*scale) then
         Q(:, i) = 0.0_wp
         cycle
      end if
      nrigid = nrigid + 1
      Q(:, i) = Q(:, i) / normq
   end do
   nvib = ndim - nrigid

   ! H averaged into symmetry first: the numerical Hessian handed over by
   ! the frequency code is only symmetric to within its finite-difference
   ! noise
   Hp = 0.5_wp * (H + transpose(H))

   ! Hp = (1 - Q Q^T) Hp (1 - Q Q^T)
   call dgemm("N", "N", ndim, 6, ndim, 1.0_wp, Hp, ndim, Q, ndim, 0.0_wp, T1, ndim)
   call dgemm("N", "T", ndim, ndim, 6, -1.0_wp, T1, ndim, Q, ndim, 1.0_wp, Hp, ndim)
   call dgemm("T", "N", 6, ndim, ndim, 1.0_wp, Q, ndim, Hp, ndim, 0.0_wp, T2, 6)
   call dgemm("N", "N", ndim, ndim, 6, -1.0_wp, Q, ndim, T2, 6, 1.0_wp, Hp, ndim)
   call symmetrise(Hp, ndim)

   ! Hp = V W V^T, overwriting Hp with V
   lwork = -1; allocate(work(1))
   call dsyev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   lwork = int(work(1)); deallocate(work); allocate(work(lwork))
   call dsyev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   if (info /= 0) then
      write(unit, "(A,I0)") "compute_compliance: DSYEV(H) info=", info
      stat = info; return
   end if

   ! C = B V D V^T B^T = Z (Z D)^T,  Z = B V,  D = diag(1/w_i)
   ! Modes of the projected Hessian that are still zero are dropped.
   tol_h = eps_svd * maxval(abs(W))
   rank_h = count(abs(W) > tol_h)
   if (rank_h /= nvib) then
      write(unit, "(A,I0,A,I0)") &
         & "  Note: Hessian rank=", rank_h, " /= 3N-rigid=", nvib
   end if
   call dgemm("N", "N", nint, ndim, ndim, 1.0_wp, B, nint, Hp, ndim, 0.0_wp, Z, nint)
   ZD = Z
   do i = 1, ndim
      if (abs(W(i)) <= tol_h) then
         ZD(:, i) = 0.0_wp
      else
         ZD(:, i) = ZD(:, i) / W(i)
      end if
   end do
   call dgemm("N", "T", nint, nint, ndim, 1.0_wp, Z, nint, ZD, nint, 0.0_wp, C, nint)
   C = 0.5_wp * (C + transpose(C))
   deallocate(Hp, Q, W, Z, ZD, T1, T2, work)

end subroutine compute_compliance


!> Average the two triangles; mirroring one of them would turn rounding
!> noise into a symmetric perturbation.
subroutine symmetrise(A, n)
   !> Order of the matrix.
   integer, intent(in) :: n
   !> Square matrix to symmetrise in place, dimension (n, n).
   real(wp), intent(inout) :: A(n, n)

   integer :: i, j
   do i = 1, n
      do j = i + 1, n
         A(j, i) = 0.5_wp * (A(i, j) + A(j, i))
         A(i, j) = A(j, i)
      end do
   end do
end subroutine symmetrise


end module xtb_compliance
