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
!> Coordinates: the redundant internal coordinates of the molecular graph,
!> built from the reference geometry.  Bonds come from a covalent-radius
!> neighbour list (TNeighbourList%generate_covalent); the graph and
!> the coordinate set are built by xtb_internals_graph and
!> xtb_internals_redundant.  The set holds bond stretches, ordinary valence
!> angles, fixed-frame pairs for bends near linearity, and defined dihedrals.
!> Dihedrals with a near-linear interior angle are omitted because their
!> torsion is undefined.  xtb_bmatrix::get_bmatrix supplies the Wilson B matrix.
!>
!> C = B H^+ B^T, H projected out of translations and rotations first;
!> handles redundant coordinate sets (symmetric tops etc.).
!>
!> Equivalent to the projected-force-constant route
!>    F = G^+ Bp Hp Bp^T G^+ ,  G = Bp Bp^T ,  C = F^+ ,
!> where Bp and Hp denote non-rigid projections.  Neither G nor F nor their
!> pseudoinverses are needed: H^+ is formed directly in the non-rigid subspace,
!> giving the same projection of each B row.  Fixed-frame linear bends need not
!> themselves annihilate rotations.
!>
!> Raw numerical Hessians (the ones handed over by the frequency code) carry
!> residual curvature along those directions, and its reciprocal would
!> otherwise dominate C.
!>
!> Works for any nint and any coordinate set (diatomic, linear, mixed,
!> general, redundant).  nint is passed explicitly -- no hardcoded 3N-6.
!>
!> Output: bonds -> angular rows -> dihedrals, with C_ii and 1/C_ii, plus local
!> mode frequencies for bonds.  Each near-linear bend has two consecutive
!> fixed-frame components.  A nonlinear near-linear molecule still has 3N-6
!> physical modes, and either projected component can be redundant; an exactly
!> linear molecule has 3N-5 modes.  An exactly zero C_ii is reported with
!> 1/C_ii = +Infinity.  The full matrix is dumped to compliance.dat, including
!> the top-20 off-diagonal couplings per coordinate sorted by |C_ij| descending.
!>
!> Units: C in Bohr^2/Hartree, 1/C in Eh/a0^2, 1 Eh/a0^2 = 15.570 N/cm.
!>
!> Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558.
!>       J. Grunenberg, Chem. Sci. 6 (2015), 4086.
!>
!> SG (with Claude), 05/26

module xtb_compliance
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_value
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   use mctc_env, only : error_type, fatal_error
   use xtb_mctc_convert, only : autoamu
   use xtb_type_molecule, only : TMolecule
   use xtb_type_neighbourlist, only : TNeighbourList, init
   use xtb_internals_graph, only : graph_type, init
   use xtb_internals_type, only : internal_coords_set_type, coord_bond, &
      & coord_angle, coord_linbend, coord_dihedral
   use xtb_internals_redundant, only : redundant_type, init
   use xtb_bmatrix, only : get_bmatrix
   use xtb_mctc_blas, only : blas_gemm
   use xtb_mctc_lapack, only : lapack_syev
   implicit none
   private

   public :: compliance_driver, compute_compliance

contains

!> Print compliance constants (bonds, angular rows, dihedrals) for the reference
!> geometry and dump the full matrix to compliance.dat.  Near-linear bends are
!> represented by two consecutive components in a fixed reference frame.
subroutine compliance_driver(unit, mol, hess, error)
   !> Molecular structure containing atom count, numbers, geometry, and masses.
   type(TMolecule), intent(in) :: mol
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Cartesian Hessian in Hartree/Bohr^2, dimension (3*mol%n, 3*mol%n).
   real(wp), intent(in) :: hess(3*mol%n, 3*mol%n)
   !> Error information.
   type(error_type), allocatable, intent(out) :: error

   type(TNeighbourList) :: neigh_list
   type(graph_type) :: graph
   type(redundant_type) :: internals
   real(wp), allocatable :: bmat(:, :), compl(:, :)

   write(unit, *)
   write(unit, *) "             ======================================="
   write(unit, *) "             |                                     |"
   write(unit, *) "             |       compliance constants          |"
   write(unit, *) "             |                                     |"
   write(unit, *) "             ======================================="
   write(unit, *)
   write(unit, *) "Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558."

   call init(neigh_list, mol%n)
   call neigh_list%generate_covalent(mol%at, mol%xyz)
   call init(graph, neigh_list)
   call init(internals, graph, mol%xyz)

   allocate(bmat(internals%ncoords, 3*mol%n), compl(internals%ncoords, internals%ncoords))
   call get_bmatrix(internals, mol%xyz, bmat)
   call compute_compliance(unit, hess, bmat, mol%xyz, mol%n, internals%ncoords, compl, error)
   if (allocated(error)) return
   call print_compl(unit, mol, internals, compl)

end subroutine compliance_driver


!> Print diagonal elements of the compliance matrix in the order of bonds,
!> angular rows, and dihedrals.  Fixed-frame linear-bend pairs are labelled
!> separately; either component can be redundant after rigid-mode projection.
!>
!> Also calls write_compliance_dat to write the full matrix to a file.
subroutine print_compl(unit, mol, internals, C)
   !> Molecular structure containing atom count, numbers, and masses.
   type(TMolecule), intent(in) :: mol
   !> Redundant coordinate set with the definitions and the reference values.
   type(redundant_type), intent(in) :: internals
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Compliance matrix, dimension (ncoords, ncoords), in Bohr^2/Hartree.
   real(wp), intent(in) :: C(internals%ncoords, internals%ncoords)

   integer :: ic, k, a1, a2
   integer :: idx_ord(internals%ncoords)
   real(wp) :: cc, invcc, mu, freq
   character(20) :: s
   ! local mode frequency: nu_a = fac * sqrt(k_a[a.u.] / mu[a.u.])
   !   k_a = 1/C_ii  in Eh/a0^2
   !   mu  = m_A*m_B/(m_A+m_B)  in electron masses
   !   fac = 1/(2*pi*c) * sqrt(Eh/(a0^2 * amu))
   !       = 5140.4869 cm^-1
   ! Convert mu to amu with autoamu before applying fac.
   real(wp), parameter :: fac = 5140.4869_wp ! cm^-1, mu converted to amu

   write(unit, *)
   write(unit, *) "units: Hartree, Bohr, radian"
   write(unit, *) "linear-bend q components are dimensionless fixed-frame projections"
   write(unit, *) "1 Eh/a0^2 (1/C = relaxed force constant) = 15.570 N/cm"
   write(unit, *) "local mode frequency nu_loc = 5140.487*sqrt(1/(mu[amu]*C[a.u.])) cm^-1"
   write(unit, *) "exact zero C gives 1/C = +Infinity (zero response)"
   write(unit, *) "Ref.: Cremer, Kraka, Zou, J. Chem. Theory Comput. 8 (2012) 2864."

   k = 0

   ! 1) bond stretches
   ! nu_a = fac * sqrt(1/(mu_AB * C_ii))   [cm^-1]
   ! mu_AB is formed from electron-mass inputs and converted to amu below.
   write(unit, "(a)") &
      "     type                atoms                   coord value" // &
      "       C       1/C    nu_loc/cm-1"
   do ic = 1, internals%ncoords
      if (internals%kind(ic) /= coord_bond) cycle
      k = k + 1
      idx_ord(k) = ic
      cc = C(ic, ic)
      invcc = reciprocal(cc)
      s = "bond stretch"
      a1 = internals%atoms(1, ic)
      a2 = internals%atoms(2, ic)
      mu = mol%atmass(a1) * mol%atmass(a2) / (mol%atmass(a1) + mol%atmass(a2)) * autoamu
      freq = fac * sqrt(invcc/mu)
      write(unit, "(i4,1x,a14,2(a2,i3,3x),16x,4f10.2)") &
         k, s, mol%sym(a1), a1, mol%sym(a2), a2, &
         internals%q(ic), cc, invcc, freq
   end do

   ! 2) angular rows
   do ic = 1, internals%ncoords
      if (internals%kind(ic) /= coord_angle .and. &
         & internals%kind(ic) /= coord_linbend) cycle
      k = k + 1
      idx_ord(k) = ic
      cc = C(ic, ic)
      invcc = reciprocal(cc)
      if (internals%kind(ic) == coord_angle) then
         s = "angle"
      else if (is_second_linbend(internals, ic)) then
         s = "linear bend 2"
      else
         s = "linear bend 1"
      end if
      write(unit, "(i4,1x,a14,3(a2,i3,3x),8x,3es16.6e3)") &
         k, s, mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
         mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic), &
         mol%sym(internals%atoms(3, ic)), internals%atoms(3, ic), &
         internals%q(ic), cc, invcc
   end do

   ! 3) dihedrals
   do ic = 1, internals%ncoords
      if (internals%kind(ic) /= coord_dihedral) cycle
      k = k + 1
      idx_ord(k) = ic
      cc = C(ic, ic)
      s = "dihedral"
      write(unit, "(i4,1x,a14,4(a2,i3,3x),3f10.4)") &
         k, s, mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
         mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic), &
         mol%sym(internals%atoms(3, ic)), internals%atoms(3, ic), &
         mol%sym(internals%atoms(4, ic)), internals%atoms(4, ic), &
         internals%q(ic), cc, reciprocal(cc)
   end do

   call write_compliance_dat(unit, mol, internals, C, idx_ord, k)

end subroutine print_compl


!> Write the full compliance matrix to a file named "compliance.dat".
!>
!> For each coordinate it outputs C_ii, its reciprocal (positive infinity when
!> C_ii is exactly zero), and the top NCOUP off-diagonal couplings |C_ij| sorted
!> in descending order.
subroutine write_compliance_dat(unit, mol, internals, C, idx_ord, ncoord)
   !> Molecular structure containing atom count and atomic numbers.
   type(TMolecule), intent(in) :: mol
   !> Internal-coordinate definitions.
   class(internal_coords_set_type), intent(in) :: internals
   !> Formatted output unit.
   integer, intent(in) :: unit
   !> Number of internal coordinates.
   integer, intent(in) :: ncoord
   !> Compliance matrix, dimension (ncoords, ncoords), in Bohr^2/Hartree.
   real(wp), intent(in) :: C(internals%ncoords, internals%ncoords)
   !> Coordinate order used for the printed table, dimension (ncoord).
   integer, intent(in) :: idx_ord(ncoord)

   integer, parameter :: NCOUP = 20
   integer :: iunit, i, j, p, q, ii, jj, nc_act
   integer :: jsort(ncoord)
   real(wp) :: aval(ncoord), tmp_r
   integer :: tmp_i
   character(20) :: lbl(ncoord)

   ! Build labels in the same order as the printed coordinates.
   call build_labels(mol, internals, ncoord, lbl)

   open(newunit=iunit, file="compliance.dat", status="replace")

   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# compliance.dat"
   write(iunit, "(a)") "#"
   write(iunit, "(a)") "# units: C   in Bohr^2/Hartree (a0^2/Eh)"
   write(iunit, "(a)") "#        1/C in Eh/a0^2  (relaxed force constant)"
   write(iunit, "(a)") "#        conversion: 1 Eh/a0^2 = 15.570 N/cm"
   write(iunit, "(a)") "#        exact zero C is reported as 1/C = +Infinity"
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
      write(iunit, "(a,i4,2x,a20,a,es22.12e3,a,es22.12e3)") &
         "# coord ", i, lbl(i), &
         "   C_ii=", C(ii, ii), "   1/C_ii=", reciprocal(C(ii, ii))

      ! diagonal entry
      write(iunit, "(2x,i4,2x,a20,2es22.12e3,a)") &
         i, lbl(i), C(ii, ii), reciprocal(C(ii, ii)), "  (diagonal)"

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
         tmp_r = aval(p)
         tmp_i = jsort(p)
         q = p - 1
         do while (q >= 1)
            if (aval(q) >= tmp_r) exit
            aval(q+1) = aval(q)
            jsort(q+1) = jsort(q)
            q = q - 1
         end do
         aval(q+1) = tmp_r
         jsort(q+1) = tmp_i
      end do

      ! write top NCOUP couplings
      do p = 1, min(NCOUP, nc_act)
         j = jsort(p)
         jj = idx_ord(j)
         if (aval(p) < 1.0e-12_wp) exit
         write(iunit, "(2x,i4,2x,a20,es22.12e3,a,i4,2x,a20)") &
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


!> Build human-readable labels for the internal coordinates, in the order
!> bonds, angular rows, dihedrals.  Fixed-frame bend pairs use lb1/lb2.
subroutine build_labels(mol, internals, ncoord, lbl)
   !> Molecular structure containing atomic numbers.
   type(TMolecule), intent(in) :: mol
   !> Internal-coordinate definitions.
   class(internal_coords_set_type), intent(in) :: internals
   !> Number of internal coordinates.
   integer, intent(in) :: ncoord
   !> Coordinate labels, dimension (ncoord), truncated to 20 characters.
   character(20), intent(out) :: lbl(ncoord)

   integer :: ic, k
   character(80) :: buf

   k = 0
   lbl = "??"

   do ic = 1, internals%ncoords ! bonds
      if (internals%kind(ic) /= coord_bond) cycle
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0)") &
         "bond ", mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
         "-", mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic)
      lbl(k) = buf(1:20)
   end do
   do ic = 1, internals%ncoords ! angular rows
      if (internals%kind(ic) /= coord_angle .and. &
         & internals%kind(ic) /= coord_linbend) cycle
      k = k + 1
      if (internals%kind(ic) == coord_angle) then
         write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0)") &
            "ang ", mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
            "-", mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic), &
            "-", mol%sym(internals%atoms(3, ic)), internals%atoms(3, ic)
      else
         if (is_second_linbend(internals, ic)) then
            buf = "lb2 "
         else
            buf = "lb1 "
         end if
         write(buf(5:), "(a2,i0,a,a2,i0,a,a2,i0)") &
            mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
            "-", mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic), &
            "-", mol%sym(internals%atoms(3, ic)), internals%atoms(3, ic)
      end if
      lbl(k) = buf(1:20)
   end do
   do ic = 1, internals%ncoords ! dihedrals
      if (internals%kind(ic) /= coord_dihedral) cycle
      k = k + 1
      write(buf, "(a,a2,i0,a,a2,i0,a,a2,i0,a,a2,i0)") &
         "dih ", mol%sym(internals%atoms(1, ic)), internals%atoms(1, ic), &
         "-", mol%sym(internals%atoms(2, ic)), internals%atoms(2, ic), &
         "-", mol%sym(internals%atoms(3, ic)), internals%atoms(3, ic), &
         "-", mol%sym(internals%atoms(4, ic)), internals%atoms(4, ic)
      lbl(k) = buf(1:20)
   end do

end subroutine build_labels

!> True for the second row of a consecutive fixed-frame linear-bend pair.
pure logical function is_second_linbend(internals, ic)
   class(internal_coords_set_type), intent(in) :: internals
   integer, intent(in) :: ic

   is_second_linbend = ic > 1
   if (is_second_linbend) then
      is_second_linbend = internals%kind(ic - 1) == coord_linbend .and. &
         & all(internals%atoms(1:3, ic - 1) == internals%atoms(1:3, ic))
   end if

end function is_second_linbend


!> Reciprocal with an explicit positive-infinity result for exact zero.
pure real(wp) function reciprocal(value)
   real(wp), intent(in) :: value

   if (value == 0.0_wp) then
      reciprocal = ieee_value(value, ieee_positive_inf)
   else
      reciprocal = 1.0_wp / value
   end if

end function reciprocal


!> Compliance constants of the reference geometry as the pseudoinverse of the
!> Hessian projected out of translations and rotations.
!>
!> C = B H^+ B^T, evaluated as Z (Z D)^T with Z = B V and D = diag(1/w_i) from
!> the eigen decomposition Hp = V W V^T of the projected, symmetrised Hessian.
!> Fixed-frame near-linear bend pairs may contain rigid-rotation components;
!> projection through H^+ removes those components, so one paired row can have
!> zero response without implying a missing physical bend.
subroutine compute_compliance(unit, H, B, xyz, nat, nint, C, error)
   !> Formatted output unit for the rank diagnostic.
   integer, intent(in) :: unit
   !> Number of atoms.
   integer, intent(in) :: nat
   !> Number of internal coordinates, passed explicitly -- no hardcoded 3N-6.
   integer, intent(in) :: nint
   !> Cartesian Hessian in Hartree/Bohr^2, dimension (3*natoms, 3*natoms).
   real(wp), intent(in) :: H(3*nat, 3*nat)
   !> Wilson B matrix, dimension (nint, 3*natoms).
   real(wp), intent(in) :: B(nint, 3*nat)
   !> Cartesian reference coordinates in Bohr, dimension (3, natoms).
   real(wp), intent(in) :: xyz(3, nat)
   !> Compliance matrix, dimension (nint, nint), in Bohr^2/Hartree.
   real(wp), intent(out) :: C(nint, nint)
   !> Error information.
   type(error_type), allocatable, intent(out) :: error

   integer :: i, j, ndim, lwork, info, nrigid, nvib, rank_h
   real(wp) :: tol_h, normq, center(3), rotation_scale
   ! Relative to the largest rotation norm, sqrt(epsilon) makes the rank test size-independent.
   real(wp), parameter :: rigid_basis_tol = sqrt(epsilon(1.0_wp))
   real(wp), allocatable :: Hp(:, :), Q(:, :), W(:), Z(:, :), ZD(:, :), &
      & T1(:, :), T2(:, :), work(:)
   character(len=64) :: message

   if (nint == 0) then
      C = 0.0_wp
      return
   end if
   ndim = 3 * nat
   allocate(Hp(ndim, ndim), Q(ndim, 6), W(ndim), Z(nint, ndim), &
      & ZD(nint, ndim), T1(ndim, 6), T2(6, ndim))

   ! orthonormal basis of the rigid (translation + rotation) space
   center = sum(xyz, dim=2) / nat
   Q = 0.0_wp
   do i = 1, 3
      Q(i::3, i) = 1.0_wp
   end do
   do j = 1, nat
      Q(3*j-2:3*j, 4) = crossProd([1.0_wp, 0.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 5) = crossProd([0.0_wp, 1.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 6) = crossProd([0.0_wp, 0.0_wp, 1.0_wp], xyz(:, j) - center)
   end do
   rotation_scale = max(norm2(Q(:, 4)), norm2(Q(:, 5)), norm2(Q(:, 6)))

   ! modified Gram-Schmidt; the rotation about the molecular axis of a
   ! linear molecule produces no displacement
   nrigid = 0
   do i = 1, 6
      do j = 1, i - 1
         normq = dot_product(Q(:, i), Q(:, j))
         Q(:, i) = Q(:, i) - normq * Q(:, j)
      end do
      normq = norm2(Q(:, i))
      if (i > 3 .and. normq <= rigid_basis_tol*rotation_scale) then
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
   call blas_gemm("N", "N", ndim, 6, ndim, 1.0_wp, Hp, ndim, Q, ndim, 0.0_wp, T1, ndim)
   call blas_gemm("N", "T", ndim, ndim, 6, -1.0_wp, T1, ndim, Q, ndim, 1.0_wp, Hp, ndim)
   call blas_gemm("T", "N", 6, ndim, ndim, 1.0_wp, Q, ndim, Hp, ndim, 0.0_wp, T2, 6)
   call blas_gemm("N", "N", ndim, ndim, 6, -1.0_wp, Q, ndim, T2, 6, 1.0_wp, Hp, ndim)
   Hp = 0.5_wp * (Hp + transpose(Hp))

   ! Hp = V W V^T, overwriting Hp with V
   lwork = -1
   allocate(work(1))
   call lapack_syev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))
   call lapack_syev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   if (info /= 0) then
      write(message, "(A,I0)") "compute_compliance: DSYEV(H) info=", info
      call fatal_error(error, trim(message))
      return
   end if

   ! Z = B V
   ! ZD = Z D, with D(i,i) = 1/W(i) for retained modes and zero otherwise
   ! C = Z (ZD)^T = B V D V^T B^T = B Hp^+ B^T
   ! Scale epsilon by matrix size and largest |w| to cover eigensolver roundoff.
   tol_h = real(ndim, wp) * epsilon(1.0_wp) * maxval(abs(W))
   rank_h = count(abs(W) > tol_h)
   if (rank_h /= nvib) then
      write(unit, "(A,I0,A,I0)") &
         & "  Note: Hessian rank=", rank_h, " /= 3N-rigid=", nvib
   end if
   call blas_gemm("N", "N", nint, ndim, ndim, 1.0_wp, B, nint, Hp, ndim, 0.0_wp, Z, nint)
   ! Apply D: divide retained modal columns by w_i and drop discarded modes.
   ZD = 0.0_wp
   do i = 1, ndim
      if (abs(W(i)) > tol_h) ZD(:, i) = Z(:, i) / W(i)
   end do
   call blas_gemm("N", "T", nint, nint, ndim, 1.0_wp, Z, nint, ZD, nint, 0.0_wp, C, nint)
   C = 0.5_wp * (C + transpose(C))
   deallocate(Hp, Q, W, Z, ZD, T1, T2, work)

end subroutine compute_compliance


end module xtb_compliance
