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

!> Z-matrix internal coordinates of a fixed reference geometry and their
!> analytic Wilson B matrix.
!>
!> Connectivity, reference atoms, coordinate kinds and the perpendicular
!> frames of linear centres are determined once from the reference geometry
!> by init and stored in a TZMatrix; values and bmatrix then evaluate the
!> coordinate vector q and its gradient at any geometry.
module xtb_zmat_type
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_math, only : crossProd
   use xtb_basic_geo, only : bangl
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, linbend_frame, &
      & bmat_torsion
   implicit none
   private

   public :: COORD_BOND, COORD_ANGLE, COORD_DIHEDRAL, COORD_LINBEND
   public :: TZMatrix, init

   !> Angle to a linear centre: i--na(i)--j counts as linear when the angle
   !> deviates from pi by less than this (radians, 5 degrees)
   real(wp), parameter :: lin_tol = 0.08726646_wp
   !> Smallest admissible reference angle for nb and nc (radians, 20 degrees)
   real(wp), parameter :: min_angle = 0.34906585_wp
   !> Window around 0 and pi in which a dihedral reference is degenerate (radians)
   real(wp), parameter :: dihedral_tol = 0.2617994_wp

   !> Atom contributes only its bond, one coordinate
   integer, parameter :: COORD_BOND = 1
   !> Atom contributes its bond and the valence angle, two coordinates
   integer, parameter :: COORD_ANGLE = 2
   !> Atom contributes its bond, the valence angle and the dihedral, three coordinates
   integer, parameter :: COORD_DIHEDRAL = 3
   !> Atom contributes its bond at a linear centre and two Decius bends, three coordinates
   integer, parameter :: COORD_LINBEND = 4

   !> Z-matrix internal coordinates of a fixed reference geometry.
   !>
   !> Atom i (i >= 2) contributes the bond i--na(i) and, depending on ctype(i),
   !> the valence angle i-na(i)-nb(i), the dihedral i-na(i)-nb(i)-nc(i), or a
   !> pair of Decius linear bends at the linear centre na(i) measured along the
   !> fixed frames e1(:,i) / e2(:,i).  Its coordinates occupy q(qoff(i)+1 ...)
   !> and the same rows of the B matrix.
   type :: TZMatrix
      !> Number of atoms, the array size of every integer component
      integer :: n = 0
      !> Number of internal coordinates, the length of q and the row count of B
      integer :: nint = 0
      !> Reference atoms of each atom: bond distance na(:), angle nb(:) and
      !> dihedral nc(:), zero where the coordinate is not used
      integer, allocatable :: na(:), nb(:), nc(:)
      !> Coordinate kind contributed by each atom, one of the COORD_* parameters
      integer, allocatable :: ctype(:)
      !> Zero-based offset of the first coordinate of each atom in q and in the
      !> rows of the B matrix
      integer, allocatable :: qoff(:)
      !> Fixed perpendicular frame e1(:), e2(:) of a linear-bend centre, zero
      !> for every other atom
      real(wp), allocatable :: e1(:, :), e2(:, :)
   contains
      !> Internal coordinate values at a geometry, q(self%nint)
      procedure :: get_coords => zmat_values
      !> Analytic Wilson B matrix at a geometry, bmat(self%nint, 3*self%n)
      procedure :: get_bmatrix => zmat_bmatrix
   end type TZMatrix

   !> Build the TZMatrix of a reference geometry
   interface init
      module procedure :: initZMatrix
   end interface init

contains

!> Build the Z-matrix of a reference geometry: connectivity, reference atoms,
!> coordinate kinds, the fixed linear-bend frames and the coordinate count.
pure subroutine initZMatrix(self, n, at, xyz)
   !> Z-matrix container to fill, its allocatable components are allocated here
   type(TZMatrix), intent(out) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Atomic numbers of the atoms, dimension (n)
   integer, intent(in) :: at(n)
   !> Cartesian reference coordinates in Bohr, dimension (3, n)
   real(wp), intent(in) :: xyz(3, n)

   logical :: bonded(n, n), inqueue(n), is_linbend
   integer :: i, j, k, q, order(n), parent_arr(n), head, tail, placed
   real(wp) :: r, rmin, angl, u(3), thr
   real(wp), parameter :: bond_f = 1.3_wp

   allocate(self%na(n), self%nb(n), self%nc(n), self%ctype(n), self%qoff(n))
   allocate(self%e1(3, n), self%e2(3, n))
   self%n = n

   ! Step 1a: bond table from covalent radii
   ! ponytail: local Alvarez radii; swap to ncoord_erf + approx_bonds when the pinned compliance
   ! reference values in test_hessian can be regenerated
   bonded = .false.
   do i = 1, n
      do j = i + 1, n
         thr = bond_f * (rcov(at(i)) + rcov(at(j)))
         r = sqrt( (xyz(1, i) - xyz(1, j))**2 &
                 + (xyz(2, i) - xyz(2, j))**2 &
                 + (xyz(3, i) - xyz(3, j))**2 )
         if (r < thr) then
            bonded(i, j) = .true.
            bonded(j, i) = .true.
         end if
      end do
   end do

   ! Step 1b: BFS on bond graph to build Z-matrix ordering + NA array.
   !   order(p) = atom at Z-matrix position p; parent(i) = NA(i).
   !   BFS guarantees every atom connects to an already-listed atom.
   !   Disconnected atoms (metals with no covalent bonds) appended last.
   order = 0; parent_arr = 0; inqueue = .false.
   head = 1; tail = 1; placed = 0
   order(1) = 1; inqueue(1) = .true.; placed = 1; tail = 2

   bfs_loop: do while (head < tail)
      i = order(head); head = head + 1
      do j = 1, n
         if (inqueue(j)) cycle
         if (.not. bonded(i, j)) cycle
         order(tail) = j; parent_arr(j) = i
         inqueue(j) = .true.; tail = tail + 1; placed = placed + 1
      end do
   end do bfs_loop

   ! Append unreached atoms (fragments / isolated metals)
   do i = 1, n
      if (inqueue(i)) cycle
      rmin = 1.0e30_wp; k = order(1)
      do j = 1, placed
         r = (xyz(1, i) - xyz(1, order(j)))**2 + (xyz(2, i) - xyz(2, order(j)))**2 &
           + (xyz(3, i) - xyz(3, order(j)))**2
         if (r < rmin) then; rmin = r; k = order(j); end if
      end do
      order(tail) = i; parent_arr(i) = k
      inqueue(i) = .true.; tail = tail + 1; placed = placed + 1
   end do

   self%na = 0
   do i = 2, n; self%na(i) = parent_arr(i); end do

   ! NB(i): best angle-reference atom (not NA(i), not i itself).
   !   Start with NA(NA(i)); validate angle > min_angle.
   !   Fall back to a best_nb() search if needed.
   ! NC(i): similarly from NA(NB(i)).
   self%nb = 0; self%nc = 0

   do i = 3, n
      k = 0
      if (self%na(i) > 0) k = self%na(self%na(i))
      if (k > 0 .and. k /= i) then
         call bangl(xyz, i, self%na(i), k, angl)
         if (angl < min_angle .or. pi - angl < lin_tol*0.5_wp) k = 0
      else
         k = 0
      end if
      if (k == 0) k = best_nb(xyz, i, self%na(i))
      self%nb(i) = k
   end do

   do i = 4, n
      k = 0
      if (self%nb(i) > 0) k = self%na(self%nb(i))
      if (k > 0 .and. k /= i .and. k /= self%na(i) .and. k /= self%nb(i)) then
         ! Validate NC: the dihedral i-NA-NB-NC is only well-defined if
         ! NC is not collinear with the NA-NB axis, i.e. the angle
         ! NC--NB--NA must not be near 0 or 180 degrees.
         call bangl(xyz, k, self%nb(i), self%na(i), angl)
         if (angl < min_angle .or. pi - angl < lin_tol) k = 0
      else
         k = 0
      end if
      if (k == 0) then
         ! Search ALL atoms (not just j<i) for a valid NC: j not collinear
         ! with the NA-NB axis.  Prefer nearest atom with good angle.
         rmin = 1.0e30_wp
         do j = 1, n
            if (j == i .or. j == self%na(i) .or. j == self%nb(i)) cycle
            call bangl(xyz, j, self%nb(i), self%na(i), angl)
            if (angl < min_angle .or. pi - angl < lin_tol) cycle
            r = (xyz(1, i) - xyz(1, j))**2 + (xyz(2, i) - xyz(2, j))**2 + (xyz(3, i) - xyz(3, j))**2
            if (r < rmin) then; rmin = r; k = j; end if
         end do
         ! Last-resort: relax angle filter, take nearest non-degenerate atom
         if (k == 0) then
            rmin = 1.0e30_wp
            do j = 1, n
               if (j == i .or. j == self%na(i) .or. j == self%nb(i)) cycle
               call bangl(xyz, j, self%nb(i), self%na(i), angl)
               if (pi - angl < lin_tol*0.5_wp) cycle ! exclude exactly linear
               r = (xyz(1, i) - xyz(1, j))**2 + (xyz(2, i) - xyz(2, j))**2 &
                 + (xyz(3, i) - xyz(3, j))**2
               if (r < rmin) then; rmin = r; k = j; end if
            end do
         end if
      end if
      ! Resolve the dihedral reference once, at the reference geometry: a
      ! reference atom that changed with the geometry would make the B row
      ! inconsistent with the coordinate it differentiates.  A reference the
      ! ladder above could not resolve at all (k = 0, e.g. an exactly linear
      ! chain) is left alone rather than re-searched with a degenerate axis.
      if (k > 0) then
         call bangl(xyz, self%na(i), self%nb(i), k, angl)
         if (angl > pi - dihedral_tol .or. angl < dihedral_tol) then
            k = find_dihedral_atom(xyz, i, self%na(i), self%nb(i), k, dihedral_tol)
         end if
      end if
      self%nc(i) = k
   end do

   self%na(1) = 0; self%nb(1) = 0; self%nc(1) = 0
   self%nb(2) = 0; self%nc(2) = 0
   if (n >= 3) self%nc(3) = 0

   ! Step 2: coordinate types and perpendicular frames.
   !
   ! LINBEND detection: atom i gets LINBEND if the chemical angle
   ! i -- na(i) -- X is near 180 deg for X = any bonded neighbour of
   ! na(i) other than i itself.  This is independent of nb(i) (which
   ! was chosen to avoid near-linear triples).
   ! For globally linear centres all bonded pairs of na(i) are ~180 deg.
   self%e1 = 0.0_wp; self%e2 = 0.0_wp
   self%ctype(1) = 0
   do i = 2, n
      self%ctype(i) = COORD_BOND
      if (i < 3) cycle

      ! Check if na(i) is a linear centre w.r.t. atom i
      ! by scanning bonded neighbours of na(i)
      is_linbend = .false.
      do j = 1, n
         if (j == i) cycle
         if (.not. bonded(self%na(i), j)) cycle
         call bangl(xyz, i, self%na(i), j, angl)
         if (pi - angl < lin_tol) then
            is_linbend = .true.
            exit
         end if
      end do

      if (is_linbend) then
         self%ctype(i) = COORD_LINBEND
         u = xyz(:, i) - xyz(:, self%na(i))
         u = u / sqrt(dot_product(u, u))
         call linbend_frame(u, self%e1(:, i), self%e2(:, i))
      else if (i < 4) then
         self%ctype(i) = COORD_ANGLE
      else
         self%ctype(i) = COORD_DIHEDRAL
      end if
   end do

   ! Step 3: offset array and nint
   q = 0; self%qoff(1) = 0
   do i = 2, n
      self%qoff(i) = q
      select case (self%ctype(i))
         case(COORD_BOND)     ; q = q + 1
         case(COORD_ANGLE)    ; q = q + 2
         case(COORD_DIHEDRAL) ; q = q + 3
         case(COORD_LINBEND)  ; q = q + 3
      end select
   end do
   self%nint = q

end subroutine initZMatrix


!> Covalent radius in Bohr for atomic number Z (Alvarez 2008)
pure function rcov(iz) result(r)
   !> Atomic number, radii are tabulated for Z = 1..86
   integer, intent(in) :: iz

   real(wp) :: r
   real(wp), parameter :: ang2bohr = 1.0_wp / 0.529177210903_wp
   ! Covalent radii in Angstrom, Z=1..86
   real(wp), parameter :: rc(86) = [ &
      0.31_wp, 0.28_wp, 1.28_wp, 0.96_wp, 0.84_wp, 0.73_wp, 0.71_wp, 0.66_wp, & !  1  H-O
      0.57_wp, 0.58_wp, 1.66_wp, 1.41_wp, 1.21_wp, 1.11_wp, 1.07_wp, 1.05_wp, & !  9  F-S
      1.02_wp, 1.06_wp, 2.03_wp, 1.76_wp, 1.70_wp, 1.60_wp, 1.53_wp, 1.39_wp, & ! 17  Cl-Cr
      1.61_wp, 1.52_wp, 1.50_wp, 1.24_wp, 1.32_wp, 1.22_wp, 1.22_wp, 1.20_wp, & ! 25  Mn-Ge
      1.19_wp, 1.20_wp, 1.20_wp, 1.16_wp, 2.20_wp, 1.95_wp, 1.90_wp, 1.75_wp, & ! 33  As-Zr
      1.64_wp, 1.54_wp, 1.47_wp, 1.46_wp, 1.42_wp, 1.39_wp, 1.45_wp, 1.44_wp, & ! 41  Nb-Cd
      1.42_wp, 1.39_wp, 1.39_wp, 1.38_wp, 1.39_wp, 1.40_wp, 2.44_wp, 2.15_wp, & ! 49  In-Ba
      2.07_wp, 2.04_wp, 2.03_wp, 2.01_wp, 1.99_wp, 1.98_wp, 1.98_wp, 1.96_wp, & ! 57  La-Gd
      1.94_wp, 1.92_wp, 1.92_wp, 1.89_wp, 1.90_wp, 1.87_wp, 1.87_wp, 1.75_wp, & ! 65  Tb-Hf
      1.70_wp, 1.62_wp, 1.51_wp, 1.44_wp, 1.41_wp, 1.36_wp, 1.36_wp, 1.32_wp, & ! 73  Ta-Hg
      1.45_wp, 1.46_wp, 1.48_wp, 1.40_wp, 1.50_wp, 1.50_wp ] ! 81  Tl-Rn
   if (iz >= 1 .and. iz <= 86) then
      r = rc(iz) * ang2bohr
   else
      r = 2.0_wp
   end if
end function rcov


!> Internal coordinate values, q(self%nint): lengths in the unit of xyz,
!> angles and dihedrals in radians.  The linear-bend frames are the ones
!> stored at build time and must not be recomputed here.
pure subroutine zmat_values(self, xyz, q)
   !> Coordinate definitions of the reference geometry
   class(TZMatrix), intent(in) :: self
   !> Cartesian coordinates, dimension (3, self%n), in the length unit of the
   !> reference geometry
   real(wp), intent(in) :: xyz(:, :)
   !> Internal coordinate values, dimension (self%nint)
   real(wp), intent(out) :: q(:)

   integer :: i, j, k, l, o
   real(wp) :: vi(3), vk(3), ri, rk

   q = 0.0_wp

   do i = 2, self%n
      j = self%na(i); k = self%nb(i); l = self%nc(i)
      o = self%qoff(i) ! 0-based offset => q(o+1) is first coord of atom i

      ! Bond (always)
      q(o+1) = dist(xyz, i, j)

      select case (self%ctype(i))

      case(COORD_BOND)
         ! nothing more

      case(COORD_ANGLE)
         call bangl(xyz, i, j, k, q(o+2))

      case(COORD_DIHEDRAL)
         call bangl(xyz, i, j, k, q(o+2))
         q(o+3) = dihedral_value(xyz(:, i), xyz(:, j), xyz(:, k), xyz(:, l))

      case(COORD_LINBEND)
         ! Two Decius linear bending coordinates at centre j=NA(i)
         ! vi = unit vector j->i,  vk = unit vector j->k=NB(i)
         vi = xyz(:, i) - xyz(:, j)
         ri = sqrt(dot_product(vi, vi)); vi = vi / ri
         vk = xyz(:, k) - xyz(:, j)
         rk = sqrt(dot_product(vk, vk)); vk = vk / rk
         q(o+2) = dot_product(self%e1(:, i), vi) + dot_product(self%e1(:, i), vk)
         q(o+3) = dot_product(self%e2(:, i), vi) + dot_product(self%e2(:, i), vk)

      end select
   end do

end subroutine zmat_values


!> Signed dihedral i-j-k-l in (-pi, pi], the coordinate whose gradient
!> bmat_torsion returns
pure function dihedral_value(xi, xj, xk, xl) result(tau)
   !> Cartesian coordinates of the four dihedral atoms i-j-k-l, dimension (3) each
   real(wp), intent(in) :: xi(3), xj(3), xk(3), xl(3)

   real(wp) :: tau, b1(3), b2(3), b3(3), n1(3), n2(3), m(3)

   b1 = xj - xi; b2 = xk - xj; b3 = xl - xk
   n1 = crossProd(b1, b2)
   n2 = crossProd(b2, b3)
   m = crossProd(n1, b2/norm2(b2))
   tau = atan2(dot_product(m, n2), dot_product(n1, n2))
end function dihedral_value


!> Analytic Wilson B matrix, bmat(self%nint, 3*self%n), assembled from the
!> xtb_bmatrix row functions.  Row qoff(i)+1 is the bond i--na(i); an ANGLE
!> adds qoff(i)+2, a DIHEDRAL qoff(i)+2 and +3, a LINBEND the two Decius
!> bends at qoff(i)+2 and +3.
pure subroutine zmat_bmatrix(self, xyz, bmat)
   !> Coordinate definitions of the reference geometry
   class(TZMatrix), intent(in) :: self
   !> Cartesian coordinates, dimension (3, self%n), in the length unit of the
   !> reference geometry
   real(wp), intent(in) :: xyz(:, :)
   !> Wilson B matrix, dimension (self%nint, 3*self%n), zeroed then filled
   real(wp), intent(out) :: bmat(:, :)

   integer :: i, j, k, l, o
   real(wp) :: vec_i(3), vec_k(3), brow(6), b9(9)
   real(wp) :: tors_xyz(3, 4), bt(3, 4)

   bmat = 0.0_wp

   do i = 2, self%n
      o = self%qoff(i) ! 0-based offset => row o+1 is the bond
      j = self%na(i)
      k = self%nb(i)
      l = self%nc(i)

      ! Bond row i--j: atoms [i, j]
      vec_i = xyz(:, i) - xyz(:, j)
      brow = bmat_bond(vec_i)
      bmat(o+1, 3*(i-1)+1:3*i) = brow(1:3)
      bmat(o+1, 3*(j-1)+1:3*j) = brow(4:6)

      select case (self%ctype(i))

      case(COORD_BOND)
         ! nothing more

      case(COORD_ANGLE)
         ! Angle i-j-k: atoms [i, j(centre), k]
         vec_k = xyz(:, k) - xyz(:, j)
         b9 = bmat_angle(vec_i, vec_k)
         bmat(o+2, 3*(i-1)+1:3*i) = b9(1:3)
         bmat(o+2, 3*(j-1)+1:3*j) = b9(4:6)
         bmat(o+2, 3*(k-1)+1:3*k) = b9(7:9)

      case(COORD_DIHEDRAL)
         ! Angle row as above, then torsion row i-j-k-l.
         vec_k = xyz(:, k) - xyz(:, j)
         b9 = bmat_angle(vec_i, vec_k)
         bmat(o+2, 3*(i-1)+1:3*i) = b9(1:3)
         bmat(o+2, 3*(j-1)+1:3*j) = b9(4:6)
         bmat(o+2, 3*(k-1)+1:3*k) = b9(7:9)
         tors_xyz(:, 1) = xyz(:, i)
         tors_xyz(:, 2) = xyz(:, j)
         tors_xyz(:, 3) = xyz(:, k)
         tors_xyz(:, 4) = xyz(:, l)
         bt = bmat_torsion(tors_xyz)
         bmat(o+3, 3*(i-1)+1:3*i) = bt(:, 1)
         bmat(o+3, 3*(j-1)+1:3*j) = bt(:, 2)
         bmat(o+3, 3*(k-1)+1:3*k) = bt(:, 3)
         bmat(o+3, 3*(l-1)+1:3*l) = bt(:, 4)

      case(COORD_LINBEND)
         ! Two Decius linear bends at centre j: atoms [i, j(centre), k]
         vec_k = xyz(:, k) - xyz(:, j)
         b9 = bmat_linbend(vec_i, vec_k, self%e1(:, i))
         bmat(o+2, 3*(i-1)+1:3*i) = b9(1:3)
         bmat(o+2, 3*(j-1)+1:3*j) = b9(4:6)
         bmat(o+2, 3*(k-1)+1:3*k) = b9(7:9)
         b9 = bmat_linbend(vec_i, vec_k, self%e2(:, i))
         bmat(o+3, 3*(i-1)+1:3*i) = b9(1:3)
         bmat(o+3, 3*(j-1)+1:3*j) = b9(4:6)
         bmat(o+3, 3*(k-1)+1:3*k) = b9(7:9)

      end select
   end do

end subroutine zmat_bmatrix


!> Distance between atoms i and j
pure function dist(xyz, i, j) result(r)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atom indices of the pair
   integer, intent(in) :: i, j

   real(wp) :: r
   r = sqrt((xyz(1, i) - xyz(1, j))**2 + (xyz(2, i) - xyz(2, j))**2 + (xyz(3, i) - xyz(3, j))**2)
end function dist


!> Replacement dihedral reference atom l for a degenerate dihedral i-j-k-l.
!>
!> Searches atoms 1..ii-1 other than j and k for the closest atom to k whose
!> angle j-k-atom stays strictly inside (tol_in, pi - tol_in); if none is
!> found, repeats the search with a 5 degree tolerance.
pure function find_dihedral_atom(xyz, ii, j, k, l_def, tol_in) result(lb)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Number of atoms, the search covers atoms 1..ii-1
   integer, intent(in) :: ii
   !> First two dihedral reference atoms, excluded from the search
   integer, intent(in) :: j, k
   !> Reference atom returned when the search finds no replacement
   integer, intent(in) :: l_def
   !> Minimum distance of the angle j-k-atom from 0 and pi (radians)
   real(wp), intent(in) :: tol_in

   integer :: lb
   integer :: i1
   real(wp) :: r, rmin, angl, tol
   tol = tol_in; rmin = 100.0_wp; lb = l_def
   do i1 = 1, ii - 1
      if (i1 == j .or. i1 == k) cycle
      r = (xyz(1, i1) - xyz(1, k))**2 + (xyz(2, i1) - xyz(2, k))**2 + (xyz(3, i1) - xyz(3, k))**2
      if (r < rmin) then
         call bangl(xyz, j, k, i1, angl)
         if (angl < pi - tol .and. angl > tol) then; rmin = r; lb = i1; end if
      end if
   end do
   if (rmin > 99.0_wp) then
      tol = 0.087266_wp; rmin = 100.0_wp
      do i1 = 1, ii - 1
         if (i1 == j .or. i1 == k) cycle
         r = (xyz(1, i1) - xyz(1, k))**2 + (xyz(2, i1) - xyz(2, k))**2 + (xyz(3, i1) - xyz(3, k))**2
         if (r < rmin) then
            call bangl(xyz, j, k, i1, angl)
            if (angl < pi - tol .and. angl > tol) then; rmin = r; lb = i1; end if
         end if
      end do
   end if
end function find_dihedral_atom


!> Atom j in 1..i-1 (j /= na_i, j /= i) that gives the largest angle
!> j-na_i-i, subject to angle > min_angle and < pi - lin_tol/2.
!> Falls back to the geometrically nearest atom when no angle qualifies;
!> returns 0 only when no candidate exists at all.
pure function best_nb(xyz, i, na_i) result(k)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atom that needs an angle-reference atom and its bond reference na(i)
   integer, intent(in) :: i, na_i

   integer :: k, j
   real(wp) :: angl, best, rmin, r
   k = 0; best = 0.0_wp
   ! First pass: largest angle in (min_angle, pi-lin_tol/2)
   do j = 1, i - 1
      if (j == i .or. j == na_i) cycle
      call bangl(xyz, i, na_i, j, angl)
      if (angl < min_angle .or. pi - angl < lin_tol*0.5_wp) cycle
      if (angl > best) then; best = angl; k = j; end if
   end do
   ! Fallback: geometrically nearest (no angle filter)
   if (k == 0) then
      rmin = 1.0e30_wp
      do j = 1, i - 1
         if (j == i .or. j == na_i) cycle
         r = (xyz(1, i) - xyz(1, j))**2 + (xyz(2, i) - xyz(2, j))**2 + (xyz(3, i) - xyz(3, j))**2
         if (r < rmin) then; rmin = r; k = j; end if
      end do
   end if
end function best_nb

end module xtb_zmat_type
