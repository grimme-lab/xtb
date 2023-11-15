! This file is part of xtb.
!
! Copyright (C) 2023 Marcel Mueller
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

!> Calculation of all required overlap matrices
module xtb_ptb_overlaps
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_integral_overlap, only: get_overlap
   use tblite_integral_dipole, only: get_dipole_integrals
   use tblite_cutoff, only: get_lattice_points
   use tblite_basis_type, only: basis_type, get_cutoff

   use xtb_ptb_vdzp, only: add_vDZP_basis, nshell, max_shell
   implicit none

   private

   public :: get_scaled_integrals

   interface get_scaled_integrals
      module procedure get_scaled_integrals_overlap_existbasis
      module procedure get_scaled_integrals_overlap_newbasis
      module procedure get_scaled_integrals_dipole_existbasis
      module procedure get_scaled_integrals_dipole_newbasis
   end interface get_scaled_integrals

contains

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_overlap_existbasis(mol, bas, overlap, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      real(wp), allocatable :: normlocal(:)
      real(wp) :: cutoff
      real(wp), allocatable :: lattr(:, :)
      integer :: i, j, ij

      allocate (normlocal(bas%nao), source=0.0_wp)
      !> Calculate cutoff and lattice points (with PTB generally turned off)
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

      !> Calculate overlap matrix
      call get_overlap(mol, lattr, cutoff, bas, overlap)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         normlocal(i) = 1.0_wp / sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j) * normlocal(i) * normlocal(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do
      if (present(norm)) then
         norm = normlocal
      end if

   end subroutine get_scaled_integrals_overlap_existbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_overlap_newbasis(mol, overlap, alpha_scal, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      !> Basis set data
      type(basis_type) :: bas

      !> Set up a new basis set with using the scaled exponents
      if (present(alpha_scal)) then
         call add_vDZP_basis(mol, alpha_scal, bas)
      else
         call add_vDZP_basis(mol, bas)
      end if

      call get_scaled_integrals(mol, bas, overlap, norm=norm)


   end subroutine get_scaled_integrals_overlap_newbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_dipole_existbasis(mol, bas, overlap, dipole, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :), dipole(:, :, :)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      real(wp), allocatable :: normlocal(:)
      real(wp) :: cutoff
      real(wp), allocatable :: lattr(:, :)
      integer :: i, j, ij

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      allocate (normlocal(bas%nao), source=0.0_wp)

      call get_dipole_integrals(mol, lattr, cutoff, bas, overlap, dipole)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         normlocal(i) = 1.0_wp / sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j) * normlocal(i) * normlocal(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do
      if (present(norm)) then
         norm = normlocal
      end if

      !> ########### AND DIPOLE!! #############

   end subroutine get_scaled_integrals_dipole_existbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_dipole_newbasis(mol, overlap, dipole, alpha_scal, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :), dipole(:, :, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      !> Basis set data
      type(basis_type) :: bas
      integer :: isp, izp, ish, nsh_id(mol%nid)
      !#####################

      !> Set up a new basis set with using the scaled exponents
      if (present(alpha_scal)) then
         call add_vDZP_basis(mol, alpha_scal, bas)
      else
         call add_vDZP_basis(mol, bas)
      end if

      !##### DEV WRITE #####
      ! nsh_id = nshell(mol%num)
      ! write (*, *) "Basis set properties:", bas%nao
      ! do isp = 1, mol%nid
      !    izp = mol%num(isp)
      !    do ish = 1, nsh_id(isp)
      !       write(*,*) bas%cgto(ish, isp)%ang
      !       write(*,*) bas%cgto(ish, isp)%nprim
      !       write(*,*) bas%cgto(ish, isp)%alpha(:)
      !       write(*,*) bas%cgto(ish, isp)%coeff(:)
      !    end do
      ! end do
      !#####################

      call get_scaled_integrals(mol, bas, overlap, dipole, norm)

   end subroutine get_scaled_integrals_dipole_newbasis
end module xtb_ptb_overlaps
