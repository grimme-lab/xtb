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

   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_param, only: nshell, max_shell
   implicit none

   private

   public :: get_scaled_integrals

   interface get_scaled_integrals
      module procedure get_scaled_integrals_overlap
      module procedure get_scaled_integrals_dipole
   end interface get_scaled_integrals

contains

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_overlap(mol, overlap, alpha_scal)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Overlap matrix as output
      real(wp), intent(out), allocatable :: overlap(:, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Basis set data
      type(basis_type) :: bas
      real(wp) :: cutoff
      real(wp), allocatable :: lattr(:, :)
      integer :: i, j, ij
      real(wp), allocatable :: norm(:)

      !> Set up a new basis set with using the scaled exponents
      if (present(alpha_scal)) then
         call add_vDZP_basis(mol, alpha_scal, bas)
      else
         call add_vDZP_basis(mol, bas)
      end if

      !> Allocate overlap matrix based on basis set dimension
      allocate (overlap(bas%nao, bas%nao), norm(bas%nao), source=0.0_wp)

      !> Calculate cutoff and lattice points (with PTB generally turned off)
      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

      !> Calculate overlap matrix
      call get_overlap(mol, lattr, cutoff, bas, overlap)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         norm(i) = 1.0_wp/sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j)*norm(i)*norm(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do

   end subroutine get_scaled_integrals_overlap

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_scaled_integrals_dipole(mol, overlap, dipole, alpha_scal)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Overlap matrix as output
      real(wp), intent(out), allocatable :: overlap(:, :), dipole(:, :, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Basis set data
      type(basis_type) :: bas
      real(wp) :: cutoff
      real(wp), allocatable :: lattr(:, :)
      integer :: i, j, ij
      real(wp), allocatable :: norm(:)
      !##### DEV WRITE #####
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

      allocate (overlap(bas%nao, bas%nao), dipole(3, bas%nao, bas%nao), &
      & source=0.0_wp)

      cutoff = get_cutoff(bas)
      call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
      allocate (norm(bas%nao), source=0.0_wp)

      call get_dipole_integrals(mol, lattr, cutoff, bas, overlap, dipole)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         norm(i) = 1.0_wp/sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j)*norm(i)*norm(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do

      !> ########### AND DIPOLE!! #############

   end subroutine get_scaled_integrals_dipole
end module xtb_ptb_overlaps
