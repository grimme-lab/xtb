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

!> Two-step SCF of the PTB method

module xtb_ptb_coulomb
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_basis_type, only: basis_type
   use tblite_coulomb_charge_effective, only: harmonic_average

   use dftd4_data_hardness, only: get_hardness

   implicit none
   private

   public :: coulomb_potential

   type :: coulomb_potential
      !> Effective Hubbard parameters
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> effective gams (chemical hardness)
      real(wp), allocatable :: gam(:, :)
      !> Coulomb matrix
      real(wp), allocatable :: cmat(:, :)
      !> Ohno-Klopman contribution
      real(wp) :: cok
   contains
      procedure :: init => init_hubbard
      procedure :: update => get_coulomb_matrix
      procedure :: get_potential
   end type coulomb_potential

contains

   subroutine get_coulomb_matrix(self, mol, bas)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> tmp variables
      integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
      real(wp) :: vec(3), r1, gam, tmp, r12

      !##  !$omp parallel do default(none) schedule(runtime) &
      !##  !$omp shared(amat, mol, nshell, offset, hubbard, gexp) &
      !##  !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, r1, r1g, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jj = bas%ish_at(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r1 = norm2(vec)
            r12 = r1**2
            do ish = 1, bas%nsh_at(iat)
               do jsh = 1, bas%nsh_at(jat)
                  gam = self%hubbard(jsh, ish, jat, iat)
                  tmp = ( self%cok / sqrt(r12 + gam**(-2)) ) + &  !> Ohno-Klopman average 
                     & ( (1.0_wp - self%cok) / ( r1 + gam**(-1) ) ) !> Mataga-Nishimoto average
                  !## !$omp atomic
                  self%cmat(jj + jsh, ii + ish) = self%cmat(jj + jsh, ii + ish) + tmp
                  !## !$omp atomic
                  self%cmat(ii + ish, jj + jsh) = self%cmat(ii + ish, jj + jsh) + tmp
               end do
            end do
         end do
         do ish = 1, bas%nsh_at(iat)
            do jsh = 1, ish - 1
               gam = self%hubbard(jsh, ish, iat, iat)
               !## !$omp atomic
               self%cmat(ii + jsh, ii + ish) = self%cmat(ii + jsh, ii + ish) + gam
               !## !$omp atomic
               self%cmat(ii + ish, ii + jsh) = self%cmat(ii + ish, ii + jsh) + gam
            end do
            !## !$omp atomic
            self%cmat(ii + ish, ii + ish) = self%cmat(ii + ish, ii + ish) + self%hubbard(ish, ish, iat, iat)
         end do
      end do

   end subroutine get_coulomb_matrix

   subroutine init_hubbard(self, mol, bas, q, gamsh, kqhubb, kok)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Shell-wise gamma scaling factors
      real(wp), intent(in) :: gamsh(:,:)
      !> Scaling factor for the charge dependent part of the Hubbard parameters
      real(wp), intent(in) :: kqhubb
      !> Ohno-Klopman contribution
      real(wp), intent(in) :: kok

      integer :: izp, iat, iid, ish, jsh, jat

      !> Initialize variables
      if (.not. allocated(self%gam)) allocate (self%gam(maxval(bas%nsh_id), mol%nat), source=0.0_wp)
      if (.not. allocated(self%hubbard)) then
         allocate (self%hubbard(maxval(bas%nsh_id), maxval(bas%nsh_id), mol%nat, mol%nat), source=0.0_wp)
      end if
      if (.not. allocated(self%cmat)) then
         allocate (self%cmat(bas%nsh, bas%nsh), source=0.0_wp)
      end if
      self%cok = kok

      !> Atom-individual chemical hardnesses per shell; Eq. 19
      do iat = 1, mol%nat
         iid = mol%id(iat)
         izp = mol%num(iid)
         do ish = 1, bas%nsh_id(iid)
            self%gam(ish, iat) = (1.0_wp + kqhubb * q(iat)) * get_hardness(izp) * &
               & gamsh(ish, iid)
         end do
      end do

      !> Effective Hubbard parameters; Eq. 16
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            self%hubbard(:, :, jat, iat) = 0.0_wp
            do ish = 1, bas%nsh_at(iat)
               do jsh = 1, bas%nsh_at(jat)
                  self%hubbard(jsh, ish, jat, iat) = &
                     & harmonic_average(self%gam(ish, iat), self%gam(jsh, jat))
               end do
            end do
         end do
      end do

   end subroutine init_hubbard

   subroutine get_potential(self, mol)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol

   end subroutine get_potential

end module xtb_ptb_coulomb
