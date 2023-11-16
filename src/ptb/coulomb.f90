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

   use xtb_ptb_data, only: TCoulombData

   implicit none
   private

   public :: coulomb_potential

   type :: coulomb_potential
      !> Effective Hubbard parameters
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> effective gams (chemical hardness)
      real(wp), allocatable :: gam(:, :)
      !> Coulomb matrix
      real(wp), allocatable :: coulomb_mat(:, :)
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
      integer :: iat, jat, izp, jzp, ii, jj, ish, jsh
      real(wp) :: vec(3), r1, r1g, gam, tmp

!     DFTB second order term J matrix
      ! ii = 0
      ! do i = 1, n
      !    ati = at(i)
      !    do ish = 1, bas_nsh(ati)
      !       ii = ii + 1
      !       gish = shell_resp(ish, ati, 2) * geff(i) ! important higher-order effect
      !       jj = 0
      !       do j = 1, n
      !          k = lin(j, i)
      !          r2 = rab(k)**2
      !          atj = at(j)
      !          do jsh = 1, bas_nsh(atj)
      !             jj = jj + 1
      !             if (jj > ii) cycle
      !             gjsh = shell_resp(jsh, atj, 2) * geff(j)
      !             xk = 2d0 / (1d0 / gish + 1d0 / gjsh)    ! harm. av.
      !             gab(jj, ii) = cok / sqrt(r2 + 1d0 / xk**2) + cmn / (rab(k) + 1d0 / xk) !Ohno-Klopman-Mataga average
      !             gab(ii, jj) = gab(jj, ii)
      !          end do
      !       end do
      !    end do
      ! end do

      !##  !$omp parallel do default(none) schedule(runtime) &
      !##  !$omp shared(amat, mol, nshell, offset, hubbard, gexp) &
      !##  !$omp private(iat, izp, ii, ish, jat, jzp, jj, jsh, gam, vec, r1, r1g, tmp)
      ! do iat = 1, mol%nat
      !    izp = mol%id(iat)
      !    ii = offset(iat)
      !    do jat = 1, iat - 1
      !       jzp = mol%id(jat)
      !       jj = offset(jat)
      !       vec = mol%xyz(:, jat) - mol%xyz(:, iat)
      !       r1 = norm2(vec)
      !       r1g = r1**gexp
      !       do ish = 1, nshell(iat)
      !          do jsh = 1, nshell(jat)
      !             gam = hubbard(jsh, ish, jzp, izp)
      !             tmp = 1.0_wp / (r1g + gam**(-gexp))**(1.0_wp / gexp)
      !             !$omp atomic
      !             amat(jj + jsh, ii + ish) = amat(jj + jsh, ii + ish) + tmp
      !             !$omp atomic
      !             amat(ii + ish, jj + jsh) = amat(ii + ish, jj + jsh) + tmp
      !          end do
      !       end do
      !    end do
      !    do ish = 1, nshell(iat)
      !       do jsh = 1, ish - 1
      !          gam = hubbard(jsh, ish, izp, izp)
      !          !$omp atomic
      !          amat(ii + jsh, ii + ish) = amat(ii + jsh, ii + ish) + gam
      !          !$omp atomic
      !          amat(ii + ish, ii + jsh) = amat(ii + ish, ii + jsh) + gam
      !       end do
      !       !$omp atomic
      !       amat(ii + ish, ii + ish) = amat(ii + ish, ii + ish) + hubbard(ish, ish, izp, izp)
      !    end do
      ! end do

   end subroutine get_coulomb_matrix

   subroutine init_hubbard(self, mol, bas, q, cdata)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Coulomb parameterization data
      type(TCoulombData), intent(in) :: cdata

      integer :: izp, iat, iid, isp, jsp, ish, jsh, jat

      if (.not. allocated(self%gam)) allocate (self%gam(maxval(bas%nsh_id), mol%nat), source=0.0_wp)
      if (.not. allocated(self%hubbard)) then
         allocate (self%hubbard(maxval(bas%nsh_id), maxval(bas%nsh_id), mol%nat, mol%nat), source=0.0_wp)
      end if
      !> Atom-individual chemical hardnesses per shell; Eq. 19
      do iat = 1, mol%nat
         iid = mol%id(iat)
         izp = mol%num(iid)
         do ish = 1, bas%nsh_id(iid)
            self%gam(ish, iat) = (1.0_wp + cdata%kQHubbard * q(iat)) * get_hardness(izp) * &
               & cdata%shellHardnessFirstIter(ish, iid)
         end do
      end do

      !> Effective Hubbard parameters; Eq. 16
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            self%hubbard(:, :, jat, iat) = 0.0_wp
            do ish = 1, bas%nsh_at(isp)
               do jsh = 1, bas%nsh_at(jsp)
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
