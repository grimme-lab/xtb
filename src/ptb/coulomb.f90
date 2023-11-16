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

   use dftd4_data_hardness, only: get_hardness

   implicit none
   private

   public :: coulomb_potential

   type :: coulomb_potential
      !> Effective Hubbard parameters
      real(wp), allocatable :: hubbard(:, :, :, :)
      !> effective gams (chemical hardness)
      real(wp), allocatable :: gam(:)
   contains
      procedure :: init => init_hubbard
      procedure :: update => get_coulomb_matrix
      procedure :: get_potential
   end type coulomb_potential

contains

   subroutine get_coulomb_matrix(self, mol)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol

      ! do i = 1, n
      !    geff(i) = (1d0 + gsc * q(i)) * gam(at(i))
      ! end do

!     ! DFTB second order term J matrix
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

   end subroutine get_coulomb_matrix

   subroutine init_hubbard(self, mol, q, gsc)
      !> Effective Hubbard parameters
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Charge dependent scaling factor
      real(wp), intent(in) :: gsc

      integer :: izp, iat, iid

      allocate (self%gam(mol%nat))
      do iat = 1, mol%nat
         iid = mol%id(iat)
         izp = mol%num(iid)
         self%gam(iat) = (1.0_wp + gsc * q(iat)) * get_hardness(izp)
         write(*,*) 'gam', self%gam(iat)
      end do

   end subroutine init_hubbard

   subroutine get_potential(self, mol)
      !> Coulomb type
      class(coulomb_potential), intent(inout) :: self
      !> Molecular structure
      type(structure_type), intent(in) :: mol

   end subroutine get_potential

end module xtb_ptb_coulomb
