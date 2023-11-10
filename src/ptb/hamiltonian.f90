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

!> Setup of the effective Hamiltonians for both iterations in PTB

module xtb_ptb_hamiltonian
   use mctc_io, only: structure_type
   use mctc_env, only: wp, error_type

   use tblite_basis_type, only: basis_type

   use xtb_ptb_data, only: THamiltonianData

   implicit none

   private

   public :: get_hamiltonian, get_selfenergy

contains

   subroutine get_hamiltonian(mol, bas, overlap, overlap_h0, overlap_xc, vecp, hamiltonian)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> (Scaled) overlap matrix
      real(wp), intent(in) :: overlap(:, :), overlap_h0(:, :), overlap_xc(:, :)
      !> Effective core potential
      real(wp), intent(in) :: vecp(:, :)
      !> Effective Hamiltonian
      real(wp), allocatable, intent(out) :: hamiltonian(:, :)

      allocate (hamiltonian(bas%nao, bas%nao))

   end subroutine get_hamiltonian

   subroutine get_selfenergy(mol, bas, hData, cn_normal, cn_star, selfenergies)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> THamiltonian data
      type(THamiltonianData), intent(in) :: hData
      !> Coordination number with fitted (normal) and literature (star) radii
      real(wp), intent(in) :: cn_normal(:), cn_star(:)
      !> Self-energies after taking into account CN dependencies
      real(wp), allocatable, intent(out) :: selfenergies(:)
      !> Loop variables
      integer :: iat, iid, ii, ish

      allocate (selfenergies(bas%nao), source=0.0_wp)

      do iat = 1, mol%nat
         iid = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_at(iid)
            selfenergies(ii + ish) = hData%selfEnergy(ish, iid)
            write(*,*) 'selfenergies', ii + ish, selfenergies(ii + ish)
         end do
      end do

   end subroutine get_selfenergy

end module xtb_ptb_hamiltonian
