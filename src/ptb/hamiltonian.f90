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

   implicit none

   private

   public :: get_hamiltonian

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

end module xtb_ptb_hamiltonian
