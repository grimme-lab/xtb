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

module xtb_ptb_potential
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use tblite_basis_type, only: basis_type
   use xtb_ptb_integral_types, only: integral_type
   use tblite_wavefunction_spin, only: magnet_to_updown
   implicit none
   private

   public :: new_potential, add_pot_to_h1

   !> Container for density dependent potential-shifts
   type, public :: potential_type
      !> Atom-resolved charge-dependent potential shift
      real(wp), allocatable :: vat(:, :)
      !> Shell-resolved charge-dependent potential shift
      real(wp), allocatable :: vsh(:, :)
      !> Orbital-resolved charge-dependent potential shift
      real(wp), allocatable :: vao(:, :)

      !> Atom-resolved dipolar potential
      real(wp), allocatable :: vdp(:, :, :)

      !> Additional H shift as potential 
      real(wp), allocatable :: vaoshift(:, :, :)
   contains
      !> Reset the density dependent potential
      procedure :: reset
   end type potential_type

contains

   !> Create a new potential object
   subroutine new_potential(self, mol, bas, nspin)
      !> Instance of the density dependent potential
      type(potential_type), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Description of the basis set
      type(basis_type), intent(in) :: bas
      !> Number of spin channels
      integer, intent(in) :: nspin

      allocate (self%vat(mol%nat, nspin))
      allocate (self%vsh(bas%nsh, nspin))
      allocate (self%vao(bas%nao, nspin))
      allocate (self%vdp(3, mol%nat, nspin))
      allocate (self%vaoshift(bas%nao, bas%nao, nspin))
   end subroutine new_potential

   subroutine reset(self)
      !> Instance of the density dependent potential
      class(potential_type), intent(inout) :: self

      self%vat(:, :) = 0.0_wp
      self%vsh(:, :) = 0.0_wp
      self%vao(:, :) = 0.0_wp
      self%vdp(:, :, :) = 0.0_wp

      self%vaoshift(:, :, :) = 0.0_wp
   end subroutine reset

   !> Add the collected potential shifts to the effective Hamiltonian
   subroutine add_pot_to_h1(bas, ints, pot, h1)
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Integral container
      type(integral_type), intent(in) :: ints
      !> Density dependent potential-shifts
      type(potential_type), intent(inout) :: pot
      !> Effective Hamiltonian
      real(wp), intent(out) :: h1(:, :, :)

      integer :: spin

      h1(:, :, 1) = ints%hamiltonian
      if (size(h1, 3) > 1) h1(:, :, 2:) = 0.0_wp

      ! call add_vat_to_vsh(bas, pot%vat, pot%vsh)
      ! call add_vsh_to_vao(bas, pot%vsh, pot%vao)
      ! call add_vao_to_h1(bas, ints%overlap, pot%vao, h1)
      ! call add_vmp_to_h1(bas, ints%dipole, pot%vdp, h1)

      h1 = h1 + pot%vaoshift

      call magnet_to_updown(h1)
   end subroutine add_pot_to_h1

end module xtb_ptb_potential
