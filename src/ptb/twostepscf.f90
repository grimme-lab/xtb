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

module xtb_ptb_scf
   use mctc_io, only: structure_type
   use mctc_env, only: wp

   use tblite_basis_type, only: basis_type

   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_param, only: kalphah0l, nshell, max_shell
   use xtb_ptb_overlaps, only: get_scaled_integrals

   implicit none
   private

   public :: twostepscf

contains

   subroutine twostepscf(mol, bas)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Array of nshells per atom ID
      integer, allocatable :: nsh_id(:)

      real(wp), allocatable :: overlap(:, :), overlap_scaled(:, :)
      real(wp), allocatable :: dipole(:, :, :)

      real(wp), allocatable :: expscal(:, :)

      integer :: i, j, isp, izp

      allocate (expscal(max_shell, mol%nid), source=0.0_wp)
      call get_scaled_integrals(mol, overlap, dipole)

      !##### DEV WRITE #####
      write (*, *) "Overlap:"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f10.6)', advance="no") overlap(i, j)
         end do
         write (*, *) ""
      end do
      !#####################

      !> Set up a new basis set with using the scaled exponents
      nsh_id = nshell(mol%num)
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = kalphah0l(:, izp)
      end do
      call get_scaled_integrals(mol, overlap_scaled, expscal)

      !##### DEV WRITE #####
      write (*, *) "Overlap H0 scaled (SS):"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f10.6)', advance="no") overlap_scaled(i, j)
         end do
         write (*, *) ""
      end do
      !#####################
      ! write (*, *) "Dipole:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f12.6)', advance="no") dipole(1, i, j)
      !    end do
      !    write (*, *) ""
      ! end do

   end subroutine twostepscf

end module xtb_ptb_scf

