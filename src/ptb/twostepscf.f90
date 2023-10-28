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
   use mctc_env, only: wp, error_type

   use xtb_mctc_lapack, only: lapack_syev

   use tblite_basis_type, only: basis_type
   use tblite_context, only: context_type
   use tblite_scf_solver, only: solver_type

   use xtb_ptb_vdzp, only: add_vDZP_basis
   use xtb_ptb_param, only: kalphah0l, klalphaxc, &
   & nshell, max_shell
   use xtb_ptb_overlaps, only: get_scaled_integrals

   implicit none
   private

   public :: twostepscf

contains

   subroutine twostepscf(ctx, mol, bas)
      !> Calculation context
      type(context_type) :: ctx
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Electronic solver
      class(solver_type), allocatable :: solver
      !> Error type
      type(error_type), allocatable :: error

      real(wp), allocatable :: overlap(:, :), overlap_h0(:, :), overlap_xc(:, :)
      real(wp), allocatable :: dipole(:, :, :)

      real(wp), allocatable :: expscal(:, :)

      integer :: i, j, isp, izp

      !> Variables for ML-pop
      real(wp), allocatable :: seig(:)
      integer :: lwork, info
      real(wp), allocatable :: aux(:)

      !> Solver for the effective Hamiltonian
      call ctx%new_solver(solver, bas%nao)

      allocate (expscal(max_shell, mol%nid), source=0.0_wp)
      call get_scaled_integrals(mol, overlap, dipole)

      !##### DEV WRITE #####
      write (*, *) "Standard overlap:"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f10.6)', advance="no") overlap(i, j)
         end do
         write (*, *) ""
      end do
      ! write (*, *) "Dipole:"
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f12.6)', advance="no") dipole(1, i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

      !> Set up a new basis set with using the scaled exponents for H0
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = kalphah0l(:, izp)
      end do
      call get_scaled_integrals(mol, overlap_h0, expscal)
      !##### DEV WRITE #####
      write (*, *) "Overlap H0 scaled (SS):"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f10.6)', advance="no") overlap_h0(i, j)
         end do
         write (*, *) ""
      end do
      !#####################

      !> Set up a new basis set with using the scaled exponents for VXC
      do isp = 1, mol%nid
         izp = mol%num(isp)
         expscal(:, isp) = klalphaxc(:, izp)
      end do
      call get_scaled_integrals(mol, overlap_xc, expscal)
      !##### DEV WRITE #####
      write (*, *) "Overlap XC scaled (SS):"
      do i = 1, bas%nao
         do j = 1, bas%nao
            write (*, '(f10.6)', advance="no") overlap_xc(i, j)
         end do
         write (*, *) ""
      end do

      allocate (aux(1))
      lwork = -1
      call lapack_syev('V', 'U', bas%nao, overlap, bas%nao, seig, aux, lwork, info)
      lwork = int(aux(1))
      deallocate (aux)
      allocate (aux(lwork))
      allocate (seig(bas%nao))
      call lapack_syev('V', 'U', bas%nao, overlap, bas%nao, seig, aux, lwork, info)
      write (*, *) "Eigenvalues:"
      do i = 1, bas%nao
         write (*, '(f10.6)') seig(i)
      end do

   end subroutine twostepscf

end module xtb_ptb_scf

