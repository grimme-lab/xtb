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

!> Pauli exchange approximation
module xtb_ptb_paulixc
   !> mctc-lib
   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type

   use tblite_blas, only: gemm
   use tblite_basis_type, only: basis_type
   use tblite_integral_multipole, only: msao

   use xtb_ptb_vdzp, only: max_shell

   implicit none
   private

   public :: calc_Vxc_pauli

   interface calc_Vxc_pauli
      module procedure calc_Vxc_pauli_shellscaling
      module procedure calc_Vxc_pauli_atomscaling
   end interface calc_Vxc_pauli

contains

   subroutine calc_Vxc_pauli_atomscaling(mol, bas, psh, Sxc, selfenergies, katom, Vxc)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> shell populations
      real(wp), intent(in) :: psh(:)
      !> Overlap matrix scaled for XC potential
      real(wp), intent(in) :: Sxc(:, :)
      !> Self-energies (levels)
      real(wp), intent(in) :: selfenergies(:)
      !> Atom-specific scaling factors
      real(wp), intent(in) :: katom(:)
      !> Vxc potential
      real(wp), allocatable, intent(out) :: Vxc(:, :)

      !> shell-specific scaling factors
      real(wp), allocatable :: kshell(:, :)
      !> Loop indices
      integer :: ii, ish

      allocate (kshell(max_shell, mol%nid), source=0.0_wp)

      !> element scaling in first iter to decouple 1. and 2. iter and
      !> to account for missing P-dependent term of 2. iter
      !> HERE: Projection of atom-wise scaling factors onto shells
      do ii = 1, mol%nid
         do ish = 1, bas%nsh_id(ii)
            kshell(ish, ii) = katom(ii)
         end do
      end do

      call calc_Vxc_pauli(mol, bas, psh, Sxc, selfenergies, kshell, Vxc)
   end subroutine calc_Vxc_pauli_atomscaling

   subroutine calc_Vxc_pauli_shellscaling(mol, bas, psh, Sxc, selfenergies, kshell, Vxc)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> shell populations
      real(wp), intent(in) :: psh(:)
      !> Overlap matrix scaled for XC potential
      real(wp), intent(in) :: Sxc(:, :)
      !> Self-energies (levels)
      real(wp), intent(in) :: selfenergies(:)
      !> Shell-specific scaling factors
      real(wp), intent(in) :: kshell(:, :)
      !> Vxc potential
      real(wp), allocatable, intent(out) :: Vxc(:, :)

      real(wp) :: f1
      real(wp), allocatable :: stmp(:, :)
      integer :: i, jat, jzp, js, jsh, jj, jao, ml

      allocate (stmp(bas%nao, bas%nao), Vxc(bas%nao, bas%nao))

      !> N^2 step
      do i = 1, bas%nao
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            do jsh = 1, bas%nsh_id(jzp) !> Iteration over core shells of atom jat
               jj = bas%iao_sh(js + jsh)
               ml = msao(bas%cgto(jsh, jat)%ang)
               f1 = psh(js + jsh) * kshell(jsh, jzp) / dble(ml) ! shell wise scaling

               !##### DEV WRITE #####
               ! if (i == 1) then
               !    write (*, '(a,i0,a,i0,a,f8.4,a,f8.4,a,f8.4,a,f8.4)') "shell: ", jsh, " atom: ", jat, " f1: ", f1, &
               ! & " psh: ", psh(js + jsh), " cnf2: ", kshell(jsh, jzp), " nl: ", dble(ml)
               ! end if
               !#####################
               do jao = 1, ml
                  stmp(jj + jao, i) = selfenergies(js + jsh) * Sxc(jj + jao, i) * f1
                  !##### DEV WRITE #####
                  ! write(*,*) "stmp: ", jj + jao, i, stmp(jj + jao, i)
                  !#####################
               end do
            end do
         end do
      end do

!     N^3 step
      call gemm(Sxc, stmp, Vxc)

   end subroutine calc_Vxc_pauli_shellscaling

end module xtb_ptb_paulixc
