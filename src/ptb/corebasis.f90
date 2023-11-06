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

!> Core basis functions relevant for approximated effective core potential in PTB

module xtb_ptb_corebasis
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use xtb_ptb_param, only: max_elem, highest_elem, &  !> General parameters
      & max_core_shell, max_core_prim, cbas_nshell, &  !> PTB core basis parameters
      & cbas_pqn, cbas_sl_exp, cbas_angshell !> PTB core basis parameters

   use tblite_basis_type, only: cgto_type, new_basis, basis_type
   use tblite_basis_slater, only: slater_to_gauss
   use tblite_integral_overlap, only: overlap_cgto, msao

   implicit none
   private

   public :: add_PTBcore_basis, core_valence_overlap

   real(wp), parameter :: cutoff = 20.0_wp

contains

   subroutine core_valence_overlap(mol, bas, cbas, cv_overlap, bas_overlap_norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(in) :: bas, cbas
      !> Core-valence overlap matrix
      real(wp), allocatable, intent(out) :: cv_overlap(:, :)
      !> Normalization factors for valence basis functions
      real(wp), intent(in) :: bas_overlap_norm(:)

      integer :: iat, jat, izp, jzp, is, js, maxl
      integer :: ish, jsh, ii, jj, iao, jao, nao
      real(wp) :: r2, vec(3), cutoff2

      integer :: i, j

      real(wp), allocatable :: stmp(:)

      cutoff2 = cutoff**2

      maxl = max(bas%maxl, cbas%maxl)

      allocate (cv_overlap(cbas%nao, bas%nao), stmp(msao(maxl)**2), source=0.0_wp)

      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            js = cbas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            if (r2 > cutoff2) cycle
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, cbas%nsh_id(jzp)
                  jj = cbas%iao_sh(js + jsh)
                  call overlap_cgto(cbas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)
                  write(*,*) "overlapelements stmp: ", stmp

                  nao = msao(cbas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, msao(cbas%cgto(jsh, jzp)%ang)
                        cv_overlap(jj + jao, ii + iao) = cv_overlap(jj + jao, ii + iao) &
                           & + stmp(jao + nao*(iao - 1))
                     end do
                  end do

               end do
            end do

         end do
      end do

      !##### DEV WRITE #####
      write (*, *) "Core-valence overlap matrix:"
      write(*,*) "core basis NAOs: ", cbas%nao
      write(*,*) "valence basis NAOs: ", bas%nao
      do i = 1, bas%nao
         do j = 1, cbas%nao
            write (*, '(f10.6)', advance="no") cv_overlap(j, i) * bas_overlap_norm(i)
         end do
         write (*, *) ""
      end do

   end subroutine core_valence_overlap

   subroutine add_PTBcore_basis(mol, cbas)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(out) :: cbas

      integer :: isp, izp, ish, stat, il
      integer, allocatable :: nsh_id(:)
      type(cgto_type), allocatable :: cgto(:, :)

      integer :: j, k

      nsh_id = cbas_nshell(mol%num)
      allocate (cgto(maxval(nsh_id), mol%nid))
      do isp = 1, mol%nid
         izp = mol%num(isp)
         !##### DEV WRITE #####
         write (*, *) "number of shells: ", nsh_id(isp)
         !#####################
         do ish = 1, nsh_id(isp)
            il = cbas_angshell(ish, izp)

            !##### DEV WRITE #####
            write (*, *) "number of primitives: ", max_core_prim
            write (*, *) "shell type: ", cbas_pqn(ish, izp), cbas_angshell(ish, izp), cbas_sl_exp(ish, izp)
            !#####################

            call slater_to_gauss(max_core_prim, cbas_pqn(ish, izp), il, &
            & cbas_sl_exp(ish, izp), cgto(ish, isp), .true., stat)

            !##### DEV WRITE #####
            write (*, *) "N_prim: ", cgto(ish, isp)%nprim
            do k = 1, cgto(ish, isp)%nprim
               write (*, *) cgto(ish, isp)%alpha(k)
            end do
            !#####################
         end do

      end do

      !##### DEV WRITE #####
      write (*, *) "---------FINAL PTB core basis-----------"
      do isp = 1, mol%nid
         write (*, *) "Atom :", mol%num(isp)
         write (*, *) "Number of shells :", nsh_id(isp)
         do j = 1, nsh_id(isp)
            write (*, *) "N_prim: ", cgto(j, isp)%nprim
            write (*, *) "Exponents:"
            write (*, *) cgto(j, isp)%alpha
            write (*, *) "Coefficients:"
            write (*, *) cgto(j, isp)%coeff
         end do
      end do
      write (*, *) "----------------------------------------"

      call new_basis(cbas, mol, nsh_id, cgto, 1.0_wp)

   end subroutine add_PTBcore_basis

end module xtb_ptb_corebasis

