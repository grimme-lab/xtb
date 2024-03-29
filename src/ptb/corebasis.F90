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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

module xtb_ptb_corebasis
#if WITH_TBLITE
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use xtb_ptb_data, only: TCorePotentialData
   use xtb_ptb_vdzp, only: new_basis

   use tblite_basis_type, only: cgto_type, basis_type
   use tblite_basis_slater, only: slater_to_gauss
   use tblite_integral_overlap, only: overlap_cgto, msao
   use tblite_blas, only: gemm

   implicit none
   private

   public :: add_core_basis, get_Vecp

   real(wp), parameter :: cutoff = 20.0_wp

contains

   subroutine get_Vecp(mol, ecpdata, bas, cbas, norm_s, vecp)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Effective core potential data
      type(TCorePotentialData), intent(in) :: ecpdata
      !> Basis set type and core basis set type
      type(basis_type), intent(in) :: bas, cbas
      !> Core-valence overlap matrix
      real(wp), intent(in) :: norm_s(:)
      !> Approximated effective core potential
      real(wp), allocatable, intent(out) :: vecp(:, :)
      !> Core-valence overlap matrix
      real(wp), allocatable :: overlap_cv(:, :)
      !> Intermediate core valence overlap matrix (scaled)
      real(wp), allocatable :: secptmp(:, :)
      !> tmp indices
      integer :: i, jat, jzp, js, jsh, jj, jao, jati, j

      logical, parameter :: debug(2) = [ .false., .false. ]

      !> Calculate the core-valence overlap matrix
      call core_valence_overlap(mol, bas, cbas, norm_s, overlap_cv)

      allocate (vecp(bas%nao, bas%nao), secptmp(cbas%nao, bas%nao), source=0.0_wp)
      !> SG: N^2 step
      do i = 1, bas%nao
         do jat = 1, mol%nat
            jzp = mol%id(jat)
            jati = mol%num(jzp)
            js = cbas%ish_at(jat)
            do jsh = 1, cbas%nsh_id(jzp) !> Iteration over core shells of atom jat
               jj = cbas%iao_sh(js + jsh)
               do jao = 1, msao(cbas%cgto(jsh, jat)%ang)

                  if (debug(1)) then !##### DEV WRITE #####
                     write (*, *) "jsh, hflev(jsh,jati), kecpepsilon(jati): ", jsh, ecpdata%hflev(jsh, jati), ecpdata%kecpepsilon(jati)
                     write (*, *) "overlap_cv(jj + jao, i): ", overlap_cv(jj + jao, i)
                     write (*, *) "jj + jao, i: ", jj + jao, i
                  endif

                  secptmp(jj + jao, i) = -ecpdata%hflev(jsh, jati) * overlap_cv(jj + jao, i) * &
                     & ecpdata%kecpepsilon(jati)
               end do
            end do
         end do
      end do

      if (debug(2)) then !##### DEV WRITE #####
         write (*, *) "Overlap_CV_scaled:"
         do i = 1, bas%nao
            do j = 1, cbas%nao
               write (*, '(f10.6)', advance="no") secptmp(j, i)
            end do
            write (*, *) ""
         end do
      endif

      !> SG: N^3 step
      call gemm(overlap_cv, secptmp, vecp, 'T', 'N')

   end subroutine get_Vecp

   subroutine core_valence_overlap(mol, bas, cbas, bas_overlap_norm, cv_overlap)
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

      logical, parameter :: debug(2) = [ .false., .false. ]

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
                  call overlap_cgto(cbas%cgto(jsh, jat), bas%cgto(ish, iat), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(cbas%cgto(jsh, jat)%ang)
                  do iao = 1, msao(bas%cgto(ish, iat)%ang)
                     do jao = 1, msao(cbas%cgto(jsh, jat)%ang)
                        cv_overlap(jj + jao, ii + iao) = cv_overlap(jj + jao, ii + iao) &
                           & + stmp(jao + nao * (iao - 1))
                     end do
                  end do

               end do
            end do

         end do
      end do

      if (debug(1)) then !##### DEV WRITE #####
         write (*, *) "Core-valence overlap matrix:"
         write (*, *) "core basis NAOs: ", cbas%nao
         write (*, *) "valence basis NAOs: ", bas%nao
      endif

      do i = 1, bas%nao
         do j = 1, cbas%nao
            cv_overlap(j, i) = cv_overlap(j, i) * bas_overlap_norm(i)
            
            if(debug(2)) & !##### DEV WRITE #####
               write (*, '(f10.6)', advance="no") cv_overlap(j, i)
         
         end do

         if(debug(2)) & !##### DEV WRITE #####
            write (*, *) ""

      end do

   end subroutine core_valence_overlap

   subroutine add_core_basis(mol, ecpdata, cbas)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Effective core potential PTB parameters
      type(TCorePotentialData), intent(in) :: ecpdata
      !> Basis set type
      type(basis_type), intent(out) :: cbas

      integer :: isp, izp, ish, stat, il, iat
      integer, allocatable :: nsh_id(:)
      type(cgto_type), allocatable :: cgto(:, :)

      !> debug mode
      logical, parameter :: debug(4) = [ .false., .false., .false., .false. ]
      integer :: j, k

      nsh_id = ecpdata%nshell(mol%num)
      allocate (cgto(maxval(nsh_id), mol%nat))
      do iat = 1, mol%nat
         isp = mol%id(iat)
         izp = mol%num(isp)

         if (debug(1)) & !##### DEV WRITE #####
            write (*, *) "number of shells: ", nsh_id(isp)
         
         do ish = 1, nsh_id(isp)
            il = ecpdata%angshell(ish, izp)

            if (debug(2)) then !##### DEV WRITE #####
               write (*, *) "shell: ", ish
               write (*, *) "number of primitives: ", ecpdata%max_prim
               write (*, *) "shell type: ", ecpdata%pqn(ish, izp), ecpdata%angshell(ish, izp), ecpdata%sl_exp(ish, izp)
            endif

            call slater_to_gauss(ecpdata%max_prim, ecpdata%pqn(ish, izp), il, &
            & ecpdata%sl_exp(ish, izp), cgto(ish, iat), .true., stat)

            if (debug(3)) then !##### DEV WRITE #####
               write (*, *) "N_prim: ", cgto(ish, iat)%nprim
               do k = 1, cgto(ish, iat)%nprim
                  write (*, *) cgto(ish, iat)%alpha(k)
               end do
            endif
         end do

      end do

      if(debug(4)) then !##### DEV WRITE #####
         write (*, *) "---------FINAL PTB core basis-----------"
         do iat = 1, mol%nat
            isp = mol%id(iat)
            write (*, *) "Atom :", mol%num(isp)
            write (*, *) "Number of shells :", nsh_id(isp)
            do j = 1, nsh_id(isp)
               write (*, *) "N_prim: ", cgto(j, iat)%nprim
               write (*, *) "Exponents:"
               write (*, *) cgto(j, iat)%alpha
               write (*, *) "Coefficients:"
               write (*, *) cgto(j, iat)%coeff
            end do
         end do
         write (*, *) "----------------------------------------"
      endif

      call new_basis(cbas, mol, nsh_id, cgto, 1.0_wp)

   end subroutine add_core_basis

#endif
end module xtb_ptb_corebasis

