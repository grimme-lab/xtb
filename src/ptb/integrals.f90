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

!> Calculation of all required overlap matrices
module xtb_ptb_integrals
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_integral_dipole, only: get_dipole_integrals
   use tblite_cutoff, only: get_lattice_points
   use tblite_basis_type, only: basis_type, get_cutoff
   use tblite_integral_overlap, only: overlap_cgto, maxl, msao
   use tblite_integral_dipole, only: dipole_cgto
   use tblite_adjlist, only: adjacency_list

   use xtb_ptb_vdzp, only: add_vDZP_basis, nshell, max_shell
   implicit none

   private

   public :: get_integrals

   interface get_integrals
      module procedure get_integrals_overlap_existbasis
      module procedure get_integrals_overlap_newbasis
      module procedure get_integrals_dipole_existbasis
      module procedure get_integrals_dipole_newbasis
   end interface get_integrals

contains

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_integrals_overlap_existbasis(mol, bas, lattr, list, overlap, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: lattr(:, :)
      !> Adjacency list
      type(adjacency_list), intent(in) :: list
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      real(wp), allocatable :: normlocal(:)
      integer :: i, j, ij

      allocate (normlocal(bas%nao), source=0.0_wp)

      !> Calculate overlap matrix
      call get_overlap(mol, lattr, list, bas, overlap)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         normlocal(i) = 1.0_wp / sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j) * normlocal(i) * normlocal(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do
      if (present(norm)) then
         norm = normlocal
      end if

   end subroutine get_integrals_overlap_existbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_integrals_overlap_newbasis(mol, lattr, list, overlap, alpha_scal, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: lattr(:, :)
      !> Adjacency list
      type(adjacency_list), intent(in) :: list
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      !> Basis set data
      type(basis_type) :: bas

      !> Set up a new basis set with using the scaled exponents
      if (present(alpha_scal)) then
         call add_vDZP_basis(mol, alpha_scal, bas)
      else
         call add_vDZP_basis(mol, bas)
      end if

      call get_integrals(mol, bas, lattr, list, overlap, norm=norm)

   end subroutine get_integrals_overlap_newbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_integrals_dipole_existbasis(mol, bas, lattr, list, overlap, dipole, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> Adjacency list
      type(adjacency_list), intent(in) :: list
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: lattr(:, :)
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :), dipole(:, :, :)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      real(wp), allocatable :: normlocal(:)
      real(wp) :: cutoff
      integer :: i, j, ij

      cutoff = get_cutoff(bas)
      allocate (normlocal(bas%nao), source=0.0_wp)

      call get_dpint(mol, lattr, list, bas, overlap, dipole)

      !> Normalize overlap
      ij = 0
      do i = 1, bas%nao
         normlocal(i) = 1.0_wp / sqrt(overlap(i, i))
      end do
      ij = 0
      do i = 1, bas%nao
         do j = 1, i
            overlap(i, j) = overlap(i, j) * normlocal(i) * normlocal(j)
            overlap(j, i) = overlap(i, j)
         end do
      end do
      if (present(norm)) then
         norm = normlocal
      end if

      !> ########### AND DIPOLE!! #############

   end subroutine get_integrals_dipole_existbasis

   !> Calculate the overlap matrix for a given scaling factor
   subroutine get_integrals_dipole_newbasis(mol, lattr, list, overlap, dipole, alpha_scal, norm)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: lattr(:, :)
      !> Adjacency list
      type(adjacency_list), intent(in) :: list
      !> Overlap matrix as output
      real(wp), intent(out) :: overlap(:, :), dipole(:, :, :)
      !> Scaling factors as input
      real(wp), intent(in), optional :: alpha_scal(max_shell, mol%nid)
      !> Normalization factors as output
      real(wp), intent(out), optional :: norm(:)
      !> Basis set data
      type(basis_type) :: bas
      integer :: isp, izp, ish, nsh_id(mol%nid)
      !#####################

      !> Set up a new basis set with using the scaled exponents
      if (present(alpha_scal)) then
         call add_vDZP_basis(mol, alpha_scal, bas)
      else
         call add_vDZP_basis(mol, bas)
      end if

      !##### DEV WRITE #####
      ! nsh_id = nshell(mol%num)
      ! write (*, *) "Basis set properties:", bas%nao
      ! do isp = 1, mol%nid
      !    izp = mol%num(isp)
      !    do ish = 1, nsh_id(isp)
      !       write(*,*) bas%cgto(ish, isp)%ang
      !       write(*,*) bas%cgto(ish, isp)%nprim
      !       write(*,*) bas%cgto(ish, isp)%alpha(:)
      !       write(*,*) bas%cgto(ish, isp)%coeff(:)
      !    end do
      ! end do
      !#####################

      call get_integrals(mol, bas, lattr, list, overlap, dipole, norm)

   end subroutine get_integrals_dipole_newbasis

   subroutine get_overlap(mol, trans, list, bas, overlap)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Overlap integral matrix
      real(wp), intent(out) :: overlap(:, :)

      integer :: iat, jat, izp, jzp, itr, img, inl
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
      real(wp) :: r2, vec(3)
      real(wp), allocatable :: stmp(:)

      overlap(:, :) = 0.0_wp

      allocate (stmp(msao(bas%maxl)**2))

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, trans, list, overlap) &
      !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, inl, img)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = list%inl(iat)
         do img = 1, list%nnl(iat)
            jat = list%nlat(img + inl)
            itr = list%nltr(img + inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js + jsh)
                  call overlap_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao * (iao - 1)
                        !$omp atomic
                        overlap(jj + jao, ii + iao) = overlap(jj + jao, ii + iao) &
                                                      + stmp(ij)

                        !$omp atomic
                        overlap(ii + iao, jj + jao) = overlap(ii + iao, jj + jao) &
                                                      + stmp(ij)
                     end do
                  end do

               end do
            end do

         end do
      end do

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, overlap) &
      !$omp private(iat, izp, itr, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         vec(:) = 0.0_wp
         r2 = 0.0_wp
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is + ish)
            do jsh = 1, bas%nsh_id(izp)
               jj = bas%iao_sh(is + jsh)
               call overlap_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

               nao = msao(bas%cgto(jsh, izp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao * (iao - 1)
                     overlap(jj + jao, ii + iao) = overlap(jj + jao, ii + iao) &
                                                   + stmp(ij)
                  end do
               end do

            end do
         end do

      end do

   end subroutine get_overlap

   subroutine get_dpint(mol, trans, list, bas, overlap, dpint)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points within a given realspace cutoff
      real(wp), intent(in) :: trans(:, :)
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Overlap integral matrix
      real(wp), intent(out) :: overlap(:, :)
      !> Dipole moment integral matrix
      real(wp), intent(out) :: dpint(:, :, :)

      integer :: iat, jat, izp, jzp, itr, img, inl, k
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
      real(wp) :: r2, vec(3), dtmpj(3)
      real(wp), allocatable :: stmp(:), dtmpi(:, :)

      overlap(:, :) = 0.0_wp
      dpint(:, :, :) = 0.0_wp

      allocate (stmp(msao(bas%maxl)**2), dtmpi(3, msao(bas%maxl)**2))

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, trans, list, overlap, dpint) &
      !$omp private(iat, jat, izp, jzp, itr, is, js, ish, jsh, ii, jj, iao, jao, nao, ij, k) &
      !$omp private(r2, vec, stmp, dtmpi, dtmpj, inl, img)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = list%inl(iat)
         do img = 1, list%nnl(iat)
            jat = list%nlat(img + inl)
            itr = list%nltr(img + inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js + jsh)
                  call dipole_cgto(bas%cgto(jsh, jzp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp, dtmpi)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao * (iao - 1)
                        call shift_operator(vec, stmp(ij), dtmpi(:, ij), &
                           & dtmpj)
                        !$omp atomic
                        overlap(jj + jao, ii + iao) = overlap(jj + jao, ii + iao) &
                           & + stmp(ij)

                        !$omp atomic
                        overlap(ii + iao, jj + jao) = overlap(ii + iao, jj + jao) &
                           & + stmp(ij)

                        do k = 1, 3
                           !$omp atomic
                           dpint(k, jj + jao, ii + iao) = dpint(k, jj + jao, ii + iao) &
                              & + dtmpi(k, ij)
                        end do
                        do k = 1, 3
                           !$omp atomic
                           dpint(k, ii + iao, jj + jao) = dpint(k, ii + iao, jj + jao) &
                              & + dtmpj(k)
                        end do
                     end do
                  end do

               end do
            end do

         end do
      end do

      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, overlap, dpint) &
      !$omp private(iat, izp, itr, is, ish, jsh, ii, jj, iao, jao, nao, ij) &
      !$omp private(r2, vec, stmp, dtmpi)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         vec(:) = 0.0_wp
         r2 = 0.0_wp
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is + ish)
            do jsh = 1, bas%nsh_id(izp)
               jj = bas%iao_sh(is + jsh)
               call overlap_cgto(bas%cgto(jsh, izp), bas%cgto(ish, izp), &
                     & r2, vec, bas%intcut, stmp)

               nao = msao(bas%cgto(jsh, izp)%ang)
               do iao = 1, msao(bas%cgto(ish, izp)%ang)
                  do jao = 1, nao
                     ij = jao + nao * (iao - 1)
                     overlap(jj + jao, ii + iao) = overlap(jj + jao, ii + iao) &
                        & + stmp(ij)
                     dpint(:, jj + jao, ii + iao) = dpint(:, jj + jao, ii + iao) &
                        & + dtmpi(:, ij)
                  end do
               end do

            end do
         end do

      end do

   end subroutine get_dpint

   !> TAKEN OVER FROM TBLITE
   !> Shift multipole operator from Ket function (center i) to Bra function (center j),
   !> the multipole operator on the Bra function can be assembled from the lower moments
   !> on the Ket function and the displacement vector using horizontal shift rules.
   pure subroutine shift_operator(vec, s, di, dj)
      !> Displacement vector of center i and j
      real(wp), intent(in) :: vec(:)
      !> Overlap integral between basis functions
      real(wp), intent(in) :: s
      !> Dipole integral with operator on Ket function (center i)
      real(wp), intent(in) :: di(:)
      !> Dipole integral with operator on Bra function (center j)
      real(wp), intent(out) :: dj(:)

      ! Create dipole operator on Bra function from Ket function and shift contribution
      ! due to monopol displacement
      dj(1) = di(1) + vec(1) * s
      dj(2) = di(2) + vec(2) * s
      dj(3) = di(3) + vec(3) * s

   end subroutine shift_operator

end module xtb_ptb_integrals
