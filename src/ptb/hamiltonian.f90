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
   use tblite_adjlist, only: adjacency_list
   use tblite_integral_multipole, only : msao

   use xtb_ptb_data, only: THamiltonianData

   implicit none

   private

   public :: get_hamiltonian, get_selfenergy

contains

   subroutine get_hamiltonian(mol, list, bas, overlap, overlap_h0, overlap_xc, &
      & vecp, selfenergies, hamiltonian)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> (Scaled) overlap matrix
      real(wp), intent(in) :: overlap(:, :), overlap_h0(:, :), overlap_xc(:, :)
      !> Effective core potential
      real(wp), intent(in) :: vecp(:, :)
      !> Self-energies
      real(wp), intent(in) :: selfenergies(:)
      !> Effective Hamiltonian
      real(wp), allocatable, intent(out) :: hamiltonian(:, :)
      !> H0 matrix
      real(wp), allocatable :: h0(:, :)

      allocate (hamiltonian(bas%nao, bas%nao), h0(bas%nao, bas%nao), source=0.0_wp)

      hamiltonian = vecp

      call get_h0(mol, list, bas, overlap_h0, selfenergies, h0)

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
      !> tmp reals
      real(wp) :: combinedcn

      allocate (selfenergies(bas%nao), source=0.0_wp)

      do iat = 1, mol%nat
         iid = mol%id(iat)
         combinedcn = (cn_normal(iat) + cn_star(iat)*hData%kcnstar(iid))

         !##### DEV WRITE #####
         ! write(*,*) 'cn_normal', cn_normal(iat)
         ! write(*,*) 'cn_star', cn_star(iat)
         ! write(*,*) 'hData%kcnstar', hData%kcnstar(iid)
         ! write(*,*) "atom: ", iat, "id: ", iid, "type: ", mol%num(iid)
         ! write(*,*) "number of shells: ", bas%nsh_id(iid)
         !#####################

         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(iid)
            selfenergies(ii + ish) = hData%selfEnergy(ish, iid) + &
               & hData%klh(ish, iid)*combinedcn + &
               & cn_star(iat)*hData%kshift(iid)

            !##### DEV WRITE #####
            ! write(*,*) 'selfenergies', ii + ish, selfenergies(ii + ish)
            !#####################

         end do
      end do

   end subroutine get_selfenergy

   subroutine get_h0(mol, list, bas, sh0, levels, h0mat)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> (Scaled) overlap matrix
      real(wp), intent(in) :: sh0(:, :)
      !> Levels
      real(wp), intent(in) :: levels(:)
      !> H0 matrix
      real(wp), intent(out) :: h0mat(:, :)
      !> Loop variables
      integer :: iat, jat, izp, jzp, itr, k, img, inl
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao, ij
      !> Interatomic vector
      real(wp) :: vec(3) = 0.0_wp
      real(wp) :: r2

      !> EXPLANATION:
      !> Loop over all neighbours (iat =/= jat) and set the Hamiltonian for off-center elements
      k=0
      do iat = 1, mol%nat
         write(*,*) "atom: ", iat, "id: ", mol%id(iat), "type: ", mol%num(mol%id(iat))
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = list%inl(iat)
         do img = 1, list%nnl(iat)
            jat = list%nlat(img + inl)
            write(*,*) "neighbour: ", jat, "id: ", mol%id(jat), "type: ", mol%num(mol%id(jat))
            itr = list%nltr(img + inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            k = k + 1
            write(*,*) "iteration: ", k, "distance: ", sqrt(r2)
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js + jsh)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        ij = jao + nao*(iao - 1)
                     end do
                  end do

               end do
            end do

         end do
      end do

      !> EXPLANATION:
      !> Loop over all atoms and set the Hamiltonian for on-center elements

   end subroutine get_h0

end module xtb_ptb_hamiltonian
