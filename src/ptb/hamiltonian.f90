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
   use tblite_integral_multipole, only: msao

   use xtb_ptb_data, only: THamiltonianData
   use xtb_ptb_param, only: ptbGlobals

   implicit none

   private

   public :: get_hamiltonian, get_selfenergy

contains

   subroutine get_hamiltonian(mol, list, bas, hData, overlap, overlap_h0, overlap_xc, &
      & vecp, selfenergies, iteration, hamiltonian)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> THamiltonian data
      type(THamiltonianData), intent(in) :: hData
      !> (Scaled) overlap matrix
      real(wp), intent(in) :: overlap(:, :), overlap_h0(:, :), overlap_xc(:, :)
      !> Effective core potential
      real(wp), intent(in) :: vecp(:, :)
      !> Self-energies
      real(wp), intent(in) :: selfenergies(:)
      !> Iteration
      integer, intent(in) :: iteration
      !> Effective Hamiltonian
      real(wp), allocatable, intent(out) :: hamiltonian(:, :)

      !> H0 matrix
      real(wp), allocatable :: h0(:, :)
      !> tmp loop variables
      integer :: i, j

      allocate (hamiltonian(bas%nao, bas%nao), h0(bas%nao, bas%nao), source=0.0_wp)

      hamiltonian = vecp

      if (iteration == 1) then
         call get_h0(mol, list, bas, hData, overlap_h0, selfenergies, h0, ptbGlobals%kpol, &
            & ptbGlobals%kitr, ptbGlobals%kitocod)
      else
         call get_h0(mol, list, bas, hData, overlap_h0, selfenergies, h0, ptbGlobals%kpol)
      end if

      !##### DEV WRITE #####
      write (*, *) "H0 ..."
      ! do i = 1, bas%nao
      !    do j = 1, bas%nao
      !       write (*, '(f10.6)', advance="no") h0(i, j)
      !    end do
      !    write (*, *) ""
      ! end do
      !#####################

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

      allocate (selfenergies(bas%nsh), source=0.0_wp)

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

   subroutine get_h0(mol, list, bas, hData, sh0, levels, h0mat, kpol, kitr, kitocod)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Neighbour list
      type(adjacency_list), intent(in) :: list
      !> Basis set data
      type(basis_type), intent(in) :: bas
      !> THamiltonian data
      type(THamiltonianData), intent(in) :: hData
      !> (Scaled) overlap matrix
      real(wp), intent(in) :: sh0(:, :)
      !> Levels
      real(wp), intent(in) :: levels(:)
      !> H0 matrix
      real(wp), intent(out) :: h0mat(:, :)
      !> Level polarization scaling in H0
      real(wp), intent(in) :: kpol
      !> Optional, iteration-dependent parameters:
      !> Radii scaling in H0
      real(wp), optional, intent(in) :: kitr
      !> One-center-off-center distance scaling in H0
      real(wp), optional, intent(in) :: kitocod
      !> Loop variables
      integer :: iat, jat, izp, jzp, itr, k, img, inl
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao
      real(wp) :: wolfsberg, polarized_levels, sum_levels, ocod_param, rscal
      real(wp) :: radii_dependence
      !> Interatomic vector
      real(wp) :: vec(3) = 0.0_wp
      real(wp) :: r2, rab

      if (present(kitocod)) then
         ocod_param = kitocod
      else
         ocod_param = 1.0_wp
      end if

      if (present(kitr)) then
         rscal = kitr
      else
         rscal = 1.0_wp
      end if

      !> Loop over all neighbours (iat =/= jat) and set the Hamiltonian for off-center elements
      k = 0
      do iat = 1, mol%nat
         write (*, *) "atom: ", iat, "id: ", mol%id(iat), "type: ", mol%num(mol%id(iat))
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = list%inl(iat)
         do img = 1, list%nnl(iat)
            jat = list%nlat(img + inl)
            write (*, *) "neighbour: ", jat, "id: ", mol%id(jat), "type: ", mol%num(mol%id(jat))
            itr = list%nltr(img + inl)
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            rab = sqrt(r2)
            k = k + 1
            write (*, *) "iteration: ", k, "distance: ", sqrt(r2)
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js + jsh)

                  nao = msao(bas%cgto(jsh, jzp)%ang)
                  do iao = 1, msao(bas%cgto(ish, izp)%ang)
                     do jao = 1, nao
                        !> Single contributions to H0
                        sum_levels = levels(js+jsh) + levels(is+ish)
                        wolfsberg =  hData%kla(ish, izp) + hData%kla(jsh, jzp)
                        polarized_levels = 1.0_wp - kpol*((levels(js + jsh) - levels(is + ish)) / &
                        & sum_levels )**2
                        radii_dependence = 1.0_wp + & 
                        & (hData%kr(izp) + hData%kr(jzp))*rscal / rab
                        !> Set H0
                        h0mat(jj + jao, ii + iao) = 0.5_wp * sh0(jj + jao, ii + iao) * &
                           & sum_levels * &
                           & wolfsberg * &
                           & polarized_levels * &
                           & radii_dependence
                        h0mat(ii + iao, jj + jao) = h0mat(jj + jao, ii + iao)
                        !##### DEV WRITE #####
                        ! write(*,'(a,i3,i3,f8.4,f8.4,f8.4,f8.4,f8.4)') "i, j, tmp, keav, pol, xk, ssh: ",&
                        ! & ii+iao, jj+jao, h0mat(jj + jao, ii + iao), wolfsberg*0.5_wp, &
                        ! & polarized_levels, radii_dependence, sh0(jj + jao, ii + iao)*sum_levels
                        !#####################
                     end do
                  end do

               end do
            end do

         end do
      end do

      !> Loop over all atoms and set the Hamiltonian for on-center elements

   end subroutine get_h0

end module xtb_ptb_hamiltonian
