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

#ifndef WITH_TBLITE
#define WITH_TBLITE 0
#endif

!> Setup of the effective Hamiltonians for both iterations in PTB

module xtb_ptb_hamiltonian
#if WITH_TBLITE
   use mctc_io, only: structure_type
   use mctc_env, only: wp, error_type

   use tblite_basis_type, only: basis_type
   use tblite_adjlist, only: adjacency_list
   use tblite_integral_multipole, only: msao

   use xtb_ptb_data, only: THamiltonianData
   use xtb_ptb_param, only: ptbGlobals

   implicit none

   private

   public :: get_hamiltonian, get_selfenergy, get_occupation

contains

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

      !> debug mode
      logical, parameter :: debug(2) = [ .false., .false. ]

      allocate (selfenergies(bas%nsh), source=0.0_wp)

      do iat = 1, mol%nat
         iid = mol%id(iat)
         combinedcn = (cn_normal(iat) + cn_star(iat) * hData%kcnstar(iid))

         if (debug(1)) then !##### DEV WRITE #####
            write(*,*) 'cn_normal', cn_normal(iat)
            write(*,*) 'cn_star', cn_star(iat)
            write(*,*) 'hData%kcnstar', hData%kcnstar(iid)
            write(*,*) "atom: ", iat, "id: ", iid, "type: ", mol%num(iid)
            write(*,*) "number of shells: ", bas%nsh_id(iid)
         endif

         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(iid)
            selfenergies(ii + ish) = hData%selfEnergy(ish, iid) + &
               & hData%klh(ish, iid) * combinedcn + &
               & cn_star(iat) * hData%kshift(iid)

            if (debug(2)) & !##### DEV WRITE #####
               write(*,*) 'selfenergies', ii + ish, selfenergies(ii + ish)

         end do
      end do

   end subroutine get_selfenergy

   subroutine get_hamiltonian(mol, list, bas, hData, wolfsberg_par, sh0, levels, h0mat, kpol, kitr, kitocod)
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
      !> Wolfsberg parameters [nsh, nid]
      real(wp), intent(in) :: wolfsberg_par(:, :)
      !> Levels
      real(wp), intent(in) :: levels(:)
      !> H0 matrix
      real(wp), intent(out) :: h0mat(:, :)
      !> Level polarization scaling in H0
      real(wp), intent(in) :: kpol
      !> Optional, iteration-dependent parameters:
      !> Radii scaling in H0
      real(wp), optional, intent(in) :: kitr
      !> One-center-off-diagonal scaling in H0
      real(wp), optional, intent(in) :: kitocod
      !> Loop variables
      integer :: iat, jat, izp, jzp, img, inl
      integer :: ish, jsh, is, js, ii, jj, iao, jao, nao
      real(wp) :: wolfsberg, polarized_levels, sum_levels, ocod_param, rscal
      real(wp) :: radii_dependence, ssquraedterm, ocodterm, sterm
      !> Interatomic vector
      real(wp) :: vec(3) = 0.0_wp
      real(wp) :: r2, rab

      !> debug mode
      logical, parameter :: debug(5) = &
               [ .false., .false., .false., .false., .false. ]

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
      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, list, h0mat, sh0, levels, hData, wolfsberg_par, kpol, rscal) &
      !$omp private(iat, jat, izp, jzp, is, js, ish, jsh, ii, jj, iao, jao, nao) &
      !$omp private(r2, vec, inl, img, wolfsberg, polarized_levels, sum_levels) &
      !$omp private(radii_dependence, ssquraedterm, ocodterm, sterm, rab)
      do iat = 1, mol%nat
         
         if (debug(1)) & !##### DEV WRITE #####
            write (*, *) "atom: ", iat, "id: ", mol%id(iat), "type: ", mol%num(mol%id(iat))
         
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         inl = list%inl(iat)
         do img = 1, list%nnl(iat)
            jat = list%nlat(img + inl)
            
            
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            rab = sqrt(r2)

            if (debug(2)) & !##### DEV WRITE #####
               write (*, *) "neighbour: ", jat, "id: ", mol%id(jat), "type: ", mol%num(mol%id(jat)), &
                           & "distance: ", sqrt(r2)
            
            do ish = 1, bas%nsh_id(izp)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_id(jzp)
                  jj = bas%iao_sh(js + jsh)

                  !> Single contributions to H0
                  sum_levels = levels(js + jsh) + levels(is + ish)
                  wolfsberg = wolfsberg_par(ish, izp) + wolfsberg_par(jsh, jzp)
                  polarized_levels = 1.0_wp - kpol * ((levels(js + jsh) - levels(is + ish)) / &
                  & sum_levels)**2
                  radii_dependence = 1.0_wp + &
                  & (hData%kr(izp) + hData%kr(jzp)) * rscal / rab

                  nao = msao(bas%cgto(jsh, jat)%ang)
                  do iao = 1, msao(bas%cgto(ish, iat)%ang)
                     do jao = 1, nao
                        !> Set H0
                        h0mat(jj + jao, ii + iao) = 0.5_wp * sh0(jj + jao, ii + iao) * &
                           & sum_levels * &
                           & wolfsberg * &
                           & polarized_levels * &
                           & radii_dependence
                        h0mat(ii + iao, jj + jao) = h0mat(jj + jao, ii + iao)
                        
                        if (debug(3)) then !##### DEV WRITE #####
                           write(*,'(a,i3,i3,f8.4,f8.4,f8.4,f8.4,f8.4)') "i, j, tmp, keav, pol, xk, ssh: ",&
                           & ii+iao, jj+jao, h0mat(jj + jao, ii + iao), wolfsberg*0.5_wp, &
                           & polarized_levels, radii_dependence, sh0(jj + jao, ii + iao)*sum_levels
                        endif
                     
                     end do
                  end do

               end do
            end do

         end do
      end do

      !> Loop over all atoms and set the Hamiltonian for one-center off(-shell)-diagonal elements
      !$omp parallel do schedule(runtime) default(none) &
      !$omp shared(mol, bas, list, h0mat, sh0, levels, hData, kpol, ocod_param, rscal) &
      !$omp private(iat, izp, is, ish, jsh, ii, jj, iao, jao, nao) &
      !$omp private(sum_levels) &
      !$omp private(ssquraedterm, ocodterm, sterm)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is + ish)
            do jsh = 1, ish - 1

               sum_levels = levels(is + jsh) + levels(is + ish)
               !> s-s', p-p', d-d' off-diagonal, li=lj because S=0 otherwise
               !> CAUTION: ksla is angular-momentum dependent, but the parameters are projected
               !> onto the shells in xtb/src/ptb/data.f90
               ocodterm = hData%ksla(ish, izp) * ocod_param

               jj = bas%iao_sh(is + jsh)
               nao = msao(bas%cgto(jsh, iat)%ang)
               do iao = 1, msao(bas%cgto(ish, iat)%ang)
                  do jao = 1, nao
                     ssquraedterm = sh0(jj + jao, ii + iao)**2 * &
                        & sum_levels * &
                        & hData%kocod(izp) * &
                        & ocodterm
                     sterm = sh0(jj + jao, ii + iao) * &
                        & sum_levels * &
                        & ocodterm
                     h0mat(jj + jao, ii + iao) = sterm + &
                        & ssquraedterm
                     h0mat(ii + iao, jj + jao) = h0mat(jj + jao, ii + iao)
                     
                     if(debug(4)) then !##### DEV WRITE #####
                        write(*,'(a,i3,i3,f8.4,f8.4,f8.4,f8.4,f8.4)') "i, j, tmp, ssquraedterm, sterm, ocodterm: ",&
                        & ii+iao, jj+jao, h0mat(jj + jao, ii + iao), ssquraedterm, sterm, ocodterm, &
                        & sh0(jj + jao, ii + iao) * sum_levels
                     endif
                  
                  end do
               end do
            end do
         end do
      end do

      !> Loop over all atoms and set the Hamiltonian for one-center diagonal elements (same AO!)
      !> CAUTION: Loop doesn't run over one-center diagonal elements of different AOs in the same shell
      !> (Overlap matrix elements are zero for these elements)
      !> NO PARALLELIZATION YET SINCE LOOP IS QUITE SMALL
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            ii = bas%iao_sh(is + ish)
            do iao = 1, msao(bas%cgto(ish, iat)%ang)
               h0mat(ii + iao, ii + iao) = 2.0_wp * levels(is + ish)

               if (debug(5)) & !##### DEV WRITE #####
                  write(*,*) "i, j, tmp: ", ii+iao, ii+iao, h0mat(ii + iao, ii + iao)
               
            end do
         end do
      end do

   end subroutine get_hamiltonian

   subroutine get_occupation(mol, bas, refocc, nocc, n0at, n0sh)
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Basis set information
      type(basis_type), intent(in) :: bas
      !> Reference occupation
      real(wp), intent(in) :: refocc(:, :)
      !> Occupation number
      real(wp), intent(out) :: nocc
      !> Reference occupation for each atom
      real(wp), intent(out) :: n0at(:)
      !> Reference occupation for each shell
      real(wp), intent(out) :: n0sh(:)

      integer :: iat, ish, izp, ii

      nocc = -mol%charge
      n0at(:) = 0.0_wp
      n0sh(:) = 0.0_wp
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ii = bas%ish_at(iat)
         do ish = 1, bas%nsh_id(izp)
            nocc = nocc + refocc(ish, izp)
            n0at(iat) = n0at(iat) + refocc(ish, izp)
            n0sh(ii + ish) = n0sh(ii + ish) + refocc(ish, izp)
         end do
      end do

   end subroutine get_occupation

#endif
end module xtb_ptb_hamiltonian
