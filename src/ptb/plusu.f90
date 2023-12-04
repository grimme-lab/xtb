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

!> DFT+U approximation within PTB
!> Consumes parameters, the density matrix and the coordination numbers;
!> returns the DFT+U potential

module xtb_ptb_plusu
   !> mctc-lib
   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type
   !> tblite-lib
   use tblite_basis_type, only: basis_type
   use tblite_integral_multipole, only: msao
   !> xtb_ptb-lib
   use xtb_ptb_data, only: TPlusU

   implicit none
   private

   public :: plusu_potential_type

   type :: plusu_potential_type
      !> DFT+U self energies for each shell
      real(wp), allocatable :: selfenergies(:, :)
      !> DFT+U interatomic damping matrix
      real(wp), allocatable :: damping(:, :)
      !##### Parameters #####
      !> Scaling of the diagonal elements of the +U nao x nao matrix
      real(wp), allocatable :: diag_scaling(:)
   contains
      procedure :: get_potential => calc_V_plusU
      procedure :: init
   end type plusu_potential_type

contains

   subroutine init(self, plusudata, mol, bas, q, cn)
      !> Initialize the data for the DFT+U potential
      !> Self class
      class(plusu_potential_type), intent(inout) :: self
      !> DFT+U parameterization data
      type(TPlusU), intent(in) :: plusudata
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Coordination numbers
      real(wp), intent(in) :: cn(:)
      integer :: iat, iid

      !> Allocate the self energies and damping matrix
      allocate (self%selfenergies(maxval(bas%nsh_id), mol%nat))
      allocate (self%damping(mol%nat, mol%nat))
      allocate (self%diag_scaling(mol%nat))

      !> Initialize the parameters with the data from the parameterization
      do iat = 1, mol%nat
         iid = mol%id(iat)
         self%diag_scaling(iat) = plusudata%cud(iid)
      end do

      call init_damping_matrix(self, mol, cn, plusudata%avcn, plusudata%ar, plusudata%arcn)
      call init_selfenergies(self, mol, bas, q, &
         & plusudata%cu1, plusudata%cu2, plusudata%cueffl)

   end subroutine init

   subroutine init_damping_matrix(self, mol, cn, average_cn, atomic_radii, atomic_radii_cn)
      !> Self type
      class(plusu_potential_type), intent(inout) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> Coordination numbers
      real(wp), intent(in) :: cn(:)
      !> Average coordination numbers
      real(wp), intent(in) :: average_cn(:)
      !> Atomic radii
      real(wp), intent(in) :: atomic_radii(:)
      !> Atomic radii for coordination numbers
      real(wp), intent(in) :: atomic_radii_cn(:)
      !> Loop variables
      integer :: iat, jat, izp, jzp
      real(wp) :: ra, rb, rij

      do iat = 1, mol%nat
         izp = mol%id(iat)
         ra = atomic_radii(izp) + (cn(iat) - average_cn(izp)) * atomic_radii_cn(izp)
         do jat = 1, iat
            jzp = mol%id(jat)
            rb = atomic_radii(jzp) + (cn(jat) - average_cn(jzp)) * atomic_radii_cn(jzp)
            rij = norm2(mol%xyz(:, iat) - (mol%xyz(:, jat)))
            self%damping(jat, iat) = 0.5_wp * (1.0_wp + erf(-1.8_wp * (rij - (ra + rb)) / (ra + rb)))
            self%damping(iat, jat) = self%damping(jat, iat)
            !##### DEV WRITE #####
            ! write (*, *) 'damping', iat, jat, self%damping(iat, jat)
            !#####################
         end do
      end do

   end subroutine init_damping_matrix

   subroutine init_selfenergies(self, mol, bas, q, q_scal, q2_scal, shell_level)
      !> Self type
      class(plusu_potential_type), intent(inout) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(in) :: bas
      !> Atomic charges
      real(wp), intent(in) :: q(:)
      !> Scaling of the first order charge
      real(wp), intent(in) :: q_scal(:)
      !> Scaling of the second order charge
      real(wp), intent(in) :: q2_scal(:)
      !> Shell level
      real(wp), intent(in) :: shell_level(:, :)
      !> Loop variables
      integer :: iat, iid, ish

      do iat = 1, mol%nat
         iid = mol%id(iat)
         do ish = 1, bas%nsh_at(iat)
            self%selfenergies(ish, iat) = shell_level(ish, iid) * &
               & (1.0_wp - (q_scal(iid) * q(iat) + q2_scal(iid) * q(iat)**2))
            !##### DEV WRITE #####
            ! write (*, *) 'selfenergies', ish, iat, self%selfenergies(ish, iat)
            !#####################
         end do
      end do
   end subroutine init_selfenergies

   subroutine calc_V_plusU(self, mol, bas, density, potential)
      !> DFT+U approximation within PTB
      !> Self class
      class(plusu_potential_type), intent(inout) :: self
      !> Structure type
      type(structure_type), intent(in) :: mol
      !> Basis set type
      type(basis_type), intent(in) :: bas
      !> Density matrix
      real(wp), intent(in) :: density(:, :)
      !> Potential matrix (nao x nao)
      real(wp), intent(inout) :: potential(:, :)
      !> Loop variables
      integer :: iat, jat, izp, jzp, ii, ish, jsh, jj
      integer :: is, js, iao, jao
      real(wp) :: sum_levels, tmp

      !> Different atom
      do iat = 1, mol%nat
         izp = mol%id(iat)
         is = bas%ish_at(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            js = bas%ish_at(jat)
            do ish = 1, bas%nsh_at(iat)
               ii = bas%iao_sh(is + ish)
               do jsh = 1, bas%nsh_at(jat)
                  jj = bas%iao_sh(js + jsh)
                  sum_levels = self%selfenergies(ish, iat) + self%selfenergies(jsh, jat)
                  do iao = 1, msao(bas%cgto(ish, iat)%ang)
                     do jao = 1, msao(bas%cgto(jsh, jat)%ang)
                        tmp = density(jj + jao, ii + iao) * sum_levels * self%damping(jat, iat)
                        potential(jj + jao, ii + iao) = potential(jj + jao, ii + iao) + tmp
                        potential(ii + iao, jj + jao) = potential(ii + iao, jj + jao) + tmp
                     end do
                  end do
               end do
            end do
         end do
         !> Same atom, different shell
         do ish = 1, bas%nsh_at(iat)
            ii = bas%iao_sh(is + ish)
            do jsh = 1, ish - 1
               jj = bas%iao_sh(is + jsh)
               sum_levels = self%selfenergies(ish, iat) + self%selfenergies(jsh, iat)
               do iao = 1, msao(bas%cgto(ish, iat)%ang)
                  do jao = 1, msao(bas%cgto(jsh, iat)%ang)
                     tmp = density(jj + jao, ii + iao) * sum_levels * self%damping(iat, iat)
                     potential(jj + jao, ii + iao) = potential(jj + jao, ii + iao) + tmp
                     potential(ii + iao, jj + jao) = potential(ii + iao, jj + jao) + tmp
                  end do
               end do
            end do
            !> Same atom, same shell, different AO
            sum_levels = 2.0_wp * self%selfenergies(ish, iat)
            do iao = 1, msao(bas%cgto(ish, iat)%ang)
               do jao = 1, iao - 1
                  tmp = density(ii + jao, ii + iao) * sum_levels * self%damping(iat, iat)
                  potential(ii + jao, ii + iao) = potential(ii + jao, ii + iao) + tmp
                  potential(ii + iao, ii + jao) = potential(ii + iao, ii + jao) + tmp
               end do
               !> Same atom, same shell, same AO
               potential(ii + iao, ii + iao) = potential(ii + iao, ii + iao) + &
                  & density(ii + iao, ii + iao) * sum_levels * self%damping(iat, iat) &
                  & * self%diag_scaling(iat)
            end do
         end do
      end do

   end subroutine calc_V_plusU

end module xtb_ptb_plusu
