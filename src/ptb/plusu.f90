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

      !> Allocate the self energies and damping matrix
      allocate(self%selfenergies(maxval(bas%nsh_id), mol%nat))
      allocate(self%damping(mol%nat, mol%nat))
      allocate(self%diag_scaling(mol%nid))

      !> Initialize the parameters with the data from the parameterization
      self%diag_scaling = plusudata%cud

      call init_damping_matrix()
      call init_selfenergies(self, mol, bas, q, &
         & plusudata%cu1, plusudata%cu2, plusudata%cueffl)

   end subroutine init

   subroutine init_damping_matrix()

      ! do i = 1, n
      !    hi = shell_cnf3(10, at(i)) + (cn(i) - avcn(at(i))) * shell_resp(10, at(i), 1)
      !    do j = 1, i
      !       k = k + 1
      !       r = hi + shell_cnf3(10, at(j)) + (cn(j) - avcn(at(j))) * shell_resp(10, at(j), 1)  ! sum of special radii
      !       t8 = (rab(k) - r) / r
      !       xab(k) = 0.5_wp * (1_wp + erf(-1.8_wp * t8)) ! parameter not important due to redundancy with shell_cnf3(10,)
      !    end do
      ! end do
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
            write(*, *) 'selfenergies', ish, iat, self%selfenergies(ish, iat)
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
      real(wp), intent(out) :: potential(:, :)
      !> Loop variables
      integer :: iat, jat, izp, jzp, ii, ish, jsh, jj
      integer :: is, js, iao, jao
      real(wp) :: sum_levels

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
                        potential(ii + iao, jj + jao) = potential(ii + iao, jj + jao) + &
                           & density(ii + iao, jj + jao) * sum_levels * self%damping(jat, iat)
                     end do
                  end do
               end do
            end do
         end do
         !> Same atom, different shell
         do ish = 1, bas%nsh_at(iat)
            ii = bas%iao_sh(is + ish)
            do jsh = 1, ish - 1
               jj = bas%iao_sh(js + jsh)
               sum_levels = self%selfenergies(ish, iat) + self%selfenergies(jsh, iat)
               do iao = 1, msao(bas%cgto(ish, iat)%ang)
                  do jao = 1, msao(bas%cgto(jsh, iat)%ang)
                     potential(ii + iao, jj + jao) = potential(ii + iao, jj + jao) + &
                        & density(ii + iao, jj + jao) * sum_levels * self%damping(iat, iat)
                  end do
               end do
            end do
            !> Same atom, same shell, different AO
            sum_levels = 2.0_wp * self%selfenergies(ish, iat)
            do iao = 1, msao(bas%cgto(ish, iat)%ang)
               do jao = 1, iao - 1
                  potential(ii + iao, ii + jao) = potential(ii + iao, ii + jao) + &
                     & density(ii + iao, ii + jao) * sum_levels * self%damping(iat, iat)
               end do
               !> Same atom, same shell, same AO
               potential(ii + iao, ii + iao) = potential(ii + iao, ii + iao) + &
                  & density(ii + iao, ii + iao) * sum_levels * self%damping(iat, iat) &
                  & * self%diag_scaling(iat)
            end do
         end do
      end do

      !> OLD CODE
      ! do i = 1, bas%nao
      !    do jat = 1, mol%nat
      !       jzp = mol%id(jat)
      !       js = bas%ish_at(jat)
      !       do jsh = 1, bas%nsh_id(jzp) !> Iteration over core shells of atom jat
      !          jj = bas%iao_sh(js + jsh)
      !          ! ia = aoat(i)
      !          ! ati = at(ia)
      !          ! ish = shell2ao(i)
      !          ! hi = shell_cnf2(ish, ati) * gq(ia)
      !          do j = 1, i - 1
      !             k = k + 1
      !             ib = aoat(j)
      !             atj = at(ib)
      !             jsh = shell2ao(j)
      !             hj = shell_cnf2(jsh, atj) * gq(ib)
      !             !                            this part is INDO two-c like
      !             Hmat(k) = Hmat(k) + P(k) * (hi + hj) * xab(lin(ib, ia))
      !          end do
      !          k = k + 1
      !          Hmat(k) = Hmat(k) + 2d0 * P(k) * shell_xi(8, ati) * hi * xab(lin(ia, ia)) ! scaled diag
      !       end do
      !    end do
      ! end do

   end subroutine calc_V_plusU

end module xtb_ptb_plusu
