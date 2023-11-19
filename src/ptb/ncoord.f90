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

!> Module for the relevant coordination numbers

module xtb_ptb_ncoord
   use mctc_env, only: wp
   use mctc_io, only: structure_type

   use tblite_data_covrad, only: get_covalent_rad

   implicit none

   private

   public :: ncoord_erf

contains

   !> @brief Calculate the coordination number from the error function
   subroutine ncoord_erf(mol, kcn, cutoff, cn, covrad)

      !> Molecular structure data
      type(structure_type), intent(in) :: mol

      !> Steepness of the counting function
      real(wp), intent(in) :: kcn

      !> Real space cutoff
      real(wp), intent(in) :: cutoff

      !> Covalent radius
      real(wp), intent(in), optional :: covrad(:)
      !> Effectively used covalent radius
      real(wp), allocatable :: rcov(:)

      !> Error function coordination number.
      real(wp), intent(out) :: cn(:)

      integer :: iat, jat, izp, jzp
      real(wp) :: r2, r1, rc, rij(3), countf, cutoff2

      cn(:) = 0.0_wp
      cutoff2 = cutoff**2

      allocate(rcov(mol%nid))
      if (present(covrad)) then
         rcov(:) = covrad
      else
         rcov(:) = get_covalent_rad(mol%num)
      end if

      !$omp parallel do default(none) reduction(+:cn) &
      !$omp shared(mol, cutoff2, rcov, kcn) &
      !$omp private(jat, izp, jzp, r2, rij, r1, rc, countf)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat
            jzp = mol%id(jat)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = erf_count(kcn, r1, rc)

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if
         end do
      end do

   end subroutine ncoord_erf

!> Error function counting function for coordination number contributions.
   pure function erf_count(k, r, r0) result(count)

      !> Steepness of the counting function.
      real(wp), intent(in) :: k

      !> Current distance.
      real(wp), intent(in) :: r

      !> Sum of covalent radii.
      real(wp), intent(in) :: r0

      real(wp) :: count

      count = 0.5_wp*(1.0_wp + erf(-k*(r - r0)/r0))

   end function erf_count

end module xtb_ptb_ncoord
