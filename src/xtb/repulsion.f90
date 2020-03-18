! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> TODO
module xtb_xtb_repulsion
   use xtb_mctc_accuracy, only : wp
   use xtb_type_molecule
   use xtb_xtb_data
   implicit none
   private

   public :: TxTBRepulsion
   public :: repulsionEnGrad


   !>
   type :: TxTBRepulsion

      !> [nId, nId]
      real(wp), allocatable :: alpha(:, :)

      !> [nId, nId]
      real(wp), allocatable :: zeff(:, :)

      !> [nId, nId]
      real(wp), allocatable :: kexp(:, :)

      !> [nId, nId]
      real(wp), allocatable :: rexp(:, :)

      !>
      real(wp) :: cutoff

   end type TxTBRepulsion


   interface init
      module procedure :: initRepulsion
   end interface init


contains


!>
subroutine initRepulsion(self)

   !>
   type(TxTBRepulsion), intent(out) :: self

end subroutine initRepulsion


!> Repulsion gradient of GFN1-xTB
subroutine repulsionEnGrad(mol, repData, energy, gradient, sigma)

   !>
   type(TMolecule), intent(in) :: mol

   !>
   type(TRepulsionData), intent(in) :: repData

   !>
   real(wp), intent(inout) :: gradient(:, :)

   !>
   real(wp), intent(inout) :: sigma(:, :)

   !>
   real(wp), intent(inout) :: energy

   integer  :: iat, jat, iZp, jZp, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, zeff, kExp, cutoff2
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   cutoff2 = repData%cutoff**2

   allocate(energies(len(mol)), source=0.0_wp)

   !$omp parallel do default(none) &
   !$omp reduction(+:energies, gradient, sigma) shared(repData, mol, cutoff2) &
   !$omp private(iat, jat, iZp, jZp, r2, rij, r1, alpha, zeff, kexp, &
   !$omp&        t16, t26, t27, dE, dG, dS)
   do iAt = 1, len(mol)
      iZp = mol%at(iAt)
      do jAt = 1, iAt-1
         jZp = mol%at(jAt)
         rij = mol%xyz(:, iAt) - mol%xyz(:, jAt)
         r2 = sum(rij**2)
         if (r2 > cutoff2) cycle
         r1 = sqrt(r2)

         alpha = sqrt(repData%alpha(iZp)*repData%alpha(jZp))
         zeff = repData%zeff(iZp)*repData%zeff(jZp)
         if (iZp > 2 .or. jZp > 2) then
            kExp = repData%kExp
         else
            kExp = repData%kExpLight
         end if
         t16 = r1**kExp
         t26 = exp(-alpha*t16)
         t27 = r1**repData%rExp
         dE = zeff * t26/t27
         dG = -(alpha*t16*kExp + repData%rExp) * dE * rij/r2
         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
         energies(iAt) = energies(iAt) + 0.5_wp * dE
         sigma = sigma + 0.5_wp * dS
         if (iAt /= jAt) then
            energies(jAt) = energies(jAt) + 0.5_wp * dE
            sigma = sigma + 0.5_wp * dS
            gradient(:, iAt) = gradient(:, iAt) + dG
            gradient(:, jAt) = gradient(:, jAt) - dG
         endif
      enddo
   enddo
   !$omp end parallel do

   energy = energy + sum(energies)

end subroutine repulsionEnGrad


!!> Repulsion gradient of GFN1-xTB
!subroutine getEnGrad(self, mol, energy, gradient, sigma)
!
!   !>
!   class(TxTBRepulsion), intent(in) :: self
!
!   !>
!   type(TMolecule), intent(in) :: mol
!
!   !>
!   real(wp), intent(inout) :: gradient(:, :)
!
!   !>
!   real(wp), intent(inout) :: sigma(:, :)
!
!   !>
!   real(wp), intent(inout) :: energy
!
!   integer  :: iat, jat, ati, atj, ij, img
!   real(wp) :: t16, t26, t27
!   real(wp) :: alpha, zeff
!   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
!   real(wp), allocatable :: energies(:)
!
!   allocate(energies(len(mol)), source=0.0_wp)
!
!   !$omp parallel do default(none) &
!   !$omp reduction(+:energies, gradient, sigma) shared(self, mol) &
!   !$omp private(ij, img, jat, ati, atj, r2, rij, r1, alpha, zeff, &
!   !$omp&        t16, t26, t27, dE, dG, dS)
!   do iAt = 1, len(mol)
!      iId = mol%id(iAt)
!      do jAt = 1, iAt-1
!         jId = mol%id(jAt)
!         rij = mol%xyz(:, iAt) - mol%xyz(:, jAt)
!         r2 = sum(rij**2)
!         if (r2 > cutoff2) cycle
!         r1 = sqrt(r2)
!
!         alpha = self%alpha(jId, iId)
!         zeff = self%zeff(jId, iId)
!         t16 = r1**self%kExp(jId, iId)
!         t26 = exp(-alpha*t16)
!         t27 = r1**self%rExp(jId, iId)
!         dE = zeff * t26/t27
!         dG = -(alpha*t16*self%kExp(jId, iId) + self%rExp(jId, iId)) * dE*rij/r2
!         dS = spread(dG, 1, 3) * spread(rij, 2, 3)
!         energies(iAt) = energies(iAt) + 0.5_wp * dE
!         sigma = sigma + 0.5_wp * dS
!         if (iAt /= jAt) then
!            energies(jAt) = energies(jAt) + 0.5_wp * dE
!            sigma = sigma + 0.5_wp * dS
!            gradient(:, iAt) = gradient(:, iAt) + dG
!            gradient(:, jAt) = gradient(:, jAt) - dG
!         endif
!      enddo
!   enddo
!   !$omp end parallel do
!
!   energy = energy + sum(energies)
!
!end subroutine repulsionEnGrad

end module xtb_xtb_repulsion
