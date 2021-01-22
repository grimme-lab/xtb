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

!> Implementation of the repulsion energy used in the xTB Hamiltonian
module xtb_xtb_repulsion
   use xtb_mctc_accuracy, only : wp
   use xtb_type_identitymap, only : TIdentityMap
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_neighbourlist, only : TNeighbourlist
   use xtb_xtb_data
   implicit none
   private

   public :: repulsionEnGrad


   interface repulsionEnGrad
      module procedure :: repulsionEnGrad_latp
      module procedure :: repulsionEnGrad_neighs
   end interface repulsionEnGrad


contains


!> Lattice point based implementation of the repulsion energy
subroutine repulsionEnGrad_latp(mol, repData, trans, cutoff, energy, gradient, &
      & sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Repulsion parametrisation
   type(TRepulsionData), intent(in) :: repData

   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   !> Repulsion energy
   real(wp), intent(inout) :: energy

   integer  :: iat, jat, iZp, jZp, itr, k, l
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, zeff, kExp, cutoff2
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   cutoff2 = cutoff**2

   allocate(energies(mol%n))

   !$acc enter data create(energies, rij, dS, dG) copyin(gradient, sigma, mol, mol%at, &
   !$acc& mol%xyz, repData, repData%alpha, repData%zeff, trans)

   !$acc kernels default(present)
   energies(:) = 0.0_wp
   !$acc end kernels

#ifdef XTB_GPU
   !$acc parallel default(present)
   !$acc loop gang collapse(2) private(iZp, jZp, alpha, zeff, kExp)
#else
   !$omp parallel do default(none) reduction(+:energies, gradient, sigma) &
   !$omp shared(repData, mol, cutoff2, trans) &
   !$omp private(iat, jat, itr, iZp, jZp, r2, rij, r1, alpha, zeff, kexp, &
   !$omp& t16, t26, t27, dE, dG, dS, k, l)
#endif
   do iAt = 1, mol%n
      do jAt = 1, mol%n
         if (jAt > iAt) cycle
         iZp = mol%at(iAt)
         jZp = mol%at(jAt)
         alpha = sqrt(repData%alpha(iZp)*repData%alpha(jZp))
         zeff = repData%zeff(iZp)*repData%zeff(jZp)
         if (iZp > 2 .or. jZp > 2) then
            kExp = repData%kExp
         else
            kExp = repData%kExpLight
         end if
         !$acc loop vector private(rij, r2, r1, t16, t26, t27, dE, dG, dS)
         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iAt) - mol%xyz(:, jAt) - trans(:, itr)
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-8_wp) cycle
            r1 = sqrt(r2)

            t16 = r1**kExp
            t26 = exp(-alpha*t16)
            t27 = r1**repData%rExp
            dE = zeff * t26/t27
            dG = -(alpha*t16*kExp + repData%rExp) * dE * rij/r2
            dS = spread(dG, 1, 3) * spread(rij, 2, 3)
            !$acc atomic
            energies(iAt) = energies(iAt) + 0.5_wp * dE
            if (iAt /= jAt) then
               !$acc atomic
               energies(jAt) = energies(jAt) + 0.5_wp * dE
               !$acc loop seq
               do k = 1, 3
                  !$acc atomic
                  gradient(k, iAt) = gradient(k, iAt) + dG(k)
                  !$acc atomic
                  gradient(k, jAt) = gradient(k, jAt) - dG(k)
               end do
               !$acc loop seq collapse(2)
               do k = 1, 3
                  do l = 1, 3
                     !$acc atomic
                     sigma(l, k) = sigma(l, k) + dS(l, k)
                  end do
               end do
            else
               !$acc loop seq collapse(2)
               do k = 1, 3
                  do l = 1, 3
                     !$acc atomic
                     sigma(l, k) = sigma(l, k) + 0.5_wp * dS(l, k)
                  end do
               end do
            endif
         enddo
      enddo
   enddo
#ifdef XTB_GPU
   !$acc end parallel

   !$acc exit data copyout(energies, gradient, sigma) delete(mol, mol%at, &
   !$acc& mol%xyz, repData, repData%alpha, repData%zeff, rij, dG, dS, trans)
#endif XTB_GPU

   energy = energy + sum(energies)

end subroutine repulsionEnGrad_latp


!> Lattice point based implementation of the repulsion energy
subroutine repulsionEnGrad_neighs(mol, repData, neighs, neighList, energy, &
      & gradient, sigma)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Repulsion parametrisation
   type(TRepulsionData), intent(in) :: repData

   !> Number of neighbours for each atom
   integer, intent(in) :: neighs(:)

   !> Neighbourlist
   class(TNeighbourList), intent(in) :: neighList

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout) :: sigma(:, :)

   !> Repulsion energy
   real(wp), intent(inout) :: energy

   integer  :: iat, jat, iZp, jZp, ij, img
   real(wp) :: t16, t26, t27
   real(wp) :: alpha, zeff, kExp, cutoff2
   real(wp) :: r1, r2, rij(3), dS(3, 3), dG(3), dE
   real(wp), allocatable :: energies(:)

   allocate(energies(len(mol)))
   energies(:) = 0.0_wp

   !$omp parallel do default(none) reduction(+:energies, gradient, sigma) &
   !$omp shared(repData, mol, neighs, neighList) &
   !$omp private(iat, jat, ij, img, iZp, jZp, r2, rij, r1, alpha, zeff, kexp, &
   !$omp& t16, t26, t27, dE, dG, dS)
   do iAt = 1, len(mol)
      iZp = mol%at(iAt)
      do ij = 1, neighs(iAt)
         img = neighList%ineigh(ij, iAt)
         jAt = neighList%image(img)
         jZp = mol%at(jAt)
         rij = neighList%coords(:, iAt) - neighList%coords(:, img)
         r2 = neighList%dist2(ij, iAt)
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

end subroutine repulsionEnGrad_neighs


end module xtb_xtb_repulsion
