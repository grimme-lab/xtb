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
module xtb_coulomb_gaussian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_coulomb_ewald, only : ewaldMatPBC3D, ewaldDerivPBC3D => ewaldDerivPBC3D_alp
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_param, only : chrg_parameter
   implicit none
   private

   public :: get_coulomb_matrix, get_coulomb_derivs


   ! √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)


   interface get_coulomb_matrix
      module procedure :: get_coulomb_matrix_0d
      module procedure :: get_coulomb_matrix_3d
   end interface get_coulomb_matrix


   interface get_coulomb_derivs
      module procedure :: get_coulomb_derivs_0d
      module procedure :: get_coulomb_derivs_3d
   end interface get_coulomb_derivs


contains


subroutine get_coulomb_matrix_0d(mol, chrgeq, amat)
   type(TMolecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(out) :: amat(:, :)
   integer :: iat, jat
   real(wp) :: r1, gamij
   Amat = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) shared(mol,chrgeq, amat) &
   !$omp private(jat, r1, gamij)
   do iat = 1, len(mol)
      ! EN of atom i
      do jat = 1, iat-1
         r1 = norm2(mol%xyz(:,jat) - mol%xyz(:,iat))
         gamij = 1.0_wp/sqrt(chrgeq%alpha(iat)**2+chrgeq%alpha(jat)**2)
         Amat(jat, iat) = erf(gamij*r1)/r1
         Amat(iat, jat) = Amat(jat,iat)
      enddo
      Amat(iat, iat) = chrgeq%gam(iat) + sqrt2pi/chrgeq%alpha(iat)
   enddo
   !$omp end parallel do
end subroutine get_coulomb_matrix_0d

subroutine get_coulomb_matrix_3d(mol, chrgeq, rTrans, gTrans, cf, amat)
   use xtb_type_molecule
   use xtb_type_param
   type(TMolecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: rTrans(:, :)
   real(wp), intent(in) :: gTrans(:, :)
   real(wp), intent(in) :: cf
   real(wp), intent(out) :: amat(:, :)
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp
   integer :: iat, jat, wscAt, iG1, iG2, iG3, iRp
   real(wp) :: gamii, gamij, riw(3), gVec(3)

   amat = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) reduction(+:amat) &
   !$omp shared(mol, chrgeq, cf, gTrans, rTrans) &
   !$omp private(iat, jat, wscAt, gamii, gamij, riw)
   do iat = 1, len(mol)
      gamii = 1.0_wp/(sqrt(2.0_wp)*chrgeq%alpha(iat))
      amat(iat, iat) = chrgeq%gam(iat) + sqrt2pi/chrgeq%alpha(iat) &
         ! reciprocal for 0th atom
         + ewaldMatPBC3D(zero, gTrans, 0.0_wp, mol%volume, cf, 1.0_wp) &
         ! direct for 0th atom
         + eeq_ewald_3d_dir(zero, rTrans, gamii, cf, 1.0_wp)
      do jat = 1, iat-1
         gamij = 1.0_wp/sqrt(chrgeq%alpha(iat)**2+chrgeq%alpha(jat)**2)
         do wscAt = 1, mol%wsc%itbl(jat,iat)
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            amat(iat,jat) = Amat(iat,jat) &
               ! reciprocal lattice sums
               + ewaldMatPBC3D(riw, gTrans, 0.0_wp, mol%volume, cf, &
                  & mol%wsc%w(jat,iat)) &
               ! direct lattice sums
               + eeq_ewald_3d_dir(riw, rTrans, gamij, cf, mol%wsc%w(jat,iat))
         end do
         amat(jat,iat) = amat(iat,jat)
      end do
   end do
   !$omp end parallel do
end subroutine get_coulomb_matrix_3d

subroutine get_coulomb_derivs_0d(mol, chrgeq, qvec, amatdr, atrace)
   use xtb_type_molecule
   use xtb_type_param
   type(TMolecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: amatdr(:, :, :)
   real(wp), intent(out) :: atrace(:, :)
   integer :: iat, jat
   real(wp) :: rij(3), r2, gamij2, arg2
   real(wp) :: dE, dG(3)
   amatdr = 0.0_wp
   atrace = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) reduction(+:atrace, amatdr) &
   !$omp shared(mol,chrgeq, qvec) private(jat, rij, r2, gamij2, arg2, dE, dG)
   do iat = 1, len(mol)
      ! EN of atom i
      do jat = 1, iat-1
         rij = mol%xyz(:,iat) - mol%xyz(:,jat)
         r2 = sum(rij**2)
         gamij2 = 1.0_wp/(chrgeq%alpha(iat)**2+chrgeq%alpha(jat)**2)
         arg2 = gamij2*r2
         dE = erf(sqrt(arg2))/sqrt(r2)
         dG = (2*sqrt(gamij2)*exp(-arg2)/sqrtpi - dE) * rij/r2
         amatdr(:, iat, jat) = amatdr(:, iat, jat) + dG*qvec(iat)
         amatdr(:, jat, iat) = amatdr(:, jat, iat) - dG*qvec(jat)
         atrace(:, iat) = atrace(:, iat) + dG*qvec(jat)
         atrace(:, jat) = atrace(:, jat) - dG*qvec(iat)
      enddo
   enddo
   !$omp end parallel do
end subroutine get_coulomb_derivs_0d

subroutine get_coulomb_derivs_3d(mol, chrgeq, qvec, rTrans, gTrans, cf, &
      & amatdr, amatdL, atrace)
   use xtb_type_molecule
   use xtb_type_param
   type(TMolecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: cf
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(in) :: gTrans(:, :)
   real(wp), intent(in) :: rTrans(:, :)
   real(wp), intent(out) :: amatdr(:, :, :)
   real(wp), intent(out) :: amatdL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp
   integer :: iat, jat, wscAt, ii
   real(wp) :: gamij, dG(3), dS(3, 3), riw(3)
   integer :: iG1, iG2, iG3, iRp
   real(wp) :: gVec(3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   amatdr = 0.0_wp
   amatdL = 0.0_wp
   atrace = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace,amatdr,amatdL) &
   !$omp shared(mol, chrgeq, qvec, cf, gTrans, rTrans) &
   !$omp private(iat, jat, ii, wscAt, riw, gamij, dG, dS)
   do iat = 1, len(mol)
      do jat = 1, iat-1
         ! over WSC partner
         gamij = 1.0_wp/sqrt(chrgeq%alpha(iat)**2 + chrgeq%alpha(jat)**2)
         do wscAt = 1, mol%wsc%itbl(jat,iat)
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            call ewaldDerivPBC3D(riw, gTrans, 0.0_wp, mol%volume, cf, &
               & mol%wsc%w(jat,iat), dG, dS)
            amatdr(:, iat, jat) = amatdr(:, iat, jat) + dG*qvec(iat)
            amatdr(:, jat, iat) = amatdr(:, jat, iat) - dG*qvec(jat)
            atrace(:, iat) = atrace(:, iat) + dG*qvec(jat)
            atrace(:, jat) = atrace(:, jat) - dG*qvec(iat)
            amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(jat)
            amatdL(:, :, jat) = amatdL(:, :, jat) + dS*qvec(iat)
            call eeq_ewald_dx_3d_dir(riw, rTrans, gamij, cf, mol%wsc%w(jat,iat), &
               &                     dG,dS)
            amatdr(:, iat, jat) = amatdr(:, iat, jat) + dG*qvec(iat)
            amatdr(:, jat, iat) = amatdr(:, jat, iat) - dG*qvec(jat)
            atrace(:, iat) = atrace(:, iat) + dG*qvec(jat)
            atrace(:, jat) = atrace(:, jat) - dG*qvec(iat)
            amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(jat)
            amatdL(:, :, jat) = amatdL(:, :, jat) + dS*qvec(iat)
         enddo  ! k WSC partner
      enddo     ! jat

      gamij = 1.0_wp/(sqrt(2.0_wp)*chrgeq%alpha(iat))
      call ewaldDerivPBC3D(zero, gTrans, 0.0_wp, mol%volume, cf, 1.0_wp, dG, dS)
      amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(iat)
      call eeq_ewald_dx_3d_dir(zero, rTrans, gamij, cf, 1.0_wp, dG, dS)
      amatdL(:, :, iat) = amatdL(:, :, iat) + (dS+unity*cf/sqrtpi/3.0_wp)*qvec(iat)
   enddo
   !$omp end parallel do
end subroutine get_coulomb_derivs_3d


pure function eeq_ewald_3d_dir(riw,rTrans,gamij,cf,scale) result(Amat)
   use xtb_mctc_constants
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   real(wp),intent(in) :: rTrans(:,:) !< direct lattice
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(in) :: scale
   real(wp) :: Amat                 !< element of interaction matrix
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz,itr
   real(wp) :: distiw,rij(3)
   real(wp) :: t(3)
   Amat = 0.0_wp
   do itr = 1, size(rTrans, dim=2)
      rij = riw + rTrans(:, itr)
      distiw = norm2(rij)
      ! self-interaction correction
      if(distiw < eps) then
         Amat = Amat - cf/sqrtpi
      else
         Amat = Amat - erf(   cf*distiw)/distiw &
            &        + erf(gamij*distiw)/distiw
      end if
   end do
   Amat = Amat * scale
end function eeq_ewald_3d_dir

pure subroutine eeq_ewald_dx_3d_dir(riw,rTrans,gamij,cf,scale,dAmat,sigma)
   use xtb_mctc_constants
   use xtb_pbc_tools
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   real(wp),intent(in) :: rTrans(:,:)
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(in) :: scale
   real(wp),intent(out) :: dAmat(3) !< element of interaction matrix
   real(wp),intent(out) :: sigma(3,3)
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz,i,itr
   real(wp) :: distiw,rij(3),arga,argb
   real(wp) :: t(3),dtmp,stmp(3)
   dAmat = 0.0_wp
   sigma = 0.0_wp
   do itr = 1, size(rTrans, dim=2)
      ! real contributions
      rij = riw + rTrans(:, itr)
      distiw = norm2(rij)
      if(distiw < eps) cycle
      arga = cf**2   *distiw**2
      stmp = + exp(-arga)/sqrtpi * cf * 2.0_wp / 3.0_wp
      do i = 1, 3
         sigma(i,i) = sigma(i,i) + stmp(i)! * rij(i)**2
      enddo
      argb = gamij**2*distiw**2
      dtmp = - 2*cf*exp(-arga)/(sqrtpi*distiw**2) &
         &   + erf(cf*distiw)/(distiw**3)           &
         &   + 2*gamij*exp(-argb)/(sqrtpi*distiw**2) &
         &   - erf(gamij*distiw)/(distiw**3)
      dAmat = dAmat + rij*dtmp
      sigma = sigma + dtmp*outer_prod_3x3(rij,rij)
   enddo
   dAmat = dAmat * scale
   sigma = sigma * scale

end subroutine eeq_ewald_dx_3d_dir


end module xtb_coulomb_gaussian
