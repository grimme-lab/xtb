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
module xtb_coulomb_klopmanohno
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_mctc_constants, only : pi, sqrtpi
   use xtb_mctc_math, only : matInv3x3, matDet3x3
   use xtb_coulomb_ewald
   use xtb_type_coulomb, only : TCoulomb
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_latticepoint, only : TLatticePoint, init
   use xtb_type_wignerseitzcell, only : TWignerSeitzCell, init
   implicit none
   private

   public :: TKlopmanOhno, init
   public :: gamAverage


   !> Possible averaging functions for hardnesses
   type :: TGamAverageEnum

      !> Harmonic average used in GFN1-xTB
      integer :: harmonic = 1

      !> Arithmetic average used in GFN2-xTB
      integer :: arithmetic = 2

      !> Geometric average
      integer :: geometric = 3

   end type TGamAverageEnum

   !> Actual enumerator for hardnesses
   type(TGamAverageEnum), parameter :: gamAverage = TGamAverageEnum()


   type, extends(TCoulomb) :: TKlopmanOhno

      !> Shell/Atomic hardnesses for each species
      real(wp), allocatable :: hardness(:, :)

      !> Generalized exponent
      real(wp) :: gExp

      !> Averaging function for the hardnesses
      procedure(funcAverage), nopass, pointer :: gamAverage => null()

   contains

      !> Returns full Coulomb matrix
      procedure :: getCoulombMatrix

      !> Returns derivatives of Coulomb matrix
      procedure :: getCoulombDerivs

   end type TKlopmanOhno


   interface init
      module procedure :: initFromMolecule
      module procedure :: initKlopmanOhno
   end interface init


   abstract interface
   pure function funcAverage(gi, gj) result(gam)
      import :: wp

      !> Hardness of shell i
      real(wp), intent(in) :: gi

      !> Hardness of shell j
      real(wp), intent(in) :: gj

      !> Averaged hardness
      real(wp) :: gam

   end function funcAverage
   end interface


contains


subroutine initFromMolecule(self, env, mol, gav, hardness, gExp, nshell, alpha, &
      & tolerance)

   !> Source of the generated error
   character(len=*), parameter :: source = 'type_coulomb_initFromMolecule'

   !> Instance of the Coulomb evaluator
   type(TKlopmanOhno), intent(out) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Averaging function
   integer, intent(in) :: gav

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: hardness(:, :)

   !> Generalized exponent
   real(wp), intent(in) :: gExp

   !> Number of shell for each species
   integer, intent(in), optional :: nshell(:)

   !> Convergence factor
   real(wp), intent(in), optional :: alpha

   !> Tolerance for the Ewald sum
   real(wp), intent(in), optional :: tolerance

   logical :: exitRun

   call init(self, env, mol%id, mol%lattice, mol%boundaryCondition, gav, &
      & hardness, gExp, nshell, alpha, tolerance)

   call env%check(exitRun)
   if (exitRun) return

   call self%update(env, mol)
   call env%check(exitRun)
   if (exitRun) then
      call env%error("Initializing internal state of evaluator failed", source)
   end if

end subroutine initFromMolecule


subroutine initKlopmanOhno(self, env, id, lattice, boundaryCond, gav, hardness, &
      & gExp, nshell, alpha, tolerance)

   !> Source of the generated error
   character(len=*), parameter :: source = 'type_coulomb_initCoulomb'

   !> Instance of the Coulomb evaluator
   type(TKlopmanOhno), intent(out) :: self

   !> Identity of each atom
   integer, intent(in) :: id(:)

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Lattice parameters
   real(wp), intent(in) :: lattice(:, :)

   !> Boundary conditions for this evaluator
   integer, intent(in) :: boundaryCond

   !> Averaging function
   integer, intent(in) :: gav

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: hardness(:, :)

   !> Generalized exponent
   real(wp), intent(in) :: gExp

   !> Number of shell for each species
   integer, intent(in), optional :: nshell(:)

   !> Convergence factor
   real(wp), intent(in), optional :: alpha

   !> Tolerance for the Ewald sum
   real(wp), intent(in), optional :: tolerance

   logical :: exitRun
   integer :: natom, ndim
   integer :: ind, iat, ish
   real(wp) :: volume, recLat(3, 3)

   self%boundaryCondition = boundaryCond
   self%hardness = hardness
   self%gExp = gExp

   natom = size(id, dim=1)
   allocate(self%itbl(2, natom))
   if (present(nshell)) then
      ind = 0
      do iat = 1, natom
         ish = nshell(id(iat))
         self%itbl(:, iat) = [ind, ish]
         ind = ind + ish
      end do
   else
      do iat = 1, natom
         self%itbl(:, iat) = [iat-1, 1]
      end do
   end if

   select case(gav)
   case default
      call env%error("Unknown averaging function specified", source)
      return
   case(gamAverage%harmonic)
      self%gamAverage => harmonicAverage
   case(gamAverage%arithmetic)
      self%gamAverage => arithmeticAverage
   case(gamAverage%geometric)
      self%gamAverage => geometricAverage
   end select

   select case(self%boundaryCondition)
   case default
      call env%error("Boundary condition not supported", source)
   case(boundaryCondition%cluster)
      ! nothing to do
   case(boundaryCondition%pbc3d)
      volume = abs(matDet3x3(lattice))
      recLat(:, :) = 2*pi*transpose(matInv3x3(lattice))

      if (present(tolerance)) then
         self%tolerance = tolerance
      else
         self%tolerance = 1.0e-8_wp
      end if

      if (present(alpha)) then
         self%alpha = alpha
      else
         call getOptimalAlpha(env, lattice, recLat, volume, self%tolerance, &
            & self%alpha)
      end if

      call getMaxR(env, self%alpha, self%tolerance, self%rCutoff)
      call getMaxG(env, self%alpha, volume, self%tolerance, self%gCutoff)

      call env%check(exitRun)
      if (exitRun) then
         call env%error("Could not setup numerical thresholds", source)
         return
      end if

      call init(self%rLatPoint, env, lattice, boundaryCond, self%rCutoff)
      call init(self%gLatPoint, env, recLat, boundaryCond, self%gCutoff)
      call init(self%wsCell, natom)
   end select

end subroutine initKlopmanOhno


subroutine getCoulombMatrix(self, mol, jmat)

   !> Instance of the Coulomb evaluator
   class(TKlopmanOhno), intent(inout) :: self

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Coulomb matrix
   real(wp), intent(out) :: jmat(:, :)

   select case(self%boundaryCondition)
   case(boundaryCondition%cluster)
      call getCoulombMatrixCluster(mol, self%itbl, self%gamAverage, self%gExp, &
         & self%hardness, jmat)
   case(boundaryCondition%pbc3d)
      call getCoulombMatrixPBC3D(self%wsCell, self%itbl, self%gamAverage, &
         & self%gExp, self%hardness, self%alpha, mol%volume, self%rTrans, &
         & self%gTrans(:, 2:), jmat)
   end select

end subroutine getCoulombMatrix


subroutine getCoulombMatrixCluster(mol, itbl, gamAverage, gExp, hardness, jmat)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Index table
   integer, intent(in) :: itbl(:, :)

   !> Averaging function for the hardnesses
   procedure(funcAverage) :: gamAverage

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: hardness(:, :)

   !> Generalized exponent
   real(wp), intent(in) :: gExp

   !> Coulomb matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, jat, ish, jsh, ii, jj, iid, jid
   real(wp) :: r1, rterm, gij

   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) shared(mol, itbl, hardness, gExp, jmat) &
   !$omp private(iat, jat, ish, jsh, ii, jj, iid, jid, r1, rterm, gij)
   do iat = 1, len(mol)
      ii = itbl(1, iat)
      iid = mol%id(iat)
      do jat = 1, iat-1
         jj = itbl(1, jat)
         jid = mol%id(jat)
         r1 = norm2(mol%xyz(:, jat) - mol%xyz(:, iat))
         rterm = 1.0_wp/r1
         do ish = 1, itbl(2, iat)
            do jsh = 1, itbl(2, jat)
               gij = gamAverage(hardness(ish, iid), hardness(jsh, jid))
               rterm = 1.0_wp/(r1**gExp + gij**(-gExp))**(1.0_wp/gExp)
               jmat(jj+jsh, ii+ish) = rterm
               jmat(ii+ish, jj+jsh) = rterm
            end do
         end do
      end do
      do ish = 1, itbl(2, iat)
         do jsh = 1, ish-1
            gij = gamAverage(hardness(ish, iid), hardness(jsh, iid))
            jmat(ii+jsh, ii+ish) = gij
            jmat(ii+ish, ii+jsh) = gij
         end do
         jmat(ii+ish, ii+ish) = hardness(ish, iid)
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombMatrixCluster


subroutine getCoulombMatrixPBC3D(wsCell, itbl, gamAverage, gExp, hardness, &
      & alpha, volume, rTrans,  gTrans, jmat)

   !> Wigner-Seitz cell
   type(TWignerSeitzCell), intent(in) :: wsCell

   !> Index table
   integer, intent(in) :: itbl(:, :)

   !> Averaging function for the hardnesses
   procedure(funcAverage) :: gamAverage

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: hardness(:, :)

   !> Generalized exponent
   real(wp), intent(in) :: gExp

   !> Cell volume
   real(wp), intent(in) :: volume

   !> Convergence factor
   real(wp), intent(in) :: alpha

   !> Real space lattice translations
   real(wp), intent(in) :: rTrans(:, :)

   !> Reciprocal space lattice translations
   real(wp), intent(in) :: gTrans(:, :)

   !> Coulomb matrix
   real(wp), intent(out) :: jmat(:, :)

   integer :: iat, ineigh, img, jat, ish, jsh, ii, jj
   real(wp) :: vec(3), rterm, gterm, weight, gij
   real(wp), parameter :: zero(3) = 0.0_wp

   jmat(:, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:jmat) &
   !$omp shared(wsCell, itbl, alpha, volume, gTrans, rTrans, gExp, hardness) &
   !$omp private(iat, ineigh, img, jat, ish, jsh, ii, jj, rterm, gterm, vec, &
   !$omp& weight, gij)
   do iat = 1, size(wsCell%neighs)
      ii = itbl(1, iat)
      do ineigh = 1, wsCell%neighs(iat)
         img = wsCell%ineigh(ineigh, iat)
         jat = wsCell%image(img)
         jj = itbl(1, jat)
         weight = wsCell%weight(ineigh, iat)
         vec(:) = wsCell%coords(:, img) - wsCell%coords(:, iat)
         gterm = ewaldMatPBC3D(vec, gTrans, 0.0_wp, volume, alpha, weight) &
            &  - pi / (volume * alpha**2) * weight
         if (iat /= jat) then
            do ish = 1, itbl(2, iat)
               do jsh = 1, itbl(2, jat)
                  gij = gamAverage(hardness(ish, iat), hardness(jsh, jat))
                  rterm = gterm + getRTerm(vec, gij, gExp, rTrans, alpha, weight)
                  jmat(jj+jsh, ii+ish) = jmat(jj+jsh, ii+ish) + rterm
                  jmat(ii+ish, jj+jsh) = jmat(ii+ish, jj+jsh) + rterm
               end do
            end do
         else
            do ish = 1, itbl(2, iat)
               do jsh = 1, ish-1
                  gij = gamAverage(hardness(ish, iat), hardness(jsh, iat))
                  rterm = gterm + getRTerm(vec, gij, gExp, rTrans, alpha, weight)
                  jmat(ii+jsh, ii+ish) = jmat(ii+jsh, ii+ish) + rterm
                  jmat(ii+ish, ii+jsh) = jmat(ii+ish, ii+jsh) + rterm
               end do
               gij = hardness(ish, iat)
               rterm = gterm + getRTerm(vec, gij, gExp, rTrans, alpha, weight)
               jmat(ii+ish, ii+ish) = jmat(ii+ish, ii+ish) + rterm
            end do
         end if
      end do
      do ish = 1, itbl(2, iat)
         do jsh = 1, ish-1
            gij = gamAverage(hardness(ish, iat), hardness(jsh, iat))
            jmat(ii+jsh, ii+ish) = jmat(ii+jsh, ii+ish) + gij
            jmat(ii+ish, ii+jsh) = jmat(ii+ish, ii+jsh) + gij
         end do
         jmat(ii+ish, ii+ish) = jmat(ii+ish, ii+ish) + hardness(ish, iat)
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombMatrixPBC3D


pure function getRTerm(vec, gam, gExp, rTrans, alpha, scale) result(rTerm)
   real(wp),intent(in) :: vec(3)
   real(wp),intent(in) :: gam
   real(wp),intent(in) :: gExp
   real(wp),intent(in) :: rTrans(:,:)
   real(wp),intent(in) :: alpha
   real(wp),intent(in) :: scale
   real(wp) :: rTerm
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: itr
   real(wp) :: r1, rij(3)
   rTerm = 0.0_wp
   do itr = 1, size(rTrans, dim=2)
      rij = vec + rTrans(:, itr)
      r1 = norm2(rij)
      ! self-interaction correction
      if(r1 < eps) then
         rTerm = rTerm - 2.0_wp*alpha/sqrtpi
      else
         rterm = 1.0_wp/(r1**gExp + gam**(-gExp))**(1.0_wp/gExp) &
            & - erf(alpha*r1)/r1
      end if
   end do
   rTerm = rTerm * scale
end function getRTerm


subroutine getCoulombDerivs(self, mol, qvec, djdr, djdtr, djdL)

   !> Instance of the Coulomb evaluator
   class(TKlopmanOhno), intent(inout) :: self

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Charges
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Coulomb matrix
   real(wp), intent(out) :: djdtr(:, :)

   !> Derivative of Coulomb matrix w.r.t. strain deformations
   real(wp), intent(out) :: djdL(:, :, :)

   select case(self%boundaryCondition)
   case(boundaryCondition%cluster)
      call getCoulombDerivsCluster(mol, self%itbl, self%gamAverage, self%gExp, &
         & self%hardness, qvec, djdr, djdtr, djdL)
   case(boundaryCondition%pbc3d)
      call getCoulombDerivsPBC3D(self%wsCell, self%itbl, self%alpha, mol%volume, &
         & self%rTrans, self%gTrans(:, 2:), qvec, djdr, djdtr, djdL)
   end select

end subroutine getCoulombDerivs


subroutine getCoulombDerivsCluster(mol, itbl, gamAverage, gExp, hardness, &
      & qvec, djdr, djdtr, djdL)

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Index table
   integer, intent(in) :: itbl(:, :)

   !> Averaging function for the hardnesses
   procedure(funcAverage) :: gamAverage

   !> Shell/Atomic hardnesses for each species
   real(wp), intent(in) :: hardness(:, :)

   !> Generalized exponent
   real(wp), intent(in) :: gExp

   !> Charges
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Coulomb matrix
   real(wp), intent(out) :: djdtr(:, :)

   !> Derivative of Coulomb matrix w.r.t. strain deformations
   real(wp), intent(out) :: djdL(:, :, :)

   integer :: iat, jat, ish, jsh, ii, jj, iid, jid
   real(wp) :: r1, g1, gij, vec(3), dG(3), dS(3, 3)

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp
   djdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) reduction(+:djdr, djdtr, djdL) &
   !$omp shared(mol, itbl, qvec, gExp, hardness) &
   !$omp private(iat, jat, ish, jsh, ii, jj, iid, jid, r1, g1, gij, vec, dG, dS)
   do iat = 1, len(mol)
      ii = itbl(1, iat)
      iid = mol%id(iat)
      do jat = 1, iat-1
         jj = itbl(1, jat)
         jid = mol%id(jat)
         vec(:) = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         do ish = 1, itbl(2, iat)
            do jsh = 1, itbl(2, jat)
               gij = gamAverage(hardness(ish, iid), hardness(jsh, jid))
               g1 = 1.0_wp / (r1**gExp + gij**(-gExp))
               dG(:) = -vec*r1**(gExp-2.0_wp) * g1 * g1**(1.0_wp/gExp)
               dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
               djdr(:, iat, jj+jsh) = djdr(:, iat, jj+jsh) - dG*qvec(ii+ish)
               djdr(:, jat, ii+ish) = djdr(:, jat, ii+ish) + dG*qvec(jj+jsh)
               djdtr(:, jj+jsh) = djdtr(:, jj+jsh) + dG*qvec(ii+ish)
               djdtr(:, ii+ish) = djdtr(:, ii+ish) - dG*qvec(jj+jsh)
               djdL(:, :, jj+jsh) = djdL(:, :, jj+jsh) - dS*qvec(ii+ish)
               djdL(:, :, ii+ish) = djdL(:, :, ii+ish) - dS*qvec(jj+jsh)
            end do
         end do
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombDerivsCluster


subroutine getCoulombDerivsPBC3D(wsCell, itbl, alpha, volume, rTrans, gTrans, &
      & qvec, djdr, djdtr, djdL)

   !> Wigner-Seitz cell
   type(TWignerSeitzCell), intent(in) :: wsCell

   !> Cell volume
   real(wp), intent(in) :: volume

   !> Convergence factor
   real(wp), intent(in) :: alpha

   !> Real space lattice translations
   real(wp), intent(in) :: rTrans(:, :)

   !> Reciprocal space lattice translations
   real(wp), intent(in) :: gTrans(:, :)

   !> Index table
   integer, intent(in) :: itbl(:, :)

   !> Charges
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. Cartesian coordinates
   real(wp), intent(out) :: djdr(:, :, :)

   !> Trace derivative of Coulomb matrix
   real(wp), intent(out) :: djdtr(:, :)

   !> Derivative of Coulomb matrix w.r.t. strain deformations
   real(wp), intent(out) :: djdL(:, :, :)

   integer :: iat, jat, ish, jsh, ii, jj, ineigh, img
   real(wp) :: weight, vec(3), dG(3), dS(3, 3), dGg(3), dGr(3), dSg(3, 3), dSr(3, 3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   djdr(:, :, :) = 0.0_wp
   djdtr(:, :) = 0.0_wp
   djdL(:, :, :) = 0.0_wp

   dGr = 0.0_wp
   dSr = 0.0_wp
   !$omp parallel do default(none) reduction(+:djdr, djdtr, djdL) &
   !$omp shared(wsCell, itbl, qvec, alpha, volume, gTrans, rTrans) &
   !$omp private(iat, ineigh, img, jat, ish, jsh, ii, jj, vec, weight, dG, dS,&
   !$omp& dGg, dGr, dSg, dSr)
   do iat = 1, size(wsCell%neighs)
      ii = itbl(1, iat)
      do ineigh = 1, wsCell%neighs(iat)
         img = wsCell%ineigh(ineigh, iat)
         jat = wsCell%image(img)
         jj = itbl(1, jat)
         weight = wsCell%weight(ineigh, iat)
         vec(:) = wsCell%coords(:, img) - wsCell%coords(:, iat)
         call ewaldDerivPBC3D(vec, gTrans, 0.0_wp, volume, alpha, weight, dGg, dSg)
         call getRDeriv(vec, rTrans, alpha, weight, dGr, dSr)
         dG(:) = dGg + dGr
         dS(:, :) = dSg + dSr + pi / (volume * alpha**2) * weight * unity
         if (iat /= jat) then
            do ish = 1, itbl(2, iat)
               do jsh = 1, itbl(2, jat)
                  djdr(:, iat, jj+jsh) = djdr(:, iat, jj+jsh) - dG*qvec(ii+ish)
                  djdr(:, jat, ii+ish) = djdr(:, jat, ii+ish) + dG*qvec(jj+jsh)
                  djdtr(:, jj+jsh) = djdtr(:, jj+jsh) + dG*qvec(ii+ish)
                  djdtr(:, ii+ish) = djdtr(:, ii+ish) - dG*qvec(jj+jsh)
                  djdL(:, :, jj+jsh) = djdL(:, :, jj+jsh) - dS*qvec(ii+ish)
                  djdL(:, :, ii+ish) = djdL(:, :, ii+ish) - dS*qvec(jj+jsh)
               end do
            end do
         else
            do ish = 1, itbl(2, iat)
               do jsh = 1, ish-1
                  djdL(:, :, ii+jsh) = djdL(:, :, ii+jsh) - dS*qvec(ii+ish)
                  djdL(:, :, ii+ish) = djdL(:, :, ii+ish) - dS*qvec(ii+jsh)
               end do
               djdL(:, :, ii+ish) = djdL(:, :, ii+ish) - dS*qvec(ii+ish)
            end do
         end if
      end do
   end do
   !$omp end parallel do

end subroutine getCoulombDerivsPBC3D


pure subroutine getRDeriv(vec, rTrans, alpha, scale, dG, dS)
   real(wp), intent(in) :: vec(:)
   real(wp), intent(in) :: rTrans(:,:)
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: scale
   real(wp), intent(out) :: dG(:)
   real(wp), intent(out) :: dS(:,:)
   real(wp), parameter :: eps = 1.0e-9_wp
   integer :: itr
   real(wp) :: r1, rij(3), arg, dd
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])
   dG = 0.0_wp
   dS = 0.0_wp
   do itr = 1, size(rTrans, dim=2)
      ! real contributions
      rij = vec + rTrans(:, itr)
      r1 = norm2(rij)
      if (r1 < eps) cycle
      arg = alpha**2*r1**2
      dd = - 2*alpha*exp(-arg)/(sqrtpi*r1**2) - erfc(alpha*r1)/(r1**3)
      dG = dG + rij*dd
      dS = dS + dd*spread(rij, 1, 3)*spread(rij, 2, 3)
   enddo
   dG = dG * scale
   dS = dS * scale

end subroutine getRDeriv


!> Harmonic averaging functions for hardnesses in GFN1-xTB
pure function harmonicAverage(gi, gj) result(gam)

   !> Hardness of shell i
   real(wp), intent(in) :: gi

   !> Hardness of shell j
   real(wp), intent(in) :: gj

   !> Averaged hardness
   real(wp) :: gam

   gam = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)

end function harmonicAverage


!> Arithmetic averaging functions for hardnesses in GFN2-xTB
pure function arithmeticAverage(gi, gj) result(gam)

   !> Hardness of shell i
   real(wp), intent(in) :: gi

   !> Hardness of shell j
   real(wp), intent(in) :: gj

   !> Averaged hardness
   real(wp) :: gam

   gam = 0.5_wp*(gi+gj)

end function arithmeticAverage


!> Geometric averaging functions for hardnesses
pure function geometricAverage(gi, gj) result(gam)

   !> Hardness of shell i
   real(wp), intent(in) :: gi

   !> Hardness of shell j
   real(wp), intent(in) :: gj

   !> Averaged hardness
   real(wp) :: gam

   gam = sqrt(gi*gj)

end function geometricAverage


end module xtb_coulomb_klopmanohno
