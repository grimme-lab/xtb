! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_eeq
   use, intrinsic :: iso_fortran_env, only : istdout => output_unit
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_blas, only : mctc_gemv, mctc_gemm, mctc_symv, mctc_copy, mctc_dot
   use xtb_mctc_lapack, only : lapack_sytrf, lapack_sytri
   use xtb_type_environment, only : TEnvironment
   use xtb_gfn0param, alp => alpg
   use xtb_coulomb_ewald, only : ewaldMatPBC3D, ewaldDerivPBC3D => ewaldDerivPBC3D_alp
   use xtb_type_molecule, only : TMolecule, len
   use xtb_type_param
   implicit none

   public :: eeq_chrgeq
   public :: goedecker_chrgeq
   private

   ! √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)

   interface eeq_chrgeq
      module procedure :: eeq_chrgeq_core
      module procedure :: eeq_chrgeq_gbsa
      module procedure :: eeq_chrgeq_qonly
   end interface eeq_chrgeq

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
         r1 = sqrt(sum((mol%xyz(:,jat) - mol%xyz(:,iat))**2))
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
   !$omp parallel do default(none) reduction(+:amat) &
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
            ! PGI 20.7 fails with `ICE: Errors in Lowering' for
            ! - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               & - (mol%lattice(:,1) * mol%wsc%lattr(1,wscAt,jat,iat) &
               &  + mol%lattice(:,2) * mol%wsc%lattr(2,wscAt,jat,iat) &
               &  + mol%lattice(:,3) * mol%wsc%lattr(3,wscAt,jat,iat))
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
            ! PGI 20.7 fails with `ICE: Errors in Lowering' for
            ! - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               & - (mol%lattice(:,1) * mol%wsc%lattr(1,wscAt,jat,iat) &
               &  + mol%lattice(:,2) * mol%wsc%lattr(2,wscAt,jat,iat) &
               &  + mol%lattice(:,3) * mol%wsc%lattr(3,wscAt,jat,iat))
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
      distiw = sqrt(sum(rij**2))
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
      distiw = sqrt(sum(rij**2))
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


!! ========================================================================
!  Purpose:
!! ------------------------------------------------------------------------
!  calculate energies and properties from the Goedecker charge model
!  using the parametrisation from the GFN0-xTB.
!
!  Input:
!! ------------------------------------------------------------------------
!  n                 - number of atoms
!  at(n)             - atom type/ordinal number
!  xyz(3,n)          - molecular geometry
!  chrg              - total charge of the system
!  cn(n)             - coordination number (usually erf-CN)
!  dcndr(3,n,n)      - derivative of the coordination number
!                      (switched indices are assumed here!)
!
!  Output:
!! ------------------------------------------------------------------------
!  q(n)              - partial charges
!  dqdr(3,n,n+1)     - derivative of the partial charges and
!                      Lagrange multiplier (n+1 instead of n!)
!
!  Input/Output:
!! ------------------------------------------------------------------------
!  energy            - isotropic electrostatic energy
!  gradient(3,n)     - molecular gradient
!
!  Citation:
!! ------------------------------------------------------------------------
!  Original work: S. Alireza Ghasemi, Albert Hofstetter, Santanu Saha, and
!                 Stefan Goedecker, PHYSICAL REVIEW B 92, 045131 (2015).
!  This work:     S. Ehlert, E. Caldeweyher, S. Spicher and S. Grimme,
!                 to be published.
!
!  Note: For computational effiency this routine does NOT work with the
!        energy expression but with the Lagrangian (n+1 instead of n),
!        therefore, all intermediates live in the space of the Lagrangian.
!
!  Implemented by SAW in 2018, see focusing lab course report for details.
!! ========================================================================
subroutine goedecker_chrgeq(n,at,xyz,chrg,cn,dcndr,q,dqdr,energy,gradient,&
                            lverbose,lgrad,lcpq)
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n             ! number of atoms
   integer, intent(in)    :: at(n)         ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)      ! geometry
   real(wp),intent(in)    :: chrg          ! total charge
   real(wp),intent(in)    :: cn(n)         ! erf-CN
   real(wp),intent(in)    :: dcndr(3,n,n)  ! derivative of erf-CN
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(n)          ! partial charges
   real(wp),intent(out)   :: dqdr(3,n,n) ! derivative of partial charges
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,n) ! molecular gradient of IES
!
!! ------------------------------------------------------------------------
!  charge model
!! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmat(:,:,:)
   real(wp),allocatable :: dXvec(:,:,:)

!! ------------------------------------------------------------------------
!  local variables
!! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
   real(wp),allocatable :: alpha(:)
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)
   real(wp),allocatable :: Xfac(:)
   real(wp),allocatable :: Afac(:,:)

!! ------------------------------------------------------------------------
!  Lapack work variables
!! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)

!! ------------------------------------------------------------------------
!  initial banner if verbose
!! ------------------------------------------------------------------------
   if (lverbose) then
      write(istdout,'(72("="),/,1x,a,/,72("="),/)') &
         "GOEDECKER CHARGE MODEL USING GFN0-xTB PARAMETRIZATION"
   endif

!! ------------------------------------------------------------------------
!  initizialization
!! ------------------------------------------------------------------------
   m    = n+1
   q    = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(m), alpha(n), source = 0.0_wp )

!! ------------------------------------------------------------------------
!  set up the A matrix and X vector
!! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Setup of the A matrix and the RHS X vector"
!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(n,at,cn,xi,cnfak,alp,gamm) &
!$omp private(i,tmp) &
!$omp shared(Xvec,Xfac,alpha)
!$omp do schedule(runtime)
   do i = 1, n
      tmp = cnfak(at(i))/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -xi(at(i)) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
      alpha(i) = alp(at(i))**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(n,at,xyz,gamm,alpha) &
!$omp private(i,j,r,gamij) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(runtime)
   ! prepare A matrix
   do i = 1, n
      ! EN of atom i
      do j = 1, i-1
         r = sqrt(sum((xyz(:,j) - xyz(:,i))**2))
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = gamm(at(i)) + sqrt2pi/sqrt(alpha(i))
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges"
   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = chrg
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u', m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) call raise('E','(goedecker_solve) DSYSV failed')

   q = Xtmp(:n)
   if(abs(sum(q)-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error')
   !print'(3f20.14)',Xtmp

   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif

!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Isotropic electrostatic (IES) energy calculation"
   call mctc_copy(Xvec, work)
   call mctc_symv(Amat, Xtmp, work, alpha=0.5_wp, beta=-1.0_wp)
   es = mctc_dot(Xtmp, work)
   energy = es + energy
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "energy",es
   endif

!! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
!! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
!! ------------------------------------------------------------------------
do_molecular_gradient: if (lgrad .or. lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmat(3,n,m), dXvec(3,n,m), Afac(3,n), source = 0.0_wp )
   !allocate( dAmat(3,n,m), Afac(3,n), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(n,dcndr,xyz,alpha,Amat,Xfac,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp shared(dXvec,dAmat) &
!$omp reduction(+:Afac)
!$omp do schedule(runtime)
   do i = 1, n
      dXvec(:,i,i) = +dcndr(:,i,i)*Xfac(i) ! merge dX and dA for speedup
      do j = 1, i-1
         dXvec(:,j,i) = dcndr(:,i,j)*Xfac(i)
         dXvec(:,i,j) = dcndr(:,j,i)*Xfac(j)
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         gamij = 1.0_wp/sqrt(alpha(i) + alpha(j))
         arg = gamij**2*r2
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
         Afac(:,i) = +dtmp*rij*Xtmp(j) + Afac(:,i)
         Afac(:,j) = -dtmp*rij*Xtmp(i) + Afac(:,j)
         dAmat(:,i,j) = +dtmp*rij*Xtmp(i) !+ dcndr(:,i,j)*Xfac(i)
         dAmat(:,j,i) = -dtmp*rij*Xtmp(j) !+ dcndr(:,j,i)*Xfac(j)
      enddo
      !dAmat(:,i,i) = +dcndr(:,i,i)*Xfac(i)
   enddo
!$omp enddo
!$omp endparallel
endif do_molecular_gradient

   if (lgrad) then
   call mctc_gemv(dAmat,Xtmp,gradient,alpha=+1.0_wp,beta=1.0_wp)
   call mctc_gemv(dXvec,Xtmp,gradient,alpha=-1.0_wp,beta=1.0_wp)
   endif

!! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
do_partial_charge_derivative: if (lcpq) then
   dqdr = 0.0_wp
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "A matrix inversion by Bunch-Kaufman factorization"
   allocate( Ainv(m,m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(goedecker_inversion) DSYTRF failed')
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(goedecker_inversion) DSYTRI failed')
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

!! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges"
   do i = 1, n
      dAmat(:,i,i) = Afac(:,i) + dAmat(:,i,i)
   enddo
   call mctc_gemm(dAmat, Ainv, dqdr, alpha=-1.0_wp, beta=1.0_wp)
   call mctc_gemm(dXvec, Ainv, dqdr, alpha=+1.0_wp, beta=1.0_wp)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

!! ------------------------------------------------------------------------
!  Clean up
!! ------------------------------------------------------------------------
   if (allocated(alpha)) deallocate(alpha)
   if (allocated(Amat))  deallocate(Amat)
   if (allocated(dAmat)) deallocate(dAmat)
   if (allocated(Afac))  deallocate(Afac)
   if (allocated(Xvec))  deallocate(Xvec)
   if (allocated(Xfac))  deallocate(Xfac)
   if (allocated(Xtmp))  deallocate(Xtmp)
   if (allocated(Atmp))  deallocate(Atmp)
   if (allocated(temp))  deallocate(temp)
   if (allocated(work))  deallocate(work)
   if (allocated(ipiv))  deallocate(ipiv)

end subroutine goedecker_chrgeq

subroutine eeq_chrgeq_qonly(mol,env,chrgeq,cn,q,lverbose)

   use xtb_type_molecule
   use xtb_type_param

   implicit none

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(TMolecule),intent(in) :: mol     !< molecular structure information
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   type(chrg_parameter),intent(in) :: chrgeq
   real(wp),allocatable   :: dcndr(:,:,:)
   real(wp),allocatable   :: dcndL(:,:,:)
   logical, intent(in), optional :: lverbose
   logical :: verbose
   type(TEnvironment), intent(inout) :: env

! ------------------------------------------------------------------------
!  Output
! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(mol%n)
   real(wp),allocatable   :: dqdr(:,:,:)
   real(wp),allocatable   :: dqdL(:,:,:)
   real(wp)               :: energy
   real(wp),allocatable   :: gradient(:,:)
   real(wp),allocatable   :: sigma(:,:)

   if (present(lverbose)) then
      verbose = lverbose
   else
      verbose = .false.
   endif

   call eeq_chrgeq_core(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &                 energy,gradient,sigma,verbose,.false.,.false.)

end subroutine eeq_chrgeq_qonly

! ======================================================================
!  Modified Version of eeq_chrgeq routine that reads also the chrgeq construct
! ======================================================================
subroutine eeq_chrgeq_core(mol,env,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &                    energy,gradient,sigma,lverbose,lgrad,lcpq)

   use xtb_mctc_convert

   use xtb_type_molecule
   use xtb_type_param

   use xtb_disp_ncoord

   use xtb_pbc_tools

   implicit none

   character(len=*), parameter :: source = 'eeq_chrgeq'

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(TMolecule),intent(in) :: mol     !< molecular structure information
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   real(wp),intent(in)    :: dcndr(3,mol%n,mol%n) ! derivative of erf-CN
   real(wp),intent(in)    :: dcndL(3,3,mol%n) ! derivative of erf-CN
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
   type(chrg_parameter),intent(in) :: chrgeq
   type(TEnvironment), intent(inout) :: env

! ------------------------------------------------------------------------
!  Output
! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(mol%n)          ! partial charges
   real(wp),intent(out)   :: dqdr(3,mol%n,mol%n) ! derivative of partial charges
   real(wp),intent(out)   :: dqdL(3,3,mol%n) ! derivative of partial charges
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,mol%n) ! molecular gradient of IES
   real(wp),intent(inout) :: sigma(3,3) ! molecular gradient of IES
!
! ------------------------------------------------------------------------
!  charge model
! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmatdr(:,:,:)
   real(wp),allocatable :: dAmatdL(:,:,:)
   real(wp),allocatable :: dXvecdr(:,:,:)
   real(wp),allocatable :: dXvecdL(:,:,:)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es

! ------------------------------------------------------------------------
!  for periodic case
! ------------------------------------------------------------------------
   integer  :: wscAt
   real(wp) :: gamii,riw(3),dAtmp(3),stmp(3,3),vec(3)
   real(wp) :: cf ! convergence factor
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), allocatable :: gTrans(:, :)
   real(wp), allocatable :: rTrans(:, :)
   integer :: iG1, iG2, iG3, iT1, iT2, iT3, iRp

! ------------------------------------------------------------------------
!  scratch variables
! ------------------------------------------------------------------------
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)
   real(wp),allocatable :: Xfac(:)
   real(wp),allocatable :: Afac(:,:)

! ------------------------------------------------------------------------
!  Lapack work variables
! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)

! ------------------------------------------------------------------------
!  initizialization
! ------------------------------------------------------------------------
   m    = mol%n+1
   q    = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(m), source = 0.0_wp )

! ------------------------------------------------------------------------
!  set up the A matrix and X vector
! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Setup of the A matrix and the RHS X vector"

   ! prepare some arrays
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol,cn,chrgeq) private(i,tmp) shared(Xvec,Xfac)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
   enddo
   !$omp end parallel do

   if (mol%npbc > 0) then

      iRp = 0
      allocate(gTrans(3, product(2*ewaldCutR+1)-1))
      do iG1 = -ewaldCutR(1), ewaldCutR(1)
         do iG2 = -ewaldCutR(2), ewaldCutR(2)
            do iG3 = -ewaldCutR(3), ewaldCutR(3)
               if (iG1 == 0 .and. iG2 == 0 .and. iG3 == 0) cycle
               iRp = iRp + 1
               vec(:) = [iG1, iG2, iG3]
               gTrans(:, iRp) = matmul(mol%rec_lat, vec)
            end do
         end do
      end do

      iRp = 0
      allocate(rTrans(3, product(2*ewaldCutD+1)))
      do iT1 = -ewaldCutD(1), ewaldCutD(1)
         do iT2 = -ewaldCutD(2), ewaldCutD(2)
            do iT3 = -ewaldCutD(3), ewaldCutD(3)
               iRp = iRp + 1
               vec(:) = [iT1, iT2, iT3]
               rTrans(:, iRp) = matmul(mol%lattice, vec)
            end do
         end do
      end do

      cf = sqrtpi/mol%volume**(1.0_wp/3.0_wp)
      ! build Ewald matrix
      call get_coulomb_matrix(mol, chrgeq, rTrans, gTrans, cf, amat)
   else
      call get_coulomb_matrix(mol, chrgeq, amat)
   endif

! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges"
   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = mol%chrg
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u', m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) then
      call env%error("Coulomb matrix is singular, cannot solve lin. eq.", source)
      return
   end if

   q = Xtmp(:mol%n)
   if(abs(sum(q)-mol%chrg) > 1.e-6_wp) then
      call env%error("Charge constraint is not satisfied", source)
      return
   end if
   !print'(3f20.14)',Xtmp

   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif

! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Isotropic electrostatic (IES) energy calculation"
   call mctc_copy(Xvec, work)
   call mctc_symv(Amat, Xtmp, work, alpha=0.5_wp, beta=-1.0_wp)
   es = mctc_dot(Xtmp, work)
   energy = es + energy
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "energy",es
   endif

! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
! ------------------------------------------------------------------------
do_molecular_gradient: if (lgrad .or. lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmatdr(3,mol%n,m), dXvecdr(3,mol%n,m), Afac(3,mol%n), source = 0.0_wp )
   if (mol%npbc > 0) then
      allocate( dAmatdL(3,3,m), dXvecdL(3,3,m), source = 0.0_wp )
      call get_coulomb_derivs(mol, chrgeq, Xtmp, rTrans, gTrans, cf, &
         & dAmatdr, dAmatdL, Afac)
      do i = 1, mol%n
         dXvecdr(:,:,i) = -dcndr(:,:,i)*Xfac(i)
         dXvecdL(:,:,i) = -dcndL(:,:,i)*Xfac(i)
      enddo
   else
      call get_coulomb_derivs(mol, chrgeq, Xtmp, dAmatdr, Afac)
      do i = 1, mol%n
         dXvecdr(:,:,i) = -dcndr(:,:,i)*Xfac(i) ! merge dX and dA for speedup
      enddo
   endif
endif do_molecular_gradient

   if (lgrad) then
   call mctc_gemv(dAmatdr, Xtmp, gradient, alpha=+1.0_wp, beta=1.0_wp)
   call mctc_gemv(dXvecdr, Xtmp, gradient, alpha=+1.0_wp, beta=1.0_wp)
   if (mol%npbc > 0) then
      call mctc_gemv(dAmatdL, Xtmp, sigma, alpha=+0.5_wp, beta=1.0_wp)
      call mctc_gemv(dXvecdL, Xtmp, sigma, alpha=+1.0_wp, beta=1.0_wp)
   endif
   endif

! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
! ------------------------------------------------------------------------
do_partial_charge_derivative: if (lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "A matrix inversion by Bunch-Kaufman factorization"
   allocate( Ainv(m,m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call env%error("Could not factorize Coulomb matrix", source)
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call env%error("Coulomb Matrix is singular, cannot invert", source)
      return
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges"
   dqdr = 0.0_wp
   if (mol%npbc > 0) then
      dqdL = 0.0_wp
   endif
   do i = 1, mol%n
      dAmatdr(:,i,i) = Afac(:,i) + dAmatdr(:,i,i)
   enddo
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmatdr,3*n,1.0_wp,dqdr,3*n)
   call mctc_gemm(dAmatdr, Ainv, dqdr, alpha=-1.0_wp, beta=1.0_wp)
   call mctc_gemm(dXvecdr, Ainv, dqdr, alpha=-1.0_wp, beta=1.0_wp)
   if (mol%npbc > 0) then
      call mctc_gemm(dAmatdL, Ainv, dqdL, alpha=-1.0_wp, beta=1.0_wp)
      call mctc_gemm(dXvecdL, Ainv, dqdL, alpha=-1.0_wp, beta=1.0_wp)
   endif
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

! ------------------------------------------------------------------------
!  Clean up
! ------------------------------------------------------------------------
   if (allocated(Amat))  deallocate(Amat)
   if (allocated(dAmatdr)) deallocate(dAmatdr)
   if (allocated(Afac))  deallocate(Afac)
   if (allocated(Xvec))  deallocate(Xvec)
   if (allocated(Xfac))  deallocate(Xfac)
   if (allocated(Xtmp))  deallocate(Xtmp)
   if (allocated(Atmp))  deallocate(Atmp)
   if (allocated(temp))  deallocate(temp)
   if (allocated(work))  deallocate(work)
   if (allocated(ipiv))  deallocate(ipiv)

end subroutine eeq_chrgeq_core


! ======================================================================
!  Modified Version of eeq_chrgeq routine that reads also the chrgeq construct
! ======================================================================
subroutine eeq_chrgeq_gbsa(mol,env,chrgeq,gbsa,cn,dcndr,q,dqdr, &
      &                    energy,gsolv,gradient,lverbose,lgrad,lcpq)

   use xtb_mctc_convert

   use xtb_type_molecule
   use xtb_type_param

   use xtb_disp_ncoord
   use xtb_solv_gbsa, only : TBorn

   use xtb_pbc_tools

   implicit none

   character(len=*), parameter :: source = "eeq_chrgeq_gbsa"

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(TMolecule),intent(in) :: mol     !< molecular structure information
   type(TBorn), intent(in) :: gbsa
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   real(wp),intent(in)    :: dcndr(3,mol%n,mol%n) ! derivative of erf-CN
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
   type(chrg_parameter),intent(in) :: chrgeq
   type(TEnvironment), intent(inout) :: env

! ------------------------------------------------------------------------
!  Output
! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(mol%n)          ! partial charges
   real(wp),intent(out)   :: dqdr(3,mol%n,mol%n) ! derivative of partial charges
   real(wp),intent(inout) :: gsolv        ! electrostatic energy
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,mol%n) ! molecular gradient of IES
!
! ------------------------------------------------------------------------
!  charge model
! ------------------------------------------------------------------------
   integer  :: m ! dimension of the Lagrangian
   real(wp),allocatable :: Amat(:,:)
   real(wp),allocatable :: Xvec(:)
   real(wp),allocatable :: Ainv(:,:)
   real(wp),allocatable :: dAmatdr(:,:,:)
   real(wp),allocatable :: dAmatdL(:,:,:)
   real(wp),allocatable :: dXvecdr(:,:,:)
   real(wp),allocatable :: dXvecdL(:,:,:)

! ------------------------------------------------------------------------
!  local variables
! ------------------------------------------------------------------------
   integer  :: i,j,k,l
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es

! ------------------------------------------------------------------------
!  for periodic case
! ------------------------------------------------------------------------
   integer  :: wscAt
   real(wp) :: gamii,riw(3),dAtmp(3),stmp(3,3)
   real(wp) :: cf ! convergence factor
   real(wp),parameter   :: zero(3) = [0.0_wp,0.0_wp,0.0_wp]
   integer              :: rep_cn(3)
   integer, parameter   :: ewaldCutD(3) = [2,2,2]
   integer, parameter   :: ewaldCutR(3) = [2,2,2]

! ------------------------------------------------------------------------
!  scratch variables
! ------------------------------------------------------------------------
   real(wp),allocatable :: xtmp(:)
   real(wp),allocatable :: atmp(:,:)
   real(wp),allocatable :: Xfac(:)
   real(wp),allocatable :: Afac(:,:)

! ------------------------------------------------------------------------
!  Lapack work variables
! ------------------------------------------------------------------------
   integer, allocatable :: ipiv(:)
   real(wp),allocatable :: temp(:)
   real(wp),allocatable :: work(:)
   integer  :: lwork
   integer  :: info
   real(wp) :: test(1)

   real(wp) :: gborn, ghb, r3, erft

   ! quick return if possible
   if (mol%n.eq.1) then
      q = mol%chrg
      ! should always hold, even for extendend systems
      es = 0.0_wp
      return
   endif

   gborn = 0.0_wp
   ghb   = 0.0_wp

! ------------------------------------------------------------------------
!  initizialization
! ------------------------------------------------------------------------
   m    = mol%n+1
   q    = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(m), source = 0.0_wp )

! ------------------------------------------------------------------------
!  set up the A matrix and X vector
! ------------------------------------------------------------------------
!  αi -> alpha(i), ENi -> xi(i), κi -> kappa(i), Jii -> gam(i)
!  γij = 1/√(αi+αj)
!  Xi  = -ENi + κi·√CNi
!  Aii = Jii + 2/√π·γii
!  Aij = erf(γij·Rij)/Rij = 2/√π·F0(γ²ij·R²ij)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Setup of the A matrix and the RHS X vector"

   ! prepare some arrays
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol,cn,chrgeq) private(i,tmp) shared(Xvec,Xfac)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
   enddo
   !$omp end parallel do

   call get_coulomb_matrix(mol, chrgeq, Amat)

   Amat(:mol%n, :mol%n) = Amat(:mol%n, :mol%n) + gbsa%bornMat

! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges"
   Amat(m,1:m) = 1.0_wp
   Amat(1:m,m) = 1.0_wp
   Amat(m,m  ) = 0.0_wp
   Xvec(m)     = mol%chrg
   ! generate temporary copy
   allocate( Atmp(m,m), source = Amat )
   allocate( Xtmp(m),   source = Xvec )

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('u', m, 1, Atmp, m, ipiv, Xtmp, m, test, -1, info)
   lwork = int(test(1))
   allocate( work(lwork), source = 0.0_wp )

   call dsysv('u',m,1,Atmp,m,ipiv,Xtmp,m,work,lwork,info)
   if(info > 0) then
      call env%error("Coulomb matrix is singular, cannot solve lin. eq.", source)
      return
   end if

   q = Xtmp(:mol%n)
   if(abs(sum(q)-mol%chrg) > 1.e-6_wp) then
      call env%error("Charge constraint is not satisfied", source)
      return
   end if
   !print'(3f20.14)',Xtmp

   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif

! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Isotropic electrostatic (IES) energy calculation"
   call mctc_copy(Xvec, work)
   call mctc_symv(Amat, Xtmp, work, alpha=0.5_wp, beta=-1.0_wp)
   es = mctc_dot(Xtmp, work)
   energy = es + energy
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "energy",es
   endif

! ------------------------------------------------------------------------
!  calculate molecular gradient of the IES energy
! ------------------------------------------------------------------------
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ½∂Aij/∂Rj -> dAmat(:,j,i)
!  dE/dR = (½∂A/∂R·q - ∂X/∂R)·q
!  ∂Aij/∂Rj = ∂Aij/∂Ri
! ------------------------------------------------------------------------
do_molecular_gradient: if (lgrad .or. lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmatdr(3,mol%n,m), dXvecdr(3,mol%n,m), Afac(3,mol%n), source = 0.0_wp )

   call get_coulomb_derivs(mol, chrgeq, Xtmp, dAmatdr, Afac)
   do i = 1, mol%n
      dAmatdr(:,:,i) = dAmatdr(:,:,i) - dcndr(:,:,i)*Xfac(i)
   enddo
   call gbsa%addBornDeriv(Xtmp, gborn, ghb, dAmatdr, Afac)
   gsolv = gsolv + gborn + ghb + gbsa%gshift
endif do_molecular_gradient

   if (lgrad) then
   call mctc_gemv(dAmatdr, Xtmp, gradient, alpha=1.0_wp, beta=1.0_wp)
   endif

! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
! ------------------------------------------------------------------------
do_partial_charge_derivative: if (lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "A matrix inversion by Bunch-Kaufman factorization"
   allocate( Ainv(m,m), source = Amat )

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call env%error("Could not factorize Coulomb matrix", source)
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call env%error("Coulomb Matrix is singular, cannot invert", source)
      return
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do i = 1, m
      do j = i+1, m
         Ainv(i,j)=Ainv(j,i)
      enddo
   enddo

! ------------------------------------------------------------------------
!  calculate gradient of the partial charge w.r.t. the nuclear coordinates
! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges"
   dqdr = 0.0_wp
   do i = 1, mol%n
      dAmatdr(:,i,i) = Afac(:,i) + dAmatdr(:,i,i)
   enddo
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmatdr,3*n,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*mol%n,mol%n,m,-1.0_wp,dAmatdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

! ------------------------------------------------------------------------
!  Clean up
! ------------------------------------------------------------------------
   if (allocated(Amat))  deallocate(Amat)
   if (allocated(dAmatdr)) deallocate(dAmatdr)
   if (allocated(Afac))  deallocate(Afac)
   if (allocated(Xvec))  deallocate(Xvec)
   if (allocated(Xfac))  deallocate(Xfac)
   if (allocated(Xtmp))  deallocate(Xtmp)
   if (allocated(Atmp))  deallocate(Atmp)
   if (allocated(temp))  deallocate(temp)
   if (allocated(work))  deallocate(work)
   if (allocated(ipiv))  deallocate(ipiv)

end subroutine eeq_chrgeq_gbsa

end module xtb_eeq
