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

module eeq_model
   use iso_fortran_env, only : wp => real64
   use mctc_constants
   use gfn0param, alp => alpg
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
   use tbdef_molecule
   use tbdef_param
   type(tb_molecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(out) :: amat(:, :)
   integer :: iat, jat
   real(wp) :: r1, gamij
   Amat = 0.0_wp
   !$omp parallel do default(none) shared(mol,chrgeq, amat) &
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

subroutine get_coulomb_matrix_3d(mol, chrgeq, cf, amat)
   use tbdef_molecule
   use tbdef_param
   type(tb_molecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: cf
   real(wp), intent(out) :: amat(:, :)
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp
   integer :: iat, jat, wscAt
   real(wp) :: gamii, gamij, riw(3)
   amat = 0.0_wp
   !$omp parallel do default(none) reduction(+:amat) &
   !$omp shared(mol, chrgeq, cf) private(jat, wscAt, gamii, gamij, riw)
   do iat = 1, len(mol)
      gamii = 1.0_wp/(sqrt(2.0_wp)*chrgeq%alpha(iat))
      amat(iat, iat) = chrgeq%gam(iat) + sqrt2pi/chrgeq%alpha(iat) &
         ! reciprocal for 0th atom
         + eeq_ewald_3d_rec(zero,ewaldCutR,mol%rec_lat,mol%volume,cf) &
         ! direct for 0th atom
         + eeq_ewald_3d_dir(zero,ewaldCutD,mol%lattice,gamii,cf)
      do jat = 1, iat-1
         gamij = 1.0_wp/sqrt(chrgeq%alpha(iat)**2+chrgeq%alpha(jat)**2)
         do wscAt = 1, mol%wsc%itbl(jat,iat)
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            amat(iat,jat) = Amat(iat,jat) + mol%wsc%w(jat,iat) * ( &
               ! reciprocal lattice sums
               + eeq_ewald_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf) &
               ! direct lattice sums
               + eeq_ewald_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf))
         end do
         amat(jat,iat) = amat(iat,jat)
      end do
   end do
   !$omp end parallel do
end subroutine get_coulomb_matrix_3d

subroutine get_coulomb_derivs_0d(mol, chrgeq, qvec, amatdr, atrace)
   use tbdef_molecule
   use tbdef_param
   type(tb_molecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: amatdr(:, :, :)
   real(wp), intent(out) :: atrace(:, :)
   integer :: iat, jat
   real(wp) :: rij(3), r2, gamij2, arg2
   real(wp) :: dE, dG(3)
   amatdr = 0.0_wp
   atrace = 0.0_wp
   !$omp parallel do default(none) reduction(+:atrace, amatdr) &
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

subroutine get_coulomb_derivs_3d(mol, chrgeq, qvec, cf, amatdr, amatdL, atrace)
   use tbdef_molecule
   use tbdef_param
   type(tb_molecule), intent(in) :: mol
   type(chrg_parameter), intent(in) :: chrgeq
   real(wp), intent(in) :: cf
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: amatdr(:, :, :)
   real(wp), intent(out) :: amatdL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: zero(3) = 0.0_wp
   integer :: iat, jat, wscAt, ii
   real(wp) :: gamij, dG(3), dS(3, 3), riw(3)
   amatdr = 0.0_wp
   amatdL = 0.0_wp
   atrace = 0.0_wp
   !$omp parallel do default(none) &
   !$omp reduction(+:atrace,amatdr,amatdL) shared(mol, chrgeq, qvec, cf) &
   !$omp private(iat, jat, ii, wscAt, riw, gamij, dG, dS)
   do iat = 1, len(mol)
      do jat = 1, iat-1
         ! over WSC partner
         gamij = 1.0_wp/sqrt(chrgeq%alpha(iat)**2 + chrgeq%alpha(jat)**2)
         do wscAt = 1, mol%wsc%itbl(jat,iat)
            riw = mol%xyz(:,iat) - mol%xyz(:,jat) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,jat,iat))
            call eeq_ewald_dx_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf, &
               &                     dG,dS)
            dG = dG*mol%wsc%w(jat,iat)
            dS = dS*mol%wsc%w(jat,iat)
            amatdr(:, iat, jat) = amatdr(:, iat, jat) + dG*qvec(iat)
            amatdr(:, jat, iat) = amatdr(:, jat, iat) - dG*qvec(jat)
            atrace(:, iat) = atrace(:, iat) + dG*qvec(jat)
            atrace(:, jat) = atrace(:, jat) - dG*qvec(iat)
            amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(jat)
            amatdL(:, :, jat) = amatdL(:, :, jat) + dS*qvec(iat)
            call eeq_ewald_dx_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf, &
               &                     dG,dS)
            dG = dG*mol%wsc%w(jat, iat)
            dS = dS*mol%wsc%w(jat, iat)
            amatdr(:, iat, jat) = amatdr(:, iat, jat) + dG*qvec(iat)
            amatdr(:, jat, iat) = amatdr(:, jat, iat) - dG*qvec(jat)
            atrace(:, iat) = atrace(:, iat) + dG*qvec(jat)
            atrace(:, jat) = atrace(:, jat) - dG*qvec(iat)
            amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(jat)
            amatdL(:, :, jat) = amatdL(:, :, jat) + dS*qvec(iat)
         enddo  ! k WSC partner
      enddo     ! jat

      gamij = 1.0_wp/(sqrt(2.0_wp)*chrgeq%alpha(iat))
      call eeq_ewald_dx_3d_rec(zero, ewaldCutR, mol%rec_lat, mol%volume, cf, &
         &                     dG, dS)
      amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(iat)
      call eeq_ewald_dx_3d_dir(zero, ewaldCutD, mol%lattice, gamij, cf, &
         &                     dG, dS)
      amatdL(:, :, iat) = amatdL(:, :, iat) + dS*qvec(iat)
      do ii = 1, 3
         amatdL(ii, ii, iat) = amatdL(ii, ii, iat) + cf/sqrtpi/3.0_wp*qvec(iat)
      enddo
   enddo
   !$omp end parallel do
end subroutine get_coulomb_derivs_3d

subroutine print_chrgeq(iunit,n,at,xyz,q,cn,dipm)
   use mctc_constants
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: q(n)
   real(wp),intent(in) :: cn(n)
   real(wp),intent(in),optional :: dipm(3,n)
   character(len=2),external :: asym
   real(wp) :: mmom(3),dmom(3)
   integer  :: i
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   write(iunit,'(a)')
   write(iunit,'(7x,"   #   Z   ")',advance='no')
   write(iunit,'("         q")',advance='no')
   if (present(dipm)) &
   write(iunit,'("      mux      muy      muz")',advance='no')
   write(iunit,'("        CN")',advance='no')
   write(iunit,'("        EN")',advance='no')
   write(iunit,'("       Aii")',advance='no')
   write(iunit,'(a)')
   do i=1,n
      write(iunit,'(i11,1x,i3,1x,a2)',advance='no') &
      &     i,at(i),asym(at(i))
      write(iunit,'(f10.3)',advance='no')q(i)
      if (present(dipm)) &
      write(iunit,'(3f9.3)',advance='no')dipm(:,i)
      write(iunit,'(f10.3)',advance='no')cn(i)
      write(iunit,'(f10.3)',advance='no')xi(at(i)) - cnfak(at(i))*sqrt(cn(i))
      write(iunit,'(f10.3)',advance='no')gamm(at(i))+sqrt2pi/alp(at(i))
      write(iunit,'(a)')
   enddo
   mmom = 0.0_wp
   do i = 1, n
      mmom = mmom + q(i)*xyz(:,i)
   enddo
   if (present(dipm)) then
      dmom = 0.0_wp
      do i = 1, n
         dmom = dmom + dipm(:,i)
      enddo
   endif
   write(iunit,'(a)')
   write(iunit,'(7x,a)')'dipole moment:'
   write(iunit,'(18x,a)')'x           y           z       tot (au)'
   if (.not.present(dipm)) then
      write(iunit,'(7x,4f12.3)')  mmom,norm2(mmom)
   else
      write(iunit,'(1x,"q only",3f12.3)')  mmom
      write(iunit,'(3x,  "full",4f12.3)')  mmom+dmom,norm2(mmom+dmom)
   endif
   write(iunit,'(a)')

end subroutine print_chrgeq

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
   use iso_fortran_env, wp => real64, istdout => output_unit
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
!$omp do
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
!$omp do
   ! prepare A matrix
   do i = 1, n
      ! EN of atom i
      do j = 1, i-1
         r = norm2(xyz(:,j) - xyz(:,i))
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
   if(info > 0) call raise('E','(goedecker_solve) DSYSV failed',1)

   q = Xtmp(:n)
   if(abs(sum(q)-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error',1)
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
   work(:m) = Xvec
   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
   es = dot_product(Xtmp,work(:m))
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
!$omp do
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
   call dgemv('n',3*n,m,+1.0_wp,dAmat,3*n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*n,m,-1.0_wp,dXvec,3*n,Xtmp,1,1.0_wp,gradient,1)
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
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(goedecker_inversion) DSYTRF failed',1)
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(goedecker_inversion) DSYTRI failed',1)
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
   !call dsymm('r','l',3*n,m,-1.0_wp,Ainv,m,dAmat,3*n,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,n,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,n,m,+1.0_wp,dXvec,3*n,Ainv,m,1.0_wp,dqdr,3*n)
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

subroutine eeq_chrgeq_qonly(mol,chrgeq,cn,q,lverbose)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use tbdef_molecule
   use tbdef_param

   implicit none

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(tb_molecule),intent(in) :: mol     !< molecular structure information
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   type(chrg_parameter),intent(in) :: chrgeq
   real(wp),allocatable   :: dcndr(:,:,:)
   real(wp),allocatable   :: dcndL(:,:,:)
   logical, intent(in), optional :: lverbose
   logical :: verbose

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

   call eeq_chrgeq_core(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &                 energy,gradient,sigma,verbose,.false.,.false.)

end subroutine eeq_chrgeq_qonly

! ======================================================================
!  Modified Version of eeq_chrgeq routine that reads also the chrgeq construct
! ======================================================================
subroutine eeq_chrgeq_core(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL, &
      &                    energy,gradient,sigma,lverbose,lgrad,lcpq)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use mctc_econv

   use tbdef_molecule
   use tbdef_param

   use ncoord

   use pbc_tools

   implicit none

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(tb_molecule),intent(in) :: mol     !< molecular structure information
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   real(wp),intent(in)    :: dcndr(3,mol%n,mol%n) ! derivative of erf-CN
   real(wp),intent(in)    :: dcndL(3,3,mol%n) ! derivative of erf-CN
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
   type(chrg_parameter),intent(in) :: chrgeq

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
   real(wp) :: gamii,riw(3),dAtmp(3),stmp(3,3)
   real(wp) :: cf ! convergence factor

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

   ! quick return if possible
   if (mol%n.eq.1) then
      q = mol%chrg
      ! should always hold, even for extendend systems
      es = 0.0_wp
      return
   endif

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
   !$omp parallel do default(none) &
   !$omp shared(mol,cn,chrgeq) private(i,tmp) shared(Xvec,Xfac)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
   enddo
   !$omp end parallel do

   if (mol%npbc > 0) then
      cf = sqrtpi/mol%volume**(1.0_wp/3.0_wp)
      ! build Ewald matrix
      call get_coulomb_matrix_3d(mol, chrgeq, cf, amat)
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
   if(info > 0) call raise('E','(eeq_solve) DSYSV failed',1)

   q = Xtmp(:mol%n)
   if(abs(sum(q)-mol%chrg) > 1.e-6_wp) &
      call raise('E','(eeq_solve) charge constrain error',1)
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
   work(:m) = Xvec
   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
   es = dot_product(Xtmp,work(:m))
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
      call get_coulomb_derivs(mol, chrgeq, Xtmp, cf, dAmatdr, dAmatdL, Afac)
      do i = 1, mol%n
         dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i)
         dXvecdL(:,:,i) = +dcndL(:,:,i)*Xfac(i)
      enddo
   else
      call get_coulomb_derivs(mol, chrgeq, Xtmp, dAmatdr, Afac)
      do i = 1, mol%n
         dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i) ! merge dX and dA for speedup
      enddo
   endif
endif do_molecular_gradient

   if (lgrad) then
   call dgemv('n',3*mol%n,m,+1.0_wp,dAmatdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*mol%n,m,+1.0_wp,dXvecdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
   if (mol%npbc > 0) then
      call dgemv('n',3*3,m,+0.5_wp,dAmatdL,3*3,Xtmp,1,1.0_wp,sigma,1)
      call dgemv('n',3*3,m,-1.0_wp,dXvecdL,3*3,Xtmp,1,1.0_wp,sigma,1)
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
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(eeq_inversion) DSYTRF failed',1)
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(eeq_inversion) DSYTRI failed',1)
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
   call dgemm('n','n',3*mol%n,mol%n,m,-1.0_wp,dAmatdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   call dgemm('n','n',3*mol%n,mol%n,m,-1.0_wp,dXvecdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   if (mol%npbc > 0) then
      call dgemm('n','n',9,mol%n,m,-1.0_wp,dAmatdL,9,Ainv,m,1.0_wp,dqdL,3*3)
      call dgemm('n','n',9,mol%n,m,+1.0_wp,dXvecdL,9,Ainv,m,1.0_wp,dqdL,3*3)
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

pure function eeq_ewald_3d_dir(riw,rep,dlat,gamij,cf) result(Amat)
   use iso_fortran_env, wp => real64
   use mctc_constants
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz
   real(wp) :: distiw,rij(3)
   real(wp) :: t(3)
   Amat = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
      distiw = norm2(rij)
      ! self-interaction correction
      if(distiw < eps) then
         Amat = Amat - cf/sqrtpi
      else
         Amat = Amat - erf(   cf*distiw)/distiw &
            &        + erf(gamij*distiw)/distiw
      end if
   end do
end function eeq_ewald_3d_dir

pure subroutine eeq_ewald_dx_3d_dir(riw,rep,dlat,gamij,cf,dAmat,sigma)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use pbc_tools
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: dlat(3,3) !< direct lattice
   real(wp),intent(in) :: gamij     !< interaction radius
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dAmat(3) !< element of interaction matrix
   real(wp),intent(out) :: sigma(3,3)
   real(wp),parameter :: eps = 1.0e-9_wp
   integer  :: dx,dy,dz,i
   real(wp) :: distiw,rij(3),arga,argb
   real(wp) :: t(3),dtmp,stmp(3)
   dAmat = 0.0_wp
   sigma = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      ! real contributions
      t = [dx,dy,dz]
      rij = riw + matmul(dlat,t)
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

end subroutine eeq_ewald_dx_3d_dir

pure function eeq_ewald_3d_rec(riw,rep,rlat,vol,cf) result(Amat)
   use iso_fortran_env, wp => real64
   use mctc_constants
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp) :: Amat                 !< element of interaction matrix
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3)
   Amat = 0.0_wp
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      if (dx==0 .and. dy==0 .and. dz==0) cycle
      t = [dx,dy,dz]
      rik = matmul(rlat,t)
      rik2 = dot_product(rik,rik)
      Amat=Amat + cos(dot_product(rik,riw) ) * 4.0_wp*pi/vol &
         * exp(-rik2/(4.0_wp*cf**2))/rik2
   end do
end function eeq_ewald_3d_rec

pure subroutine eeq_ewald_dx_3d_rec(riw,rep,rlat,vol,cf,dAmat,sigma)
   use iso_fortran_env, wp => real64
   use mctc_constants
   use pbc_tools
   implicit none
   real(wp),intent(in) :: riw(3)    !< distance from i to WSC atom
   integer, intent(in) :: rep(3)    !< images to consider
   real(wp),intent(in) :: rlat(3,3) !< reciprocal lattice
   real(wp),intent(in) :: vol       !< direct cell volume
   real(wp),intent(in) :: cf        !< convergence factor
   real(wp),intent(out) :: dAmat(3) !< element of interaction matrix
   real(wp),intent(out) :: sigma(3,3)
   integer  :: dx,dy,dz
   real(wp) :: rik2,rik(3)
   real(wp) :: t(3),dtmp,fpivol
   real(wp) :: expterm,arg
   integer  :: i
   dAmat = 0.0_wp
   sigma = 0.0_wp
   fpivol = 4.0_wp*pi/vol
   do concurrent(dx = -rep(1):rep(1), &
         &       dy = -rep(2):rep(2), &
         &       dz = -rep(3):rep(3))
      if (dx==0 .and. dy==0 .and. dz==0) cycle
      t = [dx,dy,dz]
      rik = matmul(rlat,t)
      rik2 = dot_product(rik,rik)
      expterm = exp(-rik2/(4.0_wp*cf**2))/rik2
      arg = dot_product(rik,riw)
      ! d/dx (sin**2 + cos**2) = -2*sin*cos - 2*cos*sin
      dtmp = -sin(arg) * expterm * fpivol
      dAmat = dAmat + rik*dtmp
      sigma = sigma + fpivol * expterm * cos(arg) * (&
         & reshape([-1.0_wp, 0.0_wp, 0.0_wp, &
         &           0.0_wp,-1.0_wp, 0.0_wp, &
         &           0.0_wp, 0.0_wp,-1.0_wp],shape(sigma)) &
         & * (1.0_wp + rik2/(4.0_wp*cf**2)*2.0_wp/3.0_wp) &
         & + (2.0_wp/rik2 + 0.5_wp/cf**2) * outer_prod_3x3(rik,rik))
   end do
end subroutine eeq_ewald_dx_3d_rec

! ======================================================================
!  Modified Version of eeq_chrgeq routine that reads also the chrgeq construct
! ======================================================================
subroutine eeq_chrgeq_gbsa(mol,chrgeq,gbsa,cn,dcndr,q,dqdr, &
      &                    energy,gsolv,gradient,lverbose,lgrad,lcpq)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use mctc_econv

   use tbdef_molecule
   use tbdef_param

   use ncoord
   use gbobc

   use pbc_tools

   implicit none

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(tb_molecule),intent(in) :: mol     !< molecular structure information
   type(tb_solvent), intent(in) :: gbsa
   real(wp),intent(in)    :: cn(mol%n)     ! erf-CN
   real(wp),intent(in)    :: dcndr(3,mol%n,mol%n) ! derivative of erf-CN
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
   type(chrg_parameter),intent(in) :: chrgeq

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
   !$omp parallel do default(none) &
   !$omp shared(mol,cn,chrgeq) private(i,tmp) shared(Xvec,Xfac)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
   enddo
   !$omp end parallel do

   call get_coulomb_matrix(mol, chrgeq, Amat)

   call compute_amat(gbsa,Amat)

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
   if(info > 0) call raise('E','(eeq_solve) DSYSV failed',1)

   q = Xtmp(:mol%n)
   if(abs(sum(q)-mol%chrg) > 1.e-6_wp) &
      call raise('E','(eeq_solve) charge constrain error',1)
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
   work(:m) = Xvec
   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
   es = dot_product(Xtmp,work(:m))
   energy = es + energy + gshift
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
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + dcndr(:,:,i)*Xfac(i)
   enddo
   call compute_gb_damat(gbsa,Xtmp,gborn,ghb,dAmatdr,Afac,lverbose)
   gsolv = gsolv + gborn + ghb + gshift
endif do_molecular_gradient

   if (lgrad) then
   call dgemv('n',3*mol%n,m,+1.0_wp,dAmatdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
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
   call dsytrf('L',m,Ainv,m,ipiv,test,-1,info)
   if (int(test(1)) > lwork) then
      deallocate(work)
      lwork=int(test(1))
      allocate( work(lwork), source = 0.0_wp )
   endif

   ! Bunch-Kaufman factorization A = L*D*L**T
   call dsytrf('L',m,Ainv,m,ipiv,work,lwork,info)
   if(info > 0)then
      call raise('E', '(eeq_inversion) DSYTRF failed',1)
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
   if (info > 0) then
      call raise('E', '(eeq_inversion) DSYTRI failed',1)
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

end module eeq_model
