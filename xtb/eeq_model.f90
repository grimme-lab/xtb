! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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
use gfn0param, alp => alpg
implicit none

interface eeq_chrgeq
   module procedure :: eeq_chrgeq_core
   module procedure :: eeq_chrgeq_gbsa
   module procedure :: eeq_chrgeq_qonly
end interface eeq_chrgeq

contains

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
!   real(wp),intent(in)    :: par_screen    ! dielectric screening
!   real(wp),intent(in)    :: par_gscale    ! scaling for hardnesses
!   real(wp),intent(in)    :: par_rscale    ! scaling for radii
   logical, intent(in)    :: lverbose      ! toggles printout
   logical, intent(in)    :: lgrad         ! flag for gradient calculation
   logical, intent(in)    :: lcpq          ! do partial charge derivative
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(n)          ! partial charges
   real(wp),intent(out)   :: dqdr(3,n,n+1) ! derivative of partial charges
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,n) ! molecular gradient of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
!$omp do schedule(dynamic)
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
!$omp do schedule(dynamic)
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
!$omp do schedule(dynamic)
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
   call dgemm('n','n',3*n,m,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,m,m,+1.0_wp,dXvec,3*n,Ainv,m,1.0_wp,dqdr,3*n)
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

subroutine goedecker_hessian(n,at,xyz,chrg,cn,dcndr,q,dqdr,energy,gradient,hessian,&
                             lverbose,lgrad,lcpq,lhess)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use ncoord, only : rcov, kn
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                ! number of atoms
   integer, intent(in)    :: at(n)            ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)         ! geometry
   real(wp),intent(in)    :: chrg             ! total charge
   real(wp),intent(in)    :: cn(n)            ! erf-CN
   real(wp),intent(in)    :: dcndr(3,n,n)     ! derivative of erf-CN
!   real(wp),intent(in)    :: par_screen       ! dielectric screening
!   real(wp),intent(in)    :: par_gscale       ! scaling for hardnesses
!   real(wp),intent(in)    :: par_rscale       ! scaling for radii
   logical, intent(in)    :: lverbose         ! toggles printout
   logical, intent(in)    :: lgrad            ! flag for gradient calculation
   logical, intent(in)    :: lcpq             ! do partial charge derivative
   logical, intent(in)    :: lhess            ! calculate molecular hessian
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(n)             ! partial charges
   real(wp),intent(out)   :: dqdr(3,n,n+1)    ! derivative of partial charges
   real(wp),intent(inout) :: energy           ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,n)    ! molecular gradient of IES
   real(wp),intent(inout) :: hessian(3,n,3,n) ! molecular hessian of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
   real(wp) :: arg,arg2,tmp,dtmp
   real(wp) :: lambda
   real(wp) :: es,expterm,erfterm
   real(wp) :: htmp,rxr(3,3)
   real(wp) :: rcovij,rr

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
   dqdr = 0.0_wp
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
!$omp do schedule(dynamic)
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
!$omp do schedule(dynamic)
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
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmat(3,n,m), dXvec(3,n,m), Afac(3,n), source = 0.0_wp )
   !allocate( dAmat(3,n,m), Afac(3,n), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(n,dcndr,xyz,alpha,Amat,Xfac,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp shared(dXvec,dAmat) &
!$omp reduction(+:Afac)
!$omp do schedule(dynamic)
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
   call dgemv('n',3*n,m,+1.0_wp,dAmat,3*n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*n,m,-1.0_wp,dXvec,3*n,Xtmp,1,1.0_wp,gradient,1)

!! ------------------------------------------------------------------------
!  invert the A matrix using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
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
   call dgemm('n','n',3*n,m,m,-1.0_wp,dAmat,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   call dgemm('n','n',3*n,m,m,+1.0_wp,dXvec,3*n,Ainv,m,1.0_wp,dqdr,3*n)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

!! ------------------------------------------------------------------------
!  molecular Hessian calculation
!! ------------------------------------------------------------------------
do_molecular_hessian: if (lhess) then
   do i = 1, n
      do j = 1, i-1
         rij = xyz(:,j) - xyz(:,i)
         r2 = sum(rij**2)
         r = sqrt(r2)
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = q(i)*q(j)*erf(arg)/r
         expterm = q(i)*q(j)*2*gamij*exp(-arg2)/sqrtpi
         ! ∂²(qAq)/(∂Ri∂Rj):
         ! ∂²(qAq)/(∂Xi∂Xi) = (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - (R²ij-3X²ij) erf[γij·Rij]/R⁵ij 
         ! ∂²(qAq)/(∂Xi∂Xj) = (R²ij-3X²ij) erf[γij·Rij]/R⁵ij 
         !                  - (1-3X²ij/R²ij-2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yi) = 3X²ij erf[γij·Rij]/R⁵ij
         !                  - (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         ! ∂²(qAq)/(∂Xi∂Yj) = (3X²ij/R²ij+2γ²ijX²ij) 2γij/√π exp[-γ²ij·R²ij]/R²ij
         !                  - 3X²ij erf[γij·Rij]/R⁵ij 
         rxr(1,1) = erfterm * ( 3*rij(1)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(1)**2/r2**2 + 2*gamij2*rij(1)**2/r2 - 1/r2 )
         rxr(2,2) = erfterm * ( 3*rij(2)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(2)**2/r2**2 + 2*gamij2*rij(2)**2/r2 - 1/r2 )
         rxr(3,3) = erfterm * ( 3*rij(3)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(3)**2/r2**2 + 2*gamij2*rij(3)**2/r2 - 1/r2 )
         rxr(2,1) = erfterm * 3*rij(2)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(2)*rij(1)/r2**2 + 2*gamij2*rij(2)*rij(1)/r2 )
         rxr(3,1) = erfterm * 3*rij(3)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(1)/r2**2 + 2*gamij2*rij(3)*rij(1)/r2 )
         rxr(3,2) = erfterm * 3*rij(3)*rij(2)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(2)/r2**2 + 2*gamij2*rij(3)*rij(2)/r2 )
         ! ∂CNj/∂Zi(Yij/R²ij)
         rcovij = rcov(at(j)) + rcov(at(i))
         rr = r/rcovij
         do k = 1, n
            !rij = xyz(:,k) - xyz(:,i)
            !r2 = sum(rij**2)
            !r = sqrt(r2)
            tmp = 0.5_wp*Xfac(k)/(cn(k)+1e-14_wp)*Xtmp(k)
            rxr(1,1) = rxr(1,1) + tmp*dcndr(1,k,i)*dcndr(1,k,j)
            rxr(2,2) = rxr(2,2) + tmp*dcndr(2,k,i)*dcndr(2,k,j)
            rxr(3,3) = rxr(3,3) + tmp*dcndr(3,k,i)*dcndr(3,k,j)
            rxr(2,1) = rxr(2,1) + tmp*dcndr(1,k,i)*dcndr(2,k,j)
            rxr(3,1) = rxr(3,1) + tmp*dcndr(1,k,i)*dcndr(3,k,j)
            rxr(3,2) = rxr(3,2) + tmp*dcndr(2,k,i)*dcndr(3,k,j)
         enddo
         !tmp = -0.5_wp*Xfac(i)/(cn(i)+1e-14_wp)*Xtmp(i)
         !dtmp = 2.0_wp*kn**2*(rr-1.0_wp)*rr
         !rxr(1,1) = rxr(1,1) + tmp*(dcndr(1,j,i)*rij(1)*dtmp/r2 &
         !                    - dcndr(2,j,i)*rij(2)/r2 - dcndr(3,j,i)*rij(3)/r2)
         !rxr(2,2) = rxr(2,2) + tmp*(dcndr(2,j,i)*rij(2)*dtmp/r2 &
         !                    - dcndr(1,j,i)*rij(1)/r2 - dcndr(3,j,i)*rij(3)/r2)
         !rxr(3,3) = rxr(3,3) + tmp*(dcndr(3,j,i)*rij(3)*dtmp/r2 &
         !                    - dcndr(1,j,i)*rij(1)/r2 - dcndr(2,j,i)*rij(2)/r2)
         !rxr(2,1) = rxr(2,1) + tmp*dcndr(1,j,i)*rij(2)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(2,j,i)*rij(1)*(1.0_wp+dtmp)/r2
         !rxr(3,1) = rxr(3,1) + tmp*dcndr(1,j,i)*rij(3)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(3,j,i)*rij(1)*(1.0_wp+dtmp)/r2
         !rxr(3,2) = rxr(3,2) + tmp*dcndr(2,j,i)*rij(3)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(3,j,i)*rij(2)*(1.0_wp+dtmp)/r2
         !tmp = -0.5_wp*Xfac(j)/(cn(j)+1e-14_wp)*Xtmp(j)
         !rxr(1,1) = rxr(1,1) + tmp*(dcndr(1,j,i)*rij(1)*dtmp/r2 &
         !                    - dcndr(2,j,i)*rij(2)/r2 - dcndr(3,j,i)*rij(3)/r2)
         !rxr(2,2) = rxr(2,2) + tmp*(dcndr(2,j,i)*rij(2)*dtmp/r2 &
         !                    - dcndr(1,j,i)*rij(1)/r2 - dcndr(3,j,i)*rij(3)/r2)
         !rxr(3,3) = rxr(3,3) + tmp*(dcndr(3,j,i)*rij(3)*dtmp/r2 &
         !                    - dcndr(1,j,i)*rij(1)/r2 - dcndr(2,j,i)*rij(2)/r2)
         !rxr(2,1) = rxr(2,1) + tmp*dcndr(1,j,i)*rij(2)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(2,j,i)*rij(1)*(1.0_wp+dtmp)/r2
         !rxr(3,1) = rxr(3,1) + tmp*dcndr(1,j,i)*rij(3)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(3,j,i)*rij(1)*(1.0_wp+dtmp)/r2
         !rxr(3,2) = rxr(3,2) + tmp*dcndr(2,j,i)*rij(3)*(1.0_wp+dtmp)/r2 &
         !                    + tmp*dcndr(3,j,i)*rij(2)*(1.0_wp+dtmp)/r2
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)
         hessian(:,i,:,i) = hessian(:,i,:,i) + rxr
         hessian(:,j,:,j) = hessian(:,j,:,j) + rxr
         hessian(:,i,:,j) = hessian(:,i,:,j) - rxr
         hessian(:,j,:,i) = hessian(:,j,:,i) - rxr
      enddo
   enddo

   ! ∂²(qA)/(∂Ri∂q)·∂q/∂Rj
   hessian = hessian + reshape(matmul(reshape(dqdr,(/3*n,m/)),&
      transpose(reshape(dAmat,(/3*n,m/)))),(/3,n,3,n/))
   ! ∂²X/(∂Ri∂q)·∂q/∂Rj
   hessian = hessian - reshape(matmul(reshape(dqdr,(/3*n,m/)),&
      transpose(reshape(dXvec,(/3*n,m/)))),(/3,n,3,n/))
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dqdr,3*n,dAmat,3*n,1.0_wp,hessian,3*n)
   !call dgemm('n','t',3*n,m,3*n,+1.0_wp,dAmat,3*n,dqdr,3*n,1.0_wp,hessian,3*n)
   !call dgemm('n','t',3*n,m,3*n,-1.0_wp,dqdr,3*n,dXvec,m,1.0_wp,hessian,3*n)

endif do_molecular_hessian

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

end subroutine goedecker_hessian

subroutine goedecker_multieq2(n,at,xyz,chrgeq,chrg,cn,dcndr,q,dipm,energy,&
                              gradient,dqdr,lverbose,lgrad,lcpq,lqonly)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use ncoord, only : rcov, kn
   use tbdef_param
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                ! number of atoms
   integer, intent(in)    :: at(n)            ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)         ! geometry
   real(wp),intent(in)    :: chrg             ! total charge
   real(wp),intent(in)    :: cn(n)            ! erf-CN
   real(wp),intent(in)    :: dcndr(3,n,n)     ! derivative of erf-CN
   logical, intent(in)    :: lverbose         ! toggles printout
   logical, intent(in)    :: lgrad
   logical, intent(in)    :: lcpq
   logical, intent(in)    :: lqonly
   type(chrg_parameter),intent(in) :: chrgeq
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(n)             ! partial charges
   real(wp),intent(out)   :: dipm(3,n)        ! atomic dipole moments
   real(wp),intent(inout) :: energy           ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,n)
   real(wp),intent(out)   :: dqdr(3,n,4*n+1)

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
   character(len=2),external :: asym
   integer  :: i,j,k,l,ii,jj
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,arg2,tmp,dtmp
   real(wp) :: lambda,mmom(3),dmom(3)
   real(wp) :: es,expterm,erfterm
   real(wp) :: htmp,rxr(3,3),r3(3,3,3)
   real(wp) :: rcovij,rr

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
   real(wp),allocatable :: alpha(:)
   real(wp),allocatable :: beta(:)
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
   m    = n+3*n+1
   q    = 0.0_wp
   dipm = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(n), alpha(n), beta(n), source = 0.0_wp )

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
!$omp shared(n,at,cn,chrgeq) &
!$omp private(i,ii,tmp) &
!$omp shared(Xvec,Xfac,alpha,beta)
!$omp do schedule(dynamic)
   do i = 1, n
      ii = 4*(i-1)+1
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(ii) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
      alpha(i) = chrgeq%alpha(i)**2
      beta(i) = chrgeq%alpha(i)**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(n,m,at,xyz,chrgeq,alpha,beta) &
!$omp private(i,j,ii,jj,rij,r2,r,rxr,tmp,dtmp) &
!$omp private(arg,arg2,gamij,gamij2,expterm,erfterm) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, n
      ii = 4*(i-1)+1
      ! EN of atom i
      do j = 1, i-1
         jj = 4*(j-1)+1
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         r = sqrt(r2)
         ! charge-charge
         ! ⊕ ~~~ ⊕ -> repulsive
         ! ⊖ ~~~ ⊖ -> repulsive
         ! ⊕ ~~~ ⊖ -> attractive
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         tmp = erf(gamij*r)/r
         Amat(jj,ii) = tmp
         Amat(ii,jj) = Amat(jj,ii)
         ! charge-dipole
         ! ⊖⊕ ~~~ ⊕ -> repulsive
         ! ⊕⊖ ~~~ ⊕ -> attractive
         ! ⊖⊕ ~~~ ⊖ -> attractive
         ! ⊕⊖ ~~~ ⊖ -> repulsive
         gamij = 1.0_wp/sqrt(beta(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         tmp = erf(arg)/r
         dtmp = 2.0_wp*gamij*exp(-arg2)/(sqrtpi*r2)-tmp/r2
         Amat(jj,ii+1) = dtmp*rij(1)
         Amat(jj,ii+2) = dtmp*rij(2)
         Amat(jj,ii+3) = dtmp*rij(3)
         Amat(ii+1,jj) = Amat(jj,ii+1)
         Amat(ii+2,jj) = Amat(jj,ii+2)
         Amat(ii+3,jj) = Amat(jj,ii+3)
         gamij = 1.0_wp/sqrt(alpha(i)+beta(j))
         tmp = erf(gamij*r)/r
         dtmp = 2.0_wp*gamij*exp(-arg2)/(sqrtpi*r2)-tmp/r2
         Amat(ii,jj+1) = -dtmp*rij(1)
         Amat(ii,jj+2) = -dtmp*rij(2)
         Amat(ii,jj+3) = -dtmp*rij(3)
         Amat(jj+1,ii) = Amat(ii,jj+1)
         Amat(jj+2,ii) = Amat(ii,jj+2)
         Amat(jj+3,ii) = Amat(ii,jj+3)
         ! dipole-dipole
         gamij = 1.0_wp/sqrt(beta(i)+beta(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = erf(arg)/r
         expterm = 2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = dAdxdx(rij(1),r2,gamij2,erfterm,expterm)
         rxr(2,2) = dAdxdx(rij(2),r2,gamij2,erfterm,expterm)
         rxr(3,3) = dAdxdx(rij(3),r2,gamij2,erfterm,expterm)
         rxr(2,1) = dAdxdy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,1) = dAdxdy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,2) = dAdxdy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)
         Amat(jj+1,ii+1) = -rxr(1,1)
         Amat(jj+2,ii+1) = -rxr(2,1)
         Amat(jj+3,ii+1) = -rxr(3,1)
         Amat(jj+1,ii+2) = -rxr(1,2)
         Amat(jj+2,ii+2) = -rxr(2,2)
         Amat(jj+3,ii+2) = -rxr(3,2)
         Amat(jj+1,ii+3) = -rxr(1,3)
         Amat(jj+2,ii+3) = -rxr(2,3)
         Amat(jj+3,ii+3) = -rxr(3,3)
         Amat(ii+1,jj+1) = Amat(jj+1,ii+1)
         Amat(ii+2,jj+1) = Amat(jj+2,ii+1)
         Amat(ii+3,jj+1) = Amat(jj+3,ii+1)
         Amat(ii+1,jj+2) = Amat(jj+1,ii+2)
         Amat(ii+2,jj+2) = Amat(jj+2,ii+2)
         Amat(ii+3,jj+2) = Amat(jj+3,ii+2)
         Amat(ii+1,jj+3) = Amat(jj+1,ii+3)
         Amat(ii+2,jj+3) = Amat(jj+2,ii+3)
         Amat(ii+3,jj+3) = Amat(jj+3,ii+3)
      enddo
      ! isotropic exchange correlation
      Amat(ii,ii) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha(i))
      ! anisotropic exchange correlation
      Amat(ii+1,ii+1) = chrgeq%dpol(i)
      Amat(ii+2,ii+2) = chrgeq%dpol(i)
      Amat(ii+3,ii+3) = chrgeq%dpol(i)
      ! total charge constraint
      Amat(ii,m) = 1.0_wp
      Amat(m,ii) = 1.0_wp
      ! unconstrainted dipole moments
      Amat(ii+1,m) = 0.0_wp
      Amat(ii+2,m) = 0.0_wp
      Amat(ii+3,m) = 0.0_wp
      Amat(m,ii+1) = 0.0_wp
      Amat(m,ii+2) = 0.0_wp
      Amat(m,ii+3) = 0.0_wp
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges and moments" 
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

   do i = 1, n
      ii = 4*(i-1)+1
      q(i) = Xtmp(ii)
      dipm(1,i) = Xtmp(ii+1)
      dipm(2,i) = Xtmp(ii+2)
      dipm(3,i) = Xtmp(ii+3)
   enddo
   !print'(3f20.14)',Xtmp

   !call prmat(istdout,Amat,m,m,'A matrix')
   !call prmat(istdout,Xtmp,4,n,'solutions')

   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
     write(istdout, '(1x,"    #       q        mux      muy      muz")')
      do i = 1, n
         write(istdout,'(i6,a3,f9.3,3f9.3)') i, asym(at(i)), &
            q(i), dipm(1,i), dipm(2,i), dipm(3,i)
      enddo
      mmom = 0.0_wp
      dmom = 0.0_wp
      do i = 1, n
         mmom = mmom + q(i)*xyz(:,i)
         dmom = dmom + dipm(:,i)
      enddo
      write(istdout,'(a)')'molecular dipole:'
      write(istdout,'(17x,a)')'x           y           z       tot (Debye)'
      write(istdout,'(a,3f12.3)') ' q only: ',mmom
      write(istdout,'(a,4f12.3)') '   full: ',dmom+mmom,norm2(dmom+mmom)*2.5418_wp
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif
   if(abs(sum(q)-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error',1)
  
!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Anisotropic electrostatic (AES) energy calculation"
   work(:m) = Xvec
   call dsymv('u',m,0.5_wp,Amat,m,Xtmp,1,-1.0_wp,work,1)
   es = dot_product(Xtmp,work(:m))
   energy = es + energy
   if (lverbose) then
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "energy",es
   endif

!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  calculate molecular gradient of the AES energy
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  This part usually hits people hard, so watch out carefully.
!  So here is what's gonna to happen:
!  dE/dR = q·∂A/∂R·q - q·∂X/∂R + q·∂B/∂R·µ + q·∂B/∂R·µ + µ·∂C/∂R·µ
!  Note that ∂A/∂R is a three dimensional tensor, which I have to
!  contract twice with the partial charge array, while ∂X/∂R is
!  two dimensional tensor which has to be contracted once with the
!  charges. ∂B/∂R is three dimensional and will be contracted once
!  with the partial charges and once with the atomic dipole moments.
!  Analogously ∂C/∂R will be contracted twice with the dipole moments.
!  Keeping this in mind I find that there are two kind of terms,
!  if I shuffle them around and group them together I can express them
!  as linear maps from charge to cartesian space and from dipole to
!  cartesian space:
!  dE/dR = q·(∂A/∂R·q - ∂X/∂R + ∂B/∂R·µ) + (q·∂B/∂R + µ·∂C/∂R)·µ
!  So the stategy here is to not calculate the three dimensional tensors
!  but to contract them directly with the partial charges or dipole moments
!  from one side. In contrast to that I will not contract the two
!  dimensional tensor, just to keep things interesting.
!  There is one thing to note here:  B = ∂A/∂R and C = ∂B/∂R, which
!  is something I will use extensively below. So in principle I only
!  need to calculate the interaction matrix D = ∂³A/∂³R, which is the
!  third derivative of the A matrix. So I will switch the complicated
!  first derivative of the energy for the less complicated third
!  derivative of the interaction martix of the monopole problem.
!  Having said that, you will notice that I lied to you about not
!  calculating the three dimensional tensors, already did it.
!
!  If you believe me, than the following implementation is a juggling trick
!  for you. Otherwise you will experience some magic between the lines.
!! ------------------------------------------------------------------------
!  Here is my nomenclature by the way, to keep things interesting, I use
!  some relative memory access rules and interlace partial charges and
!  atomic dipole moments. So here you go:
!  dE/dRj -> g(:,j), ∂Xi/∂Rj -> -dcn(:,i,j), ∂(Aij·qj)/∂Rj -> dAmat(:,j,ii)
!! ------------------------------------------------------------------------
do_molecular_gradient: if (lgrad .or. lcpq) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "molecular gradient calculation"
   allocate( dAmat(3,n,m), dXvec(3,n,m), Afac(3,m), source = 0.0_wp )
   !allocate( dAmat(3,n,m), Afac(3,n), source = 0.0_wp )
!$omp parallel default(none) &
!$omp shared(n,dcndr,xyz,alpha,Amat,Xfac,Xtmp,beta) &
!$omp private(i,j,ii,jj,rij,r,r2,rxr,r3,gamij,gamij2,arg,arg2,tmp,dtmp) &
!$omp private(erfterm,expterm) &
!$omp shared(dXvec) &
!$omp reduction(+:Afac,dAmat)
!$omp do schedule(dynamic)
   do i = 1, n
      ii = 4*(i-1)+1
      dXvec(:,ii,ii) = +dcndr(:,i,i)*Xfac(i) ! merge dX and dA for speedup
      do j = 1, i-1
         jj = 4*(j-1)+1
         ! ------------------------- !
         ! one side charge : dX/dR   !
         ! ------------------------- !
         dXvec(:,j,ii) = dcndr(:,i,j)*Xfac(i)
         dXvec(:,i,jj) = dcndr(:,j,i)*Xfac(j)
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         r = sqrt(r2)

         ! ------------------------- !
         ! charge-charge :  d(Aq)/dR !
         ! ------------------------- !
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         tmp = erf(arg)/r
         dtmp = 2.0_wp*gamij*exp(-arg2)/(sqrtpi*r2)-tmp/r2
         !dAmat(:,j,jj) = -dtmp*rij*Xtmp(ii) + dAmat(:,j,jj)
         !dAmat(:,i,ii) = +dtmp*rij*Xtmp(jj) + dAmat(:,i,ii)
         dAmat(:,i,jj) = +dtmp*rij*Xtmp(ii) + dAmat(:,i,jj)
         dAmat(:,j,ii) = -dtmp*rij*Xtmp(jj) + dAmat(:,j,ii)

         ! ------------------------- !
         ! charge-dipole : d(Bµ)/dR  !
         ! ------------------------- !
         gamij = 1.0_wp/sqrt(alpha(i)+beta(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = -erf(arg)/r
         expterm = -2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = dAdxdx(rij(1),r2,gamij2,erfterm,expterm)
         rxr(2,2) = dAdxdx(rij(2),r2,gamij2,erfterm,expterm)
         rxr(3,3) = dAdxdx(rij(3),r2,gamij2,erfterm,expterm)
         rxr(2,1) = dAdxdy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,1) = dAdxdy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,2) = dAdxdy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         !dAmat(1,j,jj) = +sum(rxr(:,1)*Xtmp(ii+1:ii+3)) + dAmat(1,j,jj)
         !dAmat(2,j,jj) = +sum(rxr(:,2)*Xtmp(ii+1:ii+3)) + dAmat(2,j,jj)
         !dAmat(3,j,jj) = +sum(rxr(:,3)*Xtmp(ii+1:ii+3)) + dAmat(3,j,jj)

         dAmat(1,i,jj) = -sum(rxr(:,1)*Xtmp(ii+1:ii+3)) + dAmat(1,i,jj)
         dAmat(2,i,jj) = -sum(rxr(:,2)*Xtmp(ii+1:ii+3)) + dAmat(2,i,jj)
         dAmat(3,i,jj) = -sum(rxr(:,3)*Xtmp(ii+1:ii+3)) + dAmat(3,i,jj)

         gamij = 1.0_wp/sqrt(beta(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = -erf(arg)/r
         expterm = -2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = dAdxdx(rij(1),r2,gamij2,erfterm,expterm)
         rxr(2,2) = dAdxdx(rij(2),r2,gamij2,erfterm,expterm)
         rxr(3,3) = dAdxdx(rij(3),r2,gamij2,erfterm,expterm)
         rxr(2,1) = dAdxdy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,1) = dAdxdy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,2) = dAdxdy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         !dAmat(1,i,ii) = +sum(rxr(:,1)*Xtmp(jj+1:jj+3)) + dAmat(1,i,ii)
         !dAmat(2,i,ii) = +sum(rxr(:,2)*Xtmp(jj+1:jj+3)) + dAmat(2,i,ii)
         !dAmat(3,i,ii) = +sum(rxr(:,3)*Xtmp(jj+1:jj+3)) + dAmat(3,i,ii)

         dAmat(1,j,ii) = -sum(rxr(:,1)*Xtmp(jj+1:jj+3)) + dAmat(1,j,ii)
         dAmat(2,j,ii) = -sum(rxr(:,2)*Xtmp(jj+1:jj+3)) + dAmat(2,j,ii)
         dAmat(3,j,ii) = -sum(rxr(:,3)*Xtmp(jj+1:jj+3)) + dAmat(3,j,ii)

         ! ------------------------- !
         ! dipole-charge : d(Bq)/dR  !
         ! ------------------------- !
         gamij = 1.0_wp/sqrt(beta(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = erf(arg)/r
         expterm = 2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = dAdxdx(rij(1),r2,gamij2,erfterm,expterm)
         rxr(2,2) = dAdxdx(rij(2),r2,gamij2,erfterm,expterm)
         rxr(3,3) = dAdxdx(rij(3),r2,gamij2,erfterm,expterm)
         rxr(2,1) = dAdxdy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,1) = dAdxdy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,2) = dAdxdy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         !dAmat(:,j,jj+1) = +rxr(:,1)*Xtmp(ii) + dAmat(:,j,jj+1)
         !dAmat(:,j,jj+2) = +rxr(:,2)*Xtmp(ii) + dAmat(:,j,jj+2)
         !dAmat(:,j,jj+3) = +rxr(:,3)*Xtmp(ii) + dAmat(:,j,jj+3)

         dAmat(:,i,jj+1) = -rxr(:,1)*Xtmp(ii) + dAmat(:,i,jj+1)
         dAmat(:,i,jj+2) = -rxr(:,2)*Xtmp(ii) + dAmat(:,i,jj+2)
         dAmat(:,i,jj+3) = -rxr(:,3)*Xtmp(ii) + dAmat(:,i,jj+3)

         gamij = 1.0_wp/sqrt(alpha(i)+beta(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = erf(arg)/r
         expterm = 2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = dAdxdx(rij(1),r2,gamij2,erfterm,expterm)
         rxr(2,2) = dAdxdx(rij(2),r2,gamij2,erfterm,expterm)
         rxr(3,3) = dAdxdx(rij(3),r2,gamij2,erfterm,expterm)
         rxr(2,1) = dAdxdy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,1) = dAdxdy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         rxr(3,2) = dAdxdy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)

         !dAmat(:,i,ii+1) = +rxr(:,1)*Xtmp(jj) + dAmat(:,i,ii+1)
         !dAmat(:,i,ii+2) = +rxr(:,2)*Xtmp(jj) + dAmat(:,i,ii+2)
         !dAmat(:,i,ii+3) = +rxr(:,3)*Xtmp(jj) + dAmat(:,i,ii+3)

         dAmat(:,j,ii+1) = -rxr(:,1)*Xtmp(jj) + dAmat(:,j,ii+1)
         dAmat(:,j,ii+2) = -rxr(:,2)*Xtmp(jj) + dAmat(:,j,ii+2)
         dAmat(:,j,ii+3) = -rxr(:,3)*Xtmp(jj) + dAmat(:,j,ii+3)

         ! ------------------------- !
         ! dipole-dipole : d(Cµ)/dR  !
         ! ------------------------- !
         gamij = 1.0_wp/sqrt(beta(i)+beta(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = -erf(arg)/r
         expterm = -2*gamij*exp(-arg2)/sqrtpi
         r3(1,1,1) = dAdxdxdx(rij(1),r2,gamij2,erfterm,expterm)
         r3(2,1,1) = dAdxdydy(rij(2),rij(1),r2,gamij2,erfterm,expterm)
         r3(3,1,1) = dAdxdydy(rij(3),rij(1),r2,gamij2,erfterm,expterm)
         r3(1,2,1) = r3(2,1,1)
         r3(2,2,1) = dAdxdydy(rij(1),rij(2),r2,gamij2,erfterm,expterm)
         r3(3,2,1) = -dAdxdydz(rij,r2,gamij2,erfterm,expterm)
         r3(1,3,1) = r3(3,1,1)
         r3(2,3,1) = r3(3,2,1)
         r3(3,3,1) = dAdxdydy(rij(1),rij(3),r2,gamij2,erfterm,expterm)
         r3(1,1,2) = r3(2,1,1)
         r3(2,1,2) = r3(2,2,1)
         r3(3,1,2) = r3(3,2,1)
         r3(1,2,2) = r3(2,2,1)
         r3(2,2,2) = dAdxdxdx(rij(2),r2,gamij2,erfterm,expterm)
         r3(3,2,2) = dAdxdydy(rij(3),rij(2),r2,gamij2,erfterm,expterm)
         r3(1,3,2) = r3(3,2,1)
         r3(2,3,2) = r3(3,2,2)
         r3(3,3,2) = dAdxdydy(rij(2),rij(3),r2,gamij2,erfterm,expterm)
         r3(1,1,3) = r3(3,1,1)
         r3(2,1,3) = r3(3,2,1)
         r3(3,1,3) = r3(3,3,1)
         r3(1,2,3) = r3(3,2,1)
         r3(2,2,3) = r3(3,2,2)
         r3(3,2,3) = r3(3,3,2)
         r3(1,3,3) = r3(3,3,1)
         r3(2,3,3) = r3(3,3,2)
         r3(3,3,3) = dAdxdxdx(rij(3),r2,gamij2,erfterm,expterm)
         dAmat(:,i,jj+1) = + r3(:,1,1)*Xtmp(ii+1) + r3(:,2,1)*Xtmp(ii+2) &
                           + r3(:,3,1)*Xtmp(ii+3) + dAmat(:,i,jj+1)
         dAmat(:,i,jj+2) = + r3(:,1,2)*Xtmp(ii+1) + r3(:,2,2)*Xtmp(ii+2) &
                           + r3(:,3,2)*Xtmp(ii+3) + dAmat(:,i,jj+2)
         dAmat(:,i,jj+3) = + r3(:,1,3)*Xtmp(ii+1) + r3(:,2,3)*Xtmp(ii+2) &
                           + r3(:,3,3)*Xtmp(ii+3) + dAmat(:,i,jj+3)
         dAmat(:,j,ii+1) = - r3(:,1,1)*Xtmp(jj+1) - r3(:,2,1)*Xtmp(jj+2) &
                           - r3(:,3,1)*Xtmp(jj+3) + dAmat(:,j,ii+1)
         dAmat(:,j,ii+2) = - r3(:,1,2)*Xtmp(jj+1) - r3(:,2,2)*Xtmp(jj+2) &
                           - r3(:,3,2)*Xtmp(jj+3) + dAmat(:,j,ii+2)
         dAmat(:,j,ii+3) = - r3(:,1,3)*Xtmp(jj+1) - r3(:,2,3)*Xtmp(jj+2) &
                           - r3(:,3,3)*Xtmp(jj+3) + dAmat(:,j,ii+3)
      enddo
   enddo
!$omp enddo
!$omp endparallel
endif do_molecular_gradient

   if (lgrad) then
   call dgemv('n',3*n,m,+1.0_wp,dAmat,3*n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*n,m,-1.0_wp,dXvec,3*n,Xtmp,1,1.0_wp,gradient,1)
   endif

!! ------------------------------------------------------------------------
!  You made it through, want more? Let's solve some coupled perturbed
!  equations and invert the supermatrix!
!! ------------------------------------------------------------------------
!  invert the supermatrix by using a Bunch-Kaufman factorization
!  A⁻¹ = (L·D·L^T)⁻¹ = L^T·D⁻¹·L
!! ------------------------------------------------------------------------
do_coupled_perturbed_charges: if (lcpq) then
   if (lqonly) then
      dqdr(:,:,:n+1) = 0.0_wp
   else
      dqdr = 0.0_wp
   endif
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
!  This is the most expensive part of this calculation, so we don't
!  want to waste any time with unneeded computations.
!  An easy way to reduce the effort is to only solve for the partial
!  charge derivative instead for the full system. This cuts down the
!  computational cost by a factor of 64.
!! ------------------------------------------------------------------------
   if (lqonly) then
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges"
   else
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "calculating the derivative of the partial charges and atomic dipole moments"
   endif
endif do_coupled_perturbed_charges

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

contains

pure function dAdxdx(x,r2,gam2,erfterm,expterm) result(dA)
   implicit none
   real(wp),intent(in) :: x
   real(wp),intent(in) :: r2
   real(wp),intent(in) :: gam2
   real(wp),intent(in) :: erfterm
   real(wp),intent(in) :: expterm
   real(wp) :: dA

   dA = erfterm * ( 3*x**2/r2**2 - 1.0_wp/r2 ) &
      - expterm * ( 3*x**2/r2**2 + 2*gam2*x**2/r2 - 1.0_wp/r2 )

end function dAdxdx

pure function dAdxdy(x,y,r2,gam2,erfterm,expterm) result(dA)
   implicit none
   real(wp),intent(in) :: x
   real(wp),intent(in) :: y
   real(wp),intent(in) :: r2
   real(wp),intent(in) :: gam2
   real(wp),intent(in) :: erfterm
   real(wp),intent(in) :: expterm
   real(wp) :: dA

   dA = erfterm * 3*x*y/r2**2 &
      - expterm * ( 3*x*y/r2**2 + 2*gam2*x*y/r2 )

end function dAdxdy

pure function dAdxdydz(r,r2,gam2,erfterm,expterm) result(dA)
   implicit none
   real(wp),intent(in) :: r(3)
   real(wp),intent(in) :: r2
   real(wp),intent(in) :: gam2
   real(wp),intent(in) :: erfterm
   real(wp),intent(in) :: expterm
   real(wp) :: dA
   real(wp) :: rxyz

   rxyz = product(r)

   dA = - erfterm * ( 15*rxyz/r2**3 ) &
        + expterm * ( 3*rxyz/r2 + 2*gam2*rxyz*(2*gam2*r2+3.0) &
                      - 4*gam2*rxyz )/r2**2

end function dAdxdydz

pure function dAdxdydy(x,y,r2,gam2,erfterm,expterm) result(dA)
   implicit none
   real(wp),intent(in) :: x
   real(wp),intent(in) :: y
   real(wp),intent(in) :: r2
   real(wp),intent(in) :: gam2
   real(wp),intent(in) :: erfterm
   real(wp),intent(in) :: expterm
   real(wp) :: dA

   dA = - erfterm * ( ( -2*x*r2 + 5*x*(r2-3*y**2) )/r2**3 ) &
        + expterm * ( x*(x**2*(1-2*gam2*y**2)) - x*(r2-3*y**2)/r2 &
                      - (4/r2+2*gam2)*x*(r2 - 3*y**2 - 2*gam2*y**2*r2)/r2 &
                      )/r2**2 

end function dAdxdydy

pure function dAdxdxdx(x,r2,gam2,erfterm,expterm) result(dA)
   implicit none
   real(wp),intent(in) :: x
   real(wp),intent(in) :: r2
   real(wp),intent(in) :: gam2
   real(wp),intent(in) :: erfterm
   real(wp),intent(in) :: expterm
   real(wp) :: dA

   dA = - erfterm * ( ( -4*x*r2 + 5*x*(r2-3*x**2) )/r2**3 ) &
        + expterm * ( -4*x*(gam2*x**2 + (gam2*(r2-x**2)+1)) - x*(r2-3*x**2) &
                      - (4/r2+2*gam2)*x*(r2 - 3*x**2 - 2*gam2*x**2*r2)/r2 &
                      )/r2**2

end function dAdxdxdx

end subroutine goedecker_multieq2


subroutine goedecker_multieq(n,at,xyz,chrg,cn,dcndr,q,dipm,energy,&
                             lverbose)
   use iso_fortran_env, wp => real64, istdout => output_unit
   use ncoord, only : rcov, kn
   implicit none

!! ------------------------------------------------------------------------
!  Input
!! ------------------------------------------------------------------------
   integer, intent(in)    :: n                ! number of atoms
   integer, intent(in)    :: at(n)            ! ordinal numbers
   real(wp),intent(in)    :: xyz(3,n)         ! geometry
   real(wp),intent(in)    :: chrg             ! total charge
   real(wp),intent(in)    :: cn(n)            ! erf-CN
   real(wp),intent(in)    :: dcndr(3,n,n)     ! derivative of erf-CN
!   real(wp),intent(in)    :: par_screen       ! dielectric screening
!   real(wp),intent(in)    :: par_gscale       ! scaling for hardnesses
!   real(wp),intent(in)    :: par_rscale       ! scaling for radii
   logical, intent(in)    :: lverbose         ! toggles printout
!! ------------------------------------------------------------------------
!  Output
!! ------------------------------------------------------------------------
   real(wp),intent(out)   :: q(n)             ! partial charges
   real(wp),intent(out)   :: dipm(3,n)        ! atomic dipole moments
   real(wp),intent(inout) :: energy           ! electrostatic energy

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
   character(len=2),external :: asym
   integer  :: i,j,k,l,ii,jj
   real(wp) :: r,rij(3),r2
   real(wp) :: gamij,gamij2
   real(wp) :: arg,arg2,tmp,dtmp
   real(wp) :: lambda,mmom(3),dmom(3)
   real(wp) :: es,expterm,erfterm
   real(wp) :: htmp,rxr(3,3)
   real(wp) :: rcovij,rr

!! ------------------------------------------------------------------------
!  scratch variables
!! ------------------------------------------------------------------------
   real(wp),allocatable :: alpha(:)
   real(wp),allocatable :: beta(:)
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
   m    = n+3*n+1
   q    = 0.0_wp
   dipm = 0.0_wp
   allocate( ipiv(m), source = 0 )
   allocate( Amat(m,m), Xvec(m), Xfac(n), alpha(n), beta(n), source = 0.0_wp )

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
!$omp private(i,ii,tmp) &
!$omp shared(Xvec,Xfac,alpha,beta)
!$omp do schedule(dynamic)
   do i = 1, n
      ii = 4*(i-1)+1
      tmp = cnfak(at(i))/(sqrt(cn(i))+1e-14_wp)
      Xvec(ii) = -xi(at(i)) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
      alpha(i) = alp(at(i))**2
      beta(i) = 2.5_wp*alp(at(i))**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(n,m,at,xyz,gamm,alpha,beta) &
!$omp private(i,j,ii,jj,rij,r2,r,rxr,tmp,dtmp) &
!$omp private(arg,arg2,gamij,gamij2,expterm,erfterm) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, n
      ii = 4*(i-1)+1
      ! EN of atom i
      do j = 1, i-1
         jj = 4*(j-1)+1
         rij = xyz(:,i) - xyz(:,j)
         r2 = sum(rij**2)
         r = sqrt(r2)
         ! charge-charge
         ! ⊕ ~~~ ⊕ -> repulsive
         ! ⊖ ~~~ ⊖ -> repulsive
         ! ⊕ ~~~ ⊖ -> attractive
         gamij = 1.0_wp/sqrt(alpha(i)+alpha(j))
         tmp = erf(gamij*r)/r
         Amat(jj,ii) = tmp
         Amat(ii,jj) = Amat(jj,ii)
         ! charge-dipole
         ! ⊖⊕ ~~~ ⊕ -> repulsive
         ! ⊕⊖ ~~~ ⊕ -> attractive
         ! ⊖⊕ ~~~ ⊖ -> attractive
         ! ⊕⊖ ~~~ ⊖ -> repulsive
         gamij = 1.0_wp/sqrt(beta(i)+alpha(j))
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         tmp = erf(arg)/r
         dtmp = 2.0_wp*gamij*exp(-arg2)/(sqrtpi*r2)-tmp/r2
         Amat(jj,ii+1) = dtmp*rij(1)
         Amat(jj,ii+2) = dtmp*rij(2)
         Amat(jj,ii+3) = dtmp*rij(3)
         Amat(ii+1,jj) = Amat(jj,ii+1)
         Amat(ii+2,jj) = Amat(jj,ii+2)
         Amat(ii+3,jj) = Amat(jj,ii+3)
         gamij = 1.0_wp/sqrt(alpha(i)+beta(j))
         tmp = erf(gamij*r)/r
         dtmp = 2.0_wp*gamij*exp(-arg2)/(sqrtpi*r2)-tmp/r2
         Amat(ii,jj+1) = -dtmp*rij(1)
         Amat(ii,jj+2) = -dtmp*rij(2)
         Amat(ii,jj+3) = -dtmp*rij(3)
         Amat(jj+1,ii) = Amat(ii,jj+1)
         Amat(jj+2,ii) = Amat(ii,jj+2)
         Amat(jj+3,ii) = Amat(ii,jj+3)
         ! dipole-dipole
         gamij = 1.0_wp/sqrt(beta(i)+beta(j))/4.0_wp
         gamij2 = gamij**2
         arg2 = gamij2 * r2
         arg = sqrt(arg2)
         erfterm = erf(arg)/r
         expterm = 2*gamij*exp(-arg2)/sqrtpi
         rxr(1,1) = erfterm * ( 3*rij(1)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(1)**2/r2**2 + 2*gamij2*rij(1)**2/r2 - 1/r2 )
         rxr(2,2) = erfterm * ( 3*rij(2)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(2)**2/r2**2 + 2*gamij2*rij(2)**2/r2 - 1/r2 )
         rxr(3,3) = erfterm * ( 3*rij(3)**2/r2**2 - 1.0_wp/r2 ) &
                  - expterm * ( 3*rij(3)**2/r2**2 + 2*gamij2*rij(3)**2/r2 - 1/r2 )
         rxr(2,1) = erfterm * 3*rij(2)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(2)*rij(1)/r2**2 + 2*gamij2*rij(2)*rij(1)/r2 )
         rxr(3,1) = erfterm * 3*rij(3)*rij(1)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(1)/r2**2 + 2*gamij2*rij(3)*rij(1)/r2 )
         rxr(3,2) = erfterm * 3*rij(3)*rij(2)/r2**2 &
                  - expterm * ( 3*rij(3)*rij(2)/r2**2 + 2*gamij2*rij(3)*rij(2)/r2 )
         ! symmetrize
         rxr(1,2) = rxr(2,1)
         rxr(1,3) = rxr(3,1)
         rxr(2,3) = rxr(3,2)
         Amat(jj+1,ii+1) = -rxr(1,1)
         Amat(jj+2,ii+1) = -rxr(2,1)
         Amat(jj+3,ii+1) = -rxr(3,1)
         Amat(jj+1,ii+2) = -rxr(1,2)
         Amat(jj+2,ii+2) = -rxr(2,2)
         Amat(jj+3,ii+2) = -rxr(3,2)
         Amat(jj+1,ii+3) = -rxr(1,3)
         Amat(jj+2,ii+3) = -rxr(2,3)
         Amat(jj+3,ii+3) = -rxr(3,3)
         Amat(ii+1,jj+1) = Amat(jj+1,ii+1)
         Amat(ii+2,jj+1) = Amat(jj+2,ii+1)
         Amat(ii+3,jj+1) = Amat(jj+3,ii+1)
         Amat(ii+1,jj+2) = Amat(jj+1,ii+2)
         Amat(ii+2,jj+2) = Amat(jj+2,ii+2)
         Amat(ii+3,jj+2) = Amat(jj+3,ii+2)
         Amat(ii+1,jj+3) = Amat(jj+1,ii+3)
         Amat(ii+2,jj+3) = Amat(jj+2,ii+3)
         Amat(ii+3,jj+3) = Amat(jj+3,ii+3)
      enddo
      ! isotropic exchange correlation
      Amat(ii,ii) = gamm(at(i)) + sqrt2pi/sqrt(alpha(i))
      ! anisotropic exchange correlation
      Amat(ii+1,ii+1) = sqrt2pi/sqrt(alpha(i)**4/beta(i)) ! FIXME
      Amat(ii+2,ii+2) = sqrt2pi/sqrt(alpha(i)**4/beta(i)) ! FIXME
      Amat(ii+3,ii+3) = sqrt2pi/sqrt(alpha(i)**4/beta(i)) ! FIXME
      ! total charge constraint
      Amat(ii,m) = 1.0_wp
      Amat(m,ii) = 1.0_wp
      ! unconstrainted dipole moments
      Amat(ii+1,m) = 0.0_wp
      Amat(ii+2,m) = 0.0_wp
      Amat(ii+3,m) = 0.0_wp
      Amat(m,ii+1) = 0.0_wp
      Amat(m,ii+2) = 0.0_wp
      Amat(m,ii+3) = 0.0_wp
   enddo
!$omp enddo
!$omp endparallel

!! ------------------------------------------------------------------------
!  solve the linear equations to obtain partial charges
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Solve the linear equations to obtain partial charges and moments" 
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

   do i = 1, n
      ii = 4*(i-1)+1
      q(i) = Xtmp(ii)
      dipm(1,i) = Xtmp(ii+1)
      dipm(2,i) = Xtmp(ii+2)
      dipm(3,i) = Xtmp(ii+3)
   enddo
   !print'(3f20.14)',Xtmp

   !call prmat(istdout,Amat,m,m,'A matrix')
   !call prmat(istdout,Xtmp,4,n,'solutions')

   lambda = Xtmp(m)
   if (lverbose) then
      write(istdout,'(72("-"))')
     write(istdout, '(1x,"    #       q        mux      muy      muz")')
      do i = 1, n
         write(istdout,'(i6,a3,f9.3,3f9.3)') i, asym(at(i)), &
            q(i), dipm(1,i), dipm(2,i), dipm(3,i)
      enddo
      mmom = 0.0_wp
      dmom = 0.0_wp
      do i = 1, n
         mmom = mmom + q(i)*xyz(:,i)
         dmom = dmom + dipm(:,i)
      enddo
      write(istdout,'(a)')'molecular dipole:'
      write(istdout,'(17x,a)')'x           y           z       tot (Debye)'
      write(istdout,'(a,3f12.3)') ' q only: ',mmom
      write(istdout,'(a,4f12.3)') '   full: ',dmom+mmom,norm2(dmom+mmom)*2.5418_wp
      write(istdout,'(72("-"))')
      write(istdout,'(1x,a,1x,f20.14)') &
         "lambda       :",lambda,&
         "total charge :",sum(q)
   endif
   if(abs(sum(q)-chrg) > 1.e-6_wp) &
      call raise('E','(goedecker_solve) charge constrain error',1)
  
!! ------------------------------------------------------------------------
!  calculate isotropic electrostatic (IES) energy
!! ------------------------------------------------------------------------
!  E = ∑i (ENi - κi·√CNi)·qi + ∑i (Jii + 2/√π·γii)·q²i
!      + ½ ∑i ∑j,j≠i qi·qj·2/√π·F0(γ²ij·R²ij)
!    = q·(½A·q - X)
!! ------------------------------------------------------------------------
   if (lverbose) write(istdout,'(72("="),/,1x,a)') &
      "Anisotropic electrostatic (AES) energy calculation"
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

end subroutine goedecker_multieq

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
   real(wp),intent(out)   :: dqdr(3,mol%n,mol%n+1) ! derivative of partial charges
   real(wp),intent(out)   :: dqdL(3,3,mol%n+1) ! derivative of partial charges
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,mol%n) ! molecular gradient of IES
   real(wp),intent(inout) :: sigma(3,3) ! molecular gradient of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
   real(wp),allocatable :: alpha2(:)
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
   allocate( Amat(m,m), Xvec(m), Xfac(m), alpha2(mol%n), &
      &      source = 0.0_wp )

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

!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(mol,cn,chrgeq) &
!$omp private(i,tmp) &
!$omp shared(Xvec,Xfac,alpha2)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
      alpha2(i) = chrgeq%alpha(i)**2
   enddo
!$omp enddo
!$omp endparallel

   if (mol%npbc > 0) then
   cf = sqrtpi/mol%volume**(1.0_wp/3.0_wp)
   ! build Ewald matrix
!$omp parallel default(none) &
!$omp private(i,j,wscAt,riw,gamii,gamij) &
!$omp shared(mol,chrgeq,alpha2,cf) &
!$omp reduction(+:Amat)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      gamii = 1.0/sqrt(alpha2(i)+alpha2(i))
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha2(i)) &
         ! reciprocal for 0th atom
         + eeq_ewald_3d_rec(zero,ewaldCutR,mol%rec_lat,mol%volume,cf) &
         ! direct for 0th atom
         + eeq_ewald_3d_dir(zero,ewaldCutD,mol%lattice,gamii,cf) &
         ! neutralizing background contribution
         - pi/mol%volume/cf**2
      do j = 1, i-1
         gamij = 1.0_wp/sqrt(alpha2(i)+alpha2(j))
         do wscAt = 1, mol%wsc%itbl(j,i)
            riw = mol%xyz(:,i) - mol%xyz(:,j) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,j,i))
            Amat(i,j) = Amat(i,j) + mol%wsc%w(j,i) * ( &
               ! reciprocal lattice sums
               + eeq_ewald_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf) &
               ! direct lattice sums
               + eeq_ewald_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf) &
               ! neutralizing background contribution
               - pi/mol%volume/cf**2 )
         end do
         Amat(j,i) = Amat(i,j)
      end do
   end do
!$omp enddo
!$omp endparallel
   else
!$omp parallel default(none) &
!$omp shared(mol,chrgeq,alpha2) &
!$omp private(i,j,r,gamij) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(dynamic)
   ! prepare A matrix
   do i = 1, mol%n
      do j = 1, i-1
         r = norm2(mol%xyz(:,j) - mol%xyz(:,i))
         gamij = 1.0_wp/sqrt(alpha2(i)+alpha2(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha2(i)) ! important
   enddo
!$omp enddo
!$omp endparallel
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
!$omp parallel default(none) &
!$omp shared(mol,dcndr,dcndL,alpha2,Xfac,Xtmp,cf) &
!$omp private(i,j,wscAt,riw,gamij,dAtmp,stmp) &
!$omp shared(dXvecdr,dXvecdL) &
!$omp reduction(+:Afac,dAmatdr,dAmatdL)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i)
      dXvecdL(:,:,i) = +dcndL(:,:,i)*Xfac(i)
      do j = 1, i-1
         if(mol%wsc%at(j,i).eq.0) cycle
         ! over WSC partner
         gamij = 1.0_wp/sqrt(alpha2(i) + alpha2(j))
         do wscAt = 1, mol%wsc%itbl(j,i)
            riw = mol%xyz(:,i) - mol%xyz(:,j) &
               &  - matmul(mol%lattice,mol%wsc%lattr(:,wscAt,j,i))
            call eeq_ewald_dx_3d_rec(riw,ewaldCutR,mol%rec_lat,mol%volume,cf, &
               &                     dAtmp,stmp)
            dAtmp = dAtmp*mol%wsc%w(j,i)
            stmp = stmp*mol%wsc%w(j,i)
            dAmatdr(:,i,j) = dAmatdr(:,i,j) + dAtmp(:)*Xtmp(i)
            dAmatdr(:,j,i) = dAmatdr(:,j,i) - dAtmp(:)*Xtmp(j)
            Afac(:,i)    =  Afac(:,i)   + dAtmp(:)*Xtmp(j)
            Afac(:,j)    =  Afac(:,j)   - dAtmp(:)*Xtmp(i)
            dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(j)
            dAmatdL(:,:,j) = dAmatdL(:,:,j) + stmp(:,:)*Xtmp(i)
            call eeq_ewald_dx_3d_dir(riw,ewaldCutD,mol%lattice,gamij,cf, &
               &                     dAtmp,stmp)
            dAtmp = dAtmp*mol%wsc%w(j,i)
            stmp = stmp*mol%wsc%w(j,i)
            dAmatdr(:,i,j) = dAmatdr(:,i,j) + dAtmp(:)*Xtmp(i)
            dAmatdr(:,j,i) = dAmatdr(:,j,i) - dAtmp(:)*Xtmp(j)
            Afac(:,i)    =  Afac(:,i)   + dAtmp(:)*Xtmp(j)
            Afac(:,j)    =  Afac(:,j)   - dAtmp(:)*Xtmp(i)
            dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(j)
            dAmatdL(:,:,j) = dAmatdL(:,:,j) + stmp(:,:)*Xtmp(i)
         enddo  ! k WSC partner
      enddo     ! j

      gamij = 1.0_wp/sqrt(alpha2(i) + alpha2(i))
      call eeq_ewald_dx_3d_rec(zero,ewaldCutR,mol%rec_lat,mol%volume,cf, &
         &                     dAtmp,stmp)
      !call coord_trafo(3,mol%lattice,stmp) ! WRONG
      dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(i) ! WRONG
      call eeq_ewald_dx_3d_dir(zero,ewaldCutD,mol%lattice,gamij,cf, &
         &                     dAtmp,stmp)
      dAmatdL(:,:,i) = dAmatdL(:,:,i) + stmp(:,:)*Xtmp(i)
      do j = 1,3
         dAmatdL(j,j,i) = dAmatdL(j,j,i) + cf/sqrtpi/3.0_wp*Xtmp(i)
      enddo
   enddo
!$omp enddo
!$omp endparallel
   else
!$omp parallel default(none) &
!$omp shared(mol,dcndr,alpha2,Amat,Xfac,Xtmp) &
!$omp private(i,j,rij,r2,gamij,arg,dtmp) &
!$omp shared(dXvecdr,dAmatdr) &
!$omp reduction(+:Afac)
!$omp do schedule(dynamic)
   do i = 1, mol%n
      dXvecdr(:,:,i) = +dcndr(:,:,i)*Xfac(i) ! merge dX and dA for speedup
      do j = 1, i-1
         rij = mol%xyz(:,i) - mol%xyz(:,j)
         r2 = sum(rij**2)
         gamij = 1.0_wp/sqrt(alpha2(i) + alpha2(j))
         arg = gamij**2*r2
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2)-Amat(j,i)/r2
         Afac(:,i) = +dtmp*rij*Xtmp(j) + Afac(:,i)
         Afac(:,j) = -dtmp*rij*Xtmp(i) + Afac(:,j)
         dAmatdr(:,i,j) = +dtmp*rij*Xtmp(i)
         dAmatdr(:,j,i) = -dtmp*rij*Xtmp(j)
      enddo
   enddo
!$omp enddo
!$omp endparallel
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
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dAmatdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dXvecdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   if (mol%npbc > 0) then
      call dgemm('n','n',3*3,m,m,-1.0_wp,dAmatdL,3*3,Ainv,m,1.0_wp,dqdL,3*3)
      call dgemm('n','n',3*3,m,m,+1.0_wp,dXvecdL,3*3,Ainv,m,1.0_wp,dqdL,3*3)
   endif
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

! ------------------------------------------------------------------------
!  Clean up
! ------------------------------------------------------------------------
   if (allocated(alpha2)) deallocate(alpha2)
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
   real(wp),intent(out)   :: dqdr(3,mol%n,mol%n+1) ! derivative of partial charges
   real(wp),intent(inout) :: gsolv        ! electrostatic energy
   real(wp),intent(inout) :: energy        ! electrostatic energy
   real(wp),intent(inout) :: gradient(3,mol%n) ! molecular gradient of IES

!  π itself
   real(wp),parameter :: pi = 3.1415926535897932384626433832795029_wp
!  √π
   real(wp),parameter :: sqrtpi = sqrt(pi)
!  √(2/π)
   real(wp),parameter :: sqrt2pi = sqrt(2.0_wp/pi)
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
   real(wp),allocatable :: alpha2(:)
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
   allocate( Amat(m,m), Xvec(m), Xfac(m), alpha2(mol%n), source = 0.0_wp )

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

!  prepare some arrays
!$omp parallel default(none) &
!$omp shared(mol,cn,chrgeq) &
!$omp private(i,tmp) &
!$omp shared(Xvec,Xfac,alpha2)
!$omp do schedule(runtime)
   do i = 1, mol%n
      tmp = chrgeq%kappa(i)/(sqrt(cn(i))+1e-14_wp)
      Xvec(i) = -chrgeq%en(i) + tmp*cn(i)
      Xfac(i) = 0.5_wp*tmp
      alpha2(i) = chrgeq%alpha(i)**2
   enddo
!$omp enddo
!$omp endparallel

!$omp parallel default(none) &
!$omp shared(mol,chrgeq,alpha2) &
!$omp private(i,j,r,gamij) &
!$omp shared(Xvec,Xfac,Amat)
!$omp do schedule(runtime)
   ! prepare A matrix
   do i = 1, mol%n
      do j = 1, i-1
         r = norm2(mol%xyz(:,j) - mol%xyz(:,i))
         gamij = 1.0_wp/sqrt(alpha2(i)+alpha2(j))
         Amat(j,i) = erf(gamij*r)/r
         Amat(i,j) = Amat(j,i)
      enddo
      Amat(i,i) = chrgeq%gam(i) + sqrt2pi/sqrt(alpha2(i))
   enddo
!$omp enddo
!$omp endparallel

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

!$omp parallel default(none) &
!$omp shared(mol,dcndr,alpha2,Xfac,Xtmp) &
!$omp private(i,j,rij,r,r2,gamij,arg,erft,dtmp) &
!$omp shared(dXvecdr,dAmatdr) &
!$omp reduction(+:Afac)
!$omp do schedule(runtime)
   do i = 1, mol%n
      dXvecdr(:,:,i) = dcndr(:,:,i)*Xfac(i) ! merge dX and dA for speedup
      do j = 1, i-1
         rij = mol%xyz(:,i) - mol%xyz(:,j)
         r2 = sum(rij**2)
         r  = sqrt(r2)
         gamij = 1.0_wp/sqrt(alpha2(i) + alpha2(j))
         arg = gamij**2*r2
         erft = erf(gamij*r)/r
         dtmp = 2.0_wp*gamij*exp(-arg)/(sqrtpi*r2) - erft/r2
         dAmatdr(:,i,j) = +dtmp*rij*Xtmp(i)
         dAmatdr(:,j,i) = -dtmp*rij*Xtmp(j)
         Afac(:,i) = +dtmp*rij*Xtmp(j) + Afac(:,i)
         Afac(:,j) = -dtmp*rij*Xtmp(i) + Afac(:,j)
      enddo
   enddo
!$omp enddo
!$omp endparallel
   call compute_gb_damat(gbsa,Xtmp,gborn,ghb,dAmatdr,Afac,lverbose)
   gsolv = gsolv + gborn + ghb + gshift
endif do_molecular_gradient

   if (lgrad) then
   call dgemv('n',3*mol%n,m,+1.0_wp,dAmatdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
   call dgemv('n',3*mol%n,m,+1.0_wp,dXvecdr,3*mol%n,Xtmp,1,1.0_wp,gradient,1)
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
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dAmatdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   call dgemm('n','n',3*mol%n,m,m,-1.0_wp,dXvecdr,3*mol%n,Ainv,m, &
      &       1.0_wp,dqdr,3*mol%n)
   !print'(/,"analytical gradient")'
   !print'(3f20.14)',dqdr(:,:,:n)

endif do_partial_charge_derivative

! ------------------------------------------------------------------------
!  Clean up
! ------------------------------------------------------------------------
   if (allocated(alpha2)) deallocate(alpha2)
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
