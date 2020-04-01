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
   use xtb_type_environment, only : TEnvironment
   use xtb_gfn0param, alp => alpg
   use xtb_coulomb_gaussian
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

contains


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
      call get_coulomb_derivs(mol, chrgeq, Xtmp, rTrans, gTrans, cf, &
         & dAmatdr, dAmatdL, Afac)
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
      call env%error("Could not factorize Coulomb matrix", source)
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
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


! ======================================================================
!  Modified Version of eeq_chrgeq routine that reads also the chrgeq construct
! ======================================================================
subroutine eeq_chrgeq_gbsa(mol,env,chrgeq,gbsa,cn,dcndr,q,dqdr, &
      &                    energy,gsolv,gradient,lverbose,lgrad,lcpq)

   use xtb_mctc_convert

   use xtb_type_molecule
   use xtb_type_param

   use xtb_disp_ncoord
   use xtb_solv_gbobc

   use xtb_pbc_tools

   implicit none

   character(len=*), parameter :: source = "eeq_chrgeq_gbsa"

! ------------------------------------------------------------------------
!  Input
! ------------------------------------------------------------------------
   type(TMolecule),intent(in) :: mol     !< molecular structure information
   type(TSolvent), intent(in) :: gbsa
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
      call env%error("Could not factorize Coulomb matrix", source)
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in Ainv matrix
   ! Ainv matrix is overwritten with lower triangular part of A⁻¹
   call dsytri('L',m,Ainv,m,ipiv,work,info)
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
