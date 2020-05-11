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

!> Implementation of the electronegativity equilibration model
module xtb_xtb_eeq
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_blas, only : mctc_dot, mctc_symv
   use xtb_mctc_lapack, only : lapack_sytrf, lapack_sytri
   use xtb_mctc_la, only : contract
   use xtb_type_coulomb, only : TCoulomb
   use xtb_type_environment, only : TEnvironment
   use xtb_type_molecule, only : TMolecule
   implicit none
   private

   public :: TENEquilibration, init


   !> Electronegativity equilibration model
   type :: TENEquilibration

      !> Electronegativity of each species
      real(wp), allocatable :: chi(:, :)

      !> Coordination number dependence of each species
      real(wp), allocatable :: kcn(:, :)

      !> Chemical hardness of each species
      real(wp), allocatable :: gam(:, :)

      !> Solver for the linear system
      procedure(solve), nopass, pointer :: solve_lineq => null()

      !> Inversion of the linear system (for response)
      procedure(solve), nopass, pointer :: invert_lineq => null()

   contains

      !> Evaluate charge equilibration model
      procedure :: chargeEquilibration

   end type TENEquilibration


   !> Initialize equilibration model
   interface init
      module procedure :: initENEquilibration
      module procedure :: initENEquilibrationAtomic
   end interface init


   abstract interface
   !> Interface to a solver for linear systems
   subroutine solve(env, mat, rhs, vec)
   import :: wp, TEnvironment

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Matrix containing the linear system, might contain inverse on exit
   real(wp), intent(inout) :: mat(:, :)

   !> Right hand side of the linear system
   real(wp), intent(in) :: rhs(:)

   !> Solution vector of the linear system
   real(wp), intent(out), target :: vec(:)

   end subroutine solve
   end interface


contains


!> Initialize equilibration model form parameters
subroutine initENEquilibrationAtomic(self, env, chi, kcn, gam, num)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_initENEquilibrationAtomic'

   !> Instance of the equilibration model
   type(TENEquilibration), intent(out) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Electronegativity of each species
   real(wp), intent(in) :: chi(:)

   !> Coordination number dependence of each species
   real(wp), intent(in), optional :: kcn(:)

   !> Chemical hardness of each species
   real(wp), intent(in), optional :: gam(:)

   !> Atomic number for each id
   integer, intent(in), optional :: num(:)

   integer :: ii

   if (present(num)) then
      allocate(self%chi(1, size(num)))
      do ii = 1, size(num)
         self%chi(1, ii) = chi(num(ii))
      end do
   else
      self%chi = reshape(chi, [1, size(chi)])
   end if

   self%solve_lineq => solve_sysv
   self%invert_lineq => solve_sytri

   if (present(kcn)) then
      if (any(shape(kcn) /= shape(chi))) then
         call env%error("Cannot initialize CN dependence of EN", source)
         allocate(self%kcn(1, size(self%chi)))
         self%kcn(:, :) = 0.0_wp
      else
         if (present(num)) then
            allocate(self%kcn(1, size(num)))
            do ii = 1, size(num)
               self%kcn(1, ii) = kcn(num(ii))
            end do
         else
            self%kcn = reshape(kcn, [1, size(kcn)])
         end if
      end if
   end if

   if (present(gam)) then
      if (any(shape(gam) /= shape(chi))) then
         call env%error("Cannot initialize chemical hardness", source)
         allocate(self%gam(1, size(self%chi)))
         self%gam(:, :) = 0.0_wp
      else
         if (present(num)) then
            allocate(self%gam(1, size(num)))
            do ii = 1, size(num)
               self%gam(1, ii) = gam(num(ii))
            end do
         else
            self%gam = reshape(gam, [1, size(gam)])
         end if
      end if
   end if

end subroutine initENEquilibrationAtomic


!> Initialize equilibration model form parameters
subroutine initENEquilibration(self, env, chi, kcn, gam, num)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_initENEquilibration'

   !> Instance of the equilibration model
   type(TENEquilibration), intent(out) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Electronegativity of each species
   real(wp), intent(in) :: chi(:, :)

   !> Coordination number dependence of each species
   real(wp), intent(in), optional :: kcn(:, :)

   !> Chemical hardness of each species
   real(wp), intent(in), optional :: gam(:, :)

   !> Atomic number for each id
   integer, intent(in), optional :: num(:)

   integer :: ii

   if (present(num)) then
      allocate(self%chi(size(chi, dim=1), size(num)))
      do ii = 1, size(num)
         self%chi(:, ii) = chi(:, num(ii))
      end do
   else
      self%chi = chi
   end if

   self%solve_lineq => solve_sysv
   self%invert_lineq => solve_sytri

   if (present(kcn)) then
      if (any(shape(kcn) /= shape(chi))) then
         call env%error("Cannot initialize CN dependence of EN", source)
         allocate(self%kcn(size(self%chi, 1), size(self%chi, 2)))
         self%kcn(:, :) = 0.0_wp
      else
         if (present(num)) then
            allocate(self%kcn(size(kcn, dim=1), size(num)))
            do ii = 1, size(num)
               self%kcn(:, ii) = kcn(:, num(ii))
            end do
         else
            self%kcn = kcn
         end if
      end if
   end if

   if (present(gam)) then
      if (any(shape(gam) /= shape(chi))) then
         call env%error("Cannot initialize chemical hardness", source)
         allocate(self%gam(size(self%chi, 1), size(self%chi, 2)))
         self%gam(:, :) = 0.0_wp
      else
         if (present(num)) then
            allocate(self%gam(size(gam, dim=1), size(num)))
            do ii = 1, size(num)
               self%gam(:, ii) = gam(:, num(ii))
            end do
         else
            self%gam = gam
         end if
      end if
   end if

end subroutine initENEquilibration


subroutine chargeEquilibration(self, env, mol, coulomb, cn, dcndr, dcndL, &
      & energy, gradient, sigma, qat, qsh, dqdr, dqdL)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_chargeEquilibration'

   !> Instance of the equilibration model
   class(TENEquilibration), intent(inout) :: self

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   !> Instance of the Coulomb evaluator
   class(TCoulomb), intent(inout) :: coulomb

   !> Coordination number
   real(wp), intent(in) :: cn(:)

   !> Derivative of the coordination number w.r.t. Cartesian coordinates
   real(wp), intent(in) :: dcndr(:, :, :)

   !> Derivative of the coordination number w.r.t. strain derivatives
   real(wp), intent(in) :: dcndL(:, :, :)

   !> Electrostatic energy
   real(wp), intent(inout), optional :: energy

   !> Molecular gradient
   real(wp), intent(inout), optional :: gradient(:, :)

   !> Strain derivatives
   real(wp), intent(inout), optional :: sigma(:, :)

   !> Atomic partial charges
   real(wp), intent(out), optional :: qat(:)

   !> Shell-resolved partial charges
   real(wp), intent(out), optional :: qsh(:)

   !> Derivative of the partial charges w.r.t. Cartesian coordinates
   real(wp), intent(out), optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain derivatives
   real(wp), intent(out), optional :: dqdL(:, :, :)

   real(wp), allocatable :: jmat(:, :), inv(:, :)
   real(wp), allocatable :: djdr(:, :, :), djdtr(:, :), djdL(:, :, :)
   real(wp), allocatable :: xvec(:), qvec(:), dxdcn(:), shift(:)
   real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)

   logical :: deriv, response, exitRun
   integer :: nsh, ndim
   integer :: iat, ish, ii, iid
   real(wp) :: tmp

   deriv = present(gradient) .or. present(sigma)
   response = present(dqdr) .or. present(dqdL)
   nsh = sum(coulomb%itbl(2, :))
   ndim = nsh + 1
   allocate(jmat(ndim, ndim))
   allocate(xvec(ndim))
   allocate(qvec(ndim))
   allocate(dxdcn(nsh))

   call coulomb%getCoulombMatrix(mol, jmat)

   do iat = 1, mol%n
      ii = coulomb%itbl(1, iat)
      iid = mol%id(iat)
      do ish = 1, coulomb%itbl(2, iat)
         tmp = self%kcn(ish, iid) / (sqrt(cn(iat)) + 1.0e-14_wp)
         xvec(ii+ish) = -self%chi(ish, iid) + tmp*cn(iat)
         dxdcn(ii+ish) = 0.5_wp*tmp
         jmat(ii+ish, ii+ish) = jmat(ii+ish, ii+ish) + self%gam(ish, iid)
         jmat(ii+ish, ndim) = 1.0_wp
         jmat(ndim, ii+ish) = 1.0_wp
      end do
   end do
   xvec(ndim) = mol%chrg
   jmat(ndim, ndim) = 0.0_wp

   inv = jmat
   if (response) then
      call self%invert_lineq(env, inv, xvec, qvec)
   else
      call self%solve_lineq(env, inv, xvec, qvec)
   end if

   call env%check(exitRun)
   if (exitRun) then
      call env%error("Solving linear equations failed", source)
      return
   end if

   if (present(qat)) then
      qat(:) = 0.0_wp
      do iat = 1, mol%n
         ii = coulomb%itbl(1, iat)
         do ish = 1, coulomb%itbl(2, iat)
            qat(iat) = qat(iat) + qvec(ii+ish)
         end do
      end do
   end if

   if (present(qsh)) then
      qsh(:) = qvec(:nsh)
   end if

   if (present(energy)) then
      shift = xvec
      call mctc_symv(jmat, qvec, shift, alpha=0.5_wp, beta=-1.0_wp)
      energy = energy + mctc_dot(qvec, shift)
   end if

   if (deriv .or. response) then
      allocate(djdr(3, mol%n, ndim))
      allocate(djdtr(3, ndim))
      allocate(djdL(3, 3, ndim))
      allocate(dxdr(3, mol%n, ndim))
      allocate(dxdL(3, 3, ndim))
      call coulomb%getCoulombDerivs(mol, qvec, djdr, djdtr, djdL)
      do iat = 1, mol%n
         ii = coulomb%itbl(1, iat)
         do ish = 1, coulomb%itbl(2, iat)
            dxdr(:, :, ii+ish) = -dxdcn(ii+ish)*dcndr(:, :, iat)
            dxdL(:, :, ii+ish) = -dxdcn(ii+ish)*dcndL(:, :, iat)
         end do
      end do
      dxdr(:, :, ndim) = 0.0_wp
      dxdL(:, :, ndim) = 0.0_wp

      if (present(gradient)) then
         call contract(djdr, qvec, gradient, beta=1.0_wp)
         call contract(dxdr, qvec, gradient, beta=1.0_wp)
      end if

      if (present(sigma)) then
         call contract(djdL, qvec, sigma, beta=1.0_wp)
         call contract(dxdL, qvec, sigma, beta=1.0_wp)
      end if

      if (response) then
         do iat = 1, mol%n
            ii = coulomb%itbl(1, iat)
            do ish = 1, coulomb%itbl(2, iat)
               djdr(:, iat, ii+ish) = djdr(:, iat, ii+ish) + djdtr(:, ii+ish)
            end do
         end do

         if (present(dqdr)) then
            dqdr(:, :, :) = 0.0_wp
            call contract(djdr, inv(:, :nsh), dqdr, alpha=-1.0_wp, beta=0.0_wp)
            call contract(dxdr, inv(:, :nsh), dqdr, alpha=-1.0_wp, beta=1.0_wp)
         end if

         if (present(dqdL)) then
            dqdL(:, :, :) = 0.0_wp
            call contract(djdL, inv(:, :nsh), dqdL, alpha=-2.0_wp, beta=0.0_wp)
            call contract(dxdL, inv(:, :nsh), dqdL, alpha=-1.0_wp, beta=1.0_wp)
         end if
      end if
   end if

end subroutine chargeEquilibration


subroutine solve_sysv(env, mat, rhs, vec)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_sysv'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   real(wp), intent(inout) :: mat(:, :)

   real(wp), intent(in) :: rhs(:)

   real(wp), intent(out), target :: vec(:)

   integer :: lwork, ndim, info
   real(wp), pointer :: ptr(:, :)
   real(wp) :: test(1)
   real(wp), allocatable :: work(:)
   integer, allocatable :: ipiv(:)

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim .or. size(vec) < ndim .or. size(rhs) < ndim) then
      call env%error("A not so carefully crafted algorithm did mess up", source)
      return
   end if

   vec(:ndim) = rhs(:ndim)
   allocate(ipiv(ndim))
   ptr(1:ndim, 1:1) => vec(1:ndim)

   ! assume work space query, set best value to test after first dsysv call
   call dsysv('l', ndim, 1, mat, ndim, ipiv, ptr, ndim, test, -1, info)
   lwork = int(test(1))
   allocate(work(lwork))

   call dsysv('l', ndim, 1, mat, ndim, ipiv, ptr, ndim, work, lwork, info)

   if (info > 0) then
      call env%error("LAPACK linear equation solver failed", source)
      return
   end if

end subroutine solve_sysv


subroutine solve_sytri(env, mat, rhs, vec)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_sytri'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   real(wp), intent(inout) :: mat(:, :)

   real(wp), intent(in) :: rhs(:)

   real(wp), intent(out), target :: vec(:)

   integer :: ndim
   logical :: exitRun

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim .or. size(vec) < ndim .or. size(rhs) < ndim) then
      call env%error("A not so carefully crafted algorithm did mess up", source)
      return
   end if

   call invert_sytri(env, mat)

   call env%check(exitRun)
   if (exitRun) return

   call mctc_symv(mat, rhs, vec)

end subroutine solve_sytri


subroutine invert_sytri(env, mat)

   !> Source for error generation
   character(len=*), parameter :: source = 'xtb_eeq_sytri'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   real(wp), intent(inout) :: mat(:, :)

   integer :: ii, jj
   integer :: lwork, ndim, info
   real(wp) :: test(1)
   real(wp), allocatable :: work(:)
   integer, allocatable :: ipiv(:)

   ndim = size(mat, dim=1)
   if (size(mat, dim=2) < ndim) then
      call env%error("A not so carefully crafted algorithm did mess up", source)
      return
   end if

   allocate(ipiv(ndim))

   ! assume work space query, set best value to test after first dsytrf call
   call lapack_sytrf('L', ndim, mat, ndim, ipiv, test, -1, info)
   lwork = int(test(1))
   allocate(work(lwork))

   ! Bunch-Kaufman factorization A = L*D*L**T
   call lapack_sytrf('L', ndim, mat, ndim, ipiv, work, lwork, info)
   if(info > 0)then
      call env%error('Bunch-Kaufman factorization failed', source)
      return
   endif

   ! A⁻¹ from factorized L matrix, save lower part of A⁻¹ in matrix
   ! matrix is overwritten with lower triangular part of A⁻¹
   call lapack_sytri('L', ndim, mat, ndim, ipiv, work, info)
   if (info > 0) then
      call env%error('Bunch-Kaufman inversion failed', source)
      return
   endif

   ! symmetrizes A⁻¹ matrix from lower triangular part of inverse matrix
   do ii = 1, ndim
      do jj = ii+1, ndim
         mat(ii,jj) = mat(jj,ii)
      enddo
   enddo

end subroutine invert_sytri


end module xtb_xtb_eeq
