!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  COPYRIGHT (C) 2015 by Filippo Lipparini, Benjamin Stamm, Paolo Gatto        !
!  Eric Cancès, Yvon Maday, Jean-Philip Piquemal, Louis Lagardère and          !
!  Benedetta Mennucci.                                                         !
!                             ALL RIGHT RESERVED.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! A modular implementation of COSMO using a domain decomposition linear scaling
! strategy.
!
! This code is governed by the LGPL license and abiding by the rules of
! distribution of free software.
! This program is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
! or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Lesser General Public License for more details.

module xtb_solv_ddcosmo_solver
   use xtb_solv_ddcosmo_core, only : TDomainDecomposition, hsnorm, calcv, &
      & intrhs, prtsph, adjrhs
   implicit none
   private

   public :: jacobi_diis, lx, lstarx, ldm1x, hnorm

   integer, parameter :: wp = selected_real_kind(15)


contains


!> Jacobi/DIIS solver
subroutine jacobi_diis(ddCosmo, n, lprint, diis_max, norm, tol, rhs, x, &
      & n_iter, ok, matvec, dm1vec, u_norm)

   !  Error source
   character(len=*), parameter :: source = 'ddcosmo_solver::jacobi_diis'

   !> Instance of the COSMO model
   type(TDomainDecomposition), intent(in) :: ddCosmo

   !> Integer, input, size of the matrix
   integer, intent(in) :: n

   !> Integer, input, number of points to be used for diis extrapolation
   !  if diis_max = 0, this is just a Jacobi solver.
   integer, intent(in) :: diis_max

   !> Integer, input, norm to be used to evaluate convergence
   !  1: max |x_new - x|
   !  2: rms (x_new - x)
   !  3: rms (x_new - x) and max |x_new - x|
   !  4: norm computed by the user-provided function u_norm(n, x)
   integer, intent(in) :: norm

   !> Integer, input, printing flag.
   integer, intent(in) :: lprint

   !> Real, input, convergence criterion. if norm = 3, convergence is
   !  achieved when rms (x_new - x) < tol and max |x_new - x| < 10*tol.
   real(wp), intent(in) :: tol

   !> Real, dimension(n), input, right-hand side of the linear system
   real(wp), intent(in) :: rhs(:, :)

   !> Real, dimension(n). In input, a guess of the solution (can be zero).
   !  In output, the solution
   real(wp), intent(inout), contiguous, target :: x(:, :)

   !> Integer, in input, the maximum number of iterations. In output,
   !  the number of iterations needed to converge.
   integer, intent(inout) :: n_iter

   !> Logical, output, T if the solver converged, false otherwise.
   logical, intent(inout) :: ok

   !> External, subroutine to compute the required matrix-vector multiplication
   !  format: subroutine matvec(n, x, y)
   procedure(lx) :: matvec

   !> External, subroutine to apply the inverse diagonal matrix to a vector.
   !  format: subroutine dm1vec(n, x, y)
   procedure(ldm1x) :: dm1vec

   !> External, optional function to compute the norm of a vector.
   !> Format: real(wp) function u_norm(n, x)
   procedure(hnorm), optional :: u_norm

   integer :: it, nmat, istatus, lenb
   real(wp) :: rms_norm, max_norm, tol_max
   logical :: dodiis

   real(wp), allocatable :: y(:, :), x_diis(:, :), e_diis(:, :), bmat(:, :)
   real(wp), allocatable, target :: x_new(:, :)
   real(wp), pointer :: xptr(:)

   character(len=*), parameter :: f100 = "(t3, 'iter=', i4, ' residual norm (rms, max): ', 2d14.4)"
   character(len=*), parameter :: f110 = "(t3, 'iter=', i4, ' residual norm (', a, '): ', d14.4)"
   character(len=*), parameter :: f120 = "(t3, 'iter=', i4, ' residual norm: ', d14.4)"

   ! check inputs
   if ((norm == 4) .and. (.not.present(u_norm))) then
      error stop 'must provide a function norm(n, x) to evaluate the norm of the increment'
      return
   end if

   ! DIIS extrapolation flag
   dodiis =  (diis_max /= 0)

   ! set tolerance
   tol_max = 10.0_wp * tol

   ! extrapolation required
   if (dodiis) then

      ! allocate workspaces
      lenb = diis_max + 1
      allocate(x_diis(n, diis_max), e_diis(n, diis_max), bmat(lenb, lenb))

      ! an enigmatic constant
      nmat = 1
   end if

   ! allocate workspaces
   allocate(x_new, mold=rhs)
   allocate(y, mold=rhs)

   ! Jacobi iterations
   do it = 1, n_iter

      ! y = rhs - O x
      call matvec(ddCosmo, n, x, y)
      y = rhs - y

      ! x_new = D^-1 y
      call dm1vec(ddCosmo, n, y, x_new)

      ! DIIS extrapolation
      if (dodiis) then
         x_diis(:, nmat) = reshape(x_new, [size(x_new)])
         e_diis(:, nmat) = reshape(x_new - x, [size(x_new)])

         xptr(1:size(x_new)) => x_new
         call diis(n, nmat, diis_max, x_diis, e_diis, bmat, xptr)
      end if

      ! increment
      x = x_new - x

      ! rms/max norm of increment
      if (norm <= 3) then

         ! compute norm
         xptr(1:size(x)) => x
         call rmsvec(n, xptr, rms_norm, max_norm)

         ! check norm
         if (norm == 1) then
            ok = (rms_norm < tol)
         elseif (norm == 2) then
            ok = (max_norm < tol)
         else
            ok = (rms_norm < tol) .and. (max_norm < tol_max)
         end if

         ! user-provided norm of increment
      elseif (norm == 4) then

         ! just a placeholder for printing
         max_norm = -1.0_wp

         ! compute norm
         rms_norm = u_norm(ddCosmo, n, x)

         ! check norm
         ok = (rms_norm < tol)

      end if

      ! printing
      if (lprint > 0) then
         if (norm == 1) then
            write(*, f110) it, 'max', max_norm
         else if (norm == 2) then
            write(*, f110) it, 'rms', rms_norm
         else if (norm == 3) then
            write(*, f100) it, rms_norm, max_norm
         else if (norm == 4) then
            write(*, f120) it, rms_norm
         end if
      end if

      ! update
      x = x_new

      ! EXIT Jacobi loop here

      if (ok) exit

   end do

   ! record number of Jacobi iterations
   n_iter = it

endsubroutine jacobi_diis


subroutine diis(n, nmat, ndiis, x, e, b, xnew)
   integer, intent(in) :: n, ndiis
   integer, intent(inout) :: nmat
   real(wp), intent(inout) :: x(:, :) !< [n, ndiis]
   real(wp), intent(inout) :: e(:, :) !< [n, ndiis]
   real(wp), intent(inout) :: b(:, :) !< [ndiis+1, ndiis+1]
   real(wp), intent(inout) :: xnew(:) !< [n]

   integer :: nmat1, i, istatus
   integer :: j, k
   logical :: ok

   real(wp), allocatable :: bloc(:, :), cex(:, :)

   if (nmat >= ndiis) then
      do j = 2, nmat - 10
         do k = 2, nmat - 10
            b(j, k) = b(j+10, k+10)
         end do
      end do
      do j = 1, nmat - 10
         x(:, j) = x(:, j+10)
         e(:, j) = e(:, j+10)
      end do
      nmat = nmat - 10
   end if
   nmat1 = nmat + 1
   allocate(bloc(nmat1, nmat1), cex(nmat1, 1), stat=istatus)
   if (istatus /= 0) then
      nmat = 1
      return
   end if

   call makeb(n, nmat, ndiis, e, b)
   bloc = b(1:nmat1, 1:nmat1)
   cex(:, :) = 0.0_wp
   cex(1, 1) = 1.0_wp
   call gjinv(nmat1, 1, bloc, cex, ok)
   if (.not. ok) then
      nmat = 1
      error stop "Upps!"
      return
   end if
   xnew = 0.0_wp
   do i = 1, nmat
      xnew = xnew + cex(i+1, 1)*x(:, i)
   end do
   nmat = nmat + 1

end subroutine diis


pure subroutine makeb(n, nmat, ndiis, e, b)
   integer, intent(in) :: n, nmat, ndiis
   real(wp), intent(in) :: e(n, ndiis)
   real(wp), intent(inout) :: b(ndiis+1, ndiis+1)

   integer :: i
   real(wp) :: bij

   ! 1st built
   if (nmat == 1) then
      !       [ 0 |  1  ]
      !   b = [ --+---- ]
      !       [ 1 | e*e ]
      b(1, 1) = 0.0_wp
      b(1, 2) = 1.0_wp
      b(2, 1) = 1.0_wp
      b(2, 2) = dot_product(e(:, 1), e(:, 1))
   else
      ! subsequent builts
      ! first, update the lagrangian line:
      b(nmat+1, 1) = 1.0_wp
      b(1, nmat+1) = 1.0_wp
      ! now, compute the new matrix elements:
      do i = 1, nmat - 1
         bij = dot_product(e(:, i), e(:, nmat))
         b(nmat+1, i+1) = bij
         b(i+1, nmat+1) = bij
      end do
      b(nmat+1, nmat+1) = dot_product(e(:, nmat), e(:, nmat))
   end if

end subroutine makeb


pure subroutine gjinv(n, nrhs, a, b, ok)
   integer, intent(in) :: n, nrhs
   logical, intent(inout) :: ok
   real(wp), intent(inout) :: a(:, :) ! [n, n]
   real(wp), intent(inout) :: b(:, :) ! [n, nrhs]

   integer :: i, j, k, irow, icol, istatus
   real(wp) :: big, dum, pinv

   integer, allocatable :: indxc(:), indxr(:), piv(:)
   real(wp), allocatable :: scr(:)

   ok = .false.

   allocate(indxc(n), indxr(n), piv(n), stat=istatus)
   if (istatus /= 0) then
      return
   end if
   allocate (scr(n), stat=istatus)
   if (istatus /= 0) then
      return
   end if

   piv = 0

   irow = 0
   icol = 0
   do i = 1, n
      big = 0.0_wp
      do j = 1, n
         if (piv(j) /= 1) then
            do k = 1, n
               if (piv(k) == 0) then
                  if (abs(a(j, k)) > big) then
                     big  = abs(a(j, k))
                     irow = j
                     icol = k
                  end if
               end if
            end do
         end if
      end do

      piv(icol) = piv(icol) + 1
      if (piv(icol) > 1) then
         !call warning('singular matrix', source)
         return
      end if
      if (irow /= icol) then
         scr         = a(irow, :)
         a(irow, :)   = a(icol, :)
         a(icol, :)   = scr
         scr(1:nrhs) = b(irow, :)
         b(irow, :)   = b(icol, :)
         b(icol, :)   = scr(1:nrhs)
      end if

      indxr(i) = irow
      indxc(i) = icol

      if (a(icol, icol) == 0.0_wp) then
         !call warning('singular matrix', source)
         return
      end if

      pinv = 1.0_wp/a(icol, icol)
      a(icol, icol) = 1.0_wp
      a(icol, :) = a(icol, :)*pinv
      b(icol, :) = b(icol, :)*pinv

      do j = 1, n
         if (j /= icol) then
            dum       = a(j, icol)
            a(j, icol) = 0.0_wp
            a(j, :)    = a(j, :) - a(icol, :)*dum
            b(j, :)    = b(j, :) - b(icol, :)*dum
         end if
      end do
   end do

   do j = n, 1, -1
      if (indxr(j) /= indxc(j)) then
         scr           = a(:, indxr(j))
         a(:, indxr(j)) = a(:, indxc(j))
         a(:, indxc(j)) = scr
      end if
   end do

   ok = .true.

end subroutine gjinv


!> Compute root-mean-square and max norm
pure subroutine rmsvec(n, v, vrms, vmax)
   integer, intent(in) :: n
   real(wp), intent(in) :: v(:) !< [n]
   real(wp), intent(inout) :: vrms, vmax

   integer :: i

   ! initialize
   vrms = 0.0_wp
   vmax = 0.0_wp

   ! loop over entries
   do i = 1, n

      ! max norm
      vmax = max(vmax, abs(v(i)))

      ! rms norm
      vrms = vrms + v(i)*v(i)

   end do

   ! the much neglected square root
   vrms = sqrt(vrms/dble(n))

end subroutine rmsvec


!> Given a vector x, compute y = Lx, where L is the ddCOSMO matrix
!  (off-diagonal blocks only).
subroutine lx(ddCosmo, n, x, y)

   type(TDomainDecomposition), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
   real(wp), intent(inout) :: y(:, :) ! [ddCosmo%nylm, ddCosmo%nat]

   integer :: iat, istatus
   real(wp), allocatable :: pot(:), vplm(:), basloc(:), vcos(:), vsin(:)

   ! allocate workspaces
   allocate(pot(ddCosmo%ngrid), vplm(ddCosmo%nylm), basloc(ddCosmo%nylm), &
      & vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))

   if (ddCosmo%iprint >= 5) call prtsph(ddCosmo, 'X', ddCosmo%nat, 0, x)

   ! initialize
   y(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) shared(ddCosmo, y, x) &
   !$omp private(iat, pot, basloc, vplm, vcos, vsin)
   ! loop over spheres
   do iat = 1, ddCosmo%nat

      ! compute NEGATIVE action of off-digonal blocks
      call calcv(ddCosmo, .false., iat, pot, x, basloc, vplm, vcos, vsin)
      call intrhs(ddCosmo, iat, pot, y(:, iat))

      ! action of off-diagonal blocks
      y(:, iat) = - y(:, iat)

   end do

   if (ddCosmo%iprint >= 5) call prtsph(ddCosmo, 'LX (off diagonal)', ddCosmo%nat, 0, y)

end subroutine lx


!> Given a vector x, compute y = L*x, where L* is the adjoint ddCOSMO matrix.
!  if dodiag is set to .true., L includes the diagonal blocks, otherwise
!  L only includes the off-diagonal ones.
subroutine lstarx(ddCosmo, n, x, y)

   type(TDomainDecomposition), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
   real(wp), intent(inout) :: y(:, :) ! [ddCosmo%nylm, ddCosmo%nat]

   integer :: iat, ig, istatus
   real(wp), allocatable :: xi(:, :), vplm(:), basloc(:), vcos(:), vsin(:)

   ! allocate workspaces
   allocate(xi(ddCosmo%ngrid, ddCosmo%nat), vplm(ddCosmo%nylm), &
      & basloc(ddCosmo%nylm), vcos(ddCosmo%lmax+1), vsin(ddCosmo%lmax+1))

   if (ddCosmo%iprint >= 5) call prtsph(ddCosmo, 'X', ddCosmo%nat, 0, x)

   ! initilize
   y(:, :) = 0.0_wp

   ! expand x over spherical harmonics
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(ddCosmo, xi, x) private(iat, ig)
   ! loop over spheres
   do iat = 1, ddCosmo%nat
      ! loop over griwpoints
      do ig = 1, ddCosmo%ngrid
         xi(ig, iat) = dot_product(x(:, iat), ddCosmo%basis(:, ig))
      end do
   end do

   ! compute action
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(ddCosmo, xi, y) private(iat, basloc, vplm, vcos, vsin)
   ! loop over spheres
   do iat = 1, ddCosmo%nat

      ! compute NEGATIVE action of off-digonal blocks
      call adjrhs(ddCosmo, iat, xi, y(:, iat), basloc, vplm, vcos, vsin)

      ! action of off-diagonal blocks
      y(:, iat) = - y(:, iat)

   end do

   if (ddCosmo%iprint >= 5) call prtsph(ddCosmo, 'L*X (off-diagonal)', ddCosmo%nat, 0, y)

end subroutine lstarx


!> Given a vector x, apply the inverse diagonal (block) of the L matrix:
pure subroutine ldm1x(ddCosmo, n, x, y)
   type(TDomainDecomposition), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
   real(wp), intent(inout) :: y(:, :) ! [ddCosmo%nylm, ddCosmo%nat]

   integer :: iat

   ! loop over spheres
   do iat = 1, ddCosmo%nat
      ! apply inverse
      y(:, iat) = ddCosmo%facl*x(:, iat)
   end do

end subroutine ldm1x


!> Compute the h^-1/2 norm of the increment on each sphere, then take the
!  rms value.
real(wp) function hnorm(ddCosmo, n, x)
   type(TDomainDecomposition), intent(in) :: ddCosmo
   integer, intent(in) :: n
   real(wp), intent(in) :: x(:, :) ! [ddCosmo%nylm, ddCosmo%nat]
   integer :: iat, istatus
   real(wp) :: vrms, vmax
   real(wp), allocatable :: u(:)

   ! allocate workspace
   allocate(u(ddCosmo%nat))

   ! loop over spheres
   do iat = 1, ddCosmo%nat
      ! compute norm contribution
      call hsnorm(ddCosmo, x(:, iat), u(iat))
   end do

   ! compute rms of norms
   call rmsvec(ddCosmo%nat, u, vrms, vmax)

   ! return value
   hnorm = vrms

end function hnorm


end module xtb_solv_ddcosmo_solver
