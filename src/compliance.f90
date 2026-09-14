module xtb_compliance
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : crossProd
   !---------------------------------------------------------------------------
   ! Compliance constants as the pseudoinverse of the Hessian in internal
   ! coordinates:   C = B H^+ B^T
   !
   ! Equivalent to the projected-force-constant route
   !    F = G^+ B H B^T G^+ ,  G = B B^T ,  C = F^+ ,
   ! but neither G nor F nor their pseudoinverses are needed: H^+ is formed
   ! directly in the non-rigid subspace, which is where the rows of B live.
   !
   ! H is projected out of the space of translations and rotations first.
   ! Raw numerical Hessians (the ones handed over by the frequency code)
   ! carry residual curvature along those directions, and its reciprocal
   ! would otherwise dominate C.
   !
   ! Works for any nint and any coordinate set (diatomic, linear, mixed,
   ! general, redundant).  nint is passed explicitly -- no hardcoded 3N-6.
   !
   ! Ref.: K. Brandhorst, J. Grunenberg, Chem. Soc. Rev. 37 (2008), 1558.
   !---------------------------------------------------------------------------
   implicit none
   private
   public :: compute_compliance

contains

subroutine compute_compliance(unit, H, B, xyz, natoms, nint, C, stat)
   integer, intent(in) :: unit
   integer, intent(in) :: natoms, nint
   real(wp), intent(in) :: H(3*natoms, 3*natoms)
   real(wp), intent(in) :: B(nint, 3*natoms)
   real(wp), intent(in) :: xyz(3, natoms)
   real(wp), intent(out) :: C(nint, nint)
   integer, intent(out) :: stat

   integer :: i, j, ndim, lwork, info, nrigid, nvib, rank_h
   real(wp) :: tol_h, normq, center(3), scale
   real(wp), parameter :: eps_svd = 1.0e-10_wp
   real(wp), allocatable :: Hp(:, :), Q(:, :), W(:), Z(:, :), ZD(:, :), &
      & T1(:, :), T2(:, :), work(:)

   stat = 0; ndim = 3 * natoms
   allocate(Hp(ndim, ndim), Q(ndim, 6), W(ndim), Z(nint, ndim), &
      & ZD(nint, ndim), T1(ndim, 6), T2(6, ndim))

   ! orthonormal basis of the rigid (translation + rotation) space
   center = sum(xyz, dim=2) / natoms
   scale = max(1.0_wp, maxval(abs(xyz)))
   Q = 0.0_wp
   do i = 1, 3
      Q(i::3, i) = 1.0_wp
   end do
   do j = 1, natoms
      Q(3*j-2:3*j, 4) = crossProd([1.0_wp, 0.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 5) = crossProd([0.0_wp, 1.0_wp, 0.0_wp], xyz(:, j) - center)
      Q(3*j-2:3*j, 6) = crossProd([0.0_wp, 0.0_wp, 1.0_wp], xyz(:, j) - center)
   end do

   ! modified Gram-Schmidt; the rotation about the molecular axis of a
   ! linear molecule produces no displacement
   nrigid = 0
   do i = 1, 6
      do j = 1, i - 1
         normq = dot_product(Q(:, i), Q(:, j))
         Q(:, i) = Q(:, i) - normq * Q(:, j)
      end do
      normq = norm2(Q(:, i))
      if (normq < 1.0e-8_wp*scale) then
         Q(:, i) = 0.0_wp
         cycle
      end if
      nrigid = nrigid + 1
      Q(:, i) = Q(:, i) / normq
   end do
   nvib = ndim - nrigid

   ! H averaged into symmetry first: the numerical Hessian handed over by
   ! the frequency code is only symmetric to within its finite-difference
   ! noise
   Hp = 0.5_wp * (H + transpose(H))

   ! Hp = (1 - Q Q^T) Hp (1 - Q Q^T)
   call dgemm("N", "N", ndim, 6, ndim, 1.0_wp, Hp, ndim, Q, ndim, 0.0_wp, T1, ndim)
   call dgemm("N", "T", ndim, ndim, 6, -1.0_wp, T1, ndim, Q, ndim, 1.0_wp, Hp, ndim)
   call dgemm("T", "N", 6, ndim, ndim, 1.0_wp, Q, ndim, Hp, ndim, 0.0_wp, T2, 6)
   call dgemm("N", "N", ndim, ndim, 6, -1.0_wp, Q, ndim, T2, 6, 1.0_wp, Hp, ndim)
   call symmetrise(Hp, ndim)

   ! Hp = V W V^T, overwriting Hp with V
   lwork = -1; allocate(work(1))
   call dsyev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   lwork = int(work(1)); deallocate(work); allocate(work(lwork))
   call dsyev("V", "U", ndim, Hp, ndim, W, work, lwork, info)
   if (info /= 0) then
      write(unit, "(A,I0)") "compute_compliance: DSYEV(H) info=", info
      stat = info; return
   end if

   ! C = B V D V^T B^T = Z (Z D)^T,  Z = B V,  D = diag(1/w_i)
   ! Modes of the projected Hessian that are still zero are dropped.
   tol_h = eps_svd * maxval(abs(W))
   rank_h = count(abs(W) > tol_h)
   if (rank_h /= nvib) then
      write(unit, "(A,I0,A,I0)") &
         & "  Note: Hessian rank=", rank_h, " /= 3N-rigid=", nvib
   end if
   call dgemm("N", "N", nint, ndim, ndim, 1.0_wp, B, nint, Hp, ndim, 0.0_wp, Z, nint)
   ZD = Z
   do i = 1, ndim
      if (abs(W(i)) <= tol_h) then
         ZD(:, i) = 0.0_wp
      else
         ZD(:, i) = ZD(:, i) / W(i)
      end if
   end do
   call dgemm("N", "T", nint, nint, ndim, 1.0_wp, Z, nint, ZD, nint, 0.0_wp, C, nint)
   C = 0.5_wp * (C + transpose(C))
   deallocate(Hp, Q, W, Z, ZD, T1, T2, work)

end subroutine compute_compliance

!> Average the two triangles; mirroring one of them would turn rounding
!> noise into a symmetric perturbation.
subroutine symmetrise(A, n)
   integer, intent(in) :: n
   real(wp), intent(inout) :: A(n, n)

   integer :: i, j
   do i = 1, n
      do j = i + 1, n
         A(j, i) = 0.5_wp * (A(i, j) + A(j, i))
         A(i, j) = A(j, i)
      end do
   end do
end subroutine symmetrise

end module xtb_compliance
