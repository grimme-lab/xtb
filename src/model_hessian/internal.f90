! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
! Copyright (C) 2026 Leopold M. Seidler
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

!> Shared internal-coordinate traversal for Swart and Lindh Hessians
module xtb_modelhessian_internal
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, bmat_torsion, &
      & bmat_outofplane, oop_angle, bmat_accum_packed, &
      & bmat_accum_pairblock_packed
   use xtb_chargemodel, only : new_charge_model_2019
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : pi
   use xtb_mctc_math, only : crossProd
   use xtb_modelhessian_eeq, only : add_eeq_hessian
   use xtb_param_model_hessian, only : d2_c6, d2_vander
   use xtb_modelhessian_shared, only : rcutoff, getvdw_hess
   use xtb_modelhessian_type, only : TModelHessian
   use xtb_type_param, only : chrg_parameter
   use xtb_type_environment, only : TEnvironment
   implicit none(type, external)
   private

   !> Base for model Hessians defined by redundant internal coordinates
   type, public, abstract, extends(TModelHessian) :: TInternalModelHessianBase
   contains
      procedure :: stretch
      procedure :: bend
      procedure :: torsion
      procedure :: outofplane
      procedure :: add_charge
      procedure(pair_factor_interface), deferred :: pair_factor
   end type TInternalModelHessianBase

   abstract interface
      pure function pair_factor_interface(self, at_i, at_j, r2, dispersion_scale, &
            & outofplane) result(factor)
         import :: TInternalModelHessianBase, wp
         implicit none(type, external)
         class(TInternalModelHessianBase), intent(in) :: self
         integer, intent(in) :: at_i, at_j
         real(wp), intent(in) :: r2, dispersion_scale
         logical, intent(in) :: outofplane

         real(wp) :: factor
      end function pair_factor_interface
   end interface

contains

!> Add bond-stretching and D2 contributions
pure subroutine stretch(self, xyz, n, hess, at, kr, kd, s6, lcutoff, rcut)
   class(TInternalModelHessianBase), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: kr, kd, s6
   logical, intent(inout) :: lcutoff(n, n)
   real(wp), intent(in) :: rcut

   integer :: i, j
   real(wp) :: vec(3), r2, gmm, c6ij, rv, vdw(3, 3)

   do i = 1, n
      do j = 1, i - 1
         lcutoff(i, j) = rcutoff(xyz, i, j, rcut)
         lcutoff(j, i) = lcutoff(i, j)
         vec = xyz(:, i) - xyz(:, j)
         r2 = dot_product(vec, vec)
         c6ij = sqrt(d2_c6(at(i))*d2_c6(at(j)))
         rv = d2_vander(at(i)) + d2_vander(at(j))
         call getvdw_hess(vec, c6ij, s6, rv, vdw)
         gmm = kr * self%pair_factor(at(i), at(j), r2, kd, .false.)
         call bmat_accum_packed(n, hess, [i, j], bmat_bond(vec), gmm)
         call bmat_accum_pairblock_packed(n, hess, i, j, vdw)
      end do
   end do
end subroutine stretch

!> Add angle-bending contributions
pure subroutine bend(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TInternalModelHessianBase), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: force_constant, kd
   logical, intent(in) :: lcutoff(n, n)

   real(wp), parameter :: rzero = 1.0e-10_wp
   integer :: i, j, m, ii
   real(wp) :: vec_ij(3), vec_mi(3), vec_mj(3), cross_vec(3)
   real(wp) :: rmi2, rmi, rmj2, rmj, rij2, rrij, gij, rl2, rl
   real(wp) :: sinphi, cosphi, bmat9(9), linear_bmat(2, 9)

   do m = 1, n
      do i = 1, n
         if (i == m .or. lcutoff(i, m)) cycle
         vec_mi = xyz(:, i) - xyz(:, m)
         rmi2 = dot_product(vec_mi, vec_mi)
         rmi = sqrt(rmi2)
         do j = 1, i - 1
            if (j == m) cycle
            if (lcutoff(j, i) .or. lcutoff(j, m)) cycle
            vec_mj = xyz(:, j) - xyz(:, m)
            rmj2 = dot_product(vec_mj, vec_mj)
            rmj = sqrt(rmj2)
            cosphi = dot_product(vec_mi, vec_mj) / (rmi*rmj)
            if (abs(cosphi - 1.0_wp) < 1.0e-12_wp) cycle
            vec_ij = xyz(:, j) - xyz(:, i)
            rij2 = dot_product(vec_ij, vec_ij)
            rrij = sqrt(rij2)
            gij = force_constant &
               * self%pair_factor(at(m), at(i), rmi2, 0.5_wp*kd, .false.) &
               * self%pair_factor(at(m), at(j), rmj2, 0.5_wp*kd, .false.)
            cross_vec = crossProd(vec_mi, vec_mj)
            rl2 = dot_product(cross_vec, cross_vec)
            if (rl2 < 1.0e-14_wp) then
               rl = 0.0_wp
            else
               rl = sqrt(rl2)
            end if
            if (rmj <= rzero .or. rmi <= rzero .or. rrij <= rzero) cycle
            sinphi = rl / (rmj*rmi)
            if (sinphi > rzero) then
               bmat9 = bmat_angle(vec_mi, vec_mj)
               call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
            else
               linear_bmat = bmat_linbend(vec_mi, vec_mj)
               do ii = 1, 2
                  call bmat_accum_packed(n, hess, [i, m, j], &
                     & linear_bmat(ii, :), gij)
               end do
            end if
         end do
      end do
   end do
end subroutine bend

!> Add torsional contributions using one reversal-safe orientation
pure subroutine torsion(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TInternalModelHessianBase), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: force_constant, kd
   logical, intent(in) :: lcutoff(n, n)

   real(wp), parameter :: a35 = (35.0_wp/180.0_wp) * pi
   real(wp), parameter :: cosfi_max = cos(a35)
   integer :: i, j, k, l, ij, kl
   real(wp) :: txyz(3, 4), c(3, 4), brow12(12)
   real(wp) :: rij(3), rjk(3), rkl(3), rij2, rjk2, rkl2
   real(wp) :: cosfi2, cosfi3, gij, gjk, gkl, tij

   do j = 1, n
      txyz(:, 2) = xyz(:, j)
      do k = 1, n
         if (k == j .or. lcutoff(k, j)) cycle
         txyz(:, 3) = xyz(:, k)
         do i = 1, n
            ij = n * (j - 1) + i
            if (i == j .or. i == k) cycle
            if (lcutoff(i, k) .or. lcutoff(i, j)) cycle
            txyz(:, 1) = xyz(:, i)
            do l = 1, n
               kl = n * (k - 1) + l
               if (ij <= kl) cycle
               if (l == i .or. l == j .or. l == k) cycle
               if (lcutoff(l, i) .or. lcutoff(l, k) .or. lcutoff(l, j)) cycle
               txyz(:, 4) = xyz(:, l)
               rij = xyz(:, i) - xyz(:, j)
               rjk = xyz(:, j) - xyz(:, k)
               rkl = xyz(:, k) - xyz(:, l)
               rij2 = dot_product(rij, rij)
               rjk2 = dot_product(rjk, rjk)
               rkl2 = dot_product(rkl, rkl)
               cosfi2 = dot_product(rij, rjk) / sqrt(rij2*rjk2)
               if (abs(cosfi2) > cosfi_max) cycle
               cosfi3 = dot_product(rkl, rjk) / sqrt(rkl2*rjk2)
               if (abs(cosfi3) > cosfi_max) cycle
               gij = self%pair_factor(at(i), at(j), rij2, 0.5_wp*kd, .false.)
               gjk = self%pair_factor(at(j), at(k), rjk2, 0.5_wp*kd, .false.)
               gkl = self%pair_factor(at(k), at(l), rkl2, 0.5_wp*kd, .false.)
               tij = force_constant * gij * gjk * gkl
               c = bmat_torsion(txyz)
               brow12 = [c(:, 1), c(:, 2), c(:, 3), c(:, 4)]
               call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)
            end do
         end do
      end do
   end do
end subroutine torsion

!> Add out-of-plane contributions
pure subroutine outofplane(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   class(TInternalModelHessianBase), intent(in) :: self
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: force_constant, kd
   logical, intent(in) :: lcutoff(n, n)

   integer :: i, j, k, l
   real(wp) :: txyz(3, 4), c(3, 4), brow12(12)
   real(wp) :: rij(3), rik(3), ril(3), rij2, rik2, ril2
   real(wp) :: cosfi2, cosfi3, cosfi4, gij, gik, gil, tij, tau

   do i = 1, n
      txyz(:, 4) = xyz(:, i)
      do j = 1, n
         if (j == i .or. lcutoff(j, i)) cycle
         txyz(:, 1) = xyz(:, j)
         do k = 1, n
            if (k == i .or. k == j) cycle
            if (lcutoff(k, i) .or. lcutoff(k, j)) cycle
            txyz(:, 2) = xyz(:, k)
            do l = 1, n
               if (l == i .or. l == j .or. l == k) cycle
               if (lcutoff(l, i) .or. lcutoff(l, k) .or. lcutoff(l, j)) cycle
               txyz(:, 3) = xyz(:, l)
               rij = xyz(:, i) - xyz(:, j)
               rik = xyz(:, i) - xyz(:, k)
               ril = xyz(:, i) - xyz(:, l)
               rij2 = dot_product(rij, rij)
               rik2 = dot_product(rik, rik)
               ril2 = dot_product(ril, ril)
               cosfi2 = dot_product(rij, rik) / sqrt(rij2*rik2)
               if (abs(abs(cosfi2) - 1.0_wp) < 1.0e-1_wp) cycle
               cosfi3 = dot_product(rij, ril) / sqrt(rij2*ril2)
               if (abs(abs(cosfi3) - 1.0_wp) < 1.0e-1_wp) cycle
               cosfi4 = dot_product(rik, ril) / sqrt(rik2*ril2)
               if (abs(abs(cosfi4) - 1.0_wp) < 1.0e-1_wp) cycle
               gij = self%pair_factor(at(i), at(j), rij2, 0.5_wp*kd, .true.)
               gik = self%pair_factor(at(i), at(k), rik2, 0.5_wp*kd, .true.)
               gil = self%pair_factor(at(i), at(l), ril2, 0.5_wp*kd, .true.)
               tij = force_constant * gij * gik * gil
               tau = oop_angle(txyz)
               if (abs(tau) > 45.0_wp*(pi/180.0_wp)) cycle
               c = bmat_outofplane(txyz)
               brow12 = [c(:, 4), c(:, 1), c(:, 2), c(:, 3)]
               call bmat_accum_packed(n, hess, [i, j, k, l], brow12, tij)
            end do
         end do
      end do
   end do
end subroutine outofplane

!> Add the neutral-molecule EEQ response contribution
subroutine add_charge(self, env, xyz, n, hess, at, kq)
   class(TInternalModelHessianBase), intent(in) :: self
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: n
   real(wp), intent(in) :: xyz(3, n)
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: kq

   type(chrg_parameter) :: chrgeq

   call new_charge_model_2019(chrgeq, n, at)
   call add_eeq_hessian(env, n, at, xyz, 0.0_wp, chrgeq, kq, hess)
end subroutine add_charge

end module xtb_modelhessian_internal
