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
   use xtb_bmatrix, only : bmat_bond, bmat_angle, bmat_linbend, linbend_frame, &
      & bmat_torsion, bmat_outofplane, oop_angle, bmat_accum_packed, &
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
   use xtb_type_setvar, only : modhess_setvar
   implicit none(type, external)
   private

   !> Base for model Hessians defined by redundant internal coordinates
   type, public, abstract, extends(TModelHessian) :: TInternalModelHessianBase
      private
      !> Bond-stretching force constant
      real(wp) :: kr
      !> Angle-bending force constant
      real(wp) :: kf
      !> Torsional force constant
      real(wp) :: kt
      !> Out-of-plane force constant
      real(wp) :: ko
      !> Distance-dependent scaling factor
      real(wp) :: kd
      !> Charge-dependent force constant
      real(wp) :: kq
      !> Pair distance cutoff
      real(wp) :: rcut
      !> Dispersion scaling factor
      real(wp) :: s6
   contains
      !> Initialize model-Hessian configuration
      procedure, public :: init => init_internal_model_hessian
      !> Compute packed model Hessian
      procedure, public :: compute_packed
      !> Add bond-stretching contributions
      procedure, private :: stretch
      !> Add angle-bending contributions
      procedure, private :: bend
      !> Add torsional contributions
      procedure, private :: torsion
      !> Add out-of-plane contributions
      procedure, private :: outofplane
      !> Add charge-response contribution
      procedure, private :: add_charge
      !> Evaluate pair-distance decay factor
      procedure(pair_factor_interface), deferred, public :: pair_factor
   end type TInternalModelHessianBase

   abstract interface
      !> Pair-distance decay factor for the internal-coordinate force constants
      pure function pair_factor_interface(self, at_i, at_j, r2, dispersion_scale, &
            & outofplane) result(factor)
         import :: TInternalModelHessianBase, wp
         implicit none(type, external)
         !> Model Hessian implementation
         class(TInternalModelHessianBase), intent(in) :: self
         !> Atomic numbers of the pair
         integer, intent(in) :: at_i, at_j
         !> Squared pair distance and dispersion scaling factor
         real(wp), intent(in) :: r2, dispersion_scale
         !> Out-of-plane term indicator
         logical, intent(in) :: outofplane

         real(wp) :: factor
      end function pair_factor_interface
   end interface

contains

!> Copy model-Hessian configuration into an internal-coordinate model
subroutine init_internal_model_hessian(self, modh)
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(inout) :: self
   !> Model Hessian configuration
   type(modhess_setvar), intent(in) :: modh

   self%kr = modh%kr
   self%kf = modh%kf
   self%kt = modh%kt
   self%ko = modh%ko
   self%kd = modh%kd
   self%kq = modh%kq
   self%rcut = modh%rcut
   self%s6 = modh%s6
end subroutine init_internal_model_hessian

!> Compute a packed internal-coordinate model Hessian
subroutine compute_packed(self, env, xyz, n, hess, at)
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed lower-triangle Hessian
   real(wp), intent(out) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)

   real(wp) :: kd
   logical, allocatable :: lcutoff(:, :)

   hess = 0.0_wp
   allocate(lcutoff(n, n), source=.false.)

   kd = self%kd / self%kr
   call self%stretch(xyz, n, hess, at, self%kr, kd, self%s6, lcutoff, self%rcut)
   if (self%kf /= 0.0_wp) then
      call self%bend(xyz, n, hess, at, self%kf, kd, lcutoff)
   end if
   if (self%kt /= 0.0_wp) then
      call self%torsion(xyz, n, hess, at, self%kt, kd, lcutoff)
   end if
   if (self%ko /= 0.0_wp) then
      call self%outofplane(xyz, n, hess, at, self%ko, kd, lcutoff)
   end if
   if (self%kq /= 0.0_wp) then
      call self%add_charge(env, xyz, n, hess, at, self%kq)
   end if
end subroutine compute_packed

!> Add bond-stretching and D2 contributions
pure subroutine stretch(self, xyz, n, hess, at, kr, kd, s6, lcutoff, rcut)
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Bond-stretching force constant
   real(wp), intent(in) :: kr
   !> Distance-dependent scaling factor
   real(wp), intent(in) :: kd
   !> Dispersion scaling factor
   real(wp), intent(in) :: s6
   !> Pair cutoff mask updated in place
   logical, intent(inout) :: lcutoff(n, n)
   !> Distance cutoff
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
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Internal-coordinate force constant and distance-dependent scaling factor
   real(wp), intent(in) :: force_constant, kd
   !> Pair cutoff mask
   logical, intent(in) :: lcutoff(n, n)

   ! Minimum arm and outer-pair length (Bohr).
   real(wp), parameter :: arm_pair_length_tol = 1.0e-10_wp
   ! Linear-bend sine threshold (dimensionless).
   real(wp), parameter :: linear_sine_tol = 1.0e-10_wp
   ! Same-ray cosine proximity threshold (dimensionless).
   real(wp), parameter :: same_ray_cosine_tol = 1.0e-12_wp
   integer :: i, j, m
   real(wp) :: vec_ij(3), vec_mi(3), vec_mj(3), cross_vec(3)
   real(wp) :: rmi2, rmi, rmj2, rmj, rij2, rrij, gij
   real(wp) :: sinphi, cosphi, bmat9(9), evec1(3), evec2(3)

   do m = 1, n
      do i = 1, n
         if (i == m .or. lcutoff(i, m)) cycle
         vec_mi = xyz(:, i) - xyz(:, m)
         rmi2 = dot_product(vec_mi, vec_mi)
         rmi = sqrt(rmi2)
         if (rmi <= arm_pair_length_tol) cycle
         do j = 1, i - 1
            if (j == m) cycle
            if (lcutoff(j, i) .or. lcutoff(j, m)) cycle
            vec_mj = xyz(:, j) - xyz(:, m)
            rmj2 = dot_product(vec_mj, vec_mj)
            rmj = sqrt(rmj2)
            if (rmj <= arm_pair_length_tol) cycle
            vec_ij = xyz(:, j) - xyz(:, i)
            rij2 = dot_product(vec_ij, vec_ij)
            rrij = sqrt(rij2)
            if (rrij <= arm_pair_length_tol) cycle
            cosphi = dot_product(vec_mi, vec_mj) / (rmi*rmj)
            if (abs(cosphi - 1.0_wp) < same_ray_cosine_tol) cycle
            gij = force_constant &
               * self%pair_factor(at(m), at(i), rmi2, 0.5_wp*kd, .false.) &
               * self%pair_factor(at(m), at(j), rmj2, 0.5_wp*kd, .false.)
            cross_vec = crossProd(vec_mi/rmi, vec_mj/rmj)
            sinphi = norm2(cross_vec)
            if (sinphi > linear_sine_tol) then
               bmat9 = bmat_angle(vec_mi, vec_mj)
               call bmat_accum_packed(n, hess, [i, m, j], bmat9, gij)
            else
               ! linear centre: the two Decius bends along a fixed frame
               ! perpendicular to the axis, k B^T B summed over both
               call linbend_frame(vec_mi/rmi, evec1, evec2)
               call bmat_accum_packed(n, hess, [i, m, j], &
                  & bmat_linbend(vec_mi, vec_mj, evec1), gij)
               call bmat_accum_packed(n, hess, [i, m, j], &
                  & bmat_linbend(vec_mi, vec_mj, evec2), gij)
            end if
         end do
      end do
   end do
end subroutine bend

!> Add torsional contributions using one reversal-safe orientation
pure subroutine torsion(self, xyz, n, hess, at, force_constant, kd, lcutoff)
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Internal-coordinate force constant and distance-dependent scaling factor
   real(wp), intent(in) :: force_constant, kd
   !> Pair cutoff mask
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
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Internal-coordinate force constant and distance-dependent scaling factor
   real(wp), intent(in) :: force_constant, kd
   !> Pair cutoff mask
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
   !> Model Hessian implementation
   class(TInternalModelHessianBase), intent(in) :: self
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Packed Hessian updated in place
   real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Charge-dependent force constant
   real(wp), intent(in) :: kq

   type(chrg_parameter) :: chrgeq

   call new_charge_model_2019(chrgeq, n, at)
   call add_eeq_hessian(env, n, at, xyz, 0.0_wp, chrgeq, kq, hess)
end subroutine add_charge

end module xtb_modelhessian_internal
