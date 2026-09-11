! This file is part of xtb.
!
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

!> Common interface for model Hessian implementations
module xtb_modelhessian_type
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   use xtb_type_setvar, only : modhess_setvar
   implicit none(type, external)
   private

   !> Abstract model Hessian implementation
   type, public, abstract :: TModelHessian
   contains
      procedure(model_hessian_stretch), deferred :: stretch
      procedure(model_hessian_mode), deferred :: bend
      procedure(model_hessian_mode), deferred :: torsion
      procedure(model_hessian_mode), deferred :: outofplane
      procedure(model_hessian_charge), deferred :: add_charge
      !> Compute Hessian in packed lower-triangle storage
      procedure :: compute_packed => compute_packed
      !> Compute dense symmetric Hessian
      procedure :: compute_dense => compute_dense
      generic :: compute => compute_packed, compute_dense
   end type TModelHessian

   abstract interface
      !> Add bond-stretching contributions to a packed Hessian
      subroutine model_hessian_stretch(self, xyz, n, hess, at, kr, kd, s6, &
            & lcutoff, rcut)
         import :: TModelHessian, wp
         implicit none(type, external)
         !> Model Hessian implementation
         class(TModelHessian), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: n
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(3, n)
         !> Packed Hessian updated in place
         real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
         !> Atomic numbers
         integer, intent(in) :: at(n)
         !> Stretching force constant
         real(wp), intent(in) :: kr
         !> Distance-dependent scaling factor
         real(wp), intent(in) :: kd
         !> Dispersion scaling factor
         real(wp), intent(in) :: s6
         !> Pair cutoff mask updated in place
         logical, intent(inout) :: lcutoff(n, n)
         !> Distance cutoff
         real(wp), intent(in) :: rcut
      end subroutine model_hessian_stretch

      !> Add an internal-coordinate contribution to a packed Hessian
      subroutine model_hessian_mode(self, xyz, n, hess, at, force_constant, &
            & kd, lcutoff)
         import :: TModelHessian, wp
         implicit none(type, external)
         !> Model Hessian implementation
         class(TModelHessian), intent(in) :: self
         !> Number of atoms
         integer, intent(in) :: n
         !> Cartesian coordinates
         real(wp), intent(in) :: xyz(3, n)
         !> Packed Hessian updated in place
         real(wp), intent(inout) :: hess((3*n)*(3*n + 1)/2)
         !> Atomic numbers
         integer, intent(in) :: at(n)
         !> Internal-coordinate force constant
         real(wp), intent(in) :: force_constant
         !> Distance-dependent scaling factor
         real(wp), intent(in) :: kd
         !> Pair cutoff mask
         logical, intent(in) :: lcutoff(n, n)
      end subroutine model_hessian_mode

      !> Add charge-dependent contributions to a packed Hessian
      subroutine model_hessian_charge(self, env, xyz, n, hess, at, kq)
         import :: TModelHessian, TEnvironment, wp
         implicit none(type, external)
         !> Model Hessian implementation
         class(TModelHessian), intent(in) :: self
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
      end subroutine model_hessian_charge
   end interface

contains

!> Compute Hessian in packed lower-triangle storage
subroutine compute_packed(self, env, xyz, n, hess, at, modh)
   !> Model Hessian implementation
   class(TModelHessian), intent(in) :: self
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
   !> Model Hessian parameters
   type(modhess_setvar), intent(in) :: modh

   real(wp) :: kd
   logical, allocatable :: lcutoff(:, :)

   hess = 0.0_wp
   allocate(lcutoff(n, n), source=.false.)

   kd = modh%kd / modh%kr
   call self%stretch(xyz, n, hess, at, modh%kr, kd, modh%s6, lcutoff, modh%rcut)
   if (modh%kf /= 0.0_wp) then
      call self%bend(xyz, n, hess, at, modh%kf, kd, lcutoff)
   end if
   if (modh%kt /= 0.0_wp) then
      call self%torsion(xyz, n, hess, at, modh%kt, kd, lcutoff)
   end if
   if (modh%ko /= 0.0_wp) then
      call self%outofplane(xyz, n, hess, at, modh%ko, kd, lcutoff)
   end if
   if (modh%kq /= 0.0_wp) then
      call self%add_charge(env, xyz, n, hess, at, modh%kq)
   end if
end subroutine compute_packed

!> Compute dense symmetric Hessian from packed implementation
subroutine compute_dense(self, env, xyz, n, hess, at, modh)
   !> Model Hessian implementation
   class(TModelHessian), intent(in) :: self
   !> Calculation environment
   type(TEnvironment), intent(inout) :: env
   !> Number of atoms
   integer, intent(in) :: n
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(3, n)
   !> Dense symmetric Hessian
   real(wp), intent(out) :: hess(3*n, 3*n)
   !> Atomic numbers
   integer, intent(in) :: at(n)
   !> Model Hessian parameters
   type(modhess_setvar), intent(in) :: modh

   integer :: i, j, ij
   real(wp), allocatable :: packed(:)

   allocate(packed((3*n)*(3*n + 1)/2))
   call self%compute_packed(env, xyz, n, packed, at, modh)

   ij = 0
   do i = 1, 3 * n
      do j = 1, i
         ij = ij + 1
         hess(j, i) = packed(ij)
         hess(i, j) = packed(ij)
      end do
   end do
end subroutine compute_dense

end module xtb_modelhessian_type
