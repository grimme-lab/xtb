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
   ! TODO:
   ! Could cut modelhessian module down further by using bmatrix and internal
   ! coordinate infra more extensively. Coordinate traversal is different from
   ! generic redundant internals and model specific, so needs special constructors.
   ! Coordinates may then be assigned model-specific force constants,
   ! and passed to xtb_bmatrix for shared Wilson B-row construction and
   ! accumulation of k B^T B. Pairwise Cartesian terms such as dispersion and
   ! charge contributions remain separate.
   use xtb_mctc_accuracy, only : wp
   use xtb_type_environment, only : TEnvironment
   implicit none(type, external)
   private

   !> Abstract model Hessian implementation
   type, public, abstract :: TModelHessian
   contains
      !> Compute Hessian in packed lower-triangle storage
      procedure(model_hessian_packed), deferred, public :: compute_packed
      !> Compute dense symmetric Hessian
      procedure, private :: compute_dense
      generic, public :: compute => compute_packed, compute_dense
   end type TModelHessian

   abstract interface
      !> Compute a packed lower-triangle Hessian
      subroutine model_hessian_packed(self, env, xyz, n, hess, at)
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
         !> Packed lower-triangle Hessian
         real(wp), intent(out) :: hess((3*n)*(3*n + 1)/2)
         !> Atomic numbers
         integer, intent(in) :: at(n)
      end subroutine model_hessian_packed
   end interface

contains

!> Compute dense symmetric Hessian from packed implementation
subroutine compute_dense(self, env, xyz, n, hess, at)
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

   integer :: i, j, ij
   real(wp), allocatable :: packed(:)

   allocate(packed((3*n)*(3*n + 1)/2))
   call self%compute_packed(env, xyz, n, packed, at)

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
