! This file is part of ancopt.
! SPDX-Identifier: LGPL-3.0-or-later
!
! ancopt is free software: you can redistribute it and/or modify it under
! the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! ancopt is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public License
! along with ancopt.  If not, see <https://www.gnu.org/licenses/>.

!> Proxy module for creating optimizer objects
module xtb_pbc_optimizer_optimizer
   !use mctc_env, only : error_type, fatal_error
   !use tblite_context_type, only : context_type
   use xtb_type_environment, only : TEnvironment
   !use ancopt_filter_type, only : filter_type
   use  xtb_pbc_optimizer_filter_type, only : filter_type
   !use ancopt_optimizer_bfgs, only : bfgs_optimizer, bfgs_input, new_bfgs_optimizer
   !use ancopt_optimizer_fire, only : fire_optimizer, fire_input, new_fire_optimizer
   use xtb_pbc_optimizer_lbfgs, only : optimizer_type, lbfgs_optimizer, lbfgs_input, new_lbfgs_optimizer
   !use ancopt_optimizer_rf, only : rf_optimizer, rf_input, new_rf_optimizer
   !use xtb_pbc_optimizer_driver, only : optimizer_type
   implicit none
   private

   public :: optimizer_type, optimizer_input, new_optimizer


   type :: optimizer_input
      character(len=:), allocatable :: method
      !type(bfgs_input) :: bfgs
      !type(fire_input) :: fire
      type(lbfgs_input) :: lbfgs
      !type(rf_input) :: rf
   end type


contains

subroutine new_optimizer(self, env, input, filter)
  !> Instance of the optimization driver
  class(optimizer_type), allocatable, intent(out) :: self
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
   !> Input for constructing the optimizer
   type(optimizer_input), intent(in) :: input
   !> Transformation filter
   class(filter_type), intent(in) :: filter

   !type(error_type), allocatable :: error

!   if (.not.allocated(input%method)) then
      call new_lbfgs_optimizer_(self, env, input%lbfgs, filter)
!      return
!   end if

   ! select case(input%method)
   ! case default
   !    call fatal_error(error, "Unknown optimizer '"//input%method//"' requested")
   !    call ctx%set_error(error)
   ! case("bfgs")
   !    call new_bfgs_optimizer_(self, ctx, input%bfgs, filter)
   ! case("fire")
   !    call new_fire_optimizer_(self, ctx, input%fire, filter)
   ! case("lbfgs")
   !    call new_lbfgs_optimizer_(self, ctx, input%lbfgs, filter)
   ! case("rf")
   !    call new_rf_optimizer_(self, ctx, input%rf, filter)
   ! end select
end subroutine new_optimizer

subroutine new_lbfgs_optimizer_(self, env, input, filter)
   !> Instance of the optimizer
   class(optimizer_type), allocatable, intent(out) :: self
  !> Calculation environment
  type(TEnvironment), intent(inout) :: env
   !> Input for constructing the optimizer
   type(lbfgs_input), intent(in) :: input
   !> Transformation filter
   class(filter_type), intent(in) :: filter

   type(lbfgs_optimizer), allocatable :: new

   allocate(new)
   call new_lbfgs_optimizer(new, env, input, filter)
   call move_alloc(new, self)
end subroutine new_lbfgs_optimizer_

!subroutine new_bfgs_optimizer_(self, ctx, input, filter)
!   !> Instance of the optimizer
!   class(optimizer_type), allocatable, intent(out) :: self
!   !> Calculation environment context
!   type(context_type), intent(inout) :: ctx
!   !> Input for constructing the optimizer
!   type(bfgs_input), intent(in) :: input
!   !> Transformation filter
!   class(filter_type), intent(in) :: filter

!   type(bfgs_optimizer), allocatable :: new

!   allocate(new)
!   call new_bfgs_optimizer(new, ctx, input, filter)
!   call move_alloc(new, self)
!end subroutine new_bfgs_optimizer_

!subroutine new_fire_optimizer_(self, ctx, input, filter)
!   !> Instance of the optimizer
!   class(optimizer_type), allocatable, intent(out) :: self
!   !> Calculation environment context
!   type(context_type), intent(inout) :: ctx
!   !> Input for constructing the optimizer
!   type(fire_input), intent(in) :: input
!   !> Transformation filter
!   class(filter_type), intent(in) :: filter

!   type(fire_optimizer), allocatable :: new

!   allocate(new)
!   call new_fire_optimizer(new, ctx, input, filter)
!   call move_alloc(new, self)
!end subroutine new_fire_optimizer_

!subroutine new_rf_optimizer_(self, ctx, input, filter)
!   !> Instance of the optimizer
!   class(optimizer_type), allocatable, intent(out) :: self
!   !> Calculation environment context
!   type(context_type), intent(inout) :: ctx
!   !> Input for constructing the optimizer
!   type(rf_input), intent(in) :: input
!   !> Transformation filter
!   class(filter_type), intent(in) :: filter

!   type(rf_optimizer), allocatable :: new

!   allocate(new)
!   call new_rf_optimizer(new, ctx, input, filter)
!   call move_alloc(new, self)
!end subroutine new_rf_optimizer_

end module xtb_pbc_optimizer_optimizer
