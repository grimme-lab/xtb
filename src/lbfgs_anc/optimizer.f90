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
   use xtb_type_environment, only : TEnvironment
   use  xtb_pbc_optimizer_filter_type, only : filter_type
   use xtb_pbc_optimizer_lbfgs, only : optimizer_type, lbfgs_optimizer, lbfgs_input, new_lbfgs_optimizer
   implicit none
   private

   public :: optimizer_type, optimizer_input, new_optimizer


   type :: optimizer_input
      character(len=:), allocatable :: method
      type(lbfgs_input) :: lbfgs
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

   call new_lbfgs_optimizer_(self, env, input%lbfgs, filter)

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


end module xtb_pbc_optimizer_optimizer
