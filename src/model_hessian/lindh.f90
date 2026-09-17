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

!> Lindh model Hessians with the 1995 and 2007 parameterizations
module xtb_modelhessian_lindh
   use xtb_mctc_accuracy, only : wp
   use xtb_modelhessian_internal, only : TInternalModelHessianBase
   use xtb_param_model_hessian, only : TLindhParameters, &
      & lindh_d2_parameters, lindh_parameters
   use xtb_modelhessian_shared, only : itabrow, fk_vdw
   implicit none(type, external)
   private

   type, abstract, extends(TInternalModelHessianBase) :: TLindhModelHessianBase
   contains
      procedure :: pair_factor
   end type TLindhModelHessianBase

   type, public, extends(TLindhModelHessianBase) :: TLindhModelHessian
   end type TLindhModelHessian

   type, public, extends(TLindhModelHessianBase) :: TLindhD2ModelHessian
   end type TLindhD2ModelHessian

contains

!> Evaluate the selected Lindh pair-distance decay factor
pure function pair_factor(self, at_i, at_j, r2, dispersion_scale, outofplane) &
      & result(factor)
   class(TLindhModelHessianBase), intent(in) :: self
   integer, intent(in) :: at_i, at_j
   real(wp), intent(in) :: r2, dispersion_scale
   logical, intent(in) :: outofplane

   real(wp) :: factor

   integer :: ir, jr
   type(TLindhParameters) :: parameters

   select type (self)
   class is (TLindhD2ModelHessian)
      parameters = lindh_d2_parameters
   class default
      parameters = lindh_parameters
   end select
   ir = itabrow(at_i)
   jr = itabrow(at_j)
   factor = exp(parameters%aav(ir, jr)*(parameters%rav(ir, jr)**2 - r2))
   if (outofplane) then
      factor = factor + dispersion_scale * parameters%outofplane_dispersion_scale &
         * fk_vdw(4.0_wp, parameters%dav(ir, jr), r2)
   else
      factor = factor + dispersion_scale &
         * fk_vdw(4.0_wp, parameters%dav(ir, jr), r2)
   end if
end function pair_factor

end module xtb_modelhessian_lindh
