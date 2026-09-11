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

!> Swart model Hessian augmented with D2 dispersion
module xtb_modelhessian_swart
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_param, only : covalent_radius_2009
   use xtb_modelhessian_internal, only : TInternalModelHessianBase
   use xtb_modelhessian_shared, only : fk_vdw
   implicit none(type, external)
   private

   type, public, extends(TInternalModelHessianBase) :: TSwartModelHessian
   contains
      procedure :: pair_factor
   end type TSwartModelHessian

contains

!> Evaluate the Swart pair-distance decay factor
pure function pair_factor(self, at_i, at_j, r2, dispersion_scale, outofplane) &
      & result(factor)
   class(TSwartModelHessian), intent(in) :: self
   integer, intent(in) :: at_i, at_j
   real(wp), intent(in) :: r2, dispersion_scale
   logical, intent(in) :: outofplane

   real(wp) :: factor

   real(wp) :: r0, d0

   r0 = covalent_radius_2009(at_i) + covalent_radius_2009(at_j)
   d0 = covalent_radius_2009(at_i) + covalent_radius_2009(at_j)
   factor = exp(-(sqrt(r2)/r0 - 1.0_wp)) &
      + dispersion_scale * fk_vdw(5.0_wp, d0, r2)
end function pair_factor

end module xtb_modelhessian_swart
