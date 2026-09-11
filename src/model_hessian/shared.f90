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

!> Shared data and utilities for model Hessian implementations
module xtb_modelhessian_shared
   use xtb_mctc_accuracy, only : wp
   use xtb_param_model_hessian, only : d2_damping
   implicit none(type, external)

   private
   public :: rcutoff
   public :: itabrow
   public :: getvdw_hess
   public :: fk_vdw

contains

!> Test whether an atom pair lies beyond a squared-distance cutoff
pure function rcutoff(xyz, katom, latom, rcut)
   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)
   !> Atom indices
   integer, intent(in) :: katom, latom
   !> Squared-distance cutoff
   real(wp), intent(in) :: rcut

   logical :: rcutoff
   real(wp) :: rkl(3), rkl2

   rkl = xyz(:, katom) - xyz(:, latom)
   rkl2 = sum(rkl**2)
   rcutoff = rkl2 > rcut
end function rcutoff

!> Return the legacy periodic-table row used for force-constant scaling
pure elemental function itabrow(i)
   !> Atomic number
   integer, intent(in) :: i

   integer :: itabrow

   ! NOTE: rows 4-7 all return 3 — every element past Ar (Z>18) collapses
   ! to "row 3". Intentional? Periods 4,5,6,7 (K-Rn, Z=19-86) should be
   ! rows 4,5,6,7. Carried over verbatim from legacy src/lindh.f90 iTabRow.
   itabrow = 0
   if (i > 0 .and. i <= 2) then
      itabrow = 1
   else if (i > 2 .and. i <= 10) then
      itabrow = 2
   else if (i > 10) then
      itabrow = 3
   end if
end function itabrow

!> Evaluate a mixed Cartesian derivative of the D2 pair energy
pure elemental subroutine getvdwxy(rx, ry, rz, c66, s6, r0, vdw)
   !> Cartesian displacement, C6 coefficient, scale, and reference distance
   real(wp), intent(in) :: rx, ry, rz, c66, s6, r0
   !> Mixed Cartesian second derivative
   real(wp), intent(out) :: vdw

   real(wp) :: t1, t2, t3, t4, t5, t6, t7, t11, t12, t16, t17, t25, t26, t35
   real(wp) :: t40, t41, t43, t44, t56
   t1 = s6 * c66
   t2 = rx**2
   t3 = ry**2
   t4 = rz**2
   t5 = t2 + t3 + t4
   t6 = t5**2
   t7 = t6**2
   t11 = sqrt(t5)
   t12 = 1.0_wp / r0
   t16 = exp(-d2_damping*(t11*t12 - 1.0_wp))
   t17 = 1.0_wp + t16
   t25 = t17**2
   t26 = 1.0_wp / t25
   t35 = 1.0_wp / t7
   t40 = d2_damping**2
   t41 = r0**2
   t43 = t40 / t41
   t44 = t16**2
   t56 = -48.0_wp * t1 / t7 / t5 / t17 * rx * ry + 13.0_wp * t1 / t11 / &
      t7 * t26 * rx * d2_damping * t12 * ry * t16 - 2.0_wp * t1 * t35 / t25 / &
      t17 * t43 * rx * t44 * ry + t1 * t35 * t26 * t43 * rx * ry * t16
   vdw = t56
end subroutine getvdwxy

!> Evaluate a diagonal Cartesian derivative of the D2 pair energy
pure elemental subroutine getvdwxx(rx, ry, rz, c66, s6, r0, vdw)
   !> Cartesian displacement, C6 coefficient, scale, and reference distance
   real(wp), intent(in) :: rx, ry, rz, c66, s6, r0
   !> Diagonal Cartesian second derivative
   real(wp), intent(out) :: vdw

   real(wp) :: t1, t2, t3, t4, t5, t6, t7, t10, t11, t15, t16, t17, t24, t25, t29
   real(wp) :: t33, t41, t42, t44, t45, t62
   t1 = s6 * c66
   t2 = rx**2
   t3 = ry**2
   t4 = rz**2
   t5 = t2 + t3 + t4
   t6 = t5**2
   t7 = t6**2
   t10 = sqrt(t5)
   t11 = 1.0_wp / r0
   t15 = exp(-d2_damping*(t10*t11 - 1.0_wp))
   t16 = 1.0_wp + t15
   t17 = 1.0_wp / t16
   t24 = t16**2
   t25 = 1.0_wp / t24
   t29 = t11 * t15
   t33 = 1.0_wp / t7
   t41 = d2_damping**2
   t42 = r0**2
   t44 = t41 / t42
   t45 = t15**2
   t62 = -48.0_wp * t1 / t7 / t5 * t17 * t2 + 13.0_wp * t1 / t10 / t7 * &
      t25 * t2 * d2_damping * t29 + 6.0_wp * t1 * t33 * t17 - 2.0_wp * t1 * t33 &
      / t24 / t16 * t44 * t2 * t45 - t1 / t10 / t6 / t5 * t25 * d2_damping * &
      t29 + t1 * t33 * t25 * t44 * t2 * t15
   vdw = t62
end subroutine getvdwxx

!> Build the D2 dispersion Hessian block for an atom pair
pure subroutine getvdw_hess(vec, c66, s6, r0, vdw)
   !> Cartesian pair displacement
   real(wp), intent(in) :: vec(3)
   !> C6 coefficient, scale, and reference distance
   real(wp), intent(in) :: c66, s6, r0
   !> Symmetric Cartesian pair Hessian block
   real(wp), intent(out) :: vdw(3, 3)

   call getvdwxx(vec(1), vec(2), vec(3), c66, s6, r0, vdw(1, 1))
   call getvdwxy(vec(1), vec(2), vec(3), c66, s6, r0, vdw(1, 2))
   call getvdwxy(vec(1), vec(3), vec(2), c66, s6, r0, vdw(1, 3))
   call getvdwxx(vec(2), vec(1), vec(3), c66, s6, r0, vdw(2, 2))
   call getvdwxy(vec(2), vec(3), vec(1), c66, s6, r0, vdw(2, 3))
   call getvdwxx(vec(3), vec(1), vec(2), c66, s6, r0, vdw(3, 3))
   vdw(2, 1) = vdw(1, 2)
   vdw(3, 1) = vdw(1, 3)
   vdw(3, 2) = vdw(2, 3)
end subroutine getvdw_hess

!> Evaluate the van der Waals force-constant decay factor
pure elemental function fk_vdw(alpha, r0, r2) result(gmm)
   !> Decay parameter, reference distance, and squared pair distance
   real(wp), intent(in) :: alpha, r0, r2

   real(wp) :: gmm

   gmm = exp(-alpha*(r0 - sqrt(r2))**2)
end function fk_vdw

end module xtb_modelhessian_shared
