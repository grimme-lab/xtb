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

!> Static parameter data used by model Hessian implementations
module xtb_param_model_hessian
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : aatoau
   implicit none(type, external)
   private

   public :: TLindhParameters, lindh_d2_parameters, lindh_parameters
   public :: legacy_rav, legacy_aav
   public :: d2_c6, d2_vander, d2_damping
   public :: gff_c6, gff_stretch_constant, gff_bend_constant
   public :: gff_torsion_constant, gff_distance_threshold

   !> Distance-decay parameters for a Lindh model variant
   type :: TLindhParameters
      real(wp) :: rav(3, 3)
      real(wp) :: aav(3, 3)
      real(wp) :: dav(3, 3)
      real(wp) :: outofplane_dispersion_scale
   end type TLindhParameters

   real(wp), parameter :: legacy_rav(3, 3) = reshape([ &
      1.3500_wp, 2.1000_wp, 2.5300_wp, &
      2.1000_wp, 2.8700_wp, 3.4000_wp, &
      2.5300_wp, 3.4000_wp, 3.4000_wp], [3, 3])
   real(wp), parameter :: legacy_aav(3, 3) = reshape([ &
      1.0000_wp, 0.3949_wp, 0.3949_wp, &
      0.3949_wp, 0.2800_wp, 0.2800_wp, &
      0.3949_wp, 0.2800_wp, 0.2800_wp], [3, 3])

   real(wp), parameter :: d2_dav(3, 3) = reshape([ &
      0.0000_wp, 0.0000_wp, 0.0000_wp, &
      0.0000_wp, 0.0000_wp, 0.0000_wp, &
      0.0000_wp, 0.0000_wp, 0.0000_wp], [3, 3])

   real(wp), parameter :: lindh_rav(3, 3) = reshape([ &
      1.3500_wp, 2.1000_wp, 2.5300_wp, &
      2.1000_wp, 2.8700_wp, 3.8000_wp, &
      2.5300_wp, 3.8000_wp, 4.5000_wp], [3, 3])
   real(wp), parameter :: lindh_aav(3, 3) = reshape([ &
      1.0000_wp, 0.3949_wp, 0.3949_wp, &
      0.3949_wp, 0.2800_wp, 0.1200_wp, &
      0.3949_wp, 0.1200_wp, 0.0600_wp], [3, 3])
   real(wp), parameter :: lindh_dav(3, 3) = reshape([ &
      0.0000_wp, 3.6000_wp, 3.6000_wp, &
      3.6000_wp, 5.3000_wp, 5.3000_wp, &
      3.6000_wp, 5.3000_wp, 5.3000_wp], [3, 3])

   type(TLindhParameters), parameter :: lindh_d2_parameters = TLindhParameters( &
      rav=legacy_rav, aav=legacy_aav, dav=d2_dav, &
      outofplane_dispersion_scale=1.0_wp)

   type(TLindhParameters), parameter :: lindh_parameters = TLindhParameters( &
      rav=lindh_rav, aav=lindh_aav, dav=lindh_dav, &
      outofplane_dispersion_scale=0.0_wp)

   real(wp), parameter :: d2_damping = 20.0_wp
   real(wp), parameter :: d2_vander(86) = aatoau * [ &
      0.91_wp, 0.92_wp, &
      0.75_wp, 1.28_wp, 1.35_wp, 1.32_wp, 1.27_wp, 1.22_wp, 1.17_wp, 1.13_wp, &
      1.04_wp, 1.24_wp, 1.49_wp, 1.56_wp, 1.55_wp, 1.53_wp, 1.49_wp, 1.45_wp, &
      1.35_wp, 1.34_wp, &
      1.42_wp, 1.42_wp, 1.42_wp, 1.42_wp, 1.42_wp, &
      1.42_wp, 1.42_wp, 1.42_wp, 1.42_wp, 1.42_wp, &
      1.50_wp, 1.57_wp, 1.60_wp, 1.61_wp, 1.59_wp, 1.57_wp, &
      1.48_wp, 1.46_wp, &
      1.49_wp, 1.49_wp, 1.49_wp, 1.49_wp, 1.49_wp, &
      1.49_wp, 1.49_wp, 1.49_wp, 1.49_wp, 1.49_wp, &
      1.52_wp, 1.64_wp, 1.71_wp, 1.72_wp, 1.72_wp, 1.71_wp, &
      2.00_wp, 2.00_wp, &
      2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, &
      2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, &
      2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, &
      2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, &
      2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp, 2.00_wp]

   real(wp), parameter :: d2_c6(86) = [ &
      0.14_wp, 0.08_wp, &
      1.61_wp, 1.61_wp, 3.13_wp, 1.75_wp, 1.23_wp, 0.70_wp, 0.75_wp, 0.63_wp, &
      5.71_wp, 5.71_wp, 10.79_wp, 9.23_wp, 7.84_wp, 5.57_wp, 5.07_wp, 4.61_wp, &
      10.80_wp, 10.80_wp, &
      10.80_wp, 10.80_wp, 10.80_wp, 10.80_wp, 10.80_wp, &
      10.80_wp, 10.80_wp, 10.80_wp, 10.80_wp, 10.80_wp, &
      16.99_wp, 17.10_wp, 16.37_wp, 12.64_wp, 12.47_wp, 12.01_wp, &
      24.67_wp, 24.67_wp, &
      24.67_wp, 24.67_wp, 24.67_wp, 24.67_wp, 24.67_wp, &
      24.67_wp, 24.67_wp, 24.67_wp, 24.67_wp, 24.67_wp, &
      37.32_wp, 38.71_wp, 38.44_wp, 31.74_wp, 31.50_wp, 29.99_wp, &
      50.00_wp, 50.00_wp, &
      50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, &
      50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, &
      50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, &
      50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, &
      50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp, 50.00_wp]

   real(wp), parameter :: gff_stretch_constant = 0.4500_wp
   real(wp), parameter :: gff_bend_constant = 0.3000_wp
   real(wp), parameter :: gff_torsion_constant = 0.7500_wp
   real(wp), parameter :: gff_distance_threshold = 1.0e-10_wp
   real(wp), parameter :: gff_c6(55) = [d2_c6(:54), 40.00_wp]

end module xtb_param_model_hessian
