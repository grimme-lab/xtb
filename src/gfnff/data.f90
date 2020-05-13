! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Type holding the parametrisation data of the force field
module xtb_gfnff_data
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: TGFFData


   !> Parametrisation data for the force field
   type :: TGFFData

      !> repulsion scaling
      real(wp) :: repscaln
      real(wp) :: repscalb

      !> bend/tors angle damping
      real(wp) :: atcuta
      real(wp) :: atcutt

      !> bend/tors nci angle damping for HB term
      real(wp) :: atcuta_nci
      real(wp) :: atcutt_nci

      !> damping HB
      real(wp) :: hbacut
      real(wp) :: hbscut

      !> damping XB
      real(wp) :: xbacut
      real(wp) :: xbscut

      !> damping HB/XB
      real(wp) :: hbalp

      !> damping HB/XB
      real(wp) :: hblongcut
      real(wp) :: hblongcut_xb

      !> charge scaling HB/XB
      real(wp) :: hbst
      real(wp) :: hbsf
      real(wp) :: xbst
      real(wp) :: xbsf

      !> HB AH-B
      real(wp) :: xhaci_globabh

      !> HB AH-O=C
      real(wp) :: xhaci_coh

      !> acidity
      real(wp) :: xhaci_glob

      !> HB AH-B
      real(wp) :: hbabmix

      !> new parameter for neighbour angle
      real(wp) :: hbnbcut

      !> new parameter for HB NCI angle term
      real(wp) :: tors_hb

      !> new parameter for HB NCI torsion term
      real(wp) :: bend_hb

      !> new parameter for FC scaling of bonds in HB
      real(wp) :: vbond_scale

      !> max CN cut-off
      real(wp) :: cnmax

   end type TGFFData


contains


end module xtb_gfnff_data
