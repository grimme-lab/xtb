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

   public :: TGFFData, init


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

      !> D3 scaling
      real(wp) :: dispscale

      !> Constant data
      real(wp), allocatable :: en(:)
      real(wp), allocatable :: rad(:)
      real(wp), allocatable :: rcov(:)
      integer, allocatable :: metal(:)
      integer, allocatable :: group(:)
      integer, allocatable :: normcn(:)

      !> rep alpha bond
      real(wp), allocatable :: repa (:)
      real(wp), allocatable :: repan(:)

      !> prefactor (Zval), 3atm bond term
      real(wp), allocatable :: repz (:)
      real(wp), allocatable :: zb3atm(:)

      !> HB/XB
      real(wp), allocatable :: xhaci(:)
      real(wp), allocatable :: xhbas(:)
      real(wp), allocatable :: xbaci(:)

      !> EN dep. in EEQ.
      real(wp), allocatable :: chi(:)
      real(wp), allocatable :: gam(:)
      real(wp), allocatable :: cnf(:)
      real(wp), allocatable :: alp(:)

      !> Elem. bond param.
      real(wp), allocatable :: bond(:)

      !> Elem. angular param.
      real(wp), allocatable :: angl(:)

      !> Elem. angular param.
      real(wp), allocatable :: angl2(:)

      !> Elem. torsion param_alloc.
      real(wp), allocatable :: tors(:)

      !> Elem. torsion param.
      real(wp), allocatable :: tors2(:)

      !> BJ radii set in gnff_ini()
      real(wp), allocatable :: d3r0(:)

   end type TGFFData


   !> Initialize a new instance of the parametrisation data
   interface init
      module procedure :: initGFFData
   end interface init


contains


!> Initialize a new instance of the parametrisation data
subroutine initGFFData(self, ndim)

   !> Instance of the parametrisation data
   type(TGFFData), intent(out) :: self

   !> Dimension for allocating space
   integer, intent(in) :: ndim

   allocate(self%en(ndim))
   allocate(self%rad(ndim))
   allocate(self%metal(ndim))
   allocate(self%group(ndim))
   allocate(self%normcn(ndim))
   allocate(self%rcov(ndim))

   allocate(self%repa (ndim))
   allocate(self%repan(ndim))

   allocate(self%repz (ndim))
   allocate(self%zb3atm(ndim))

   allocate(self%xhaci(ndim))
   allocate(self%xhbas(ndim))
   allocate(self%xbaci(ndim))

   allocate(self%chi(ndim))
   allocate(self%gam(ndim))
   allocate(self%cnf(ndim))
   allocate(self%alp(ndim))

   allocate(self%bond(ndim))

   allocate(self%angl(ndim))

   allocate(self%angl2(ndim))

   allocate(self%tors(ndim))

   allocate(self%tors2(ndim))

   allocate(self%d3r0(ndim*(1+ndim)/2))

end subroutine initGFFData


end module xtb_gfnff_data
