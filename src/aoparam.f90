! This file is part of xtb.
!
! Copyright (C) 2017-2020 Stefan Grimme
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

module xtb_aoparam
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_param, only: pauling_en, covalent_radius_2010, chemical_hardness, &
      &                  covalent_radius_d3
   use xtb_type_param, only : TxTBParameter
   implicit none
   private

   integer, public, parameter :: max_elem = 118
   integer, public, parameter :: max_sh = 10
   type, public :: tb_parameter
      real(wp) :: en = 1.50_wp
      real(wp) :: mc = 0.0_wp
      real(wp) :: rad = 0.0_wp
      real(wp) :: gam = 0.0_wp
      real(wp) :: gam3 = 0.0_wp
      real(wp) :: alp0 = 0.0_wp
      real(wp) :: wll(max_sh) = 0.0_wp
      real(wp) :: rep(2) = 0.0_wp
      real(wp) :: polyr(4) = 0.0_wp
      real(wp) :: cxb = 0.0_wp
      real(wp) :: ao_exp(max_sh) = 0.0_wp
      real(wp) :: ao_lev(max_sh) = 0.0_wp
      real(wp) :: lpar(0:2) = 0.0_wp
      real(wp) :: kcnat(0:2) = 0.0_wp
      real(wp) :: kqat(3) = 0.0_wp
      real(wp) :: radaes = 5.0_wp
      real(wp) :: dpolc = 0.0_wp
      real(wp) :: qpolc = 0.0_wp
      integer  :: ao_pqn(max_sh) = 0
      integer  :: ao_l(max_sh) = 0
      integer  :: ao_n = 0
      integer  :: ao_typ(max_sh) = 0
      integer  :: metal = 0
      integer  :: cnval = 0
      character(len=30) :: timestp='------------------------------'
   end type tb_parameter

   character, parameter :: flag = '$'
   character, parameter :: space = ' '
   character, parameter :: equal = '='
   character, parameter :: hash = '#'
   character, parameter :: dot = '.'
   character(len=*), parameter :: flag_end = '$end'

   real(wp), public :: en(max_elem) = pauling_en
   real(wp), public :: mc(max_elem)
   real(wp), public :: rad(max_elem) = covalent_radius_2010 * autoaa
   real(wp), public :: gam(max_elem) = chemical_hardness
   real(wp), public :: gam3(max_elem)
   real(wp), public :: alp0(max_elem)
   real(wp), public :: wll(max_elem,10)
   real(wp), public :: rep(2,max_elem)
   real(wp), public :: polyr(4,max_elem)
   real(wp), public :: cxb(max_elem)
   real(wp), public :: eeqkcn(max_elem)
   real(wp), public :: kqat2(max_elem)
   real(wp), public :: ao_exp(10,max_elem)
   real(wp), public :: ao_lev(10,max_elem)
   real(wp), public :: lpar(0:2,max_elem)
   real(wp), public :: kpair(max_elem,max_elem)
   real(wp), public :: kcnat(0:2,max_elem)
   real(wp), public :: kqat(3,max_elem)
   real(wp), public :: radaes(max_elem) = covalent_radius_d3
   real(wp), public :: eeqEN(max_elem)
   real(wp), public :: dpolc(max_elem)
   real(wp), public :: qpolc(max_elem)
   integer, public  :: ao_pqn(10,max_elem)
   integer, public  :: ao_l(10,max_elem)
   integer, public  :: ao_n(max_elem)
   integer, public  :: ao_typ(10,max_elem)
   integer, public  :: metal(max_elem)
   integer, public  :: cnval(max_elem)
   character(len=30), public :: timestp(max_elem)

   data metal / &
   &  0,                                                 0, &! H-He
   &  1, 1,                               1, 0, 0, 0, 0, 0, &! Li-Ne
   &  1, 1,                               1, 0, 0, 0, 0, 0, &! Na-Ar
   &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! K-Kr
   &  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Rb-Xe
   &  1, 1, &! Cs/Ba
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &!La-Lu
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, &! Lu-Rn
   &  1, 1, &
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, &!Fr-
   &        1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0 /! -Og

      data    cnval   / & ! normal CN used in CN dep. AES damping
     & 1,                                                             1, &
     & 1, 2,                                           3, 3, 3, 2, 1, 1, &
     & 1, 2,                                           3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 1, 2, 4,14*6,     4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
     & 2*0, 14*0, 10*0, 6*0 /

contains


end module xtb_aoparam
