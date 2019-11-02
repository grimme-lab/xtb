! This file is part of xtb.
!
! Copyright (C) 2017-2019 Stefan Grimme
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

module aoparam
   use iso_fortran_env, wp => real64
   use tbdef_hamiltonian
   implicit none
   private :: wp

   type(tb_hamiltonian) :: gfn

   integer,private,parameter :: max_sh = 10
   type :: tb_parameter
      real(wp) :: en = 1.50_wp
      real(wp) :: mc = 0.0_wp
      real(wp) :: gam = 0.0_wp
      real(wp) :: gam3 = 0.0_wp
      real(wp) :: alp0 = 0.0_wp
      real(wp) :: rad = 0.0_wp
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

   integer,private,parameter :: max_elem = 94

   public

!      data    cnval   / & ! normal CN used in CN dep. AES damping
!     & 1,                                                             1, &
!     & 1, 2,                                           3, 3, 3, 2, 1, 1, &
!     & 1, 2,                                           3, 3, 3, 3, 1, 1, &
!     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
!     & 1, 2, 4,          4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
!     & 1, 2, 4,14*6,     4, 6, 6, 6, 6, 6, 4, 4,  2,   3, 3, 3, 3, 1, 1, &
!     & 8*0 /

contains

subroutine use_parameterset(name,globpar,exist)
   implicit none
   character(len=*),intent(in) :: name
   logical,intent(out)  :: exist
   real(wp),intent(out) :: globpar(25)
   exist = .false.
   select case(name)
   case('.param_gfn.xtb')
      call copy_gfn1_parameterset(globpar)
   case('.param_ipea.xtb')
      call copy_ipea_parameterset(globpar)
   case('.param_gfn2.xtb')
      call copy_gfn2_parameterset(globpar)
   case default
      return
   end select
   exist = .true.
contains
include 'paramcopy_gfn1.inc'
include 'paramcopy_ipea.inc'
include 'paramcopy_gfn2.inc'
end subroutine use_parameterset


end module aoparam
