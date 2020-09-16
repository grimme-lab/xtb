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

module xtb_grad_core
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert, only : evtoau
   use xtb_xtb_data

   use xtb_mctc_la

!! ========================================================================
!  we use the complete xtb_scc_core, so we inherit all interfaces and all
!  parameters. Therefore, we don't need to declare them again, note that
!  by including the xtb_grad_core you also inherit the xtb_scc_core!
   use xtb_scc_core

   implicit none

   integer, private, parameter :: mmm(20)=(/1,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4/)

contains

!! ========================================================================
!  derivative of S(R) enhancement factor
!! ========================================================================
pure subroutine dshellPoly(iPoly,jPoly,iRad,jRad,rab2,xyz1,xyz2,rf,dxyz)
   !$acc routine seq
   real(wp), intent(in)  :: iPoly,jPoly
   real(wp), intent(in)  :: iRad,jRad
   real(wp), intent(in)  :: rab2
   real(wp), intent(out) :: dxyz(3),rf
   real(wp), intent(in)  :: xyz1(3),xyz2(3)
   real(wp) :: rab,k1,k2,rr,r,a,dum,rf1,rf2,dx,dy,dz
   real(wp) :: t14,t15,t17,t22,t20,t23,t10,t11,t35,t13

   a=0.5            ! R^a dependence 0.5 in GFN1

   dx=xyz1(1)-xyz2(1)
   dy=xyz1(2)-xyz2(2)
   dz=xyz1(3)-xyz2(3)

   rab=sqrt(rab2)

   ! this sloppy conv. factor has been used in development, keep it
   r=iRad+jRad

   rr=rab/r

   k1=iPoly*0.01_wp
   k2=jPoly*0.01_wp

   t14 = rr**a
   t15 = k1*t14
   t17 = 1.0_wp/rab2
   t22 = rr**a
   t23 = k2*t22
   rf=(1.0_wp+t15)*(1.0_wp+k2*t22)
   dxyz(:)=(t15*(1.0_wp+t23)+(1.0_wp+t15)*k2*t22)*a*t17*[dx,dy,dz]

end subroutine dshellPoly

end module xtb_grad_core
