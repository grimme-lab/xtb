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

subroutine set_gfn1_parameter(xpar,globpar)
   use iso_fortran_env, wp => real64
   use tbdef_param
   use tbpar_dftd3
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   real(wp),intent(in) :: globpar(25)
   integer :: i,j
   call setpair(1)
   xpar%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         xpar%kmagic(i,j)=(xpar%kspd(i)+xpar%kspd(j))*0.50_wp
      enddo
   enddo
   if(xpar%kspd(5).gt.0.1)then   ! s - p
      xpar%kmagic(1,2)=xpar%kspd(5)
      xpar%kmagic(2,1)=xpar%kspd(5)
   endif
   xpar%gscal    =globpar(8)*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar(9) !p
   xpar%gam3l(2) =globpar(10)!d-pol
   xpar%gam3l(3) =globpar(11)!d-val
   xpar%kcnsh(1) =globpar(11)*0.01_wp
   xpar%kcnsh(2) =globpar(12)*0.01_wp
   xpar%kcnsh(3) =globpar(13)*0.01_wp
   xpar%kcnsh(4) = 0.005_wp
   xpar%kenscal  =globpar(23) ! kenscal in scf.f
   xpar%xbdamp   =globpar(24)
   xpar%xbrad    =globpar(25)
   xpar%disp%s9  =0.0_wp! d3atm 
   xpar%alphaj   =globpar(18)
   xpar%disp%a1  =globpar(20)
   xpar%disp%a2  =globpar(21)
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar(22)
   xpar%ipshift  =globpar(7)*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

end subroutine set_gfn1_parameter


subroutine set_gfn2_parameter(xpar,globpar,n,at)
   use iso_fortran_env, wp => real64
   use tbdef_param
   use dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   real(wp),intent(in) :: globpar(25)
   integer, intent(in) :: n,at(n)
   integer :: i,j
   call setpair(2)
   xpar%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         xpar%kmagic(i,j)=(xpar%kspd(i)+xpar%kspd(j))*0.50_wp
      enddo
   enddo
   if(xpar%kspd(5).gt.0.1)then   ! s - p
      xpar%kmagic(1,2)=xpar%kspd(5)
      xpar%kmagic(2,1)=xpar%kspd(5)
   endif
   xpar%kmagic(1,3)=xpar%kspd(4)
   xpar%kmagic(3,1)=xpar%kspd(4)
   xpar%kmagic(2,3)=xpar%kspd(6)
   xpar%kmagic(3,2)=xpar%kspd(6)
   xpar%gscal    =globpar(8)*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar(9) !p
   xpar%gam3l(2) =globpar(10)!d-pol
   xpar%gam3l(3) =globpar(11)!d-val
   xpar%cn_shift =globpar(14) ! R AES CN val offset
   xpar%cn_expo  =globpar(15) ! R AES CN steepness
   xpar%cn_rmax  =globpar(16) ! R AES CN Rmax
   xpar%kenscal  =globpar(13) ! kenscal in scf.f
   xpar%disp%s9  =globpar(23) ! d3atm
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   xpar%xbrad    =globpar(24)
   xpar%xbdamp   =globpar(25)
   xpar%alphaj   =globpar(18)
   xpar%disp%a1  =globpar(20)
   xpar%disp%a2  =globpar(21)
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar(22)
   xpar%ipshift  =globpar(7)*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   call d4init(n,at,xpar%g_a,xpar%g_c,p_refq_gfn2xtb,i)
end subroutine set_gfn2_parameter

subroutine set_gfn0_parameter(xpar,globpar,n,at)
   use iso_fortran_env, wp => real64
   use tbdef_param
   !use gfn0_module
   use aoparam
   use dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   real(wp),intent(in) :: globpar(25)
   integer, intent(in) :: n,at(n)
   integer :: i,j

   call setpair(0)
   xpar%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         xpar%kmagic(i,j)=(xpar%kspd(i)+xpar%kspd(j))*0.50_wp
      enddo
   enddo
   if(xpar%kspd(5).gt.0.1)then   ! s - p
      xpar%kmagic(1,2)=xpar%kspd(5)
      xpar%kmagic(2,1)=xpar%kspd(5)
   endif
   xpar%gscal    =globpar(8)*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar(9) !p
   xpar%gam3l(2) =globpar(10)!d-pol
   xpar%gam3l(3) =globpar(11)!d-val
   xpar%kcnsh(1) =globpar(17)         ! K 2s - 2s 
   xpar%kcnsh(2) =globpar( 9)         ! SRB shift
   xpar%kcnsh(3) =globpar(10)         ! SRB prefactor
   xpar%kcnsh(4) =globpar(11)         ! SRB steepnes 
   xpar%gscal    =globpar(12)         ! EN dep  
   xpar%kenscal  =globpar(19)         ! ken² (d3atm in old main.f)
   xpar%xbdamp   =globpar(24)         ! ken⁴
   xpar%xbrad    =globpar(25)         ! rep dEN
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   gam3          =gam3*10.0_wp        ! *10 since gam3 is read with the factor 0.1 in readparam.f
   cxb           =cxb*10.0_wp         ! *10 since cxb  is read with the factor 0.1 in readparam.f
   xpar%alphaj   =globpar(18)
   xpar%disp%a1  =globpar(20)
   xpar%disp%a2  =globpar(21)
   xpar%disp%s6  =1.0_wp
   xpar%disp%s9  =0.0_wp
   xpar%disp%s8  =globpar(22)
   xpar%ipshift  =globpar(7)*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   call d4init(n,at,xpar%g_a,xpar%g_c,p_refq_goedecker,i)
end subroutine set_gfn0_parameter
