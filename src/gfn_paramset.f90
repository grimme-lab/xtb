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

module xtb_paramset
contains

subroutine set_gfn1_parameter(xpar,globpar)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_param
   use xtb_disp_dftd3param
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   integer :: i,j
   call setpair(1)
   xpar%kspd(1:6)=[globpar%ks, globpar%kp, globpar%kd, globpar%kf, &
      & globpar%kdiffa, globpar%kdiffb]
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
   xpar%gscal    =globpar%gscal*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar%zcnf !p
   xpar%gam3l(2) =globpar%tscal!d-pol
   xpar%gam3l(3) =globpar%kcn!d-val
   xpar%kcnsh(1) =globpar%kcn*0.01_wp
   xpar%kcnsh(2) =globpar%fpol*0.01_wp
   xpar%kcnsh(3) =globpar%ken*0.01_wp
   xpar%kcnsh(4) = 0.005_wp
   xpar%kenscal  =globpar%dispatm ! kenscal in scf.f
   xpar%xbdamp   =globpar%xbdamp
   xpar%xbrad    =globpar%xbrad
   xpar%disp%s9  =0.0_wp! d3atm
   xpar%alphaj   =globpar%alphaj
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%ipshift  =globpar%wllscal*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

end subroutine set_gfn1_parameter


subroutine set_gfn2_parameter(xpar,globpar)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_param
   use xtb_disp_dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   integer :: i,j
   call setpair(2)
   xpar%kspd(1:6)=[globpar%ks, globpar%kp, globpar%kd, globpar%kf, &
      & globpar%kdiffa, globpar%kdiffb]
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
   xpar%gscal    =globpar%gscal*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar%zcnf !p
   xpar%gam3l(2) =globpar%tscal!d-pol
   xpar%gam3l(3) =globpar%kcn!d-val
   xpar%cn_shift =globpar%lshift ! R AES CN val offset
   xpar%cn_expo  =globpar%lshifta ! R AES CN steepness
   xpar%cn_rmax  =globpar%split ! R AES CN Rmax
   xpar%kenscal  =globpar%ken ! kenscal in scf.f
   xpar%disp%s9  =globpar%dispatm ! d3atm
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   xpar%xbrad    =globpar%xbdamp
   xpar%xbdamp   =globpar%xbrad
   xpar%alphaj   =globpar%alphaj
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%ipshift  =globpar%wllscal*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   call d4init(xpar%g_a,xpar%g_c,p_refq_gfn2xtb)
end subroutine set_gfn2_parameter

subroutine set_gfn0_parameter(xpar,globpar)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_param
   !use gfn0_module
   use xtb_aoparam
   use xtb_disp_dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   integer :: i,j

   call setpair(0)
   xpar%kspd(1:6)=[globpar%ks, globpar%kp, globpar%kd, globpar%kf, &
      & globpar%kdiffa, globpar%kdiffb]
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
   xpar%gscal    =globpar%gscal*0.1      ! purpose changed in scf.f !
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar%zcnf !p
   xpar%gam3l(2) =globpar%tscal!d-pol
   xpar%gam3l(3) =globpar%kcn!d-val
   xpar%kcnsh(1) =globpar%zqf         ! K 2s - 2s
   xpar%kcnsh(2) =globpar%zcnf         ! SRB shift
   xpar%kcnsh(3) =globpar%tscal         ! SRB prefactor
   xpar%kcnsh(4) =globpar%kcn         ! SRB steepnes
   xpar%gscal    =globpar%fpol         ! EN dep
   xpar%kenscal  =globpar%kexpo         ! ken² (d3atm in old main.f)
   xpar%xbdamp   =globpar%xbdamp         ! ken⁴
   xpar%xbrad    =globpar%xbrad         ! rep dEN
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   xpar%alphaj   =globpar%alphaj
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s9  =0.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%ipshift  =globpar%wllscal*0.1
   xpar%eashift  =xpar%ipshift
   xpar%ken1     =1.0
   xpar%zqf      =0
   xpar%zcnf     =0
   xpar%fpol     =0
   xpar%wllscal  =1
   call d4init(xpar%g_a,xpar%g_c,p_refq_goedecker)
end subroutine set_gfn0_parameter

end module xtb_paramset
