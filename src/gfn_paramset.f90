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
   use xtb_mctc_accuracy, only : wp
   use xtb_xtb_data
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2
   use xtb_type_param

contains

subroutine set_gfn1_parameter(xpar,globpar,xtbData)
   use xtb_disp_dftd3param
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   type(TxTBData), intent(inout) :: xtbData
   integer :: i,j
   xpar%kcnsh(1) =globpar%kcn*0.01_wp
   xpar%kcnsh(2) =globpar%fpol*0.01_wp
   xpar%kcnsh(3) =globpar%ken*0.01_wp
   xpar%kcnsh(4) = 0.005_wp
   xpar%kenscal  =globpar%dispatm ! kenscal in scf.f
   xpar%disp%s9  =0.0_wp! d3atm
   xpar%alphaj   =globpar%alphaj
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%ipshift  =globpar%ipeashift*0.1
   xpar%eashift  =xpar%ipshift
   if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

end subroutine set_gfn1_parameter


subroutine set_gfn2_parameter(xpar,globpar,xtbData)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_param
   use xtb_disp_dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   type(TxTBData), intent(inout) :: xtbData
   integer :: i,j
   xpar%gam3l(0) =1.00_wp     !s
   xpar%gam3l(1) =globpar%zcnf !p
   xpar%gam3l(2) =globpar%tscal!d-pol
   xpar%gam3l(3) =globpar%kcn!d-val
   xpar%kenscal  =globpar%ken ! kenscal in scf.f
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%disp%s9  =globpar%dispatm ! d3atm
   xpar%ipshift  =globpar%ipeashift*0.1
   xpar%eashift  =xpar%ipshift
   call d4init(xpar%g_a,xpar%g_c,p_refq_gfn2xtb)
end subroutine set_gfn2_parameter

subroutine set_gfn0_parameter(xpar,globpar,xtbData)
   use xtb_mctc_accuracy, only : wp
   use xtb_type_param
   !use gfn0_module
   use xtb_aoparam
   use xtb_disp_dftd4
   implicit none
   type(scc_parameter),intent(inout) :: xpar
   type(TxTBParameter), intent(in) :: globpar
   type(TxTBData), intent(inout) :: xtbData
   integer :: i,j

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
   xpar%g_a      =3.0_wp
   xpar%g_c      =2.0_wp
   xpar%wf       =6.0_wp
   xpar%alphaj   =globpar%alphaj
   xpar%disp%a1  =globpar%dispa
   xpar%disp%a2  =globpar%dispb
   xpar%disp%s6  =1.0_wp
   xpar%disp%s9  =0.0_wp
   xpar%disp%s8  =globpar%dispc
   xpar%ipshift  =globpar%ipeashift*0.1
   xpar%eashift  =xpar%ipshift
   call d4init(xpar%g_a,xpar%g_c,p_refq_goedecker)
end subroutine set_gfn0_parameter

subroutine use_parameterset(name,globpar,xtbData,exist)
   use xtb_aoparam, only : kpair
   implicit none
   character(len=*),intent(in) :: name
   logical,intent(out)  :: exist
   type(TxTBParameter), intent(out) :: globpar
   type(TxTBData), intent(out) :: xtbData
   exist = .false.
   select case(name)
   case('.param_gfn0.xtb')
      return
      globpar = gfn0Globals
      call initGFN0(xtbData)
   case('.param_gfn.xtb')
      globpar = gfn1Globals
      call initGFN1(xtbData)
   case('.param_gfn2.xtb')
      globpar = gfn2Globals
      call initGFN2(xtbData)
   case default
      return
   end select
   exist = .true.
end subroutine use_parameterset

! global, predefined pair parameters
subroutine setpair(gfn_method, pairParam)
   real(wp), intent(inout) :: pairParam(:, :)
   integer gfn_method
   integer i,j,ii,jj
   real*8  kp(3)
   real*8  kparam
   integer tmgroup(3)
   logical notset

   if(gfn_method.eq.1)then
      kp(1)=1.1    ! 3d
      kp(2)=1.2    ! 4d
      kp(3)=1.2    ! 5d or 4f
   elseif(gfn_method.eq.0)then
      kp(1)=1.10
      kp(2)=1.10
      kp(3)=1.10
      kparam=0.9
      tmgroup=(/29,47,79/)
      do i=1,3
         do j=1,3
            ii=tmgroup(i)
            jj=tmgroup(j)
            pairParam(ii,jj)=kparam
            pairParam(jj,ii)=kparam
         enddo
      enddo
   elseif(gfn_method.gt.1)then
      kp(1)=1.   ! 3d
      kp(2)=1.   ! 4d
      kp(3)=1.   ! 5d or 4f
      !     write(*,'(''KAB for pair M(3d)-M(3d) :'',f8.4)')kp(1)
      !     write(*,'(''KAB for pair M(4d)-M(4d) :'',f8.4)')kp(2)
      !     write(*,'(''KAB for pair M(5d)-M(5d) :'',f8.4)')kp(3)
   endif

   do i=21,79
      do j=21,i
         ii=tmmetal(i)
         jj=tmmetal(j)
         !           metal-metal interaction
         notset=abs(pairParam(i,j)-1.0d0).lt.1.d-6 .and. &
            &             abs(pairParam(j,i)-1.0d0).lt.1.d-6
         if(ii.gt.0.and.jj.gt.0.and.notset) then
            pairParam(i,j)=0.5*(kp(ii)+kp(jj))
            pairParam(j,i)=0.5*(kp(ii)+kp(jj))
         endif
      enddo
   enddo

end subroutine setpair

integer function tmmetal(i)
   integer i,j

   j=0
   if(i.gt.20.and.i.lt.30) j=1
   if(i.gt.38.and.i.lt.48) j=2
   if(i.gt.56.and.i.lt.80) j=3

   tmmetal=j

end function tmmetal

end module xtb_paramset
