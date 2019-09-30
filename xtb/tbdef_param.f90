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

module tbdef_param
   use iso_fortran_env, wp => real64
   implicit none

   public :: dftd_parameter

   public ::  scc_parameter
   public :: gfn0_param_from_globpar
   public :: gfn1_param_from_globpar
   public :: gfn2_param_from_globpar

   public :: chrg_parameter

   public :: gfn_parameter
   private

   type :: dftd_parameter
      real(wp) :: s6  = -1.0_wp
      real(wp) :: s8  = -1.0_wp
      real(wp) :: s10 =  0.0_wp
      real(wp) :: a1  = -1.0_wp
      real(wp) :: a2  = -1.0_wp
      real(wp) :: s9  =  1.0_wp
      integer  :: alp = 16
      ! for MBD@rsSCS
      real(wp) :: beta = 1.0_wp
   end type dftd_parameter

   type :: scc_parameter
      real(wp) :: kspd(6)
      real(wp) :: gscal
      real(wp) :: gam3l(0:3)
      real(wp) :: kcnsh(4)
      real(wp) :: kmagic(4,4)
      type(dftd_parameter) :: disp
      real(wp) :: cn_shift
      real(wp) :: cn_expo
      real(wp) :: cn_rmax
      real(wp) :: xbdamp
      real(wp) :: xbrad
      real(wp) :: ipshift
      real(wp) :: eashift
      real(wp) :: ken1
      real(wp) :: zqf
      real(wp) :: zcnf
      real(wp) :: fpol
      real(wp) :: wllscal
      real(wp) :: tfac
      real(wp) :: kcn
      real(wp) :: alphaj
      real(wp) :: lshift
      real(wp) :: lshifta
      real(wp) :: split
      real(wp) :: kenscal
      real(wp) :: g_a
      real(wp) :: g_c
      real(wp) :: wf
   end type scc_parameter

   integer, private, parameter :: max_elem  = 94
   integer, private, parameter :: max_shell = 10
   type gfn_parameter(nelem,nshell)
      integer, len :: nelem  = max_elem
      integer, len :: nshell = max_shell
      integer  :: method = -1
      real(wp) :: en(nelem)
      real(wp) :: mc(nelem)
      real(wp) :: gam(nelem)
      real(wp) :: gam3(nelem)
      real(wp) :: rad(nelem)
      real(wp) :: wll(nelem,nshell)
      real(wp) :: rep(2,nelem)
      real(wp) :: polyr(4,nelem)
      real(wp) :: cxb(nelem)
      real(wp) :: ao_exp(nshell,nelem)
      real(wp) :: ao_lev(nshell,nelem)
      real(wp) :: lpar(0:2,nelem)
      real(wp) :: kpair(nelem,nelem)
      real(wp) :: kcnat(0:2,nelem)
      real(wp) :: radaes(nelem)
      real(wp) :: dpolc(nelem)
      real(wp) :: qpolc(nelem)
      integer  :: ao_pqn(nshell,nelem)
      integer  :: ao_l(nshell,nelem)
      integer  :: ao_n(nelem)
      integer  :: ao_typ(nshell,nelem)
      integer  :: metal(nelem)
      integer  :: cnval(nelem)
      character(len=30) :: timestp(nelem)
   end type gfn_parameter

   type chrg_parameter
      integer  :: n
      real(wp),allocatable :: en(:)
      real(wp),allocatable :: gam(:)
      real(wp),allocatable :: kappa(:)
      real(wp),allocatable :: alpha(:)
      real(wp),allocatable :: dpol(:)
      real(wp),allocatable :: beta(:)
   contains
      procedure :: allocate => allocate_chrgeq
      procedure :: deallocate => deallocate_chrgeq
   end type chrg_parameter

contains

subroutine allocate_chrgeq(self,n,extended)
   implicit none
   class(chrg_parameter) :: self
   integer,intent(in) :: n
   logical,intent(in),optional :: extended
   logical :: multipoles
   if (present(extended)) then
      multipoles = extended
   else
      multipoles = .false.
   endif
   call self%deallocate
   allocate( self%en(n),    source = 0.0_wp )
   allocate( self%gam(n),   source = 0.0_wp )
   allocate( self%kappa(n), source = 0.0_wp )
   allocate( self%alpha(n), source = 0.0_wp )
   if (multipoles) then
   allocate( self%dpol(n),  source = 0.0_wp )
   allocate( self%beta(n),  source = 0.0_wp )
   endif
end subroutine allocate_chrgeq

subroutine deallocate_chrgeq(self)
   implicit none
   class(chrg_parameter) :: self
   if (allocated(self%en))    deallocate(self%en)
   if (allocated(self%gam))   deallocate(self%gam)
   if (allocated(self%kappa)) deallocate(self%kappa)
   if (allocated(self%alpha)) deallocate(self%alpha)
   if (allocated(self%dpol))  deallocate(self%dpol)
   if (allocated(self%beta))  deallocate(self%beta)
end subroutine deallocate_chrgeq

subroutine gfn0_param_defaults(param)
   implicit none
   type(scc_parameter)  :: param
   integer :: i,j
   param%kspd(1:6) = [2.0_wp, 2.4868_wp, 2.27_wp, 0.6_wp, 0.0_wp, -0.1_wp]
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         param%kmagic(i,j)=(param%kspd(i)+param%kspd(j))*0.50_wp
      enddo
   enddo
   param%gscal    = 0.0_wp   ! purpose changed in scf.f !
   param%kcnsh    =[0.0_wp, 0.0537_wp, -0.0129_wp, 3.4847_wp ]
   param%gscal    = 0.5097_wp         ! EN dep
   param%kenscal  =-0.2_wp         ! ken² (d3atm in old main.f)
   param%xbdamp   = 4.0_wp         ! ken⁴
   param%xbrad    =-0.09_wp         ! rep dEN
   param%g_a      = 3.0_wp
   param%g_c      = 2.0_wp
   param%wf       = 6.0_wp
   param%alphaj   = 1.1241_wp
   param%disp = dftd_parameter( s6 = 1.0_wp, s8 = 2.85_wp, s9 = 0.0_wp, &
      &                         a1 = 0.8_wp, a2 = 4.6_wp )
   param%ipshift  =0.0_wp
   param%eashift  =0.0_wp
   param%ken1     =1.0_wp
   param%zqf      =0.0_wp
   param%zcnf     =0.0_wp
   param%fpol     =0.0_wp
   param%wllscal  =1.0_wp

end subroutine gfn0_param_defaults

subroutine gfn0_param_from_globpar(param,globpar)
   implicit none
   type(scc_parameter)  :: param
   real(wp), intent(in) :: globpar(25)
   integer :: i,j
   param%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         param%kmagic(i,j)=(param%kspd(i)+param%kspd(j))*0.50_wp
      enddo
   enddo
   if(param%kspd(5).gt.0.1)then   ! s - p
      param%kmagic(1,2)=param%kspd(5)
      param%kmagic(2,1)=param%kspd(5)
   endif
   param%gscal    =globpar(8)*0.1_wp   ! purpose changed in scf.f !
   param%gam3l(0) =1.00_wp     !s
   param%gam3l(1) =globpar(9) !p
   param%gam3l(2) =globpar(10)!d-pol
   param%gam3l(3) =globpar(11)!d-val
   param%kcnsh(1) =globpar(17)         ! K 2s - 2s
   param%kcnsh(2) =globpar( 9)         ! SRB shift
   param%kcnsh(3) =globpar(10)         ! SRB prefactor
   param%kcnsh(4) =globpar(11)         ! SRB steepnes
   param%gscal    =globpar(12)         ! EN dep
   param%kenscal  =globpar(19)         ! ken² (d3atm in old main.f)
   param%xbdamp   =globpar(24)         ! ken⁴
   param%xbrad    =globpar(25)         ! rep dEN
   param%g_a      =3.0_wp
   param%g_c      =2.0_wp
   param%wf       =6.0_wp
   param%alphaj   =globpar(18)
   param%disp%a1  =globpar(20)
   param%disp%a2  =globpar(21)
   param%disp%s6  =1.0_wp
   param%disp%s9  =0.0_wp
   param%disp%s8  =globpar(22)
   param%ipshift  =globpar(7)*0.1_wp
   param%eashift  =param%ipshift
   param%ken1     =1.0_wp
   param%zqf      =0.0_wp
   param%zcnf     =0.0_wp
   param%fpol     =0.0_wp
   param%wllscal  =1.0_wp
end subroutine gfn0_param_from_globpar

subroutine gfn1_param_from_globpar(param,globpar)
   implicit none
   type(scc_parameter)  :: param
   real(wp), intent(in) :: globpar(25)
   integer :: i,j
   param%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         param%kmagic(i,j)=(param%kspd(i)+param%kspd(j))*0.50_wp
      enddo
   enddo
   if(param%kspd(5).gt.0.1_wp)then   ! s - p
      param%kmagic(1,2)=param%kspd(5)
      param%kmagic(2,1)=param%kspd(5)
   endif
   param%gscal    =globpar(8)*0.1_wp    ! purpose changed in scf.f !
   param%gam3l(0) =1.00_wp     !s
   param%gam3l(1) =globpar(9) !p
   param%gam3l(2) =globpar(10)!d-pol
   param%gam3l(3) =globpar(11)!d-val
   param%kcnsh(1) =globpar(11)*0.01_wp
   param%kcnsh(2) =globpar(12)*0.01_wp
   param%kcnsh(3) =globpar(13)*0.01_wp
   param%kcnsh(4) = 0.005_wp
   param%kenscal  =globpar(23) ! kenscal in scf.f
   param%xbdamp   =globpar(24)
   param%xbrad    =globpar(25)
   param%disp%s9  =0.0_wp! d3atm
   param%alphaj   =globpar(18)
   param%disp%a1  =globpar(20)
   param%disp%a2  =globpar(21)
   param%disp%s6  =1.0_wp
   param%disp%s8  =globpar(22)
   param%ipshift  =globpar(7)*0.1
   param%eashift  =param%ipshift
   param%ken1     =1.0_wp
   param%zqf      =0.0_wp
   param%zcnf     =0.0_wp
   param%fpol     =0.0_wp
   param%wllscal  =1.0_wp
end subroutine gfn1_param_from_globpar

subroutine gfn2_param_from_globpar(param,globpar)
   implicit none
   type(scc_parameter)  :: param
   real(wp), intent(in) :: globpar(25)
   integer :: i,j
   param%kspd(1:6)=globpar(1:6)
   ! ini prop factors for Hav(l1,l2), NO f-AO !
   do i=1,3
      do j=1,3
         param%kmagic(i,j)=(param%kspd(i)+param%kspd(j))*0.50_wp
      enddo
   enddo
   if(param%kspd(5).gt.0.1)then   ! s - p
      param%kmagic(1,2)=param%kspd(5)
      param%kmagic(2,1)=param%kspd(5)
   endif
   param%kmagic(1,3)=param%kspd(4)
   param%kmagic(3,1)=param%kspd(4)
   param%kmagic(2,3)=param%kspd(6)
   param%kmagic(3,2)=param%kspd(6)
   param%gscal    =globpar(8)*0.1_wp   ! purpose changed in scf.f !
   param%gam3l(0) =1.00_wp     !s
   param%gam3l(1) =globpar(9) !p
   param%gam3l(2) =globpar(10)!d-pol
   param%gam3l(3) =globpar(11)!d-val
   param%cn_shift =globpar(14) ! R AES CN val offset
   param%cn_expo  =globpar(15) ! R AES CN steepness
   param%cn_rmax  =globpar(16) ! R AES CN Rmax
   param%kenscal  =globpar(13) ! kenscal in scf.f
   param%disp%s9  =globpar(23) ! d3atm
   param%g_a      =3.0_wp
   param%g_c      =2.0_wp
   param%wf       =6.0_wp
   param%xbrad    =globpar(24)
   param%xbdamp   =globpar(25)
   param%alphaj   =globpar(18)
   param%disp%a1  =globpar(20)
   param%disp%a2  =globpar(21)
   param%disp%s6  =1.0_wp
   param%disp%s8  =globpar(22)
   param%ipshift  =globpar(7)*0.1_wp
   param%eashift  =param%ipshift
   param%ken1     =1.0_wp
   param%zqf      =0.0_wp
   param%zcnf     =0.0_wp
   param%fpol     =0.0_wp
   param%wllscal  =1.0_wp
end subroutine gfn2_param_from_globpar

end module tbdef_param
