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

subroutine read_gfn_param &
      (iunit,globpar,initialize)
   use iso_fortran_env, wp => real64

   use aoparam

   use readin, only : getline => strip_line

   implicit none

   logical, parameter :: debug = .false.

   integer, intent(in) :: iunit
   real(wp), intent(inout) :: globpar(25)
   logical, intent(in) :: initialize

   character(len=1), parameter :: equal = '='
   character(len=1), parameter :: space = ' '
   character(len=1), parameter :: flag = '$'
   character(len=*), parameter :: flag_end = flag//'end'
   integer, parameter :: p_str_length = 48
   integer, parameter :: p_arg_length = 24

   character(len=:), allocatable :: line

   integer :: err

   if (initialize) then
      globpar = 0.0_wp

      gfn = tb_hamiltonian()
      gfn%radaes=5.0_wp ! default atom radius
      gfn%radaes(1) =1.4_wp
      gfn%radaes(2) =3.0_wp
      gfn%radaes(6) =3.0_wp
      gfn%radaes(7) =1.9_wp
      gfn%radaes(8) =1.8_wp
      gfn%radaes(9) =2.4_wp
      gfn%radaes(14)=3.9_wp
      gfn%radaes(15)=2.1_wp
      gfn%radaes(16)=3.1_wp
      gfn%radaes(17)=2.5_wp
      gfn%radaes(34)=3.9_wp
      gfn%radaes(35)=4.0_wp

      gfn%kpair =1.0_wp
   endif

   call getline(iunit,line,err)
   if (debug) print'(">",a)',line
   readgroups: do
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('globpar')
            call read_globpar
         case('pairpar')
            call read_pairpar
         case default
            if (index(line,'Z').eq.2) then
               call read_elempar
            else
               call getline(iunit,line,err)
               if (debug) print'(">",a)',line
            endif
         end select
      else
         call getline(iunit,line,err)
         if (debug) print'(">",a)',line
      endif
      if (err.ne.0) exit readgroups
      !if (index(line,flag_end).gt.0) exit readgroups
   enddo readgroups

contains

subroutine read_globpar
   implicit none
   character(len=:), allocatable :: key, val
   integer :: ie
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      ie = index(line,space)
      if (line.eq.'') cycle ! skip empty lines
      if (ie.eq.0) cycle

      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))

      call gfn_globpar(key,val,globpar)
      
   enddo
end subroutine read_globpar

subroutine gfn0_globpar(key,val,param)
   use tbdef_param
   use readin, only : get_value
   implicit none
   character(len=*), intent(in) :: key, val
   type(scc_parameter), intent(inout) :: param
   real(wp) :: ddum
   select case(key)
   case default
      call raise('S',"Unknown key '"//key//"' for '"//flag//"globpar'",1)
   case('ks'); if (get_value(val,ddum)) param%kspd(1) = ddum
   case('kp'); if (get_value(val,ddum)) param%kspd(2) = ddum
   case('kd'); if (get_value(val,ddum)) param%kspd(3) = ddum
   case('kf'); if (get_value(val,ddum)) param%kspd(4) = ddum
   case('kdiffa'); if (get_value(val,ddum)) param%kspd(5) = ddum
   case('kdiffb'); if (get_value(val,ddum)) param%kspd(6) = ddum
   case('wllscal')
      if (get_value(val,ddum)) then
         param%ipshift = ddum * 0.1_wp
         param%eashift = ddum * 0.1_wp
      endif
   case('gscal'); if (get_value(val,ddum)) param%gscal = ddum * 0.1_wp
   case('zcnf'); if (get_value(val,ddum)) param%gam3l(1) = ddum
   case('tscal'); if (get_value(val,ddum)) param%gam3l(2) = ddum
   case('kcn'); if (get_value(val,ddum)) param%gam3l(3) = ddum
   case('fpol'); if (get_value(val,ddum)) param%gscal = ddum
   case('zqf'); if (get_value(val,ddum)) param%kcnsh(1) = ddum
   case('alphaj'); if (get_value(val,ddum)) param%alphaj = ddum
   case('kexpo'); if (get_value(val,ddum)) param%kenscal = ddum
   case('dispa'); if (get_value(val,ddum)) param%disp%a1 = ddum
   case('dispb'); if (get_value(val,ddum)) param%disp%a2 = ddum
   case('dispc'); if (get_value(val,ddum)) param%disp%s8 = ddum
   case('xbdamp'); if (get_value(val,ddum)) param%xbdamp = ddum
   case('xbrad'); if (get_value(val,ddum)) param%xbrad = ddum
   end select
end subroutine gfn0_globpar

subroutine gfn1_globpar(key,val,param)
   use tbdef_param
   use readin, only : get_value
   implicit none
   character(len=*), intent(in) :: key, val
   type(scc_parameter), intent(inout) :: param
   real(wp) :: ddum
   select case(key)
   case default
      call raise('S',"Unknown key '"//key//"' for '"//flag//"globpar'",1)
   case('ks'); if (get_value(val,ddum)) param%kspd(1) = ddum
   case('kp'); if (get_value(val,ddum)) param%kspd(2) = ddum
   case('kd'); if (get_value(val,ddum)) param%kspd(3) = ddum
   case('kf'); if (get_value(val,ddum)) param%kspd(4) = ddum
   case('kdiffa'); if (get_value(val,ddum)) param%kspd(5) = ddum
   case('kdiffb'); if (get_value(val,ddum)) param%kspd(6) = ddum
   case('wllscal')
      if (get_value(val,ddum)) then
         param%ipshift = ddum * 0.1_wp
         param%eashift = ddum * 0.1_wp
      endif
   case('gscal'); if (get_value(val,ddum)) param%gscal = ddum * 0.1_wp
   case('zcnf'); if (get_value(val,ddum)) param%gam3l(1) = ddum
   case('tscal'); if (get_value(val,ddum)) param%gam3l(2) = ddum
   case('kcn')
      if (get_value(val,ddum)) then
         param%gam3l(3) = ddum
         param%kcnsh(1) = ddum * 0.01_wp
      endif
   case('fpol'); if (get_value(val,ddum)) param%kcnsh(2) = ddum * 0.01_wp
   case('ken'); if (get_value(val,ddum)) param%kcnsh(3) = ddum * 0.01_wp
   case('alphaj'); if (get_value(val,ddum)) param%alphaj = ddum
   case('dispa'); if (get_value(val,ddum)) param%disp%a1 = ddum
   case('dispb'); if (get_value(val,ddum)) param%disp%a2 = ddum
   case('dispc'); if (get_value(val,ddum)) param%disp%s8 = ddum
   case('dispatm'); if (get_value(val,ddum)) param%kenscal = ddum
   case('xbdamp'); if (get_value(val,ddum)) param%xbdamp = ddum
   case('xbrad'); if (get_value(val,ddum)) param%xbrad = ddum
   end select
end subroutine gfn1_globpar

subroutine gfn2_globpar(key,val,param)
   use tbdef_param
   use readin, only : get_value
   implicit none
   character(len=*), intent(in) :: key, val
   type(scc_parameter), intent(inout) :: param
   real(wp) :: ddum
   select case(key)
   case default
      call raise('S',"Unknown key '"//key//"' for '"//flag//"globpar'",1)
   case('ks'); if (get_value(val,ddum)) param%kspd(1) = ddum
   case('kp'); if (get_value(val,ddum)) param%kspd(2) = ddum
   case('kd'); if (get_value(val,ddum)) param%kspd(3) = ddum
   case('kf'); if (get_value(val,ddum)) param%kspd(4) = ddum
   case('kdiffa'); if (get_value(val,ddum)) param%kspd(5) = ddum
   case('kdiffb'); if (get_value(val,ddum)) param%kspd(6) = ddum
   case('wllscal')
      if (get_value(val,ddum)) then
         param%ipshift = ddum * 0.1_wp
         param%eashift = ddum * 0.1_wp
      endif
   case('gscal'); if (get_value(val,ddum)) param%gscal = ddum * 0.1_wp
   case('zcnf'); if (get_value(val,ddum)) param%gam3l(1) = ddum
   case('tscal'); if (get_value(val,ddum)) param%gam3l(2) = ddum
   case('kcn'); if (get_value(val,ddum)) param%gam3l(3) = ddum
   case('lshift'); if (get_value(val,ddum)) param%cn_shift = ddum
   case('lshifta'); if (get_value(val,ddum)) param%cn_expo = ddum
   case('split'); if (get_value(val,ddum)) param%cn_rmax = ddum
   case('alphaj'); if (get_value(val,ddum)) param%alphaj = ddum
   case('dispa'); if (get_value(val,ddum)) param%disp%a1 = ddum
   case('dispb'); if (get_value(val,ddum)) param%disp%a2 = ddum
   case('dispc'); if (get_value(val,ddum)) param%disp%s8 = ddum
   case('dispatm'); if (get_value(val,ddum)) param%disp%s9 = ddum
   case('xbdamp'); if (get_value(val,ddum)) param%xbdamp = ddum
   case('xbrad'); if (get_value(val,ddum)) param%xbrad = ddum
   end select
end subroutine gfn2_globpar

subroutine gfn_globpar(key,val,globpar)
   use readin, only : get_value
   implicit none
   character(len=*), intent(in) :: key, val
   real(wp), intent(inout) :: globpar(25)
   real(wp) :: ddum
   select case(key)
   case default
      call raise('S',"Unknown key '"//key//"' for '"//flag//"globpar'",1)
   case('ks'); if (get_value(val,ddum)) globpar(1) = ddum
   case('kp'); if (get_value(val,ddum)) globpar(2) = ddum
   case('kd'); if (get_value(val,ddum)) globpar(3) = ddum
   case('kf'); if (get_value(val,ddum)) globpar(4) = ddum
   case('kdiffa'); if (get_value(val,ddum)) globpar(5) = ddum
   case('kdiffb'); if (get_value(val,ddum)) globpar(6) = ddum
   case('wllscal'); if (get_value(val,ddum)) globpar(7) = ddum
   case('gscal'); if (get_value(val,ddum)) globpar(8) = ddum
   case('zcnf'); if (get_value(val,ddum)) globpar(9) = ddum
   case('tscal'); if (get_value(val,ddum)) globpar(10) = ddum
   case('kcn'); if (get_value(val,ddum)) globpar(11) = ddum
   case('fpol'); if (get_value(val,ddum)) globpar(12) = ddum
   case('ken'); if (get_value(val,ddum)) globpar(13) = ddum
   case('lshift'); if (get_value(val,ddum)) globpar(14) = ddum
   case('lshifta'); if (get_value(val,ddum)) globpar(15) = ddum
   case('split'); if (get_value(val,ddum)) globpar(16) = ddum
   case('zqf'); if (get_value(val,ddum)) globpar(17) = ddum
   case('alphaj'); if (get_value(val,ddum)) globpar(18) = ddum
   case('kexpo'); if (get_value(val,ddum)) globpar(19) = ddum
   case('dispa'); if (get_value(val,ddum)) globpar(20) = ddum
   case('dispb'); if (get_value(val,ddum)) globpar(21) = ddum
   case('dispc'); if (get_value(val,ddum)) globpar(22) = ddum
   case('dispatm'); if (get_value(val,ddum)) globpar(23) = ddum
   case('xbdamp'); if (get_value(val,ddum)) globpar(24) = ddum
   case('xbrad'); if (get_value(val,ddum)) globpar(25) = ddum
   end select
end subroutine gfn_globpar

subroutine read_pairpar
   use mctc_strings
   use aoparam
   use readin, only : get_value
   implicit none
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer  :: iAt,jAt
   real(wp) :: ddum
   do
      call getline(iunit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) exit
      if (index(line,flag).gt.0) exit

      call parse(line,space,argv,narg)
      if (narg .ne. 3) cycle
      if (get_value(trim(argv(1)),iAt) .and. &
         &get_value(trim(argv(2)),jAt) .and. &
         &get_value(trim(argv(3)),ddum)) then
            gfn%kpair(iAt,jAt) = ddum
            gfn%kpair(jAt,iAt) = ddum
      endif
   enddo
end subroutine read_pairpar

subroutine read_elempar
   use mctc_strings
   use aoparam
   use readin
   implicit none
   character(len=:), allocatable :: key, val
   integer :: iz, ie
   if (get_value(line(4:5),iz)) then
      gfn%timestp(iz) = line(7:35)
      do
         call getline(iunit,line,err)
         if (debug) print'("->",a)',line
         if (err.ne.0) exit
         if (index(line,flag).gt.0) exit

         ie = index(line,equal)
         if (line.eq.'') cycle ! skip empty lines
         if (ie.eq.0) cycle

         key = lowercase(trim(line(:ie-1)))
         val = trim(adjustl(line(ie+1:)))

         call gfn_elempar(key,val,iz)

      enddo
   else
      call getline(iunit,line,err)
   endif
end subroutine read_elempar

subroutine gfn_elempar(key,val,iz)
   use mctc_strings
   use aoparam, only : gfn
   use readin
   implicit none
   character(len=*), intent(in) :: key, val
   integer, intent(in) :: iz
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   integer :: i, ii
   integer :: idum
   real(wp) :: ddum
   select case(key)
   case default
      call raise('S',"Unknown key '"//key//"' for '"//flag//"Z'",1)
   case('ao')
      !print'(a,":",a)',key,val
      if (mod(len(val),2).eq.0) then
         gfn%ao_n(iz) = len(val)/2
         do i = 1, gfn%ao_n(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (get_value(val(ii:ii),idum)) then
               gfn%ao_pqn(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); gfn%ao_l(i,iz) = 0
               case('p'); gfn%ao_l(i,iz) = 1
               case('d'); gfn%ao_l(i,iz) = 2
               case('f'); gfn%ao_l(i,iz) = 3
               case('g'); gfn%ao_l(i,iz) = 4
               case('S'); gfn%ao_l(i,iz) = 11
               end select
            endif
         enddo
      endif
   case('lev')
      call parse(val,space,argv,narg)
      if (narg .eq. gfn%ao_n(iz)) then
         do i = 1, gfn%ao_n(iz)
            if (get_value(trim(argv(i)),ddum)) gfn%ao_lev(i,iz) = ddum
         enddo
      endif
   case('exp')
      call parse(val,space,argv,narg)
      if (narg .eq. gfn%ao_n(iz)) then
         do i = 1, gfn%ao_n(iz)
            if (get_value(trim(argv(i)),ddum)) gfn%ao_exp(i,iz) = ddum
         enddo
      endif
   case('en');  if (get_value(val,ddum)) gfn%en(iz)    = ddum
   case('gam'); if (get_value(val,ddum)) gfn%gam(iz)   = ddum
   case('epr'); if (get_value(val,ddum)) gfn%mc(iz)    = ddum
   case('xi');  if (get_value(val,ddum)) gfn%dpolc(iz) = ddum
   case('alpg')
      if (get_value(val,ddum)) then
         gfn%radaes(iz) = ddum
         gfn%alp0(iz)   = ddum
      endif
   case('gam3');  if (get_value(val,ddum)) gfn%gam3(iz)    = ddum * 0.1_wp
   case('cxb');   if (get_value(val,ddum)) gfn%cxb(iz)     = ddum * 0.1_wp
   case('dpol');  if (get_value(val,ddum)) gfn%dpolc(iz)   = ddum * 0.01_wp
   case('qpol');  if (get_value(val,ddum)) gfn%qpolc(iz)   = ddum * 0.01_wp
   case('repa');  if (get_value(val,ddum)) gfn%rep(1,iz)   = ddum
   case('repb');  if (get_value(val,ddum)) gfn%rep(2,iz)   = ddum
   case('polys'); if (get_value(val,ddum)) gfn%polyr(1,iz) = ddum
   case('polyp'); if (get_value(val,ddum)) gfn%polyr(2,iz) = ddum
   case('polyd'); if (get_value(val,ddum)) gfn%polyr(3,iz) = ddum
   case('polyf'); if (get_value(val,ddum)) gfn%polyr(4,iz) = ddum
   case('lpars'); if (get_value(val,ddum)) gfn%lpar(0,iz)  = ddum * 0.1_wp
   case('lparp'); if (get_value(val,ddum)) gfn%lpar(1,iz)  = ddum * 0.1_wp
   case('lpard'); if (get_value(val,ddum)) gfn%lpar(2,iz)  = ddum * 0.1_wp
   case('kcns');  if (get_value(val,ddum)) gfn%kcnat(0,iz) = ddum * 0.1_wp
   case('kcnp');  if (get_value(val,ddum)) gfn%kcnat(1,iz) = ddum * 0.1_wp
   case('kcnd');  if (get_value(val,ddum)) gfn%kcnat(2,iz) = ddum * 0.1_wp
   case('kqs');   if (get_value(val,ddum)) gfn%kqat(1,iz)  = ddum
   case('kqp');   if (get_value(val,ddum)) gfn%kqat(2,iz)  = ddum
   case('kqd');   if (get_value(val,ddum)) gfn%kqat(3,iz)  = ddum
   end select
end subroutine gfn_elempar
   
end subroutine read_gfn_param


      logical function maingroup(i)
      integer i
      logical main_group(107)
      data main_group /& 
     &  2*.true.,&                              ! H  - He
     &  8*.true.,&                              ! Li - Ne
     &  8*.true.,&                              ! Na - Ar
     &  2*.true., 9*.false., 7*.true.,&         ! K  - Kr
     &  2*.true., 9*.false., 7*.true.,&         ! Rb - Xe
     &  2*.true.,23*.false., 7*.true.,&         ! Cs - Rn
     & 21*.true.                               / ! Fr - Tv

      maingroup = main_group(i)

      end function maingroup

! global, predefined pair parameters
      subroutine setpair(gfn_method)
      use aoparam
      implicit none 
      integer gfn_method
      integer i,j,ii,jj
      integer tmmetal
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
            gfn%kpair(ii,jj)=kparam
            gfn%kpair(jj,ii)=kparam
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
            notset=abs(gfn%kpair(i,j)-1.0d0).lt.1.d-6 .and. &
     &             abs(gfn%kpair(j,i)-1.0d0).lt.1.d-6
            if(ii.gt.0.and.jj.gt.0.and.notset) then  
               gfn%kpair(i,j)=0.5*(kp(ii)+kp(jj))
               gfn%kpair(j,i)=0.5*(kp(ii)+kp(jj))
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
