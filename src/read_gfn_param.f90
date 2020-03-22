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

module xtb_readparam
   use xtb_xtb_data
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2

contains

subroutine readParam &
      (env, iunit,globpar,xtbData,initialize)
   use xtb_mctc_accuracy, only : wp

   use xtb_aoparam

   use xtb_readin, only : getline => strip_line
   use xtb_type_environment, only : TEnvironment
   use xtb_type_param, only : TxTBParameter

   implicit none

   logical, parameter :: debug = .false.

   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: iunit
   type(TxTBParameter), intent(inout) :: globpar
   type(TxTBData), intent(out) :: xtbData
   logical, intent(in) :: initialize

   character(len=1), parameter :: equal = '='
   character(len=1), parameter :: space = ' '
   character(len=1), parameter :: flag = '$'
   character(len=*), parameter :: flag_end = flag//'end'
   integer, parameter :: p_str_length = 48
   integer, parameter :: p_arg_length = 24

   character(len=:), allocatable :: line

   integer :: version
   integer :: err
   logical :: newFormat

   if (initialize) then
      globpar = TxTBParameter()

      ao_pqn=0
      ao_l  =0
      ao_n  =0
      ao_lev=0.0_wp
      ao_exp=0.0_wp
      ao_typ=0
      polyr =0.0_wp
      cxb   =0.0_wp
      rep   =0.0_wp
      mc    =0.0_wp
      lpar  =0.0_wp
      gam3  =0.0_wp
      kqat2 =0.0_wp
      eeqkcn=0.0_wp
      eeqen =0.0_wp
      kcnat =0.0_wp

      dpolc =0.0_wp ! read values are scaled by 0.01
      qpolc =0.0_wp !  "     "     "    "    "   "
      radaes=5.0_wp ! default atom radius
      radaes(1) =1.4_wp
      radaes(2) =3.0_wp
      radaes(6) =3.0_wp
      radaes(7) =1.9_wp
      radaes(8) =1.8_wp
      radaes(9) =2.4_wp
      radaes(14)=3.9_wp
      radaes(15)=2.1_wp
      radaes(16)=3.1_wp
      radaes(17)=2.5_wp
      radaes(34)=3.9_wp
      radaes(35)=4.0_wp

      kpair =1.0_wp
   endif

   version = -1
   newFormat = .false.
   call getline(iunit,line,err)
   if (debug) print'(">",a)',line
   readgroups: do
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('level 0')
            newFormat = .true.
            version = 0
            call getline(iunit,line,err)
         case('level 1')
            newFormat = .true.
            version = 1
            call getline(iunit,line,err)
         case('level 2')
            newFormat = .true.
            version = 2
            call getline(iunit,line,err)
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

   if (.not.newFormat) then
      call env%error("Old format parameter file is not supported anymore")
   end if

   call setpair(version)

   select case(version)
   case(0)
      call initGFN0(xtbData)
   case(1)
      call initGFN1(xtbData)
   case(2)
      call initGFN2(xtbData)
   end select

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

subroutine gfn_globpar(key,val,globpar)
   use xtb_readin, only : getValue
   implicit none
   character(len=*), intent(in) :: key, val
   type(TxTBParameter), intent(inout) :: globpar
   real(wp) :: ddum
   select case(key)
   case default
      call env%warning("Unknown key '"//key//"' for '"//flag//"globpar'")
   case('ks'); if (getValue(env,val,ddum)) globpar%ks = ddum
   case('kp'); if (getValue(env,val,ddum)) globpar%kp = ddum
   case('kd'); if (getValue(env,val,ddum)) globpar%kd = ddum
   case('kf'); if (getValue(env,val,ddum)) globpar%kf = ddum
   case('kdiffa'); if (getValue(env,val,ddum)) globpar%kdiffa = ddum
   case('kdiffb'); if (getValue(env,val,ddum)) globpar%kdiffb = ddum
   case('wllscal'); if (getValue(env,val,ddum)) globpar%wllscal = ddum
   case('ipeashift'); if (getValue(env,val,ddum)) globpar%ipeashift = ddum
   case('gscal'); if (getValue(env,val,ddum)) globpar%gscal = ddum
   case('zcnf'); if (getValue(env,val,ddum)) globpar%zcnf = ddum
   case('tscal'); if (getValue(env,val,ddum)) globpar%tscal = ddum
   case('kcn'); if (getValue(env,val,ddum)) globpar%kcn = ddum
   case('fpol'); if (getValue(env,val,ddum)) globpar%fpol = ddum
   case('ken'); if (getValue(env,val,ddum)) globpar%ken = ddum
   case('lshift'); if (getValue(env,val,ddum)) globpar%lshift = ddum
   case('lshifta'); if (getValue(env,val,ddum)) globpar%lshifta = ddum
   case('split'); if (getValue(env,val,ddum)) globpar%split = ddum
   case('zqf'); if (getValue(env,val,ddum)) globpar%zqf = ddum
   case('alphaj'); if (getValue(env,val,ddum)) globpar%alphaj = ddum
   case('kexpo'); if (getValue(env,val,ddum)) globpar%kexpo = ddum
   case('dispa'); if (getValue(env,val,ddum)) globpar%dispa = ddum
   case('dispb'); if (getValue(env,val,ddum)) globpar%dispb = ddum
   case('dispc'); if (getValue(env,val,ddum)) globpar%dispc = ddum
   case('dispatm'); if (getValue(env,val,ddum)) globpar%dispatm = ddum
   case('xbdamp'); if (getValue(env,val,ddum)) globpar%xbdamp = ddum
   case('xbrad'); if (getValue(env,val,ddum)) globpar%xbrad = ddum
   case('aesdmp3'); if (getValue(env,val,ddum)) globpar%aesdmp3 = ddum
   case('aesdmp5'); if (getValue(env,val,ddum)) globpar%aesdmp5 = ddum
   case('aesshift'); if (getValue(env,val,ddum)) globpar%aesshift = ddum
   case('aesexp'); if (getValue(env,val,ddum)) globpar%aesexp = ddum
   case('aesrmax'); if (getValue(env,val,ddum)) globpar%aesrmax = ddum
   end select
end subroutine gfn_globpar

subroutine read_pairpar
   use xtb_mctc_strings
   use xtb_aoparam
   use xtb_readin, only : getValue
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
      if (getValue(env,trim(argv(1)),iAt) .and. &
         &getValue(env,trim(argv(2)),jAt) .and. &
         &getValue(env,trim(argv(3)),ddum)) then
            kpair(iAt,jAt) = ddum
            kpair(jAt,iAt) = ddum
      endif
   enddo
end subroutine read_pairpar

subroutine read_elempar
   use xtb_mctc_strings
   use xtb_aoparam
   use xtb_readin
   implicit none
   character(len=:), allocatable :: key, val
   integer :: iz, ie
   if (getValue(env,line(4:5),iz)) then
      timestp(iz) = line(7:35)
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
   use xtb_mctc_strings
   use xtb_aoparam
   use xtb_readin
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
      call env%warning("Unknown key '"//key//"' for '"//flag//"Z'")
   case('ao')
      !print'(a,":",a)',key,val
      if (mod(len(val),2).eq.0) then
         ao_n(iz) = len(val)/2
         do i = 1, ao_n(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (getValue(env,val(ii:ii),idum)) then
               ao_pqn(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); ao_l(i,iz) = 0
               case('p'); ao_l(i,iz) = 1
               case('d'); ao_l(i,iz) = 2
               case('f'); ao_l(i,iz) = 3
               case('g'); ao_l(i,iz) = 4
               case('S'); ao_l(i,iz) = 0
               end select
            endif
         enddo
      endif
   case('lev')
      call parse(val,space,argv,narg)
      if (narg .eq. ao_n(iz)) then
         do i = 1, ao_n(iz)
            if (getValue(env,trim(argv(i)),ddum)) ao_lev(i,iz) = ddum
         enddo
      endif
   case('exp')
      call parse(val,space,argv,narg)
      if (narg .eq. ao_n(iz)) then
         do i = 1, ao_n(iz)
            if (getValue(env,trim(argv(i)),ddum)) ao_exp(i,iz) = ddum
         enddo
      endif
   case('en');  if (getValue(env,val,ddum)) en(iz)    = ddum
   case('gam'); if (getValue(env,val,ddum)) gam(iz)   = ddum
   case('epr'); if (getValue(env,val,ddum)) mc(iz)    = ddum
   case('xi');  if (getValue(env,val,ddum)) eeqEN(iz) = ddum
   case('alpg')
      if (getValue(env,val,ddum)) then
         radaes(iz) = ddum
         alp0(iz)   = ddum
      endif
   case('gam3');  if (getValue(env,val,ddum)) gam3(iz)    = ddum * 0.1_wp
   case('kappa'); if (getValue(env,val,ddum)) eeqkCN(iz)  = ddum
   case('cxb');   if (getValue(env,val,ddum)) cxb(iz)     = ddum * 0.1_wp
   case('kqat2'); if (getValue(env,val,ddum)) kqat2(iz)   = ddum
   case('dpol');  if (getValue(env,val,ddum)) dpolc(iz)   = ddum * 0.01_wp
   case('qpol');  if (getValue(env,val,ddum)) qpolc(iz)   = ddum * 0.01_wp
   case('repa');  if (getValue(env,val,ddum)) rep(1,iz)   = ddum
   case('repb');  if (getValue(env,val,ddum)) rep(2,iz)   = ddum
   case('polys'); if (getValue(env,val,ddum)) polyr(1,iz) = ddum
   case('polyp'); if (getValue(env,val,ddum)) polyr(2,iz) = ddum
   case('polyd'); if (getValue(env,val,ddum)) polyr(3,iz) = ddum
   case('polyf'); if (getValue(env,val,ddum)) polyr(4,iz) = ddum
   case('lpars'); if (getValue(env,val,ddum)) lpar(0,iz)  = ddum * 0.1_wp
   case('lparp'); if (getValue(env,val,ddum)) lpar(1,iz)  = ddum * 0.1_wp
   case('lpard'); if (getValue(env,val,ddum)) lpar(2,iz)  = ddum * 0.1_wp
   case('kcns');  if (getValue(env,val,ddum)) kcnat(0,iz) = ddum * 0.1_wp
   case('kcnp');  if (getValue(env,val,ddum)) kcnat(1,iz) = ddum * 0.1_wp
   case('kcnd');  if (getValue(env,val,ddum)) kcnat(2,iz) = ddum * 0.1_wp
   case('kqs');   if (getValue(env,val,ddum)) kqat(1,iz)  = ddum
   case('kqp');   if (getValue(env,val,ddum)) kqat(2,iz)  = ddum
   case('kqd');   if (getValue(env,val,ddum)) kqat(3,iz)  = ddum
   end select
end subroutine gfn_elempar

end subroutine readParam

end module xtb_readparam


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
      use xtb_aoparam
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
            kpair(ii,jj)=kparam
            kpair(jj,ii)=kparam
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
            notset=abs(kpair(i,j)-1.0d0).lt.1.d-6 .and. &
     &             abs(kpair(j,i)-1.0d0).lt.1.d-6
            if(ii.gt.0.and.jj.gt.0.and.notset) then
               kpair(i,j)=0.5*(kp(ii)+kp(jj))
               kpair(j,i)=0.5*(kp(ii)+kp(jj))
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
