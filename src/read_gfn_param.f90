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
   use xtb_paramset

contains

subroutine readParam &
      (env, iunit,globpar,xtbData,initialize)
   use xtb_mctc_accuracy, only : wp

   use xtb_readin, only : getline => strip_line
   use xtb_type_environment, only : TEnvironment
   use xtb_type_param, only : TxTBParameter
   use xtb_param_paulingen, only : paulingEN
   use xtb_param_atomicrad, only : atomicRad
   use xtb_mctc_param, only: chemical_hardness

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

   integer, parameter :: max_elem = 118
   integer :: nShell(max_elem)
   integer :: principalQuantumNumber(10,max_elem)
   integer :: angShell(10,max_elem)
   real(wp) :: shellPoly(1:4,max_elem)
   real(wp) :: selfEnergy(10,max_elem)
   real(wp) :: slaterExponent(10,max_elem)
   real(wp) :: thirdOrderAtom(max_elem)
   real(wp) :: atomicHardness(max_elem)
   real(wp) :: shellHardness(1:4,max_elem)
   real(wp) :: electronegativity(max_elem)
   real(wp) :: repAlpha(max_elem)
   real(wp) :: repZeff(max_elem)
   real(wp) :: halogenBond(max_elem)
   real(wp) :: dipKernel(max_elem)
   real(wp) :: quadKernel(max_elem)
   real(wp) :: eeqkcn(max_elem)
   real(wp) :: chargeWidth(max_elem)
   real(wp) :: kqat2(max_elem)
   real(wp) :: kqat(3,max_elem)
   real(wp) :: kpair(max_elem,max_elem)
   real(wp) :: kcnat(0:2,max_elem)
   real(wp) :: eeqEN(max_elem)
   real(wp) :: kExpLight
   character(len=30) :: timestp(max_elem)

   character(len=:), allocatable :: line

   integer :: mShell
   integer :: version
   integer :: err
   logical :: newFormat

   if (initialize) then
      globpar = TxTBParameter()

      principalQuantumNumber=0
      angShell  =0
      nShell  =0
      selfEnergy=0.0_wp
      slaterExponent=0.0_wp
      shellPoly =0.0_wp
      halogenBond = 0.0_wp
      atomicHardness = chemical_hardness(:max_elem)
      shellHardness  =0.0_wp
      thirdOrderAtom  =0.0_wp
      kqat2 =0.0_wp
      eeqkcn=0.0_wp
      eeqen =0.0_wp
      kcnat =0.0_wp

      electronegativity = paulingEN(:max_elem)
      repAlpha = 0.0_wp
      repZeff = 0.0_wp

      dipKernel = 0.0_wp ! read values are scaled by 0.01
      quadKernel = 0.0_wp !  "     "     "    "    "   "

      kpair =1.0_wp
   endif

   kExpLight = 0.0_wp
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
            kExpLight = 1.5_wp
            call getline(iunit,line,err)
         case('level 1')
            newFormat = .true.
            version = 1
            kExpLight = 1.5_wp
            call getline(iunit,line,err)
         case('level 2')
            newFormat = .true.
            version = 2
            kExpLight = 1.0_wp
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

   call setpair(version, kpair)

   xtbData%nShell = nShell
   ! Repulsion
   call init(xtbData%repulsion, 1.5_wp, kExpLight, 1.0_wp, globpar%renscale, &
      & repAlpha, repZeff, electronegativity)
   ! Coulomb
   call init(xtbData%coulomb, nShell, atomicHardness, shellHardness, &
      & thirdOrderAtom, eeqen, eeqkcn, chargeWidth)
   ! Hamiltonian
   mShell = maxval(xtbData%nShell)
   xtbData%hamiltonian%angShell = angShell(:mShell, :)
   xtbData%hamiltonian%kScale = 0.5_wp * (spread(globpar%kshell, 1, 4) &
      & + spread(globpar%kshell, 2, 4))
   if (globpar%ksp > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,1) = globpar%ksp
      xtbData%hamiltonian%kScale(1,0) = globpar%ksp
   end if
   if (globpar%ksd > 0.0_wp) then
      xtbData%hamiltonian%kScale(0,2) = globpar%ksd
      xtbData%hamiltonian%kScale(2,0) = globpar%ksd
   end if
   if (globpar%kpd > 0.0_wp) then
      xtbData%hamiltonian%kScale(1,2) = globpar%kpd
      xtbData%hamiltonian%kScale(2,1) = globpar%kpd
   end if
   xtbData%hamiltonian%enScale = 0.005_wp * (spread(globpar%enshell, 1, 4) &
      & + spread(globpar%enshell, 2, 4))
   xtbData%hamiltonian%enScale4 = globpar%enscale4
   xtbData%hamiltonian%kDiff = globpar%kDiff
   xtbData%hamiltonian%electronegativity = electronegativity(:)
   xtbData%hamiltonian%atomicRad = atomicRad(:)
   xtbData%hamiltonian%shellPoly = shellPoly(:, :)
   xtbData%hamiltonian%pairParam = kpair(:, :)
   xtbData%hamiltonian%kCN = kcnat(:, :)
   xtbData%hamiltonian%selfEnergy = selfEnergy(:mShell, :)
   xtbData%hamiltonian%slaterExponent = slaterExponent(:mShell, :)
   xtbData%hamiltonian%principalQuantumNumber = principalQuantumNumber(:mShell, :)
   xtbData%hamiltonian%kQShell = kqat(:, :)
   xtbData%hamiltonian%kQAtom = kqat2(:)
   allocate(xtbData%hamiltonian%valenceShell(mShell, max_elem))
   call generateValenceShellData(xtbData%hamiltonian%valenceShell, &
      & xtbData%nShell, xtbData%hamiltonian%angShell)
   select case(version)
   case(0)
      ! Hamiltonian
      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)
      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN0NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

   case(1)
      ! Halogen
      allocate(xtbData%halogen)
      call init(xtbData%halogen, globpar%xbrad, globpar%xbdamp, halogenBond)
      ! Hamiltonian
      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN1ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)
      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN1NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

   case(2)
      ! Multipole
      allocate(xtbData%multipole)
      call init(xtbData%multipole, globpar%aesshift, globpar%aesexp, &
         & globpar%aesrmax, globpar%aesdmp3, globpar%aesdmp5, &
         & dipKernel, quadKernel)
      ! Hamiltonian
      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)
      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN2NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

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
   case('ks'); if (getValue(env,val,ddum)) globpar%kshell(0) = ddum
   case('kp'); if (getValue(env,val,ddum)) globpar%kshell(1) = ddum
   case('kd'); if (getValue(env,val,ddum)) globpar%kshell(2) = ddum
   case('kf'); if (getValue(env,val,ddum)) globpar%kshell(3) = ddum
   case('ksp'); if (getValue(env,val,ddum)) globpar%ksp = ddum
   case('ksd'); if (getValue(env,val,ddum)) globpar%ksd = ddum
   case('kpd'); if (getValue(env,val,ddum)) globpar%kpd = ddum
   case('kdiff'); if (getValue(env,val,ddum)) globpar%kdiff = ddum
   case('kdiffa'); if (getValue(env,val,ddum)) globpar%kdiffa = ddum
   case('kdiffb'); if (getValue(env,val,ddum)) globpar%kdiffb = ddum
   case('ens'); if (getValue(env,val,ddum)) globpar%enshell(0) = ddum
   case('enp'); if (getValue(env,val,ddum)) globpar%enshell(1) = ddum
   case('end'); if (getValue(env,val,ddum)) globpar%enshell(2) = ddum
   case('enf'); if (getValue(env,val,ddum)) globpar%enshell(3) = ddum
   case('enscale'); if (getValue(env,val,ddum)) globpar%enshell = ddum
   case('enscale4'); if (getValue(env,val,ddum)) globpar%enscale4 = ddum
   case('renscale'); if (getValue(env,val,ddum)) globpar%renscale = ddum
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
         nShell(iz) = len(val)/2
         do i = 1, nShell(iz)
            ii = 2*i-1
            !print*,i,ii,val(ii:ii),val(ii+1:ii+1)
            if (getValue(env,val(ii:ii),idum)) then
               principalQuantumNumber(i,iz) = idum
               select case(val(ii+1:ii+1))
               case('s'); angShell(i,iz) = 0
               case('p'); angShell(i,iz) = 1
               case('d'); angShell(i,iz) = 2
               case('f'); angShell(i,iz) = 3
               case('g'); angShell(i,iz) = 4
               case('S'); angShell(i,iz) = 0
               end select
            endif
         enddo
      endif
   case('lev')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) selfEnergy(i,iz) = ddum
         enddo
      endif
   case('exp')
      call parse(val,space,argv,narg)
      if (narg .eq. nShell(iz)) then
         do i = 1, nShell(iz)
            if (getValue(env,trim(argv(i)),ddum)) slaterExponent(i,iz) = ddum
         enddo
      endif
   case('en');  if (getValue(env,val,ddum)) electronegativity(iz)    = ddum
   case('gam'); if (getValue(env,val,ddum)) atomicHardness(iz)   = ddum
   case('xi');  if (getValue(env,val,ddum)) eeqEN(iz) = ddum
   case('alpg'); if (getValue(env,val,ddum)) chargeWidth(iz)   = ddum
   case('gam3');  if (getValue(env,val,ddum)) thirdOrderAtom(iz)    = ddum * 0.1_wp
   case('kappa'); if (getValue(env,val,ddum)) eeqkCN(iz)  = ddum
   case('cxb');   if (getValue(env,val,ddum)) halogenBond(iz) = ddum * 0.1_wp
   case('kqat2'); if (getValue(env,val,ddum)) kqat2(iz)   = ddum
   case('dpol');  if (getValue(env,val,ddum)) dipKernel(iz)   = ddum * 0.01_wp
   case('qpol');  if (getValue(env,val,ddum)) quadKernel(iz)   = ddum * 0.01_wp
   case('repa');  if (getValue(env,val,ddum)) repAlpha(iz)   = ddum
   case('repb');  if (getValue(env,val,ddum)) repZeff(iz)   = ddum
   case('polys'); if (getValue(env,val,ddum)) shellPoly(1,iz) = ddum
   case('polyp'); if (getValue(env,val,ddum)) shellPoly(2,iz) = ddum
   case('polyd'); if (getValue(env,val,ddum)) shellPoly(3,iz) = ddum
   case('polyf'); if (getValue(env,val,ddum)) shellPoly(4,iz) = ddum
   case('lpars'); if (getValue(env,val,ddum)) shellHardness(1,iz)  = ddum * 0.1_wp
   case('lparp'); if (getValue(env,val,ddum)) shellHardness(2,iz)  = ddum * 0.1_wp
   case('lpard'); if (getValue(env,val,ddum)) shellHardness(3,iz)  = ddum * 0.1_wp
   case('lparf'); if (getValue(env,val,ddum)) shellHardness(4,iz)  = ddum * 0.1_wp
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
