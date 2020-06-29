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
   use xtb_disp_dftd3param, only : copy_c6, reference_c6
   use xtb_disp_dftd4, only : newD4Model, p_refq_gfn2xtb, p_refq_goedecker

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
   real(wp) :: kqat(0:2,max_elem)
   real(wp) :: kpair(max_elem,max_elem)
   real(wp) :: kcnat(0:2,max_elem)
   real(wp) :: eeqEN(max_elem)
   real(wp) :: kExpLight, kExp
   character(len=30) :: timestp(max_elem)
   type(dftd_parameter) :: disp

   character(len=:), allocatable :: line

   integer :: mShell, iSh, jSh
   integer :: level
   integer :: err
   logical :: newFormat

   disp = dftd_parameter(s6=1.0_wp, s8=0.0_wp, a1=0.0_wp, a2=0.0_wp, s9=0.0_wp)
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

   kExpLight = 0.0_wp
   level = -1
   newFormat = .false.
   call getline(iunit,line,err)
   if (debug) print'(">",a)',line
   readgroups: do
      if (index(line,flag).eq.1) then
         select case(line(2:))
         case('info')
            newFormat = .true.
            call read_info
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

   call setpair(level, kpair)

   mShell = maxval(nShell)
   xtbData%level = level
   xtbData%nShell = nShell
   xtbData%ipeashift = globpar%ipeashift * 0.1_wp

   ! Repulsion
   call init(xtbData%repulsion, kExp, kExpLight, 1.0_wp, globpar%renscale, &
      & repAlpha, repZeff, electronegativity)

   ! Coulomb
   xtbData%coulomb%gExp = globpar%alphaj
   xtbData%coulomb%chemicalHardness = atomicHardness(:max_elem)
   allocate(xtbData%coulomb%shellHardness(mShell, max_elem))
   call setGFN1ShellHardness(xtbData%coulomb%shellHardness, nShell, angShell, &
      & atomicHardness, shellHardness)
   xtbData%coulomb%thirdOrderAtom = thirdOrderAtom(:max_elem)
   xtbData%coulomb%electronegativity = eeqEN(:max_elem)
   xtbData%coulomb%kCN = eeqkCN(:max_elem)
   xtbData%coulomb%chargeWidth = chargeWidth(:max_elem)

   ! Dispersion
   xtbData%dispersion%dpar = disp
   xtbData%dispersion%g_a = 3.0_wp
   xtbData%dispersion%g_c = 2.0_wp
   xtbData%dispersion%wf  = 6.0_wp

   ! Hamiltonian
   xtbData%hamiltonian%angShell = angShell(:mShell, :)

   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%kScale(jSh, iSh) = 0.5_wp * (globpar%kShell(iSh) &
            & + globpar%kShell(jSh))
      end do
   end do
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
   xtbData%hamiltonian%kDiff = globpar%kDiff

   do iSh = 0, 3
      do jSh = 0, 3
         xtbData%hamiltonian%enScale(jSh, iSh) = 0.005_wp * (globpar%enShell(iSh) &
            & + globpar%enShell(jSh))
      end do
   end do
   xtbData%hamiltonian%enScale4 = globpar%enscale4

   xtbData%hamiltonian%electronegativity = electronegativity(:)
   xtbData%hamiltonian%atomicRad = atomicRad(:)
   xtbData%hamiltonian%shellPoly = shellPoly(:, :)
   xtbData%hamiltonian%pairParam = kpair(:, :)
   xtbData%hamiltonian%selfEnergy = selfEnergy(:mShell, :)
   xtbData%hamiltonian%slaterExponent = slaterExponent(:mShell, :)
   xtbData%hamiltonian%principalQuantumNumber = principalQuantumNumber(:mShell, :)

   allocate(xtbData%hamiltonian%valenceShell(mShell, max_elem))
   call generateValenceShellData(xtbData%hamiltonian%valenceShell, &
      & xtbData%nShell, xtbData%hamiltonian%angShell)

   select case(level)
   case(0)
      ! Hamiltonian
      xtbData%hamiltonian%wExp = 1.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)

      allocate(xtbData%hamiltonian%kQShell(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kQShell, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kqat)
      xtbData%hamiltonian%kQAtom = kqat2(:)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN0NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%srb)
      xtbData%srb%shift = globpar%srbshift
      xtbData%srb%prefactor = globpar%srbpre
      xtbData%srb%steepness = globpar%srbexp
      xtbData%srb%enScale = globpar%srbken

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_goedecker)

   case(1)
      ! Halogen
      allocate(xtbData%halogen)
      call init(xtbData%halogen, globpar%xbrad, globpar%xbdamp, halogenBond)

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.0_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call setGFN1kCN(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, xtbData%hamiltonian%selfEnergy, &
         & globpar%cnshell)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN1ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN1NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      if (.not.allocated(reference_c6)) call copy_c6(reference_c6)

   case(2)
      ! Coulomb
      if (any(globpar%gam3shell > 0.0_wp)) then
         allocate(xtbData%Coulomb%thirdOrderShell(mShell, max_elem))
         call setGFN2ThirdOrderShell(xtbData%Coulomb%thirdOrderShell, &
            & xtbData%nShell, xtbData%hamiltonian%angShell, thirdOrderAtom, &
            & globpar%gam3shell)
      end if

      ! Multipole
      allocate(xtbData%multipole)
      call init(xtbData%multipole, globpar%aesshift, globpar%aesexp, &
         & globpar%aesrmax, globpar%aesdmp3, globpar%aesdmp5, &
         & dipKernel, quadKernel)

      ! Hamiltonian
      xtbData%hamiltonian%wExp = 0.5_wp

      allocate(xtbData%hamiltonian%kCN(mShell, max_elem))
      call angToShellData(xtbData%hamiltonian%kCN, xtbData%nShell, &
         & xtbData%hamiltonian%angShell, kcnat)

      allocate(xtbData%hamiltonian%referenceOcc(mShell, max_elem))
      call setGFN2ReferenceOcc(xtbData%hamiltonian, xtbData%nShell)

      allocate(xtbData%hamiltonian%numberOfPrimitives(mShell, max_elem))
      call setGFN2NumberOfPrimitives(xtbData%hamiltonian, xtbData%nShell)

      ! Dispersion
      call newD4Model(xtbData%dispersion%dispm, xtbData%dispersion%g_a, &
         & xtbData%dispersion%g_c, p_refq_gfn2xtb)

   end select

contains

subroutine read_info
   use xtb_readin, only : getValue
   character(len=:), allocatable :: key, val
   integer :: ie, idum
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

      select case(key)
      case default
         call env%warning("Unknown key '"//key//"' for '"//flag//"info'")
      case('level')
         if (getValue(env,val,idum)) level = idum
      case('name')
         xtbData%name = val
      case('doi')
         xtbData%doi = val
      end select
   enddo
end subroutine read_info

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
   case('kexp'); if (getValue(env,val,ddum)) kExp = ddum
   case('kexplight'); if (getValue(env,val,ddum)) kExpLight = ddum
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
   case('cns'); if (getValue(env,val,ddum)) globpar%cnshell(:, 0) = ddum
   case('cnp'); if (getValue(env,val,ddum)) globpar%cnshell(:, 1) = ddum
   case('cnd'); if (getValue(env,val,ddum)) globpar%cnshell(:, 2) = ddum
   case('cnf'); if (getValue(env,val,ddum)) globpar%cnshell(:, 3) = ddum
   case('cnd1'); if (getValue(env,val,ddum)) globpar%cnshell(1, 2) = ddum
   case('cnd2'); if (getValue(env,val,ddum)) globpar%cnshell(2, 2) = ddum
   case('gam3s'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 0) = ddum
   case('gam3p'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 1) = ddum
   case('gam3d'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 2) = ddum
   case('gam3f'); if (getValue(env,val,ddum)) globpar%gam3shell(:, 3) = ddum
   case('gam3d1'); if (getValue(env,val,ddum)) globpar%gam3shell(1, 2) = ddum
   case('gam3d2'); if (getValue(env,val,ddum)) globpar%gam3shell(2, 2) = ddum
   case('srbshift'); if (getValue(env,val,ddum)) globpar%srbshift = ddum
   case('srbpre'); if (getValue(env,val,ddum)) globpar%srbpre = ddum
   case('srbexp'); if (getValue(env,val,ddum)) globpar%srbexp = ddum
   case('srbken'); if (getValue(env,val,ddum)) globpar%srbken = ddum
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
   case('a1'); if (getValue(env,val,ddum)) disp%a1 = ddum
   case('a2'); if (getValue(env,val,ddum)) disp%a2 = ddum
   case('s6'); if (getValue(env,val,ddum)) disp%s6 = ddum
   case('s8'); if (getValue(env,val,ddum)) disp%s8 = ddum
   case('s9'); if (getValue(env,val,ddum)) disp%s9 = ddum
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
      timestp(iz) = line(7:len_trim(line))
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
   case('kqs');   if (getValue(env,val,ddum)) kqat(0,iz)  = ddum
   case('kqp');   if (getValue(env,val,ddum)) kqat(1,iz)  = ddum
   case('kqd');   if (getValue(env,val,ddum)) kqat(2,iz)  = ddum
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
