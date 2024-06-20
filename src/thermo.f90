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

module xtb_thermo
   use xtb_mctc_accuracy, only : wp
   implicit none

contains

subroutine getsymmetry (pr, iunit, n, iat, xyz, symthr, maxatdesy, sfsym)
   use iso_c_binding, only: c_char, c_null_char
   use xtb_mctc_accuracy, only : wp
   implicit none
   integer, intent(in) :: iunit
   integer n, iat(n), maxatdesy
   real(wp)  xyz(3,n)
   real(wp) symthr
   Character(len=*) sfsym
   logical pr
   Character(len=6) atmp
   integer,allocatable :: ictdum(:,:)

   Real(wp) :: paramar (11)  !parameter array for get_schoenflies_

   if(n.gt.maxatdesy) then
      if(pr)write(iunit,*) 'symmetry recognition skipped because # atoms >',maxatdesy
      return
   endif

   if (pr) write(iunit,'(a)')
   !parameters for symmetry recognition:
   paramar (1) = - 1         ! verbose, increase for more detailed output (to stdout)
   paramar (2) = 10          ! MaxAxisOrder
   paramar (3) = 100         ! MaxOptCycles
   paramar (4) = 0.001d0     ! ToleranceSame
   paramar (5) = 0.5d0       ! TolerancePrimary
   paramar (6) = symthr      ! ToleranceFinal, THIS IS THE IMPORTANT VALUE
   paramar (7) = 0.5d0       ! MaxOptStep
   paramar (8) = 1.0D-7      ! MinOptStep
   paramar (9) = 1.0D-7      ! GradientStep
   paramar (10) = 1.0D-8     ! OptChangeThreshold
   paramar (11) = 5          ! OptChangeHits

   atmp='    '
   Call get_schoenflies (n, iat, xyz, atmp, paramar)
   !call flush(iunit)

   !TM stuff (trafo table)
   sfsym(1:3)=atmp(1:3)
   if(sfsym(1:1).eq.'D') sfsym(1:1)='d'
   if(sfsym(1:1).eq.'C') sfsym(1:1)='c'
   if(sfsym(1:1).eq.'T') sfsym(1:1)='t'
   if(sfsym(1:1).eq.'O') sfsym(1:1)='o'
   if(sfsym(1:1).eq.'I') sfsym(1:1)='i'
   if(sfsym.eq.'dih') sfsym='d6h'
   if(sfsym.eq.'civ') sfsym='c6v'
   if(sfsym(3:3).gt.'v'.or.sfsym(3:3).lt.'a') sfsym(3:3)=' '

   if(pr) then
      write(iunit,'(a3,'' symmetry found (for desy threshold: '',e9.2,'') used in thermo'')')sfsym, symthr
   endif

End subroutine getsymmetry

subroutine thermodyn(iunit,A_rcm,B_rcm,C_rcm,avmom_si,linear,atom,sym,molmass, &
      &              vibs,nvibs,escf,T,sthr_rcm,et,ht,g,ts,zp,pr)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   implicit none
   integer, intent(in)  :: iunit       !< output_unit
   integer, intent(in)  :: nvibs       !< number of vibrational frequencies
   real(wp),intent(in)  :: escf        !< total energy from electronic structure
   real(wp),intent(in)  :: A_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: B_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: C_rcm       !< rotational constants in cm-1
   real(wp),intent(in)  :: avmom_si    !< average moment of inertia in whatever
   real(wp),intent(in)  :: sym         !< symmetry number
   real(wp),intent(in)  :: molmass     !< molecular mass in in amu
   real(wp),intent(in)  :: T           !< temperature in K
   real(wp),intent(in)  :: sthr_rcm    !< rotor cutoff in cm-1
   real(wp),intent(out) :: et          !< enthalpy in Eh
   real(wp),intent(out) :: ht          !< enthalpy in Eh
   real(wp),intent(out) :: g           !< free energy in Eh
   real(wp),intent(out) :: ts          !< entropy in Eh
   real(wp),intent(in)  :: zp          !< zero point vibrational energy in Eh
   real(wp),intent(in)  :: vibs(nvibs) !< vibrational frequencies in cm-1
   logical, intent(in)  :: linear      !< is linear
   logical, intent(in)  :: atom        !< only one atom
   logical, intent(in)  :: pr          !< clutter the screen with printout
   real(wp),parameter :: R = 1.98726D0    ! GAS CONSTANT IN CALORIES/MOLE
   real(wp),parameter :: H = 6.626176D-27 ! PLANCK'S CONSTANT IN ERG-SECONDS
   real(wp),parameter :: AK = 1.3807D-16  ! BOLTZMANN CONSTANT IN ERG/DEGREE
   real(wp),parameter :: conv3 = amutokg*1000 ! 1.6606d-24
   real(wp),parameter :: magic4 = 2.2868d0 ! unknown
   real(wp),parameter :: magic5 = 2.3135d0 ! unknown
   real(wp),parameter :: caltoj = autokj/autokcal

   integer  :: i
   real(wp) :: s_tr,s_rot,s_vib,s_int,s_tot
   real(wp) :: h_tr,h_rot,h_vib,h_int,h_tot
   real(wp) :: q_tr,q_rot,q_vib,q_int,q_tot
   real(wp) :: cptr,cprot,cpvib,cpint,cptot
   real(wp) :: beta,sthr,avmom,A,B,C
   real(wp) :: ewj,omega,e,xxmom,xmom,mu
   real(wp) :: wofrot,RT
   real(wp) :: cp_ho,cp_rr
   real(wp) :: sv_ho,sv_rr
   !*******************************************************************

   ! convert EVERYTHING to atomic units NOW and avoid horror and dispair later
   beta=1.0_wp/kB/T ! beta in 1/Eh
   !c1=h*ac/ak/T    ! beta in cm
   sthr  = sthr_rcm * rcmtoau ! sthr in Eh
   avmom = avmom_si*kgtome*aatoau**2*1.0e+20_wp ! in me·α² (au)
   A = A_rcm * rcmtoau
   B = B_rcm * rcmtoau
   C = C_rcm * rcmtoau

   !   ***   INITIALISE SOME VARIABLES   ***
   q_vib = 1.0_wp
   h_vib = 0.0_wp
   cpvib = 0.0_wp
   s_vib = 0.0_wp
   q_rot = 1.0_wp
   h_rot = 0.0_wp
   cprot = 0.0_wp
   s_rot = 0.0_wp
   ! free rotor heat capacity is constant 0.5*R
   ! we use the below to work with Eh² units 
   ! conversion to cal/mol/K is done later
   cp_rr = 0.5_wp/(beta**2) 

   ! construct the frequency dependent parts of partition function
   do i=1,nvibs
      omega=vibs(i)
      ! omega in Eh, beta in 1/Eh
      ewj=exp(-omega*beta)
      q_vib=q_vib/(1.0_wp-ewj)
      ! h_vib in Eh
      h_vib=h_vib+omega*ewj/(1.0_wp-ewj)
      ! cp_ho in Eh²
      cp_ho=omega**2 * ewj/(1.0_wp-ewj)/(1.0_wp-ewj)
      ! replace low-lying vibs for S by rotor approx.
      mu = 0.5_wp / (omega + 1.0e-14_wp)
      ! this reduced moment limits the rotational moment of inertia for
      ! this vibration to that of the total molecule rotation/3
      ! avmom and xmom are in me·α² (au)
      mu = mu*avmom/(mu+avmom)
      if(omega.gt.0)then
         ! sv is S/R which is dimensionless
         ! harm. osc. entropy
         sv_ho = vibs(i)*beta*ewj/(1.0_wp-ewj) - log(1.0_wp-ewj)
         ! free rotor entropy
         ! Cramer, page 328 for one degree of freedom or
         ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
         !sv_rr = (0.5_wp+log(sqrt(8.0_wp*pi**3*xmom*sik*t)/sih))
         sv_rr = 0.5_wp + log(sqrt(pi/beta*2.0_wp*mu))
      else
         sv_ho = 0.0_wp
         sv_rr = 0.0_wp
      endif
      ! fermi weigthing
      ! wofrot=1./(1.+exp( (omega-sthr)/20.0 ) )
      ! Head-Gordon weighting
      wofrot=1.0_wp-chg_switching(omega,sthr)
      ! heat capacity (cp_rr is a constant), all in Eh²
      cpvib = cpvib + ((1.0_wp-wofrot)*cp_ho + wofrot*cp_rr)
      ! entropy s_vib is converted to cal/mol/K... by multiplying with R
      s_vib=s_vib+R*((1.0_wp-wofrot)*sv_ho + wofrot*sv_rr)
   enddo
   !   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
   ! now unit conversion again...
   ! h_vib in Eh, beta is in 1/Eh, T is in K, R is in cal/mol/K,
   h_vib=h_vib*R*beta*T
   ! same here
   ! cpvib is in Eh², beta in 1/Eh, R in cal/mol/K
   cpvib=cpvib*R*beta**2
   !   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST LINEAR MOLECULES)
   if (.not.atom) then
      if(linear) then
         ! A is in Eh, beta is in 1/Eh, q_rot is dimensionless
         q_rot=1/(beta*A*sym)
         ! h_rot is in cal/mol
         h_rot=R*T
         ! cprot  is in cal/mol/K
         cprot=R
         ! s_rot is in cal/mol/K
         s_rot=R*((log(1.0_wp/(beta*A*sym)))+1.0_wp)
      else
         ! see above
         q_rot=sqrt(pi/(A*B*C*beta**3))/sym
         h_rot=3.0_wp*R*T/2.0_wp
         cprot=3.0_wp*R/2.0_wp
         s_rot=0.5_wp*R*(-2.0_wp*log(sym)+log(pi/(A*B*C*beta**3))+3.0_wp)
      endif
   endif
   !   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
   q_int=q_vib*q_rot
   h_int=h_vib+h_rot
   cpint=cpvib+cprot
   s_int=s_vib+s_rot
   !   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
   q_tr=(sqrt(2.0_wp*pi*molmass*t*ak*conv3)/h)**3
   ! this is 3/2rt+pv=5/2rt
   h_tr=5.0_wp*R*T/2.0_wp
   cptr=5.0_wp*R/2.0_wp
   s_tr=magic4*(5.0_wp*log10(t)+3.0_wp*log10(molmass))-magic5
   !   ***   CONSTRUCT TOTALS   ***
   cptot=cptr+cpint
   s_tot=s_tr+s_int
   h_tot=h_tr+h_int

   if(pr)then
      write(iunit,'(a)')
      write(iunit,'("   temp. (K)  partition function ", &
         & "  enthalpy   heat capacity  entropy")')
      write(iunit,'(  "                                   ", &
         & "cal/mol     cal/K/mol   cal/K/mol   J/K/mol")')
      write(iunit,'(  f7.2,"  VIB ",G10.3,10X,3F11.3)') &
         & T,q_vib,  h_vib,  cpvib,  s_vib
      write(iunit,'(7X,"  ROT ",G10.3,10X,3F11.3)') &
         & q_rot,  h_rot,  cprot,  s_rot
      write(iunit,'(7X,"  INT ",G10.3,10X,3F11.3)') &
         & q_int,h_int,cpint,s_int
      write(iunit,'(7X,"  TR  ",G10.3,10X,3F11.3)') &
         & q_tr, h_tr, cptr, s_tr
      write(iunit,'(7X,"  TOT ",21X,F11.4,3F11.4)') &
         & h_tot,cptot,s_tot,s_tot*caltoj
   endif

   ht=h_tot/1000.0_wp*kcaltoau
   et=ht+zp
   ts=s_tot*t/1000.0_wp*kcaltoau

   g=et-ts

end subroutine thermodyn

subroutine oldthermo(a,b,c,avmom,linear,atom,sym,molmass,vibs,nvibs,escf, &
      & T,sthr,et,ht,g,ts,zp,pr)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   implicit none
   integer, intent(in)  :: nvibs       !< number of vibrational frequencies
   real(wp),intent(in)  :: escf        !< total energy from electronic structure
   real(wp),intent(in)  :: A,B,C       !< rotational constants in cm-1
   real(wp),intent(in)  :: avmom       !< average moment of inertia in whatever
   real(wp),intent(in)  :: sym         !< symmetry number
   real(wp),intent(in)  :: molmass     !< molecular mass in in amu
   real(wp),intent(in)  :: T           !< temperature in K
   real(wp),intent(in)  :: sthr        !< rotor cutoff in cm-1
   real(wp),intent(out) :: et          !< enthalpy in Eh
   real(wp),intent(out) :: ht          !< enthalpy in Eh
   real(wp),intent(out) :: g           !< free energy in Eh
   real(wp),intent(out) :: ts          !< entropy in Eh
   real(wp),intent(out) :: zp          !< zero point vibrational energy in Eh
   real(wp),intent(in)  :: VIBS(nvibs) !< vibrational frequencies in cm-1
   LOGICAL, intent(in)  :: LINEAR      !< is linear
   LOGICAL, intent(in)  :: ATOM        !< only one atom
   LOGICAL, intent(in)  :: pr          !< clutter the screen with printout
   !*******************************************************************
   !*  THE FOLLOWING CONSTANTS ARE NOW DEFINED:
   !*          PI  = CIRCUMFERENCE TO DIAMETER OF A CIRCLE
   !*          R   = GAS CONSTANT IN CALORIES/MOLE
   !*          H   = PLANCK'S CONSTANT IN ERG-SECONDS
   !*          AK  = BOLTZMANN CONSTANT IN ERG/DEGREE
   !*          AC  = SPEED OF LIGHT IN CM/SEC
   !*******************************************************************
   real(wp),parameter :: R = 1.98726D0
   real(wp),parameter :: H = 6.626176D-27
   real(wp),parameter :: AK = 1.3807D-16
   real(wp),parameter :: AC = 2.99792458D+10
   real(wp),parameter :: siK = 1.38066D-23
   real(wp),parameter :: siH = 6.626076D-34
   real(wp),parameter :: siNA = 6.022137D+23
   real(wp),parameter :: conv1 = rcmtoau*autokcal*0.5_wp ! 1.4295718d-3
   real(wp),parameter :: conv2 = rcmtoau*autokJ*1000 ! 0.01196266D+3
   real(wp),parameter :: conv3 = amutokg*1000 ! 1.6606d-24
   real(wp),parameter :: conv4 = 2.2868d0 ! unknown
   real(wp),parameter :: conv5 = 2.3135d0 ! unknown
   real(wp),parameter :: caltoj = autokj/autokcal
   real(wp),parameter :: hbar = H/twopi
   real(wp),parameter :: sihbar = siH/twopi

   integer  :: i
   real(wp) :: s_tr,s_rot,s_vib,s_int,s_tot
   real(wp) :: h_tr,h_rot,h_vib,h_int,h_tot
   real(wp) :: q_tr,q_rot,q_vib,q_int,q_tot
   real(wp) :: cptr,cprot,cpvib,cpint,cptot
   real(wp) :: beta,c1
   real(wp) :: ewj,omega,e,xxmom,xmom
   real(wp) :: sv1,sv2,sv3,wofrot
   !*******************************************************************

   !   ***   INITIALISE SOME VARIABLES   ***
   beta=1.0_wp/kB/T ! beta in 1/Eh
   c1=h*ac/ak/T     ! beta in cm
   q_vib = 1.0_wp
   h_vib = 0.0_wp
   cpvib = 0.0_wp
   s_vib = 0.0_wp
   q_rot = 0.0_wp
   h_rot = 0.0_wp
   cprot = 0.0_wp
   s_rot = 0.0_wp
   zp = 0.5_wp * sum(vibs(1:nvibs))
   !   ***   CONSTRUCT THE FREQUENCY DEPENDENT PARTS OF PARTITION FUNCTION
   do i=1,nvibs
      ! frequency is in cm-1
      omega=vibs(i)*autorcm
      ! omega in cm-1, c1 in cm
      ewj=exp(-omega*c1)
      q_vib=q_vib/(1.0_wp-ewj)
      ! h_vib in cm-1
      h_vib=h_vib+omega*ewj/(1.0_wp-ewj)
      ! cpvib in cm-2
      cpvib=cpvib+omega*omega*ewj/(1.0_wp-ewj)/(1.0_wp-ewj)
      ! replace low-lying vibs for S by rotor approx.
      ! E in J (omega in cm-1)
      e=(1.d-14+omega)*rcmtoau*autokj*1000.0_wp/siNA
      ! mom of intertia corresponding to the rotor with frequency omega
      ! xxmom is in kg·m²
      xxmom=siHbar**2.0d0/(2.0d0*e)
      ! write(*,*) xxmom,avmom,xxmom*avmom/(xxmom+avmom)
      ! this reduced moment limits the rotational moment of
      ! inertia for this vib to that of the total molecule
      ! rotation/3
      ! avmom and xmom are in kg·m²
      xmom=xxmom*avmom/(xxmom+avmom)
      !              free rotor entropy
      ! Cramer, page 328 for one degree of freedom or
      ! http://cccbdb.nist.gov/thermo.asp, eq. 35, sigma=1
      !              harm. osc. entropy
      if(omega.gt.0)then
         ! sv1 is S/R which is dimensionless
         sv1=log(1.0_wp-ewj)
         ! sv2 is in cm-1
         sv2=omega*ewj/(1.0_wp-ewj)
         ! sv3 is cal/mol/K
         sv3=r*(0.5_wp+log(sqrt(8.0_wp*pi**3*xmom*sik*t)/sih))
      else
         sv1=0
         sv2=0
         sv3=0
      endif
      ! fermi weigthing
      ! wofrot=1./(1.+exp( (omega-sthr)/20.0 ) )
      ! Head-Gordon weighting
      wofrot=1.0_wp-chg_switching(omega,sthr)
      if(sthr.lt.0) wofrot=0.0_wp
      ! now convert everything to cal/mol...
      ! sv1 [dimensionless] -> multiply with R    -> cal/mol/K
      ! sv2 [cm-1]          -> multiply with R*c1 -> cal/mol/K
      ! sv3 [cal/mol/K]     ->
      s_vib=s_vib+(1.0_wp-wofrot)*(sv2*r*c1-r*sv1)+wofrot*sv3
      if(omega.lt.300.0_wp .and. pr) &
         & write(*,'("omega=",f6.2,5x,"gvib=",f10.3, &
         &       5x,"grot=",f10.3,5x,"gused=",f10.3)') &
         &  omega,-(sv2*r*c1-r*sv1)*t/1000.0_wp,-sv3*t/1000.0_wp, &
         & -t*((1.0_wp-wofrot)*(sv2*r*c1-r*sv1)+wofrot*sv3)/1000.0_wp
   enddo
   !   ***   FINISH CALCULATION OF VIBRATIONAL PARTS   ***
   ! now unit conversion again...
   ! h*ac/ak is c1*T -> cm*K, R is in cal/mol/K, h_vib in cm-1
   h_vib=h_vib*r*h*ac/ak
   ! same here
   ! cpvib is in cm-2, c1 in cm, R in cal/mol/K
   cpvib=cpvib*r*c1*c1
   !   ***   NOW CALCULATE THE ROTATIONAL PARTS  (FIRST LINEAR MOLECULES
   if (.not.atom) then
      if(linear) then
         ! A is in cm-1, q_rot is dimensionless
         q_rot=1/(c1*a*sym)
         ! h_rot is in cal/mol
         h_rot=r*t
         ! cprot  is in cal/mol/K
         cprot=r
         ! s_rot is in cal/mol/K
         s_rot=r*(log(t*ak/(h*ac*a*sym)))+r
      else
         ! see above
         q_rot=sqrt(pi/(a*b*c*c1*c1*c1))/sym
         h_rot=3.0_wp*r*t/2.0_wp
         cprot=3.0_wp*r/2.0_wp
         s_rot=0.5_wp*r*(3.0_wp*log(1.0_wp/c1)-2.0_wp*log(sym)+log(pi/(a*b*c))+3.0_wp)
      endif
   endif
   !   ***   CALCULATE INTERNAL CONTRIBUTIONS   ***
   q_int=q_vib*q_rot
   h_int=h_vib+h_rot
   cpint=cpvib+cprot
   s_int=s_vib+s_rot
   !   ***   CONSTRUCT TRANSLATION CONTRIBUTIONS   ***
   q_tr=(sqrt(2.0_wp*pi*molmass*t*ak*conv3)/h)**3
   ! this is 3/2rt+pv=5/2rt
   h_tr=5.0_wp*r*t/2.0_wp
   cptr=5.0_wp*r/2.0_wp
   s_tr=conv4*(5.0_wp*log10(t)+3.0_wp*log10(molmass))-conv5
   !   ***   CONSTRUCT TOTALS   ***
   cptot=cptr+cpint
   s_tot=s_tr+s_int
   h_tot=h_tr+h_int

   if(pr)then
      write(*,*)
      write(*,'("   temp. (K)  partition function ", &
         & "  enthalpy   heat capacity  entropy")')
      write(*,'(  "                                   ", &
         & "cal/mol     cal/K/mol   cal/K/mol   J/K/mol")')
      write(*,'(  f7.2,"  VIB ",G10.3,10X,3F11.3)') &
         & T,q_vib,  h_vib,  cpvib,  s_vib
      write(*,'(7X,"  ROT ",G10.3,10X,3F11.3)') &
         & q_rot,  h_rot,  cprot,  s_rot
      write(*,'(7X,"  INT ",G10.3,10X,3F11.3)') &
         & q_int,h_int,cpint,s_int
      write(*,'(7X,"  TR  ",G10.3,10X,3F11.3)') &
         & q_tr, h_tr, cptr, s_tr
      write(*,'(7X,"  TOT ",21X,F11.4,3F11.4)') &
         & h_tot,cptot,s_tot,s_tot*caltoj
   endif

   ht=h_tot/1000.0_wp*kcaltoau
   et=ht+zp
   ts=s_tot*t/1000.0_wp*kcaltoau

   g=et-ts

end subroutine oldthermo

pure elemental function lnqrot(temp,f,avmom) result(lnq_r)
   use xtb_mctc_constants, only : pi
   use xtb_mctc_convert
   implicit none
   real(wp), parameter :: rcmtoj = rcmtoau*autokj*1000.0_wp
   real(wp), parameter :: avogadro = 6.0221413e23_wp ! 1/mol
   real(wp), parameter :: planck = 6.62606957e-34_wp ! J*s
   real(wp), parameter :: hbar = planck/(2.0_wp*pi)
   real(wp), parameter :: kb = 1.3806488e-23_wp ! J/K
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp), intent(in) :: avmom  !< average moment of inertia in kg·m²
   real(wp):: e
   real(wp):: mu
   real(wp):: lnq_r
   real(wp):: t_rot

   ! moment of inertia corresponding to the rotor with frequency f(ifreq)
   ! convert frequency first from cm⁻¹ to J, we add also a little offset to avoid Infinities
   e = ((f*rcmtoJ)+1.0e-14_wp)/avogadro
   mu = hbar**2/(2.0_wp*e)
   ! the vibrational moment of inertia is now in SI
   ! reduce the moment relative to the total rotational moment
   mu = avmom*mu/(avmom+mu)
   ! now we need a rotational temperature of mu,
   ! since we are SI already no unit conversion needed
   t_rot = hbar**2/(2.0_wp*mu*kb)
   ! ln(q) of the free rotor partition function, we assume σ=1
   lnq_r = log(sqrt(pi*temp/t_rot))

end function lnqrot

pure elemental function lnqvib(temp,f) result(lnq_v)
   implicit none
   real(wp), parameter  :: planck = 6.62606957e-34_wp ! J*s
   real(wp), parameter  :: kb = 1.3806488e-23_wp ! J/K
   real(wp), parameter  :: speed_of_light = 299792458.0_wp ! m/s
   real(wp), parameter  :: factor = planck*100.0_wp*speed_of_light/kb
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp):: lnq_v
   real(wp):: t_vib
   ! get the vibrational temperature (which is in K, BTW)
   t_vib = factor*f
   ! modified oscillator, for sthr = 0 -> HO.
   ! ln(q) of the harmonic oscillator partition function
   lnq_v = - 0.5_wp*t_vib/temp - log(1.0_wp - exp(-t_vib/temp))

end function lnqvib

pure elemental function lnqvibmod(temp,f,sthr,avmom) result(lnq)
   implicit none
   real(wp), intent(in) :: temp   !< temperature in K
   real(wp), intent(in) :: f      !< vibrational frequency in cm⁻¹
   real(wp), intent(in) :: sthr   !< rotor cutoff in cm⁻¹
   real(wp), intent(in) :: avmom  !< average moment of inertia in kg·m²
   real(wp):: fswitch
   real(wp):: lnq_r,lnq_v,lnq
   ! ln(q) of the harmonic oscillator partition function
   lnq_v = lnqvib(temp,f)
   ! ln(q) of the free rotor partition function, we assume σ=1
   lnq_r =  lnqrot(temp,f,avmom)

   ! Chai--Head-Gordon weighting
   fswitch = 1.0_wp - chg_switching(sthr,f)

   ! now final modified vibrational partiation function
   lnq = (1.0_wp-fswitch) * lnq_v + fswitch * lnq_r

end function lnqvibmod

pure elemental function chg_switching(omega,sthr) result(f)
   real(wp),intent(in) :: omega
   real(wp),intent(in) :: sthr
   real(wp) :: f
   if(sthr.ge.0.0_wp) then
      f = 1.0_wp/(1.0_wp+(sthr/omega)**4)
   else
      f = 1.0_wp
   endif
end function chg_switching

pure elemental function chg_inverted(f,sthr) result(omega)
   real(wp),intent(in) :: f
   real(wp),intent(in) :: sthr
   real(wp) :: omega
   omega = sthr/(1.0_wp/f - 1.0_wp)**0.25_wp

end function chg_inverted

end module xtb_thermo
