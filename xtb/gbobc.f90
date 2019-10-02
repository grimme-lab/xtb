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

module gbobc
   use iso_fortran_env, wp => real64
   use mctc_constants
   use mctc_econv
   use tbdef_solvent
   implicit none

   public :: lgbsa,lhb,lsalt
   public :: init_gbsa,new_gbsa
   public :: gshift
   public :: ionst,ion_rad
   public :: tb_solvent
   public :: allocate_gbsa,deallocate_gbsa
   public :: compute_brad_sasa
   public :: compute_amat
   public :: compute_fgb
   public :: compute_gb_egrad
   public :: compute_gb_damat
   public :: compute_hb
   public :: update_nnlist_gbsa

   private
! ========================================================================
!  PARAMETER
! ------------------------------------------------------------------------
!  flag for gbsa
   logical :: lgbsa = .false.
!  cutoffs
!  Born radii
   real(wp),parameter :: lrcut_a=35._wp
!  SASA = 2*(w+maxrasasa) + srcut_add
   real(wp),parameter :: srcut_add=2._wp
! ------------------------------------------------------------------------
!  GBOBC PARAMETER
!  offset parameter (fitted)
!  real(wp) :: soset=0.09_wp*aatoau
   real(wp) :: soset
!  van der Waals to Lee-Richard's surface correction
   real(wp), parameter :: alp=1._wp
   real(wp), parameter :: bet=0.8_wp
   real(wp), parameter :: gam=4.85_wp
! ------------------------------------------------------------------------
!  Smoothing dielectric function parameters
   real(wp), parameter :: w=0.3_wp*aatoau
   real(wp), parameter :: w3=w*(w*w)
   real(wp), parameter :: ah0=0.5_wp
   real(wp), parameter :: ah1=3._wp/(4._wp*w)
   real(wp), parameter :: ah3=-1._wp/(4._wp*w3)
! ------------------------------------------------------------------------
!  Angular grid. 38 points lead rot.inv.errors for c60 of 2 kcal
   integer :: nangsa=230
!  real(wp) :: grida(4,nangsa)
!  include 'grida230.inc'
!  integer, parameter :: nangsa=86
!  include 'grida86.inc'
!  integer, parameter :: nangsa=110
!  include 'grida110.inc'
! ------------------------------------------------------------------------
!  real space cut-offs
   real(wp), parameter :: tolsesp=1.d-6
! ------------------------------------------------------------------------
!  Dielectric data
   real(wp) :: epsv
   real(wp) :: epsu
   real(wp) :: keps
! ------------------------------------------------------------------------
!  Surface tension (mN/m=dyn/cm)
   real(wp), parameter :: mNmkcal=4.0305201015221386d-4
   real(wp) :: gammas
   real(wp) :: gamscale(94)
! ------------------------------------------------------------------------
!  Solvent density (g/cm^3) and molar mass (g/mol)
   real(wp), parameter :: molcm3au=8.92388d-2
   real(wp) :: smass
   real(wp) :: rhos
! ------------------------------------------------------------------------
!  Born radii
   real(wp) :: c1
! ------------------------------------------------------------------------
!  Salt screening
   logical :: lsalt=.false.
   real(wp)  :: ionst=0._wp
   real(wp)  :: kappa_const=0.7897d-3
   real(wp)  :: ion_rad=0._wp
   real(wp)  :: kappa=0._wp
! ------------------------------------------------------------------------
!  Atomic surfaces
   real(wp) :: rprobe
   real(wp) :: sasamol
! ------------------------------------------------------------------------
!  Hydrogen bond contribution
   logical :: lhb=.true.

! ------------------------------------------------------------------------
!  Parameters:
! ------------------------------------------------------------------------
!  van der Waals radii
   real(wp) :: rvdw(94)
!  dielectric descreening parameters
   real(wp) :: sx(94)
!  solvent accesible surface radii
   real(wp) :: rasasa(94)
!  HB correction present if zero no HB correction
   integer :: at_hb(94)
!  solvent HB donor or acceptor strength
   real(wp) :: hb_mag(94)
!  Gshift (gsolv=reference vs. gsolv)
   real(wp) :: gshift

   type,private :: gbsa_parameter
!     Dielectric data
      real(wp) :: epsv
!     Solvent density (g/cm^3) and molar mass (g/mol)
      real(wp) :: smass
      real(wp) :: rhos
!     Born radii
      real(wp) :: c1
!     Atomic surfaces
      real(wp) :: rprobe
!     Gshift (gsolv=reference vs. gsolv)
      real(wp) :: gshift
!     offset parameter (fitted)
!     real(wp) :: soset=0.09_wp*aatoau
      real(wp) :: soset
      real(wp) :: dum
!     Surface tension (mN/m=dyn/cm)
      real(wp) :: gamscale(94)
!     dielectric descreening parameters
      real(wp) :: sx(94)
      real(wp) :: tmp(94)
   end type gbsa_parameter


   include 'param_gbsa_acetone.inc'
   include 'param_gbsa_acetonitrile.inc'
   include 'param_gbsa_benzene.inc'
   include 'param_gbsa_ch2cl2.inc'
   include 'param_gbsa_chcl3.inc'
   include 'param_gbsa_cs2.inc'
   include 'param_gbsa_dmso.inc'
   include 'param_gbsa_ether.inc'
   include 'param_gbsa_h2o.inc'
   include 'param_gbsa_methanol.inc'
   include 'param_gbsa_thf.inc'
   include 'param_gbsa_toluene.inc'
   include 'param_gbsa_dmf.inc'
   include 'param_gbsa_nhexan.inc'

contains

subroutine init_gbsa(iunit,n,at,sname,mode,temp,gfn_method,ngrida)
!  use setparam
   use mctc_strings
   use readin
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   character(len=*),intent(in) :: sname
   integer, intent(in) :: mode
   real(wp),intent(in) :: temp
   integer, intent(in) :: gfn_method
   integer, intent(in) :: ngrida

   integer :: i,fix,inum,ich
   real(wp) :: rad
   real(wp) :: gamma_in, rvdwscal, tmp(94), gstate, dum
   character(len=:),allocatable :: fname
   logical ex
   character(len=80) a80
   type(gbsa_parameter) :: gfn_solvent

   !     D3 cut-off radii
   rvdw(1:94)= (/ &
   &1.09155,0.86735,1.7478 ,1.5491 ,1.608  ,1.45515,1.31125,1.24085, &
   &1.1498 ,1.0687 ,1.8541 ,1.74195,2.0053 ,1.89585,1.75085,1.65535, &
   &1.5523 ,1.4574 ,2.12055,2.05175,1.94515,1.8821 ,1.86055,1.7207, &
   &1.7731 ,1.72105,1.71635,1.6731 ,1.6504 ,1.61545,1.97895,1.93095, &
   &1.83125,1.7634 ,1.6831 ,1.6048 ,2.3088 ,2.2382 ,2.1098 ,2.02985, &
   &1.9298 ,1.87715,1.7845 ,1.73115,1.69875,1.67625,1.6654 ,1.731, &
   &2.13115,2.0937 ,2.0075 ,1.94505,1.869  ,1.79445,2.52835,2.5907, &
   &2.31305,2.31005,2.2851 ,2.26355,2.2448 ,2.22575,2.2117 ,2.06215, &
   &2.12135,2.07705,2.1397 ,2.1225 ,2.1104 ,2.0993 ,2.0065 ,2.1225, &
   &2.049  ,1.99275,1.94775,1.8745 ,1.7228 ,1.67625,1.6282 ,1.67995, &
   &2.15635,2.1382 ,2.05875,2.0027 ,1.9322 ,1.8608 ,2.5398 ,2.4647, &
   &2.35215,2.2126 ,2.2297 ,2.19785,2.17695,2.21705/)

   !     hydrogen bonding parameters
   lhb=.false.

   at_hb=0
   at_hb(1)=1
   at_hb(6)=1
   at_hb(7)=1
   at_hb(8)=1
   at_hb(9)=1
   at_hb(15)=1
   at_hb(16)=1
   at_hb(17)=1
   at_hb(34)=1
   at_hb(35)=1
   at_hb(53)=1

   rvdwscal=1.0_wp

   if (gfn_method.gt.1) then
      fname = '.param_gbsa2_'//trim(sname)
   else if (gfn_method.eq.0) then
      fname = '.param_gbsa0_'//trim(sname)
   else
      fname = '.param_gbsa_'//trim(sname)
   endif
   fname = xfind(fname)
   write(iunit,*) 'Solvent             : ', trim(sname)

   inquire(file=fname,exist=ex)
   if(ex)then
      write(iunit,*) 'GBSA parameter file : ', trim(fname)
      open(newunit=ich,file=fname)
      read(ich,*)epsv
      read(ich,*)smass
      read(ich,*)rhos
      read(ich,*)c1
      read(ich,*)rprobe
      read(ich,*)gshift
      read(ich,*)soset
      read(ich,*)dum
      do i=1,94
         read(ich,*)gamscale(i),sx(i),tmp(i)
         if(abs(tmp(i)).gt.1.d-3) lhb=.true.
      enddo
      close(ich)
   else
      !call raise('W','Could not find GBSA parameters in XTBPATH,'//&
      !   ' trying internal parameters')
      if (gfn_method.gt.1.or.gfn_method == 0) then
         select case(lowercase(trim(sname)))
         case default
            call raise('E','solvent : '//trim(sname)//&
               ' not parametrized for GFN2-xTB Hamiltonian',1)
         case('acetone');      gfn_solvent = gfn2_acetone
         case('acetonitrile'); gfn_solvent = gfn2_acetonitrile
!        case('benzene');      gfn_solvent = gfn2_benzene
         case('ch2cl2','dichlormethane'); gfn_solvent = gfn2_ch2cl2
         case('chcl3','chloroform');      gfn_solvent = gfn2_chcl3
         case('cs2');          gfn_solvent = gfn2_cs2
         case('dmso');         gfn_solvent = gfn2_dmso
         case('ether');        gfn_solvent = gfn2_ether
         case('h2o','water');  gfn_solvent = gfn2_h2o
         case('methanol');     gfn_solvent = gfn2_methanol
         case('thf');          gfn_solvent = gfn2_thf
         case('toluene');      gfn_solvent = gfn2_toluene
         case('dmf');          gfn_solvent = gfn2_dmf
         case('nhexan','n-hexan','nhexane','n-hexane');
            gfn_solvent = gfn2_nhexan
         end select
      !else if (gfn_method.eq.0) then
            !call raise('E','solvent : '//trim(sname)//&
               !' not parametrized for GFN0-xTB Hamiltonian',1)
      else
         select case(lowercase(trim(sname)))
         case default
            call raise('E','solvent : '//trim(sname)//&
               ' not parametrized for GFN-xTB Hamiltonian',1)
         case('acetone');      gfn_solvent = gfn1_acetone
         case('acetonitrile'); gfn_solvent = gfn1_acetonitrile
         case('benzene');      gfn_solvent = gfn1_benzene
         case('ch2cl2','dichlormethane'); gfn_solvent = gfn1_ch2cl2
         case('chcl3','chloroform');      gfn_solvent = gfn1_chcl3
         case('cs2');          gfn_solvent = gfn1_cs2
         case('dmso');         gfn_solvent = gfn1_dmso
         case('ether');        gfn_solvent = gfn1_ether
         case('h2o','water');  gfn_solvent = gfn1_h2o
         case('methanol');     gfn_solvent = gfn1_methanol
         case('thf');          gfn_solvent = gfn1_thf
         case('toluene');      gfn_solvent = gfn1_toluene
!        case('dmf');          gfn_solvent = gfn1_dmf
!        case('nhexan');       gfn_solvent = gfn1_nhexan
         end select
      endif
      write(iunit,'(1x,"Using internal parameter file:",1x,a)') fname
      epsv   = gfn_solvent % epsv
      smass  = gfn_solvent % smass
      rhos   = gfn_solvent % rhos
      c1     = gfn_solvent % c1
      rprobe = gfn_solvent % rprobe
      gshift = gfn_solvent % gshift
      soset  = gfn_solvent % soset
      dum    = gfn_solvent % dum
      do i = 1, 94
         gamscale(i) = gfn_solvent % gamscale(i)
         sx(i)       = gfn_solvent % sx(i)
         tmp(i)      = gfn_solvent % tmp(i)
         if(abs(tmp(i)).gt.1.d-3) lhb=.true.
      enddo
   endif

   if(mode.eq.1) then ! gsolv=reference option in COSMOTHERM
      !               RT*(ln(ideal gas mol volume)+ln(rho/M))
      gstate=(temp*8.31451/1000./4.184)* &
      &      (log(24.79_wp*temp/298.15)+ &
      &       log(1000.0_wp*rhos/smass))
      gshift=(gshift+gstate)/autokcal
      write(iunit,*) 'Gsolv state corr. (kcal):',gstate
      a80='gsolv=reference [X=1]'
   elseif(mode.eq.0)then !gsolv option in COSMOTHERM to which it was fitted
      gshift=gshift/autokcal
      a80='gsolv [1 M gas/solution]'
   elseif(mode.eq.2)then ! 1 bar gas/ 1 M solution is not implemented in COSMOTHERM although its the canonical choice
      gstate=(temp*8.31451/1000./4.184)*log(24.79_wp*temp/298.15)
      gshift=(gshift+gstate)/autokcal
      write(iunit,*) 'Gsolv state corr. (kcal):',gstate
      a80='gsolv [1 bar gas/ 1 M solution]'
   endif

!  if(fit)then !penalty to avoid small sx which lead to numerical instabs
!  dum=0
!  do i=1,n
!     dum=dum+2.*(sx(at(i))-0.8)**4
!  enddo
!  gshift=gshift+dum/autokcal
!  endif

!  hydrogen bonding magnitude
   hb_mag = -(tmp**2)/autokcal

!  scaling of the van der Waals radius
   rvdw = rvdw * rvdwscal

!  add the probe radius to the molecular surface
   rasasa=rvdw+rprobe

!  surface tension scaling
   gamma_in=(1.0d-5)*autokcal/mNmkcal

!  dielectric scaling
   epsu=1.0_wp
   keps=((1.0_wp/epsv)-(1.0_wp/epsu))

!  set the salt term
   if(lsalt) then
!     convert to au
      ion_rad=ion_rad*aatoau
!     inverse Debye screening length
      kappa=sqrt(epsv*temp*kappa_const/ionst)*aatoau
      kappa=1.0_wp/kappa
   endif

!  print parameters
   write(iunit,*) 'Gsolv ref. state (COSMO-RS): ',trim(a80)
   write(iunit,*) 'temperature (mdtemp)       : ',temp
   write(iunit,*) 'dielectric constant        : ',epsv
   write(iunit,*) 'rho                        : ',rhos
   write(iunit,*) 'mass                       : ',smass
   write(iunit,*) 'surface tension            : ',gamma_in
   write(iunit,*) 'probe radius               : ',rprobe
   write(iunit,*) 'vdW radii scaling          : ',rvdwscal
   write(iunit,*) 'Gshift (Eh)                : ',gshift
   write(iunit,*) 'c1                         : ',c1
   write(iunit,*) 'soset                      : ',soset
   write(iunit,*) 'HB correction              : ',lhb
   if(lsalt) then
      write(iunit,*) 'Debye screening length     : ',1.0_wp/kappa/aatoau
   endif

   soset=0.1*soset*aatoau

   rhos=rhos*molcm3au/smass

   nangsa = ngrida

end subroutine init_gbsa

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine new_gbsa(this,n,at)
   use grid_module
   implicit none
   type(tb_solvent), intent(inout) :: this

   integer, intent(in) :: n
   integer, intent(in) :: at(n)

   integer i,j,k
   integer ierr
   real(wp) minvdwr
   real(wp) maxrasasa
   real(wp) r
   real(wp), allocatable :: xang(:),yang(:),zang(:),wang(:)

   ! get some space
   call allocate_gbsa(this,n,nangsa)

   ! initialize the vdw radii array
   this%at = at
   this%maxvdwr=0.0_wp
   minvdwr=1000.0_wp
   do i=1,this%nat
      this%vdwr(i)=rvdw(this%at(i))*aatoau
      this%rho(i)=this%vdwr(i)*sx(this%at(i))
      this%svdw(i)=this%vdwr(i)-soset
      this%maxvdwr=max(this%maxvdwr,this%vdwr(i))
      minvdwr=min(minvdwr,this%vdwr(i))
   enddo

   ! nearest-neighbor list preparation
   this%lrcut = lrcut_a*aatoau
   k=0
   do i=1,this%nat
      do j = 1,i-1
         k=k+1
         this%ppind(1,k)=i
         this%ppind(2,k)=j
      enddo
   enddo

   ! initialize solvent-accessible atomic surface area computation (SASA)
   maxrasasa=0.0_wp
   do i = 1, this%nat
      this%vdwsa(i) = rasasa(this%at(i))*aatoau
      maxrasasa=max(maxrasasa,this%vdwsa(i))
      this%trj2(1,i) = (this%vdwsa(i)-w)**2
      this%trj2(2,i) = (this%vdwsa(i)+w)**2
      r=this%vdwsa(i)+w
      this%wrp(i)=(0.25_wp/w+ &
         &            3.0_wp*ah3*(0.2_wp*r*r-0.5*r*this%vdwsa(i)+ &
         &            this%vdwsa(i)*this%vdwsa(i)/3.))*r*r*r
      r=this%vdwsa(i)-w
      this%wrp(i)=this%wrp(i)-(0.25/w+ &
         &    3.0_wp*ah3*(0.2_wp*r*r-0.5*r*this%vdwsa(i)+ &
         &            this%vdwsa(i)*this%vdwsa(i)/3.))*r*r*r
   enddo

   this%srcut = 2.0_wp*(w + maxrasasa) + srcut_add*aatoau
   gammas=(1.0d-5)
   this%sasagam=fourpi*gammas
   do i = 1, this%nat
      this%gamsasa(i)=gamscale(this%at(i))*fourpi*gammas
   enddo

   allocate(xang(this%nang),yang(this%nang),zang(this%nang),wang(this%nang))
   select case(this%nang)
   case(   6);   call ld0006(xang,yang,zang,wang,this%nang)
   case(  14);   call ld0014(xang,yang,zang,wang,this%nang)
   case(  26);   call ld0026(xang,yang,zang,wang,this%nang)
   case(  38);   call ld0038(xang,yang,zang,wang,this%nang)
   case(  50);   call ld0050(xang,yang,zang,wang,this%nang)
   case(  74);   call ld0074(xang,yang,zang,wang,this%nang)
   case(  86);   call ld0086(xang,yang,zang,wang,this%nang)
   case( 110);   call ld0110(xang,yang,zang,wang,this%nang)
   case( 146);   call ld0146(xang,yang,zang,wang,this%nang)
   case( 170);   call ld0170(xang,yang,zang,wang,this%nang)
   case( 194);   call ld0194(xang,yang,zang,wang,this%nang)
   case( 230);   call ld0230(xang,yang,zang,wang,this%nang)
   case( 266);   call ld0266(xang,yang,zang,wang,this%nang)
   case( 302);   call ld0302(xang,yang,zang,wang,this%nang)
   case( 350);   call ld0350(xang,yang,zang,wang,this%nang)
   case( 434);   call ld0434(xang,yang,zang,wang,this%nang)
   case( 590);   call ld0590(xang,yang,zang,wang,this%nang)
   case( 770);   call ld0770(xang,yang,zang,wang,this%nang)
   case( 974);   call ld0974(xang,yang,zang,wang,this%nang)
   case(1202);   call ld1202(xang,yang,zang,wang,this%nang)
   case(1454);   call ld1454(xang,yang,zang,wang,this%nang)
   case(1730);   call ld1730(xang,yang,zang,wang,this%nang)
   case(2030);   call ld2030(xang,yang,zang,wang,this%nang)
   case(2354);   call ld2354(xang,yang,zang,wang,this%nang)
   case(2702);   call ld2702(xang,yang,zang,wang,this%nang)
   case(3074);   call ld3074(xang,yang,zang,wang,this%nang)
   case(3470);   call ld3470(xang,yang,zang,wang,this%nang)
   case(3890);   call ld3890(xang,yang,zang,wang,this%nang)
   case(4334);   call ld4334(xang,yang,zang,wang,this%nang)
   case(4802);   call ld4802(xang,yang,zang,wang,this%nang)
   case(5294);   call ld5294(xang,yang,zang,wang,this%nang)
   case(5810);   call ld5810(xang,yang,zang,wang,this%nang)
   case default; call raise('E',"(gengrid) unknown grid size!",1)
   end select
   this%grida(1,:) = xang
   this%grida(2,:) = yang
   this%grida(3,:) = zang
   this%grida(4,:) = wang
   deallocate(xang,yang,zang,wang)


end subroutine new_gbsa

pure subroutine compute_amat(this,Amat)
   implicit none
   type(tb_solvent),intent(in) :: this

   real(wp),intent(inout) :: Amat(:,:)

   integer  :: i,j,nnj
   integer  :: kk
   real(wp), allocatable :: fhb(:)
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,gg,iepsu,arg,bp
   real(wp) :: dd,expd,fgb,fgb2,dfgb

   if(.not.lsalt) then

      ! compute energy and Amat direct and radii derivatives
      do kk = 1, this%ntpair
         r2 = this%ddpair(1,kk)
         r2 = r2*r2

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb = 1.0_wp/sqrt(fgb2)
         Amat(i,j) = keps*dfgb + Amat(i,j)
         Amat(j,i) = keps*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, this%nat
         bp = 1._wp/this%brad(i)
         Amat(i,i) = Amat(i,i) + keps*bp
      enddo

   else

      iepsu=1.0_wp/epsu

      ! compute energy and Amat direct and radii derivatives
      do kk = 1, this%ntpair
         r2=this%ddpair(1,kk)
         r2=r2*r2

         i=this%ppind(1,kk)
         j=this%ppind(2,kk)

         aa=this%brad(i)*this%brad(j)
         dd=a4*r2/aa
         expd=exp(-dd)
         fgb2=r2+aa*expd
         fgb=sqrt(r2+aa*expd)
         dfgb=1._wp/fgb
         gg=this%ionscr(i)+this%ionscr(j)
         Amat(i,j)=(exp(-kappa*fgb)*gg-iepsu)*dfgb + Amat(i,j)
         Amat(j,i)=(exp(-kappa*fgb)*gg-iepsu)*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, this%nat
         gg=this%ionscr(i)*2.0_wp
         Amat(i,i)= Amat(i,i) + (exp(-kappa*this%brad(i))*gg-iepsu)/this%brad(i)
      enddo


   endif

   ! compute the HB term
   if(lhb) then
      do i = 1, this%nat
         Amat(i,i) = Amat(i,i) + this%hbw(i)
      enddo
   endif

end subroutine compute_amat

pure subroutine compute_fgb(this,fgb,fhb)
   implicit none
   type(tb_solvent),intent(in) :: this

   real(wp),intent(out) :: fgb(this%nat,this%nat)
   real(wp),intent(out) :: fhb(this%nat)

   integer  :: i,j,nnj
   integer  :: kk
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa,r2,gg,iepsu
   real(wp) :: dd,expd,dfgb,hkeps

!  initialize
   fgb=0.0_wp

   hkeps=keps*autoeV

   if(lsalt) then

      iepsu=1.0_wp/epsu

!     compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r2=this%ddpair(1,kk)
         r2=r2*r2

         i=this%ppind(1,kk)
         j=this%ppind(2,kk)

         aa=this%brad(i)*this%brad(j)
         dd=a4*r2/aa
         expd=exp(-dd)
         dfgb=sqrt(r2+aa*expd)
         gg=this%ionscr(i)+this%ionscr(j)
         fgb(i,j)=autoeV*(exp(-kappa*dfgb)*gg-iepsu)/dfgb
         fgb(j,i)=fgb(i,j)
      enddo

!     self-energy part
      do i = 1, this%nat
         gg=this%ionscr(i)*2.0_wp
         fgb(i,i)=autoeV*(exp(-kappa*this%brad(i))*gg-iepsu)/this%brad(i)
      enddo

   else

!     compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r2=this%ddpair(1,kk)
         r2=r2*r2

         i=this%ppind(1,kk)
         j=this%ppind(2,kk)

         aa=this%brad(i)*this%brad(j)
         dd=a4*r2/aa
         expd=exp(-dd)
         dfgb=1.0_wp/(r2+aa*expd)
         fgb(i,j)=hkeps*sqrt(dfgb)
         fgb(j,i)=fgb(i,j)
      enddo

!     self-energy part
      do i = 1, this%nat
         fgb(i,i)=hkeps/this%brad(i)
      enddo

   endif

!  compute the HB term
   if(lhb) then
      fhb=this%hbw*autoeV
   else
      fhb=0.0_wp
   endif

end subroutine compute_fgb

pure subroutine compute_gb_damat(this,q,gborn,ghb,dAmatdr,Afac,lpr)
   implicit none
   type(tb_solvent), intent(in) :: this

   real(wp), intent(in)    :: q(this%nat)
   real(wp), intent(inout) :: dAmatdr(3,this%nat,this%nat)
   real(wp), intent(inout) :: Afac(3,this%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(out)   :: ghb
   logical,  intent(in)    :: lpr

   integer :: i,j,k,nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,fgb,fgb2,br3
   real(wp) :: qq,dd,expd,dfgb,dfgb2,dfgb3,ap,bp,qfg
   real(wp) :: gg,expa,epu,aii,egb
   real(wp) :: r0vdw,r01,r02,ar02
   real(wp) :: grddbi,grddbj
   real(wp) :: dr(3),r

   egb = 0._wp

   if(.not.lsalt) then
      ! GB energy and gradient

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r = this%ddpair(1,kk)
         r2 = r*r

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         ! dielectric scaling of the charges
         qq = q(i)*q(j)*keps
         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*keps

         egb = egb + qq*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*this%ddpair(2:4,kk)
         dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*q(i)
         dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*q(j)
         Afac(:,i) = Afac(:,i) - dr*q(j)
         Afac(:,j) = Afac(:,j) + dr*q(i)

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = this%dbrdp(i)*this%brad(j)*bp
         grddbj = this%dbrdp(j)*this%brad(i)*bp

         dAmatdr(:,:,j) = dAmatdr(:,:,j) + this%brdr(:,:,j)*grddbj*q(i)
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%brdr(:,:,i)*grddbi*q(j)

      enddo

      ! self-energy part
      do i = 1, this%nat
         bp = 1._wp/this%brad(i)
         qq = q(i)*bp
         egb = egb + 0.5_wp*q(i)*qq*keps
         grddbi = -this%dbrdp(i)*keps*bp*bp*0.5_wp
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%brdr(:,:,i)*grddbi*q(i)
      enddo

   else
      ! GB-SE energy and dAmatdr

      epu=1._wp/epsu

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r = this%ddpair(1,kk)
         r2 = r*r

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         qq = q(i)*q(j)
         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         fgb = sqrt(fgb2)
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         aa = kappa*fgb
         expa = exp(-aa)
         gg = (this%ionscr(i)+this%ionscr(j))*expa

         egb = egb + qq*dfgb*(gg-epu)

         dfgb3 = (gg*(1._wp+aa)-epu)*dfgb*dfgb2

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*this%ddpair(2:4,kk)
         dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*q(i)
         dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*q(j)
         Afac(:,i) = Afac(:,i) - dr*q(j)
         Afac(:,j) = Afac(:,j) + dr*q(i)

         qfg = dfgb*expa
         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = this%dbrdp(i)*(this%brad(j)*bp+qfg*this%discr(i))*q(j)
         grddbj = this%dbrdp(j)*(this%brad(i)*bp+qfg*this%discr(j))*q(i)

         dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%brdr(:,:,i)*grddbi
         dAmatdr(:,:,j) = dAmatdr(:,:,j) + this%brdr(:,:,j)*grddbj

      enddo

      ! self-energy part
      do i = 1, this%nat
         gg = exp(-kappa*this%brad(i))
         aa = 2._wp*this%ionscr(i)*gg-epu
         qq = q(i)/this%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-this%brad(i)*2._wp*(this%discr(i)+this%ionscr(i)*kappa)*gg
         grddbi = -this%dbrdp(i)*0.5_wp*qq*ap/this%brad(i)
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%brdr(:,:,i)*grddbi
      enddo

   endif

   gborn = egb

   if(lhb) then
      call compute_ahb(this,q,ghb,dAmatdr)
   endif

!  if(lopt.and.lpr) then
!   write(*,'(/,a)') 'Results GBOBC:'
!   write(*,*) 'At #,  Z , GBOBC (A), RVDW (A)'
!   do i = 1, this%nat
!    write(*,'(I5,2x,I2,6F12.4)') i,at(i),this%brad(i)/aatoau, &
!&   rvdw(at(i)),sx(at(i)),xyz(1:3,i)/aatoau
!   enddo
!   write(*,'(/,a)') 'Free Energy (kcal/mol):'
!   write(*,'(''G-EL  = '',F8.3)') this%gborn*autokcal
!   write(*,'(''GCAV  = '',F8.3)') this%gsasa*autokcal
!   write(*,'(''G-HB  = '',F8.3)') this%ghb*autokcal
!   write(*,'(''GSOL  = '',F8.3)') (this%gborn+this%gsasa+this%ghb)*autokcal
!  endif

end subroutine compute_gb_damat

subroutine compute_gb_egrad(this,q,gborn,ghb,gradient,lpr)
   implicit none
   type(tb_solvent), intent(in) :: this

   real(wp), intent(in)    :: q(this%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: gradient(3,this%nat)
   logical,  intent(in)    :: lpr

   integer :: i,j,k,nnj
   integer :: kk
   real(wp), parameter :: a13=1._wp/3._wp
   real(wp), parameter :: a4=0.25_wp
   real(wp) :: aa,r2,fgb,fgb2,br3
   real(wp) :: qq,dd,expd,dfgb,dfgb2,dfgb3,egb,ap,bp,qfg
   real(wp) :: gg,expa,epu
   real(wp) :: r0vdw,r01,r02,ar02
   real(wp) :: grddbi,grddbj
   real(wp) :: dr(3),r
   real(wp),allocatable :: grddb(:)

   allocate(grddb(this%nat), source = 0.0_wp )

   egb = 0._wp
   grddb = 0._wp

   if(.not.lsalt) then
      ! GB energy and gradient

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r = this%ddpair(1,kk)
         r2 = r*r

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         ! dielectric scaling of the charges
         qq = q(i)*q(j)
         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*keps

         egb = egb + qq*keps*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*this%ddpair(2:4,kk)
         gradient(:,i) = gradient(:,i) - dr*qq
         gradient(:,j) = gradient(:,j) + dr*qq

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = this%dbrdp(i)*this%brad(j)*bp
         grddbj = this%dbrdp(j)*this%brad(i)*bp
         grddb(i) = grddb(i) + grddbi*qq
         !gradient = gradient + this%brdr(:,:,i) * grddbi*qq
         grddb(j) = grddb(j) + grddbj*qq
         !gradient = gradient + this%brdr(:,:,j) * grddbj*qq

      enddo

      ! self-energy part
      do i = 1, this%nat
         bp = 1._wp/this%brad(i)
         qq = q(i)*bp
         egb = egb + 0.5_wp*q(i)*qq*keps
         grddbi = -this%dbrdp(i)*0.5_wp*keps*qq*bp
         grddb(i) = grddb(i) + grddbi*q(i)
         !gradient = gradient + this%brdr(:,:,i) * grddbi*q(i)
      enddo

   else
      ! GB-SE energy and gradient

      epu=1._wp/epsu

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r = this%ddpair(1,kk)
         r2 = r*r

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         qq = q(i)*q(j)
         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         dfgb = r2+aa*expd
         fgb = sqrt(dfgb)
         aa = kappa*fgb
         expa = exp(-aa)
         gg = (this%ionscr(i)+this%ionscr(j))*expa
         qfg = qq/fgb

         egb = egb + qfg*(gg-epu)

         dfgb3 = qfg*(gg*(1._wp+aa)-epu)/dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*this%ddpair(2:4,kk)
         gradient(:,i) = gradient(:,i) - dr
         gradient(:,j) = gradient(:,j) + dr

         qfg = qfg*expa
         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = this%dbrdp(i)*(this%brad(j)*bp+qfg*this%discr(i))
         grddbj = this%dbrdp(j)*(this%brad(i)*bp+qfg*this%discr(j))
         grddb(i) = grddb(i) + grddbi
         grddb(j) = grddb(j) + grddbj

      enddo

      ! self-energy part
      do i = 1, this%nat
         gg = exp(-kappa*this%brad(i))
         aa = 2._wp*this%ionscr(i)*gg-epu
         qq = q(i)/this%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-this%brad(i)*2._wp*(this%discr(i)+this%ionscr(i)*kappa)*gg
         grddbi = -this%dbrdp(i)*0.5_wp*qq*qq*ap
         grddb(i) = grddb(i) + grddbi
      enddo

   endif
   ! contract with the Born radii derivatives
   call dgemv('n',3*this%nat,this%nat,1.0_wp,this%brdr,3*this%nat,grddb,1, &
      &       1.0_wp,gradient,1)

   gborn = egb
   gradient = gradient + this%dsdr

   if(lhb) then
      call compute_hb(this,q,ghb,gradient)
   else
      ghb = 0.0_wp
   endif

!  if(lopt.and.lpr) then
!   write(*,'(/,a)') 'Results GBOBC:'
!   write(*,*) 'At #,  Z , GBOBC (A), RVDW (A)'
!   do i = 1, this%nat
!    write(*,'(I5,2x,I2,6F12.4)') i,at(i),this%brad(i)/aatoau, &
!&   rvdw(at(i)),sx(at(i)),xyz(1:3,i)/aatoau
!   enddo
!   write(*,'(/,a)') 'Free Energy (kcal/mol):'
!   write(*,'(''G-EL  = '',F8.3)') this%gborn*autokcal
!   write(*,'(''GCAV  = '',F8.3)') this%gsasa*autokcal
!   write(*,'(''G-HB  = '',F8.3)') this%ghb*autokcal
!   write(*,'(''GSOL  = '',F8.3)') (this%gborn+this%gsasa+this%ghb)*autokcal
!  endif

end subroutine compute_gb_egrad

pure subroutine compute_fhb(this,xyz)
   implicit none
   type(tb_solvent), intent(inout) :: this

   real(wp), intent(in) :: xyz(3,this%nat)

   integer  :: i
   integer  :: iz,nhb
   real(wp) :: hbed,dhbed
   real(wp) :: smaxd,sasad,sasaw
   real(wp) :: sfw,dsfw,w3,w2,w1
   integer  :: j
   real(wp) :: wbh,wah

   this%hbw=0.0_wp
   this%dhbdw=0.0_wp

   do i = 1, this%nat
      ! atomic Z
      iz=this%at(i)
      ! number of HB
      nhb=at_hb(iz)
      if(nhb.le.0) cycle
      ! SASA-D for HB
      smaxd=1.0_wp/(this%vdwsa(i)*this%vdwsa(i)*this%gamsasa(i))
      sasad=this%sasa(i)*smaxd
      this%hbw(i)=hb_mag(iz)*sasad
      this%dhbdw(i)=hb_mag(iz)*smaxd
   enddo

end subroutine compute_fhb

pure subroutine compute_hb(this,q,ghb,gradient)
   implicit none
   type(tb_solvent), intent(in) :: this

   real(wp), intent(in)    :: q(this%nat)
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: gradient(3,this%nat)

   integer  :: i,j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, this%nat
      qq = q(i)*q(i)
      ghb = ghb + this%hbw(i)*qq
   enddo

   do i = 1, this%nat
      dhbed=this%dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)*q(i)
      do j = 1, this%nat
         gradient(:,j) = gradient(:,j) + this%dsdrt(:,j,i)*dhbed
      enddo
   enddo

end subroutine compute_hb

pure subroutine compute_ahb(this,q,ghb,dAmatdr)
   implicit none
   type(tb_solvent), intent(in) :: this

   real(wp), intent(in)    :: q(this%nat)
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: dAmatdr(3,this%nat,this%nat)

   integer  :: i,j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, this%nat
      qq = q(i)*q(i)
      ghb = ghb + this%hbw(i)*qq
   enddo

   do i = 1, this%nat
      dhbed=this%dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%dsdrt(:,:,i)*dhbed
   enddo

end subroutine compute_ahb

subroutine compute_brad_sasa(this,xyz)
   implicit none
   type(tb_solvent), intent(inout) :: this

   real(wp), intent(in) :: xyz(3,this%nat)

   integer i,j,kk

   this%brad=0.0_wp
   this%dsdr=0.0_wp
   this%dsdrt=0.0_wp
   this%brdr=0.0_wp

   call compute_psi(this)

   call compute_bornr(this)

   ! compute solvent accessible surface and its derivatives
   call compute_numsa(this,xyz)

   ! compute the Debye-Hueckel ion exclusion term
   if (lsalt) call compute_debye_hueckel(this)

   ! compute the HB term
   if (lhb) call compute_fhb(this,xyz)

end subroutine compute_brad_sasa

!> compute the Debye-Hueckel ion exclusion term
pure subroutine compute_debye_hueckel(this)
   implicit none
   type(tb_solvent), intent(inout) :: this

   integer  :: i
   real(wp) :: aa,gg

   aa=0.5_wp/epsv
   do i = 1, this%nat
      gg=kappa*(this%brad(i)+ion_rad)
      this%ionscr(i)=aa*exp(gg)/(1.0_wp+gg)
      this%discr(i)=this%ionscr(i)*kappa*gg/(1.0_wp+gg)
   enddo

end subroutine compute_debye_hueckel

pure subroutine compute_bornr(this)
   implicit none
   type(tb_solvent), intent(inout) :: this

   integer iat
   real(wp) br
   real(wp) dpsi
   real(wp) svdwi,vdwri
   real(wp) s1,v1,s2
   real(wp) arg,arg2,th,ch
   real(wp) alpi,beti,gami

   do iat = 1, this%nat

      br = this%brad(iat)

      svdwi=this%svdw(iat)
      vdwri=this%vdwr(iat)
      s1=1.0_wp/svdwi
      v1=1.0_wp/vdwri
      s2=0.5_wp*svdwi

      br=br*s2

      arg2=br*(gam*br-bet)
      arg=br*(alp+arg2)
      arg2=2.0_wp*arg2+alp+gam*br*br

      th=tanh(arg)
      ch=cosh(arg)

      br=1.0_wp/(s1-v1*th)
      ! Include GBMV2-like scaling
      br=c1*br

      dpsi=ch*(s1-v1*th)
      dpsi=s2*v1*arg2/(dpsi*dpsi)
      dpsi=c1*dpsi

      this%brad(iat) = br
      this%dbrdp(iat) = dpsi

   enddo

end subroutine compute_bornr

pure subroutine compute_psi(this)
   implicit none
   type(tb_solvent), intent(inout) :: this

   real(wp),allocatable :: br(:),brdr(:,:,:),brdrt(:,:)
   integer  :: kk
   integer  :: i,ii,jj,nn
   real(wp) :: dr(3),r,rhoi,rhoj
   real(wp) :: gi,gj,ap,am,lnab,rhab,ab,dgi,dgj
   real(wp) :: drjj(3)
   real(wp) :: rh1,rhr1,r24,rh2,r1,aprh1,r12
   real(wp) :: rvdwi,rvdwj
   integer  :: ovij,ovji,ov

   allocate( br(this%nat),brdr(3,this%nat,this%nat),brdrt(3,this%nat), &
      &      source = 0.0_wp)

   do kk = 1, this%nnrad

      ii=this%nnlistr(1,kk)
      jj=this%nnlistr(2,kk)
      nn=this%nnlistr(3,kk)

      r=this%ddpair(1,nn)
      dr(:)=this%ddpair(2:4,nn)

      rhoi=this%rho(ii)
      rhoj=this%rho(jj)
      rvdwi=this%vdwr(ii)
      rvdwj=this%vdwr(jj)

      ovij=1
      ovji=1
      if(r.ge.(rvdwi+rhoj)) ovij=0
      if(r.ge.(rhoi+rvdwj)) ovji=0
      ov=ovij+10*ovji

      select case(ov)
      case(0) ! ij do not overlap; ji do not overlap
         ! nonoverlaping spheres
         if(abs(rhoi-rhoj).lt.1.d-8) then
            ! equal reduced radii
            r1=1.0_wp/r
            ap=r+rhoj
            am=r-rhoj
            ab=ap*am
            rhab=rhoj/ab
            lnab=0.5_wp*log(am/ap)*r1
            gi=rhab+lnab
            dgi=-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            br(ii)=br(ii)+gi
            br(jj)=br(jj)+gi
            ! accumulate psi gradient
            drjj(:)=dgi*dr(:)
            brdrt(:,ii)=brdrt(:,ii)+drjj(:)
            brdr(:,jj,ii)=brdr(:,jj,ii)-drjj(:)
            brdrt(:,jj)=brdrt(:,jj)-drjj(:)
            brdr(:,ii,jj)=brdr(:,ii,jj)+drjj(:)
         else
            ! unequal reduced radii
            ! ij contribution
            r1=1.0_wp/r
            ap=r+rhoj
            am=r-rhoj
            ab=ap*am
            rhab=rhoj/ab
            lnab=0.5_wp*log(am/ap)*r1
            gi=rhab+lnab
            dgi=-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! ji contribution
            ap=r+rhoi
            am=r-rhoi
            ab=ap*am
            rhab=rhoi/ab
            lnab=0.5_wp*log(am/ap)*r1
            gj=rhab+lnab
            dgj=-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
            ! accumulate psi
            br(ii)=br(ii)+gi
            br(jj)=br(jj)+gj
            ! accumulate psi gradient
            drjj(:)=dgi*dr(:)
            brdrt(:,ii)=brdrt(:,ii)+drjj(:)
            brdr(:,jj,ii)=brdr(:,jj,ii)-drjj(:)

            drjj(:)=dgj*dr(:)
            brdrt(:,jj)=brdrt(:,jj)-drjj(:)
            brdr(:,ii,jj)=brdr(:,ii,jj)+drjj(:)
         endif

      case(10) ! ij do not overlap; ji overlap

         ! ij contribution
         r1=1.0_wp/r
         ap=r+rhoj
         am=r-rhoj
         ab=ap*am
         rhab=rhoj/ab
         lnab=0.5_wp*log(am/ap)*r1
         gi=rhab+lnab
         dgi=-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
         ! accumulate psi
         br(ii)=br(ii)+gi
         ! accumulate psi gradient
         drjj(:)=dgi*dr(:)
         brdrt(:,ii)=brdrt(:,ii)+drjj(:)
         brdr(:,jj,ii)=brdr(:,jj,ii)-drjj(:)

         if((r+rhoi).gt.rvdwj) then
            ! ji contribution
            r1=1.0_wp/r
            r12=0.5_wp*r1
            r24=r12*r12

            ap=r+rhoi
            am=r-rhoi
            rh1=1.0_wp/rvdwj
            rhr1=1.0_wp/ap
            aprh1=ap*rh1
            lnab=log(aprh1)

            gj=rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgj=rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoi*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgj=dgj*r1
            ! accumulate psi
            br(jj)=br(jj)+gj
            ! accumulate psi gradient
            drjj(:)=dgj*dr(:)
            brdrt(:,jj)=brdrt(:,jj)-drjj(:)
            brdr(:,ii,jj)=brdr(:,ii,jj)+drjj(:)
         endif

      case(1) ! ij overlap; ji do not overlap

         if((r+rhoj).gt.rvdwi) then
            ! ij contribution
            r1=1.0_wp/r
            r12=0.5_wp*r1
            r24=r12*r12

            ap=r+rhoj
            am=r-rhoj
            rh1=1.0_wp/rvdwi
            rhr1=1.0_wp/ap
            aprh1=ap*rh1
            lnab=log(aprh1)

            gi=rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgi=rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoj*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgi=dgi*r1
            ! accumulate psi
            br(ii)=br(ii)+gi
            ! accumulate psi gradient
            drjj(:)=dgi*dr(:)
            brdrt(:,ii)=brdrt(:,ii)+drjj(:)
            brdr(:,jj,ii)=brdr(:,jj,ii)-drjj(:)
         endif

         ! ji contribution
         ap=r+rhoi
         am=r-rhoi
         ab=ap*am
         rhab=rhoi/ab
         lnab=0.5_wp*log(am/ap)*r1
         gj=rhab+lnab
         dgj=-2.0_wp*rhab/ab+(rhab-lnab)*r1*r1
         ! accumulate psi
         br(jj)=br(jj)+gj
         ! accumulate psi gradient
         drjj(:)=dgj*dr(:)
         brdrt(:,jj)=brdrt(:,jj)-drjj(:)
         brdr(:,ii,jj)=brdr(:,ii,jj)+drjj(:)

      case(11) ! ij and ji overlap
         ! overlaping spheres
         if((r+rhoj).gt.rvdwi) then
            ! ij contribution
            r1=1.0_wp/r
            r12=0.5_wp*r1
            r24=r12*r12

            ap=r+rhoj
            am=r-rhoj
            rh1=1.0_wp/rvdwi
            rhr1=1.0_wp/ap
            aprh1=ap*rh1
            lnab=log(aprh1)

            gi=rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgi=rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoj*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgi=dgi*r1
            ! accumulate psi
            br(ii)=br(ii)+gi
            ! accumulate psi gradient
            drjj(:)=dgi*dr(:)
            brdrt(:,ii)=brdrt(:,ii)+drjj(:)
            brdr(:,jj,ii)=brdr(:,jj,ii)-drjj(:)
         endif

         if((r+rhoi).gt.rvdwj) then
            ! ji contribution
            r1=1.0_wp/r
            r12=0.5_wp*r1
            r24=r12*r12

            ap=r+rhoi
            am=r-rhoi
            rh1=1.0_wp/rvdwj
            rhr1=1.0_wp/ap
            aprh1=ap*rh1
            lnab=log(aprh1)

            gj=rh1-rhr1+r12*(0.5_wp*am*(rhr1-rh1*aprh1)-lnab)

            dgj=rhr1*rhr1*(1.0_wp-0.25_wp*am*r1*(1.0_wp+aprh1*aprh1))+ &
               &         rhoi*r24*(rhr1-rh1*aprh1)+ &
               &         r12*(r1*lnab-rhr1)
            dgj=dgj*r1
            ! accumulate psi
            br(jj)=br(jj)+gj
            ! accumulate psi gradient
            drjj(:)=dgj*dr(:)
            brdrt(:,jj)=brdrt(:,jj)-drjj(:)
            brdr(:,ii,jj)=brdr(:,ii,jj)+drjj(:)
         endif

      end select

   enddo

   ! save Born radii
   this%brad = br
   ! save derivative of Born radii w.r.t. atomic positions
   this%brdr = brdr
   ! save one-center terms
   do i = 1, this%nat
      this%brdr(:,i,i) = brdrt(:,i)
   enddo

   deallocate( br,brdr,brdrt)

end subroutine compute_psi

pure subroutine compute_numsa(this,xyz)
   implicit none
   type(tb_solvent), intent(inout) :: this

   real(wp), intent(in) :: xyz(3,this%nat)

   integer iat,jat
   integer ip,jj,nnj,nnk
   real(wp) rsas,sasai
   real(wp) xyza(3),xyzp(3),sasap
   real(wp) wr,wsa,drjj(3)
   integer :: nno
   real(wp), allocatable :: grds(:,:)
   real(wp), allocatable :: grads(:,:)
   integer :: nni
   integer, allocatable :: grdi(:)

   allocate(grads(3,this%nat), source = 0.0_wp)

   do iat = 1, this%nat

      rsas = this%vdwsa(iat)
      ! allocate space for the gradient storage
      nno=this%nnsas(iat)
      allocate(grds(3,nno))
      allocate(grdi(nno))

      ! initialize storage
      grads = 0.0_wp
      sasai = 0.0_wp

      ! atomic position
      xyza(:)=xyz(:,iat)
      ! radial atomic weight
      wr=this%wrp(iat)*this%gamsasa(iat)

      ! loop over grid points
      do ip=1,this%nang
         ! grid point position
         xyzp(:) = xyza(:) + rsas*this%grida(1:3,ip)
         ! atomic surface function at the grid point
         call compute_w_sp(this%nat,this%nnlists,this%trj2,this%vdwsa, &
            &              xyz,iat,nno,xyzp,sasap,grds,nni,grdi)

         if(sasap.gt.tolsesp) then
            ! numerical quadrature weight
            wsa = this%grida(4,ip)*wr*sasap
            ! accumulate the surface area
            sasai = sasai + wsa
            ! accumulate the surface gradient
            do jj = 1, nni
               nnj = grdi(jj)
               drjj(:) = wsa*grds(:,jj)
               grads(:,iat)=grads(:,iat)+drjj(:)
               grads(:,nnj)=grads(:,nnj)-drjj(:)
            enddo
         endif
      enddo

      deallocate(grds)
      deallocate(grdi)

      this%sasa(iat) = sasai
      this%dsdrt(:,:,iat) = grads

   enddo

   ! contract surface gradient
   this%dsdr = sum(this%dsdrt, dim=3)

   this%gsasa = sum(this%sasa)

end subroutine compute_numsa

pure subroutine compute_w_sp(nat,nnlists,trj2,vdwsa,xyza,iat,nno,xyzp,sasap,grds, &
      &                      nni,grdi)
   implicit none

   integer, intent(in)  :: nat
   integer, intent(in)  :: nnlists(nat,nat)
   integer, intent(in)  :: iat
   integer, intent(in)  :: nno
   integer, intent(out) :: nni
   real(wp),intent(in)  :: xyza(3,nat)
   real(wp),intent(in)  :: xyzp(3)
   real(wp),intent(out) :: sasap
   real(wp),intent(out) :: grds(3,nno)
   integer, intent(out) :: grdi(nno)
   real(wp),intent(in)  :: trj2(2,nat)
   real(wp),intent(in)  :: vdwsa(nat)

   integer  :: i,ia
   real(wp) :: tj(3),tj2,sqtj
   real(wp) :: uj,uj3,ah3uj2
   real(wp) :: sasaij,dsasaij

   ! initialize storage
   nni=0
   sasap=1.0_wp
   do i = 1, nno
      ia = nnlists(i,iat)
      ! compute the distance to the atom
      tj(:)=xyzp(:)-xyza(:,ia)
      tj2=tj(1)*tj(1)+tj(2)*tj(2)+tj(3)*tj(3)
      ! if within the outer cut-off compute
      if(tj2.lt.trj2(2,ia)) then
         if(tj2.le.trj2(1,ia)) then
            sasap=0.0_wp
            return
         else
            sqtj = sqrt(tj2)
            uj = sqtj - vdwsa(ia)
            ah3uj2 = ah3*uj*uj
            dsasaij = ah1+3.0_wp*ah3uj2
            sasaij =  ah0+(ah1+ah3uj2)*uj

            ! accumulate the molecular surface
            sasap = sasap*sasaij
            ! compute the gradient wrt the neighbor
            dsasaij = dsasaij/(sasaij*sqtj)
            nni=nni+1
            grdi(nni)=ia
            grds(:,nni) = dsasaij*tj(:)
         endif
         ! check if the point is already buried
!        if(sasap.lt.tolsesp) return
      endif
   enddo

end subroutine compute_w_sp

!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
subroutine update_nnlist_gbsa(this,xyz,parallel)
   implicit none
   type(tb_solvent),intent(inout) :: this

   real(wp),intent(in) :: xyz(3,this%nat)
   logical, intent(in) :: parallel

   if (parallel) then
      call update_nnlist_gbsa_parallel(this,xyz)
   else
      call update_nnlist_gbsa_sequential(this,xyz)
   endif

end subroutine update_nnlist_gbsa
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
!  OMP parallel version.
subroutine update_nnlist_gbsa_parallel(this,xyz)
!$ use omp_lib
   implicit none
   type(tb_solvent),intent(inout) :: this

   real(wp),intent(in) :: xyz(3,this%nat)

   integer kk,i1,i2
   real(wp) rcutn2,lrcut2,srcut2
   real(wp) x,y,z,dr2
   integer ip,ip2,thrid,nproc
   integer, allocatable :: npid(:)
   integer, allocatable :: plisttr(:,:,:)
   integer, allocatable :: nntmp(:)
   integer, allocatable :: nnls(:,:)

   nproc=1
!$ nproc=omp_get_max_threads()

   allocate(plisttr(3,this%ntpair,nproc),nnls(this%nat,this%nat))
   allocate(nntmp(this%nat),npid(nproc))
   npid = 0

   lrcut2 = this%lrcut*this%lrcut
   srcut2 = this%srcut*this%srcut

   this%nnsas=0
   this%nnlists=0
!$omp parallel default(none) &
!$omp&         shared ( this,xyz,lrcut2,srcut2 ) &
!$omp&         private( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls ) &
!$omp&         shared ( plisttr, npid )
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   thrid=1
!$ thrid=omp_get_thread_num() + 1
!$omp do
   do kk=1,this%ntpair
      i1=this%ppind(1,kk)
      i2=this%ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      this%ddpair(2,kk)=x
      this%ddpair(3,kk)=y
      this%ddpair(4,kk)=z
      this%ddpair(1,kk)=sqrt(dr2)
      if(dr2.lt.lrcut2) then
         ip = ip + 1
         plisttr(1,ip,thrid)=i1
         plisttr(2,ip,thrid)=i2
         plisttr(3,ip,thrid)=kk
         if(dr2.lt.srcut2) then
            nntmp(i1) = nntmp(i1) + 1
            nntmp(i2) = nntmp(i2) + 1
            nnls(nntmp(i1),i1)=i2
            nnls(nntmp(i2),i2)=i1
         endif
      endif
   enddo
!$omp end do
   npid(thrid)=ip
!$omp critical
   do i1=1,this%nat
      do i2=1,nntmp(i1)
         this%nnlists(this%nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      this%nnsas(i1)=this%nnsas(i1)+nntmp(i1)
   enddo
!$omp end critical
!$omp end parallel

   this%nnrad=0
   do thrid=1,nproc
      do kk = this%nnrad+1,this%nnrad+npid(thrid)
         this%nnlistr(1:3,kk)=plisttr(1:3,kk-this%nnrad,thrid)
      enddo
      this%nnrad = this%nnrad + npid(thrid)
   enddo

   deallocate(nntmp,nnls)

end subroutine update_nnlist_gbsa_parallel
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut
!  Sequential version.
pure subroutine update_nnlist_gbsa_sequential(this,xyz)
   implicit none
   type(tb_solvent),intent(inout) :: this

   real(wp),intent(in) :: xyz(3,this%nat)

   integer kk,i1,i2
   real(wp) rcutn2,lrcut2,srcut2
   real(wp) x,y,z,dr2
   integer ip,ip2
   integer, allocatable :: nntmp(:)
   integer, allocatable :: nnls(:,:)

   allocate(nnls(this%nat,this%nat))
   allocate(nntmp(this%nat))

   lrcut2 = this%lrcut*this%lrcut
   srcut2 = this%srcut*this%srcut

   this%nnsas=0
   this%nnlists=0
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   do kk=1,this%ntpair
      i1=this%ppind(1,kk)
      i2=this%ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      this%ddpair(2,kk)=x
      this%ddpair(3,kk)=y
      this%ddpair(4,kk)=z
      this%ddpair(1,kk)=sqrt(dr2)
      if(dr2.lt.lrcut2) then
         ip = ip + 1
         this%nnlistr(1,ip)=i1
         this%nnlistr(2,ip)=i2
         this%nnlistr(3,ip)=kk
         if(dr2.lt.srcut2) then
            nntmp(i1) = nntmp(i1) + 1
            nntmp(i2) = nntmp(i2) + 1
            nnls(nntmp(i1),i1)=i2
            nnls(nntmp(i2),i2)=i1
         endif
      endif
   enddo
   this%nnrad = ip
   do i1=1,this%nat
      do i2=1,nntmp(i1)
         this%nnlists(this%nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      this%nnsas(i1)=this%nnsas(i1)+nntmp(i1)
   enddo

   deallocate(nntmp,nnls)

end subroutine update_nnlist_gbsa_sequential

pure subroutine update_dist_gbsa(this,xyz)
   implicit none
   type(tb_solvent),intent(inout) :: this

   real(wp),intent(in) :: xyz(3,this%nat)

   integer i1,i2,kk

   do kk = 1, this%ntpair
      i1=this%ppind(1,kk)
      i2=this%ppind(2,kk)
      this%ddpair(2:4,kk)=xyz(1:3,i1)-xyz(1:3,i2)
      this%ddpair(1,kk)=sqrt(this%ddpair(2,kk)**2+ &
      &                      this%ddpair(3,kk)**2+ &
      &                      this%ddpair(4,kk)**2)
   enddo

end subroutine update_dist_gbsa

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

end module gbobc
