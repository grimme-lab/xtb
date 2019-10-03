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
   public :: load_custom_parameters

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
!  real space cut-offs
   real(wp), parameter :: tolsesp=1.d-6
! ------------------------------------------------------------------------
!  Surface tension (mN/m=dyn/cm)
   real(wp), parameter :: mNmkcal=4.0305201015221386d-4
! ------------------------------------------------------------------------
!  Solvent density (g/cm^3) and molar mass (g/mol)
   real(wp), parameter :: molcm3au=8.92388d-2
! ------------------------------------------------------------------------
!  Salt screening
   logical  :: lsalt=.false.
   real(wp) :: ionst=0._wp
   real(wp) :: ion_rad=0._wp
   real(wp), parameter :: kappa_const=0.7897d-3
! ------------------------------------------------------------------------
!  Hydrogen bond contribution
   logical :: lhb=.true.
! ------------------------------------------------------------------------
!  Gshift (gsolv=reference vs. gsolv)
   real(wp) :: gshift

   real(wp), parameter :: d3_cutoff_radii(1:94) = [&
      & 1.09155_wp, 0.86735_wp, 1.74780_wp, 1.54910_wp, &
      & 1.60800_wp, 1.45515_wp, 1.31125_wp, 1.24085_wp, &
      & 1.14980_wp, 1.06870_wp, 1.85410_wp, 1.74195_wp, &
      & 2.00530_wp, 1.89585_wp, 1.75085_wp, 1.65535_wp, &
      & 1.55230_wp, 1.45740_wp, 2.12055_wp, 2.05175_wp, &
      & 1.94515_wp, 1.88210_wp, 1.86055_wp, 1.72070_wp, &
      & 1.77310_wp, 1.72105_wp, 1.71635_wp, 1.67310_wp, &
      & 1.65040_wp, 1.61545_wp, 1.97895_wp, 1.93095_wp, &
      & 1.83125_wp, 1.76340_wp, 1.68310_wp, 1.60480_wp, &
      & 2.30880_wp, 2.23820_wp, 2.10980_wp, 2.02985_wp, &
      & 1.92980_wp, 1.87715_wp, 1.78450_wp, 1.73115_wp, &
      & 1.69875_wp, 1.67625_wp, 1.66540_wp, 1.73100_wp, &
      & 2.13115_wp, 2.09370_wp, 2.00750_wp, 1.94505_wp, &
      & 1.86900_wp, 1.79445_wp, 2.52835_wp, 2.59070_wp, &
      & 2.31305_wp, 2.31005_wp, 2.28510_wp, 2.26355_wp, &
      & 2.24480_wp, 2.22575_wp, 2.21170_wp, 2.06215_wp, &
      & 2.12135_wp, 2.07705_wp, 2.13970_wp, 2.12250_wp, &
      & 2.11040_wp, 2.09930_wp, 2.00650_wp, 2.12250_wp, &
      & 2.04900_wp, 1.99275_wp, 1.94775_wp, 1.87450_wp, &
      & 1.72280_wp, 1.67625_wp, 1.62820_wp, 1.67995_wp, &
      & 2.15635_wp, 2.13820_wp, 2.05875_wp, 2.00270_wp, &
      & 1.93220_wp, 1.86080_wp, 2.53980_wp, 2.46470_wp, &
      & 2.35215_wp, 2.21260_wp, 2.22970_wp, 2.19785_wp, &
      & 2.17695_wp, 2.21705_wp]

   type, private :: gbsa_parameter
!     Dielectric data
      real(wp) :: epsv = 0.0_wp
!     Solvent density (g/cm^3) and molar mass (g/mol)
      real(wp) :: smass = 0.0_wp
      real(wp) :: rhos = 0.0_wp
!     Born radii
      real(wp) :: c1 = 0.0_wp
!     Atomic surfaces
      real(wp) :: rprobe = 0.0_wp
!     Gshift (gsolv=reference vs. gsolv)
      real(wp) :: gshift = 0.0_wp
!     offset parameter (fitted)
      real(wp) :: soset = 0.0_wp
      real(wp) :: dum = 0.0_wp
!     Surface tension (mN/m=dyn/cm)
      real(wp) :: gamscale(94) = 0.0_wp
!     dielectric descreening parameters
      real(wp) :: sx(94) = 0.0_wp
      real(wp) :: tmp(94) = 0.0_wp
   end type gbsa_parameter

   type, extends(gbsa_parameter) :: gbsa_model
      integer :: mode = -1
      real(wp) :: temp = -1.0_wp
      !> Dielectric data
      real(wp) :: epsu = 1.0_wp
      real(wp) :: keps = 0.0_wp
      !> Surface tension (mN/m=dyn/cm)
      real(wp) :: gammas = 1.0e-5_wp
      !> Atomic surfaces
      real(wp) :: sasamol = 0.0_wp
      !> Hydrogen bond contribution
      logical :: lhb = .false.
      !> van der Waals radii
      real(wp) :: rvdw(94) = 0.0_wp
      !> solvent accesible surface radii
      real(wp) :: rasasa(94) = 0.0_wp
      !> HB correction present if zero no HB correction
      integer :: at_hb(94) = 0
      !> solvent HB donor or acceptor strength
      real(wp) :: hb_mag(94) = 0.0_wp
      !> Salt screening
      logical :: lsalt = .false.
      real(wp) :: ionst = 0._wp
      real(wp) :: ion_rad = 0._wp
      real(wp) :: kappa = 0._wp
      !> SASA grid size
      integer :: nangsa = 230
   end type gbsa_model

   type(gbsa_model), private :: gbm
   type(gbsa_parameter), private :: custom_solvent

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

subroutine load_custom_parameters(epsv,smass,rhos,c1,rprobe,gshift,soset,dum, &
      &                           gamscale,sx,tmp)
   implicit none
   !> Dielectric data
   real(wp), intent(in), optional :: epsv
   !> Solvent density (g/cm^3) and molar mass (g/mol)
   real(wp), intent(in), optional :: smass
   real(wp), intent(in), optional :: rhos
   !> Born radii
   real(wp), intent(in), optional :: c1
   !> Atomic surfaces
   real(wp), intent(in), optional :: rprobe
   !> Gshift (gsolv=reference vs. gsolv)
   real(wp), intent(in), optional :: gshift
   !> offset parameter (fitted)
   real(wp), intent(in), optional :: soset
   real(wp), intent(in), optional :: dum
   !> Surface tension (mN/m=dyn/cm)
   real(wp), intent(in), optional :: gamscale(94)
   !> dielectric descreening parameters
   real(wp), intent(in), optional :: sx(94)
   real(wp), intent(in), optional :: tmp(94)

   if (present(epsv))    custom_solvent%epsv     = epsv
   if (present(smass))   custom_solvent%smass    = smass
   if (present(rhos))    custom_solvent%rhos     = rhos
   if (present(c1))      custom_solvent%c1       = c1
   if (present(rprobe))  custom_solvent%rprobe   = rprobe
   if (present(gshift))  custom_solvent%gshift   = gshift
   if (present(soset))   custom_solvent%soset    = soset
   if (present(dum))     custom_solvent%dum      = dum
   if (present(gamscale))custom_solvent%gamscale = gamscale
   if (present(sx))      custom_solvent%sx       = sx
   if (present(tmp))     custom_solvent%tmp      = tmp

end subroutine load_custom_parameters

subroutine init_gbsa(iunit,sname,mode,temp,gfn_method,ngrida)
   use mctc_strings
   use readin
   implicit none
   integer, intent(in) :: iunit
   character(len=*),intent(in) :: sname
   integer, intent(in) :: mode
   real(wp),intent(in) :: temp
   integer, intent(in) :: gfn_method
   integer, intent(in) :: ngrida

   integer :: i,fix,inum,ich
   real(wp) :: rad
   real(wp) :: gamma_in, tmp(94), gstate, dum
   character(len=:),allocatable :: fname
   logical ex
   character(len=80) a80
   type(gbsa_parameter) :: gfn_solvent

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
      call read_gbsa_parameters(ich, gfn_solvent)
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
         case('custom'); gfn_solvent = custom_solvent
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
         case('custom'); gfn_solvent = custom_solvent
         end select
      endif
      write(iunit,'(1x,"Using internal parameter file:",1x,a)') fname
   endif

   call new_gbsa_model(gbm,gfn_solvent,mode,temp,ngrida)

   lhb = gbm%lhb
   gshift = gbm%gshift

   call gbsa_info(iunit,gbm)

end subroutine init_gbsa

subroutine gbsa_info(iunit,gbm)
   implicit none
   integer, intent(in) :: iunit
   type(gbsa_model), intent(in) :: gbm

   write(iunit,'(1x,a,1x)',advance='no') 'Gsolv ref. state (COSMO-RS):'
   select case(gbm%mode)
   case(1)
      write(iunit,'(a)') 'gsolv=reference [X=1]'
   case(0)
      write(iunit,'(a)') 'gsolv [1 M gas/solution]'
   case(2)
      write(iunit,'(a)') 'gsolv [1 bar gas/ 1 M solution]'
   case default
      write(iunit,'(i0)') gbm%mode
   end select
   write(iunit,*) 'temperature (mdtemp)       : ',gbm%temp
   write(iunit,*) 'dielectric constant        : ',gbm%epsv
   write(iunit,*) 'rho                        : ',gbm%rhos / (molcm3au/gbm%smass)
   write(iunit,*) 'mass                       : ',gbm%smass
   write(iunit,*) 'surface tension            : ',(1.0d-5)*autokcal/mNmkcal
   write(iunit,*) 'probe radius               : ',gbm%rprobe
   write(iunit,*) 'Gshift (Eh)                : ',gbm%gshift
   write(iunit,*) 'c1                         : ',gbm%c1
   write(iunit,*) 'soset                      : ',gbm%soset / (0.1_wp*aatoau)
   write(iunit,*) 'HB correction              : ',gbm%lhb
   if(gbm%lsalt) then
      write(iunit,*) 'Debye screening length     : ',1.0_wp/gbm%kappa/aatoau
   endif
end subroutine gbsa_info

subroutine new_gbsa_model(gbm,solvent,mode,temp,ngrida)
   use mctc_strings
   use readin
   implicit none
   type(gbsa_model), intent(inout) :: gbm
   type(gbsa_parameter), intent(inout) :: solvent
   integer, intent(in) :: mode
   real(wp),intent(in) :: temp
   integer, intent(in) :: ngrida

   integer :: i,fix,inum,ich
   real(wp) :: rad
   real(wp) :: gamma_in, rvdwscal, tmp(94), gstate, dum
   character(len=:),allocatable :: fname
   logical ex
   character(len=80) a80
   type(gbsa_parameter) :: gfn_solvent

   gbm%lsalt = lsalt

   gbm%temp = temp
   gbm%mode = mode

   ! D3 cut-off radii
   gbm%rvdw(1:94) = d3_cutoff_radii

   ! hydrogen bonding parameters
   gbm%lhb = .false.

   gbm%at_hb = 0
   gbm%at_hb([1,6,7,8,9,15,16,17,34,35,53]) = 1

   rvdwscal=1.0_wp

   gbm%epsv   = solvent%epsv
   gbm%smass  = solvent%smass
   gbm%rhos   = solvent%rhos*molcm3au/gbm%smass
   gbm%c1     = solvent%c1
   gbm%rprobe = solvent%rprobe
   gbm%gshift = solvent%gshift
   gbm%soset  = solvent%soset*0.1_wp*aatoau
   gbm%dum    = solvent%dum
   do i = 1, 94
      gbm%gamscale(i) = solvent%gamscale(i)
      gbm%sx(i)       = solvent%sx(i)
      tmp(i)      = solvent%tmp(i)
      if(abs(tmp(i)).gt.1.d-3) gbm%lhb=.true.
   enddo

   if(mode.eq.1) then ! gsolv=reference option in COSMOTHERM
      !               RT*(ln(ideal gas mol volume)+ln(rho/M))
      gstate=(gbm%temp*8.31451/1000./4.184)* &
      &      (log(24.79_wp*gbm%temp/298.15)+ &
      &       log(1000.0_wp*gbm%rhos/gbm%smass))
      gbm%gshift=(gbm%gshift+gstate)*kcaltoau
   elseif(mode.eq.0)then !gsolv option in COSMOTHERM to which it was fitted
      gbm%gshift=gbm%gshift*kcaltoau
   elseif(mode.eq.2)then ! 1 bar gas/ 1 M solution is not implemented in COSMOTHERM although its the canonical choice
      gstate=(gbm%temp*8.31451/1000./4.184)*log(24.79_wp*gbm%temp/298.15)
      gbm%gshift=(gbm%gshift+gstate)*kcaltoau
   endif

!  if(fit)then !penalty to avoid small sx which lead to numerical instabs
!  dum=0
!  do i=1,n
!     dum=dum+2.*(gbm%sx(at(i))-0.8)**4
!  enddo
!  gbm%gshift=gbm%gshift+dum/autokcal
!  endif

!  hydrogen bonding magnitude
   gbm%hb_mag = -(tmp**2)*kcaltoau

!  scaling of the van der Waals radius
   gbm%rvdw = gbm%rvdw * rvdwscal

!  add the probe radius to the molecular surface
   gbm%rasasa=gbm%rvdw+gbm%rprobe

!  surface tension scaling
   gamma_in=(1.0d-5)*autokcal/mNmkcal

!  dielectric scaling
   gbm%epsu=1.0_wp
   gbm%keps=((1.0_wp/gbm%epsv)-(1.0_wp/gbm%epsu))

!  set the salt term
   if(gbm%lsalt) then
!     convert to au
      gbm%ion_rad=ion_rad*aatoau
!     inverse Debye screening length
      gbm%kappa=sqrt(gbm%epsv*gbm%temp*kappa_const/ionst)*aatoau
      gbm%kappa=1.0_wp/gbm%kappa
   endif

   gbm%nangsa = ngrida

end subroutine new_gbsa_model

subroutine read_gbsa_parameters(ifile, param)
   integer, intent(in) :: ifile
   type(gbsa_parameter) :: param
   integer :: i
   read(ifile,*) param%epsv
   read(ifile,*) param%smass
   read(ifile,*) param%rhos
   read(ifile,*) param%c1
   read(ifile,*) param%rprobe
   read(ifile,*) param%gshift
   read(ifile,*) param%soset
   read(ifile,*) param%dum
   do i = 1, 94
      read(ifile,*) param%gamscale(i), param%sx(i), param%tmp(i)
   enddo
end subroutine read_gbsa_parameters

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
   call allocate_gbsa(this,n,gbm%nangsa)

   ! initialize the vdw radii array
   this%at = at
   this%maxvdwr=0.0_wp
   minvdwr=1000.0_wp
   do i=1,this%nat
      this%vdwr(i)=gbm%rvdw(this%at(i))*aatoau
      this%rho(i)=this%vdwr(i)*gbm%sx(this%at(i))
      this%svdw(i)=this%vdwr(i)-gbm%soset
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
      this%vdwsa(i) = gbm%rasasa(this%at(i))*aatoau
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
   this%sasagam=fourpi*gbm%gammas
   do i = 1, this%nat
      this%gamsasa(i)=gbm%gamscale(this%at(i))*fourpi*gbm%gammas
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

   if(.not.gbm%lsalt) then

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
         Amat(i,j) = gbm%keps*dfgb + Amat(i,j)
         Amat(j,i) = gbm%keps*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, this%nat
         bp = 1._wp/this%brad(i)
         Amat(i,i) = Amat(i,i) + gbm%keps*bp
      enddo

   else

      iepsu=1.0_wp/gbm%epsu

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
         Amat(i,j)=(exp(-gbm%kappa*fgb)*gg-iepsu)*dfgb + Amat(i,j)
         Amat(j,i)=(exp(-gbm%kappa*fgb)*gg-iepsu)*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, this%nat
         gg=this%ionscr(i)*2.0_wp
         Amat(i,i)= Amat(i,i) + (exp(-gbm%kappa*this%brad(i))*gg-iepsu)/this%brad(i)
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

   hkeps=gbm%keps*autoeV

   if(gbm%lsalt) then

      iepsu=1.0_wp/gbm%epsu

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
         fgb(i,j)=autoeV*(exp(-gbm%kappa*dfgb)*gg-iepsu)/dfgb
         fgb(j,i)=fgb(i,j)
      enddo

!     self-energy part
      do i = 1, this%nat
         gg=this%ionscr(i)*2.0_wp
         fgb(i,i)=autoeV*(exp(-gbm%kappa*this%brad(i))*gg-iepsu)/this%brad(i)
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

   if(.not.gbm%lsalt) then
      ! GB energy and gradient

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, this%ntpair
         r = this%ddpair(1,kk)
         r2 = r*r

         i = this%ppind(1,kk)
         j = this%ppind(2,kk)

         ! dielectric scaling of the charges
         qq = q(i)*q(j)*gbm%keps
         aa = this%brad(i)*this%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*gbm%keps

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
         egb = egb + 0.5_wp*q(i)*qq*gbm%keps
         grddbi = -this%dbrdp(i)*gbm%keps*bp*bp*0.5_wp
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + this%brdr(:,:,i)*grddbi*q(i)
      enddo

   else
      ! GB-SE energy and dAmatdr

      epu=1._wp/gbm%epsu

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
         aa = gbm%kappa*fgb
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
         gg = exp(-gbm%kappa*this%brad(i))
         aa = 2._wp*this%ionscr(i)*gg-epu
         qq = q(i)/this%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-this%brad(i)*2._wp*(this%discr(i)+this%ionscr(i)*gbm%kappa)*gg
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

   if(.not.gbm%lsalt) then
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
         dfgb3 = dfgb2*dfgb*gbm%keps

         egb = egb + qq*gbm%keps*dfgb

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
         egb = egb + 0.5_wp*q(i)*qq*gbm%keps
         grddbi = -this%dbrdp(i)*0.5_wp*gbm%keps*qq*bp
         grddb(i) = grddb(i) + grddbi*q(i)
         !gradient = gradient + this%brdr(:,:,i) * grddbi*q(i)
      enddo

   else
      ! GB-SE energy and gradient

      epu=1._wp/gbm%epsu

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
         aa = gbm%kappa*fgb
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
         gg = exp(-gbm%kappa*this%brad(i))
         aa = 2._wp*this%ionscr(i)*gg-epu
         qq = q(i)/this%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-this%brad(i)*2._wp*(this%discr(i)+this%ionscr(i)*gbm%kappa)*gg
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
      nhb=gbm%at_hb(iz)
      if(nhb.le.0) cycle
      ! SASA-D for HB
      smaxd=1.0_wp/(this%vdwsa(i)*this%vdwsa(i)*this%gamsasa(i))
      sasad=this%sasa(i)*smaxd
      this%hbw(i)=gbm%hb_mag(iz)*sasad
      this%dhbdw(i)=gbm%hb_mag(iz)*smaxd
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
   if (gbm%lsalt) call compute_debye_hueckel(this)

   ! compute the HB term
   if (lhb) call compute_fhb(this,xyz)

end subroutine compute_brad_sasa

!> compute the Debye-Hueckel ion exclusion term
pure subroutine compute_debye_hueckel(this)
   implicit none
   type(tb_solvent), intent(inout) :: this

   integer  :: i
   real(wp) :: aa,gg

   aa=0.5_wp/gbm%epsv
   do i = 1, this%nat
      gg=gbm%kappa*(this%brad(i)+gbm%ion_rad)
      this%ionscr(i)=aa*exp(gg)/(1.0_wp+gg)
      this%discr(i)=this%ionscr(i)*gbm%kappa*gg/(1.0_wp+gg)
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
      br=gbm%c1*br

      dpsi=ch*(s1-v1*th)
      dpsi=s2*v1*arg2/(dpsi*dpsi)
      dpsi=gbm%c1*dpsi

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
