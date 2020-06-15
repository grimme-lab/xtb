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

module xtb_solv_gbobc
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants
   use xtb_mctc_convert
   use xtb_mctc_blas, only : mctc_gemv
   use xtb_solv_kernel, only : gbKernel
   use xtb_type_solvent
   implicit none

   public :: lhb,lsalt
   public :: initGBSA,new_gbsa
   public :: gshift
   public :: ionst,ion_rad
   public :: TSolvent
   public :: allocate_gbsa,deallocate_gbsa
   public :: compute_brad_sasa
   public :: compute_amat
   public :: compute_gb_egrad
   public :: compute_gb_damat
   public :: compute_hb_egrad
   public :: update_nnlist_gbsa
   public :: load_custom_parameters
   private


! ========================================================================
!  PARAMETER
! ------------------------------------------------------------------------
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
      real(wp) :: alpha = 0.0_wp
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
      !> β = 1 / ε, we save αβ because they are always used together
      real(wp) :: alpbet = 0.0_wp
      !> Interaction kernel
      integer :: kernel = gbKernel%still
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

   !> P16 zeta parameter
   real(wp), parameter :: zetaP16 = 1.028_wp

   !> P16 zeta parameter over 16
   real(wp), parameter :: zetaP16o16 = zetaP16 / 16.0_wp


contains

subroutine load_custom_parameters(epsv,smass,rhos,c1,rprobe,gshift,soset,alpha, &
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
   real(wp), intent(in), optional :: alpha
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
   if (present(alpha))   custom_solvent%alpha    = alpha
   if (present(gamscale))custom_solvent%gamscale = gamscale
   if (present(sx))      custom_solvent%sx       = sx
   if (present(tmp))     custom_solvent%tmp      = tmp

end subroutine load_custom_parameters

subroutine initGBSA(env,sname,mode,temp,gfn_method,ngrida,alpb,kernel,verbose)
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=*), parameter :: source = 'solv_gbobc_initGBSA'
   type(TEnvironment), intent(inout) :: env
   character(len=*),intent(in) :: sname
   integer, intent(in) :: mode
   real(wp),intent(in) :: temp
   integer, intent(in) :: gfn_method
   integer, intent(in) :: ngrida
   logical, intent(in) :: verbose
   logical, intent(in) :: alpb
   integer, intent(in) :: kernel

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
   if (verbose) then
      write(env%unit,*) 'Solvent             : ', trim(sname)
   end if

   inquire(file=fname,exist=ex)
   if(ex)then
      if (verbose) then
         write(env%unit,*) 'GBSA parameter file : ', trim(fname)
      end if
      open(newunit=ich,file=fname)
      call read_gbsa_parameters(ich, gfn_solvent)
      close(ich)
   else
      !call env%warning('Could not find GBSA parameters in XTBPATH,'//&
      !   ' trying internal parameters', source)
      if (gfn_method.gt.1.or.gfn_method == 0) then
         select case(lowercase(trim(sname)))
         case default
            call env%error('solvent : '//trim(sname)//&
               ' not parametrized for GFN2-xTB Hamiltonian', source)
            return
         case('acetone');      gfn_solvent = gfn2_acetone
         case('acetonitrile'); gfn_solvent = gfn2_acetonitrile
         case('benzene');      gfn_solvent = gfn2_benzene
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
            !call env%error('solvent : '//trim(sname)//&
               !' not parametrized for GFN0-xTB Hamiltonian',source)
      else
         select case(lowercase(trim(sname)))
         case default
            call env%error('solvent : '//trim(sname)//&
               ' not parametrized for GFN-xTB Hamiltonian', source)
            return
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
      if (verbose) then
         write(env%unit,'(1x,"Using internal parameter file:",1x,a)') fname
      end if
   endif

   if (alpb) gfn_solvent%alpha = 0.571412_wp

   call new_gbsa_model(gbm,gfn_solvent,mode,temp,ngrida)

   gbm%kernel = kernel

   lhb = gbm%lhb
   gshift = gbm%gshift

   if (verbose) then
      call gbsa_info(env%unit,gbm)
   end if

end subroutine initGBSA

subroutine gbsa_info(iunit,gbm)
   implicit none
   integer, intent(in) :: iunit
   type(gbsa_model), intent(in) :: gbm

   if (gbm%alpbet > 0.0_wp) then
      write(iunit,*) 'solvation model            : ALPB'
   else
      write(iunit,*) 'solvation model            : GBSA'
   end if
   select case(gbm%kernel)
   case(gbKernel%still)
      write(iunit,*) 'interacton kernel          : Still'
   case(gbKernel%p16)
      write(iunit,*) 'interacton kernel          : P16'
   end select
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
   use xtb_mctc_strings
   use xtb_readin
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
   gbm%alpha  = solvent%alpha
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
   ! ALPB asymmetric correction
   if (gbm%alpha > 0.0_wp) then
      gbm%alpbet = gbm%alpha / gbm%epsv
   else
      gbm%alpbet = 0.0_wp
   end if

   gbm%epsu=1.0_wp
   gbm%keps=((1.0_wp/gbm%epsv)-(1.0_wp/gbm%epsu)) / (1.0_wp + gbm%alpbet)

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
   read(ifile,*) param%alpha
   do i = 1, 94
      read(ifile,*) param%gamscale(i), param%sx(i), param%tmp(i)
   enddo
end subroutine read_gbsa_parameters

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine new_gbsa(self,n,at)
   use xtb_solv_lebedev
   implicit none
   type(TSolvent), intent(inout) :: self

   integer, intent(in) :: n
   integer, intent(in) :: at(n)

   integer i,j,k,iAng
   integer ierr
   real(wp) minvdwr
   real(wp) maxrasasa
   real(wp) r
   real(wp), allocatable :: xang(:),yang(:),zang(:),wang(:)

   ! get some space
   call allocate_gbsa(self,n,gbm%nangsa)

   ! initialize the vdw radii array
   self%at = at
   self%maxvdwr=0.0_wp
   minvdwr=1000.0_wp
   do i=1,self%nat
      self%vdwr(i)=gbm%rvdw(self%at(i))*aatoau
      self%rho(i)=self%vdwr(i)*gbm%sx(self%at(i))
      self%svdw(i)=self%vdwr(i)-gbm%soset
      self%maxvdwr=max(self%maxvdwr,self%vdwr(i))
      minvdwr=min(minvdwr,self%vdwr(i))
   enddo

   ! nearest-neighbor list preparation
   self%lrcut = lrcut_a*aatoau
   k=0
   do i=1,self%nat
      do j = 1,i-1
         k=k+1
         self%ppind(1,k)=i
         self%ppind(2,k)=j
      enddo
   enddo

   ! initialize solvent-accessible atomic surface area computation (SASA)
   maxrasasa=0.0_wp
   do i = 1, self%nat
      self%vdwsa(i) = gbm%rasasa(self%at(i))*aatoau
      maxrasasa=max(maxrasasa,self%vdwsa(i))
      self%trj2(1,i) = (self%vdwsa(i)-w)**2
      self%trj2(2,i) = (self%vdwsa(i)+w)**2
      r=self%vdwsa(i)+w
      self%wrp(i)=(0.25_wp/w+ &
         &            3.0_wp*ah3*(0.2_wp*r*r-0.5*r*self%vdwsa(i)+ &
         &            self%vdwsa(i)*self%vdwsa(i)/3.))*r*r*r
      r=self%vdwsa(i)-w
      self%wrp(i)=self%wrp(i)-(0.25/w+ &
         &    3.0_wp*ah3*(0.2_wp*r*r-0.5*r*self%vdwsa(i)+ &
         &            self%vdwsa(i)*self%vdwsa(i)/3.))*r*r*r
   enddo

   self%srcut = 2.0_wp*(w + maxrasasa) + srcut_add*aatoau
   self%sasagam=fourpi*gbm%gammas
   do i = 1, self%nat
      self%gamsasa(i)=gbm%gamscale(self%at(i))*fourpi*gbm%gammas
   enddo

   iAng = 0
   do i = 1, size(gridSize)
      if (self%nang == gridSize(i)) iAng = i
   end do
   call getAngGrid(iAng, self%angGrid, self%angWeight, ierr)

end subroutine new_gbsa

subroutine compute_amat(self,Amat)
   implicit none
   type(TSolvent),intent(in) :: self

   real(wp),intent(inout) :: Amat(:,:)

   integer  :: i,j,nnj
   integer  :: kk
   real(wp), allocatable :: fhb(:)
   real(wp), parameter :: a13=1.0_wp/3.0_wp
   real(wp), parameter :: a4=0.25_wp
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp) :: aa,r2,gg,iepsu,arg,bp
   real(wp) :: dd,expd,fgb,fgb2,dfgb

   select case(gbm%kernel)
   case(gbKernel%still)
      call compute_amat_still(self, Amat)
   case(gbKernel%p16)
      call compute_amat_p16(self, Amat)
   end select

   ! compute the HB term
   if(lhb) then
      do i = 1, self%nat
         Amat(i,i) = Amat(i,i) + 2*self%hbw(i)
      enddo
   endif

   ! ALPB shape dependent correction for charged systems
   if (gbm%alpbet > 0.0_wp) then
      Amat(:, :) = Amat + gbm%keps * gbm%alpbet / self%aDet
   end if

end subroutine compute_amat

subroutine compute_amat_still(self, Amat)
   implicit none
   type(TSolvent),intent(in) :: self

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
      do kk = 1, self%ntpair
         r2 = self%ddpair(1,kk)
         r2 = r2*r2

         i = self%ppind(1,kk)
         j = self%ppind(2,kk)

         aa = self%brad(i)*self%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb = 1.0_wp/sqrt(fgb2)
         Amat(i,j) = gbm%keps*dfgb + Amat(i,j)
         Amat(j,i) = gbm%keps*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, self%nat
         bp = 1._wp/self%brad(i)
         Amat(i,i) = Amat(i,i) + gbm%keps*bp
      enddo

   else

      iepsu=1.0_wp/gbm%epsu

      ! compute energy and Amat direct and radii derivatives
      do kk = 1, self%ntpair
         r2=self%ddpair(1,kk)
         r2=r2*r2

         i=self%ppind(1,kk)
         j=self%ppind(2,kk)

         aa=self%brad(i)*self%brad(j)
         dd=a4*r2/aa
         expd=exp(-dd)
         fgb2=r2+aa*expd
         fgb=sqrt(r2+aa*expd)
         dfgb=1._wp/fgb
         gg=self%ionscr(i)+self%ionscr(j)
         Amat(i,j)=(exp(-gbm%kappa*fgb)*gg-iepsu)*dfgb + Amat(i,j)
         Amat(j,i)=(exp(-gbm%kappa*fgb)*gg-iepsu)*dfgb + Amat(j,i)
      enddo

      ! self-energy part
      do i = 1, self%nat
         gg=self%ionscr(i)*2.0_wp
         Amat(i,i)= Amat(i,i) + (exp(-gbm%kappa*self%brad(i))*gg-iepsu)/self%brad(i)
      enddo


   endif

end subroutine compute_amat_still

subroutine compute_amat_p16(self, Amat)
   implicit none
   type(TSolvent),intent(in) :: self

   real(wp),intent(inout) :: Amat(:,:)

   integer :: kk
   integer :: iat, jat
   real(wp) :: r1, ab, arg, eab, fgb, dfgb, bp

   ! compute energy and Amat direct and radii derivatives
   !$omp parallel do default(none) shared(Amat, self, gbm) &
   !$omp private(kk, iat, jat, r1, ab, arg, fgb, dfgb)
   do kk = 1, self%ntpair
      r1 = self%ddpair(1,kk)

      iat = self%ppind(1,kk)
      jat = self%ppind(2,kk)

      ab = sqrt(self%brad(iat) * self%brad(jat))
      arg = ab / (ab + zetaP16o16*r1) ! ab / (1 + ζR/(16·ab))
      arg = arg * arg ! ab / (1 + ζR/(16·ab))²
      arg = arg * arg ! ab / (1 + ζR/(16·ab))⁴
      arg = arg * arg ! ab / (1 + ζR/(16·ab))⁸
      arg = arg * arg ! ab / (1 + ζR/(16·ab))¹⁶
      fgb = r1 + ab*arg
      dfgb = 1.0_wp / fgb

      Amat(iat,jat) = gbm%keps*dfgb + Amat(iat,jat)
      Amat(jat,iat) = gbm%keps*dfgb + Amat(jat,iat)
   enddo
   !$omp end parallel do

   ! self-energy part
   do iat = 1, self%nat
      bp = 1._wp/self%brad(iat)
      Amat(iat,iat) = Amat(iat,iat) + gbm%keps*bp
   enddo

end subroutine compute_amat_p16

pure subroutine compute_gb_damat(self,q,gborn,ghb,dAmatdr,Afac,lpr)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(inout) :: dAmatdr(3,self%nat,self%nat)
   real(wp), intent(inout) :: Afac(3,self%nat)
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
      do kk = 1, self%ntpair
         r = self%ddpair(1,kk)
         r2 = r*r

         i = self%ppind(1,kk)
         j = self%ppind(2,kk)

         ! dielectric scaling of the charges
         qq = q(i)*q(j)*gbm%keps
         aa = self%brad(i)*self%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*gbm%keps

         egb = egb + qq*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*self%ddpair(2:4,kk)
         dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*q(i)
         dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*q(j)
         Afac(:,i) = Afac(:,i) - dr*q(j)
         Afac(:,j) = Afac(:,j) + dr*q(i)

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = self%dbrdp(i)*self%brad(j)*bp
         grddbj = self%dbrdp(j)*self%brad(i)*bp

         dAmatdr(:,:,j) = dAmatdr(:,:,j) + self%brdr(:,:,j)*grddbj*q(i)
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + self%brdr(:,:,i)*grddbi*q(j)

      enddo

      ! self-energy part
      do i = 1, self%nat
         bp = 1._wp/self%brad(i)
         qq = q(i)*bp
         egb = egb + 0.5_wp*q(i)*qq*gbm%keps
         grddbi = -self%dbrdp(i)*gbm%keps*bp*bp*0.5_wp
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + self%brdr(:,:,i)*grddbi*q(i)
      enddo

   else
      ! GB-SE energy and dAmatdr

      epu=1._wp/gbm%epsu

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, self%ntpair
         r = self%ddpair(1,kk)
         r2 = r*r

         i = self%ppind(1,kk)
         j = self%ppind(2,kk)

         qq = q(i)*q(j)
         aa = self%brad(i)*self%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         fgb = sqrt(fgb2)
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         aa = gbm%kappa*fgb
         expa = exp(-aa)
         gg = (self%ionscr(i)+self%ionscr(j))*expa

         egb = egb + qq*dfgb*(gg-epu)

         dfgb3 = (gg*(1._wp+aa)-epu)*dfgb*dfgb2

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*self%ddpair(2:4,kk)
         dAmatdr(:,i,j) = dAmatdr(:,i,j) - dr*q(i)
         dAmatdr(:,j,i) = dAmatdr(:,j,i) + dr*q(j)
         Afac(:,i) = Afac(:,i) - dr*q(j)
         Afac(:,j) = Afac(:,j) + dr*q(i)

         qfg = dfgb*expa
         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = self%dbrdp(i)*(self%brad(j)*bp+qfg*self%discr(i))*q(j)
         grddbj = self%dbrdp(j)*(self%brad(i)*bp+qfg*self%discr(j))*q(i)

         dAmatdr(:,:,i) = dAmatdr(:,:,i) + self%brdr(:,:,i)*grddbi
         dAmatdr(:,:,j) = dAmatdr(:,:,j) + self%brdr(:,:,j)*grddbj

      enddo

      ! self-energy part
      do i = 1, self%nat
         gg = exp(-gbm%kappa*self%brad(i))
         aa = 2._wp*self%ionscr(i)*gg-epu
         qq = q(i)/self%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-self%brad(i)*2._wp*(self%discr(i)+self%ionscr(i)*gbm%kappa)*gg
         grddbi = -self%dbrdp(i)*0.5_wp*qq*ap/self%brad(i)
         dAmatdr(:,:,i) = dAmatdr(:,:,i) + self%brdr(:,:,i)*grddbi
      enddo

   endif

   gborn = egb

   if(lhb) then
      call compute_ahb(self,q,ghb,dAmatdr)
   endif

end subroutine compute_gb_damat

subroutine compute_gb_egrad(self,xyz,q,gborn,ghb,gradient,lpr)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: xyz(3,self%nat)
   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: gradient(3,self%nat)
   logical,  intent(in)    :: lpr

   select case(gbm%kernel)
   case(gbKernel%still)
      call compute_gb_egrad_still(self,q,gborn,gradient,lpr)
   case(gbKernel%p16)
      call compute_gb_egrad_p16(self,q,gborn,gradient,lpr)
   end select

   gradient = gradient + self%dsdr

   if(lhb) then
      call compute_hb_egrad(self,q,ghb,gradient)
   else
      ghb = 0.0_wp
   endif

   if (gbm%alpbet > 0.0_wp) then
      gborn = gborn + sum(q)**2 * gbm%alpbet / self%aDet * gbm%kEps
      call getADetDeriv(self%nat, xyz, self%vdwr, gbm%kEps*gbm%alpbet, q, gradient)
   end if

end subroutine compute_gb_egrad

subroutine compute_gb_egrad_still(self,q,gborn,gradient,lpr)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(inout) :: gradient(3,self%nat)
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

   allocate(grddb(self%nat), source = 0.0_wp )

   egb = 0._wp
   grddb = 0._wp

   if(.not.gbm%lsalt) then
      ! GB energy and gradient

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, self%ntpair
         r = self%ddpair(1,kk)
         r2 = r*r

         i = self%ppind(1,kk)
         j = self%ppind(2,kk)

         ! dielectric scaling of the charges
         qq = q(i)*q(j)
         aa = self%brad(i)*self%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         fgb2 = r2+aa*expd
         dfgb2 = 1._wp/fgb2
         dfgb = sqrt(dfgb2)
         dfgb3 = dfgb2*dfgb*gbm%keps

         egb = egb + qq*gbm%keps*dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*self%ddpair(2:4,kk)
         gradient(:,i) = gradient(:,i) - dr*qq
         gradient(:,j) = gradient(:,j) + dr*qq

         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = self%dbrdp(i)*self%brad(j)*bp
         grddbj = self%dbrdp(j)*self%brad(i)*bp
         grddb(i) = grddb(i) + grddbi*qq
         !gradient = gradient + self%brdr(:,:,i) * grddbi*qq
         grddb(j) = grddb(j) + grddbj*qq
         !gradient = gradient + self%brdr(:,:,j) * grddbj*qq

      enddo

      ! self-energy part
      do i = 1, self%nat
         bp = 1._wp/self%brad(i)
         qq = q(i)*bp
         egb = egb + 0.5_wp*q(i)*qq*gbm%keps
         grddbi = -self%dbrdp(i)*0.5_wp*gbm%keps*qq*bp
         grddb(i) = grddb(i) + grddbi*q(i)
         !gradient = gradient + self%brdr(:,:,i) * grddbi*q(i)
      enddo

   else
      ! GB-SE energy and gradient

      epu=1._wp/gbm%epsu

      ! compute energy and fgb direct and radii derivatives
      do kk = 1, self%ntpair
         r = self%ddpair(1,kk)
         r2 = r*r

         i = self%ppind(1,kk)
         j = self%ppind(2,kk)

         qq = q(i)*q(j)
         aa = self%brad(i)*self%brad(j)
         dd = a4*r2/aa
         expd = exp(-dd)
         dfgb = r2+aa*expd
         fgb = sqrt(dfgb)
         aa = gbm%kappa*fgb
         expa = exp(-aa)
         gg = (self%ionscr(i)+self%ionscr(j))*expa
         qfg = qq/fgb

         egb = egb + qfg*(gg-epu)

         dfgb3 = qfg*(gg*(1._wp+aa)-epu)/dfgb

         ap = (1._wp-a4*expd)*dfgb3
         dr = ap*self%ddpair(2:4,kk)
         gradient(:,i) = gradient(:,i) - dr
         gradient(:,j) = gradient(:,j) + dr

         qfg = qfg*expa
         bp = -0.5_wp*expd*(1._wp+dd)*dfgb3
         grddbi = self%dbrdp(i)*(self%brad(j)*bp+qfg*self%discr(i))
         grddbj = self%dbrdp(j)*(self%brad(i)*bp+qfg*self%discr(j))
         grddb(i) = grddb(i) + grddbi
         grddb(j) = grddb(j) + grddbj

      enddo

      ! self-energy part
      do i = 1, self%nat
         gg = exp(-gbm%kappa*self%brad(i))
         aa = 2._wp*self%ionscr(i)*gg-epu
         qq = q(i)/self%brad(i)
         egb = egb + 0.5_wp*qq*q(i)*aa
         ap = aa-self%brad(i)*2._wp*(self%discr(i)+self%ionscr(i)*gbm%kappa)*gg
         grddbi = -self%dbrdp(i)*0.5_wp*qq*qq*ap
         grddb(i) = grddb(i) + grddbi
      enddo

   endif
   ! contract with the Born radii derivatives
   call mctc_gemv(self%brdr, grddb, gradient, beta=1.0_wp)

   gborn = egb

end subroutine compute_gb_egrad_still

subroutine compute_gb_egrad_p16(self,q,gborn,gradient,lpr)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(out)   :: gborn
   real(wp), intent(inout) :: gradient(3,self%nat)
   logical,  intent(in)    :: lpr

   integer :: iat, jat, kk
   real(wp) :: vec(3), r2, r1, ab, arg1, arg16, qq, fgb, fgb2, dfgb, dfgb2, egb
   real(wp) :: dEdbri, dEdbrj, dG(3), ap, bp, dS(3, 3)
   real(wp),allocatable :: dEdbr(:)

   allocate(dEdbr(self%nat), source = 0.0_wp )

   egb = 0._wp
   dEdbr(:) = 0._wp

   ! GB energy and gradient
   !$omp parallel do default(none) reduction(+:egb, gradient, dEdbr) &
   !$omp private(iat, jat, vec, r1, r2, ab, arg1, arg16, fgb, dfgb, dfgb2, ap, &
   !$omp& bp, qq, dEdbri, dEdbrj, dG, dS) shared(self, gbm, q)
   do kk = 1, self%ntpair
      vec(:) = self%ddpair(2:4, kk)
      r1 = self%ddpair(1,kk)
      r2 = r1*r1

      iat = self%ppind(1,kk)
      jat = self%ppind(2,kk)
      qq = q(iat)*q(jat)

      ab = sqrt(self%brad(iat) * self%brad(jat))
      arg1 = ab / (ab + zetaP16o16*r1) ! 1 / (1 + ζR/(16·ab))
      arg16 = arg1 * arg1 ! 1 / (1 + ζR/(16·ab))²
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁴
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))⁸
      arg16 = arg16 * arg16 ! 1 / (1 + ζR/(16·ab))¹⁶

      fgb = r1 + ab*arg16
      dfgb = 1.0_wp / fgb
      dfgb2 = dfgb * dfgb

      egb = egb + qq*gbm%keps*dfgb

      ! (1 - ζ/(1 + Rζ/(16 ab))^17)/(R + ab/(1 + Rζ/(16 ab))¹⁶)²
      ap = (1.0_wp - zetaP16 * arg1 * arg16) * dfgb2
      dG(:) = ap * vec * gbm%keps / r1 * qq
      gradient(:, iat) = gradient(:, iat) - dG
      gradient(:, jat) = gradient(:, jat) + dG

      ! -(Rζ/(2·ab²·(1 + Rζ/(16·ab))¹⁷) + 1/(2·ab·(1 + Rζ/(16·ab))¹⁶))/(R + ab/(1 + Rζ/(16·ab))¹⁶)²
      bp = -0.5_wp*(r1 * zetaP16 / ab * arg1 + 1.0_wp) / ab * arg16 * dfgb2
      dEdbri = self%dbrdp(iat)*self%brad(jat) * bp * gbm%keps * qq
      dEdbrj = self%dbrdp(jat)*self%brad(iat) * bp * gbm%keps * qq
      dEdbr(iat) = dEdbr(iat) + dEdbri
      dEdbr(jat) = dEdbr(jat) + dEdbrj

   end do
   !$omp end parallel do

   ! self-energy part
   do iat = 1, self%nat
      bp = 1._wp/self%brad(iat)
      qq = q(iat)*bp
      egb = egb + 0.5_wp*q(iat)*qq*gbm%keps
      dEdbri = -self%dbrdp(iat)*0.5_wp*gbm%keps*qq*bp
      dEdbr(iat) = dEdbr(iat) + dEdbri*q(iat)
      !gradient = gradient + self%brdr(:,:,i) * dEdbri*q(i)
   enddo

   ! contract with the Born radii derivatives
   call mctc_gemv(self%brdr, dEdbr, gradient, beta=1.0_wp)

   gborn = egb

end subroutine compute_gb_egrad_p16

!> Compute contributions to potential for hydrogen bonding correction
pure subroutine compute_fhb(self)
   implicit none
   type(TSolvent), intent(inout) :: self

   integer  :: i
   integer  :: iz,nhb
   real(wp) :: hbed,dhbed
   real(wp) :: smaxd,sasad,sasaw
   real(wp) :: sfw,dsfw,w3,w2,w1
   integer  :: j
   real(wp) :: wbh,wah

   self%hbw=0.0_wp
   self%dhbdw=0.0_wp

   do i = 1, self%nat
      ! atomic Z
      iz=self%at(i)
      ! number of HB
      nhb=gbm%at_hb(iz)
      if(nhb.le.0) cycle
      ! SASA-D for HB
      smaxd=1.0_wp/(self%vdwsa(i)*self%vdwsa(i)*self%gamsasa(i))
      sasad=self%sasa(i)*smaxd
      self%hbw(i)=gbm%hb_mag(iz)*sasad
      self%dhbdw(i)=gbm%hb_mag(iz)*smaxd
   enddo

end subroutine compute_fhb

!> Compute contributions to energy and gradient for hydrogen bonding correction
pure subroutine compute_hb_egrad(self,q,ghb,gradient)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: gradient(3,self%nat)

   integer  :: i,j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, self%nat
      qq = q(i)*q(i)
      ghb = ghb + self%hbw(i)*qq
   enddo

   do i = 1, self%nat
      dhbed=self%dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)*q(i)
      do j = 1, self%nat
         gradient(:,j) = gradient(:,j) + self%dsdrt(:,j,i)*dhbed
      enddo
   enddo

end subroutine compute_hb_egrad

pure subroutine compute_ahb(self,q,ghb,dAmatdr)
   implicit none
   type(TSolvent), intent(in) :: self

   real(wp), intent(in)    :: q(self%nat)
   real(wp), intent(out)   :: ghb
   real(wp), intent(inout) :: dAmatdr(3,self%nat,self%nat)

   integer  :: i,j
   real(wp) :: dhbed
   real(wp) :: qq

   ghb=0.0_wp
   do i = 1, self%nat
      qq = q(i)*q(i)
      ghb = ghb + self%hbw(i)*qq
   enddo

   do i = 1, self%nat
      dhbed=self%dhbdw(i)
      if(abs(dhbed).le.0.0_wp) cycle
      dhbed=dhbed*q(i)
      dAmatdr(:,:,i) = dAmatdr(:,:,i) + self%dsdrt(:,:,i)*dhbed
   enddo

end subroutine compute_ahb

subroutine compute_brad_sasa(self,xyz)
   implicit none
   type(TSolvent), intent(inout) :: self

   real(wp), intent(in) :: xyz(3,self%nat)

   integer i,j,kk

   self%brad=0.0_wp
   self%dsdr=0.0_wp
   self%dsdrt=0.0_wp
   self%brdr=0.0_wp

   call compute_psi(self)

   call compute_bornr(self)

   ! compute solvent accessible surface and its derivatives
   call compute_numsa(self,xyz)

   ! compute the Debye-Hueckel ion exclusion term
   if (gbm%lsalt) call compute_debye_hueckel(self)

   ! compute the HB term
   if (lhb) call compute_fhb(self)

   if (gbm%alpbet > 0.0_wp) then
      call getADet(self%nat, xyz, self%vdwr, self%aDet)
   end if

end subroutine compute_brad_sasa

!> compute the Debye-Hueckel ion exclusion term
pure subroutine compute_debye_hueckel(self)
   implicit none
   type(TSolvent), intent(inout) :: self

   integer  :: i
   real(wp) :: aa,gg

   aa=0.5_wp/gbm%epsv
   do i = 1, self%nat
      gg=gbm%kappa*(self%brad(i)+gbm%ion_rad)
      self%ionscr(i)=aa*exp(gg)/(1.0_wp+gg)
      self%discr(i)=self%ionscr(i)*gbm%kappa*gg/(1.0_wp+gg)
   enddo

end subroutine compute_debye_hueckel

pure subroutine compute_bornr(self)
   implicit none
   type(TSolvent), intent(inout) :: self

   integer iat
   real(wp) br
   real(wp) dpsi
   real(wp) svdwi,vdwri
   real(wp) s1,v1,s2
   real(wp) arg,arg2,th,ch
   real(wp) alpi,beti,gami

   do iat = 1, self%nat

      br = self%brad(iat)

      svdwi=self%svdw(iat)
      vdwri=self%vdwr(iat)
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

      self%brad(iat) = br
      self%dbrdp(iat) = dpsi

   enddo

end subroutine compute_bornr

pure subroutine compute_psi(self)
   implicit none
   type(TSolvent), intent(inout) :: self

   real(wp),allocatable :: br(:),brdr(:,:,:),brdrt(:,:)
   integer  :: kk
   integer  :: i,ii,jj,nn
   real(wp) :: dr(3),r,rhoi,rhoj
   real(wp) :: gi,gj,ap,am,lnab,rhab,ab,dgi,dgj
   real(wp) :: drjj(3)
   real(wp) :: rh1,rhr1,r24,rh2,r1,aprh1,r12
   real(wp) :: rvdwi,rvdwj
   integer  :: ovij,ovji,ov

   allocate( br(self%nat),brdr(3,self%nat,self%nat),brdrt(3,self%nat), &
      &      source = 0.0_wp)

   do kk = 1, self%nnrad

      ii=self%nnlistr(1,kk)
      jj=self%nnlistr(2,kk)
      nn=self%nnlistr(3,kk)

      r=self%ddpair(1,nn)
      dr(:)=self%ddpair(2:4,nn)

      rhoi=self%rho(ii)
      rhoj=self%rho(jj)
      rvdwi=self%vdwr(ii)
      rvdwj=self%vdwr(jj)

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
   self%brad = br
   ! save derivative of Born radii w.r.t. atomic positions
   self%brdr = brdr
   ! save one-center terms
   do i = 1, self%nat
      self%brdr(:,i,i) = brdrt(:,i)
   enddo

   deallocate( br,brdr,brdrt)

end subroutine compute_psi

pure subroutine compute_numsa(self,xyz)
   implicit none
   type(TSolvent), intent(inout) :: self

   real(wp), intent(in) :: xyz(3,self%nat)

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

   allocate(grads(3,self%nat), source = 0.0_wp)

   do iat = 1, self%nat

      rsas = self%vdwsa(iat)
      ! allocate space for the gradient storage
      nno=self%nnsas(iat)
      allocate(grds(3,nno))
      allocate(grdi(nno))

      ! initialize storage
      grads = 0.0_wp
      sasai = 0.0_wp

      ! atomic position
      xyza(:)=xyz(:,iat)
      ! radial atomic weight
      wr=self%wrp(iat)*self%gamsasa(iat)

      ! loop over grid points
      do ip=1,self%nang
         ! grid point position
         xyzp(:) = xyza(:) + rsas*self%angGrid(1:3,ip)
         ! atomic surface function at the grid point
         call compute_w_sp(self%nat,self%nnlists,self%trj2,self%vdwsa, &
            &              xyz,iat,nno,xyzp,sasap,grds,nni,grdi)

         if(sasap.gt.tolsesp) then
            ! numerical quadrature weight
            wsa = self%angWeight(ip)*wr*sasap
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

      self%sasa(iat) = sasai
      self%dsdrt(:,:,iat) = grads

   enddo

   ! contract surface gradient
   self%dsdr = sum(self%dsdrt, dim=3)

   self%gsasa = sum(self%sasa)

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
subroutine update_nnlist_gbsa(self,xyz,parallel)
   implicit none
   type(TSolvent),intent(inout) :: self

   real(wp),intent(in) :: xyz(3,self%nat)
   logical, intent(in) :: parallel

   if (parallel) then
      call update_nnlist_gbsa_parallel(self,xyz)
   else
      call update_nnlist_gbsa_sequential(self,xyz)
   endif

end subroutine update_nnlist_gbsa
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut.
!  OMP parallel version.
subroutine update_nnlist_gbsa_parallel(self,xyz)
!$ use omp_lib
   implicit none
   type(TSolvent),intent(inout) :: self

   real(wp),intent(in) :: xyz(3,self%nat)

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

   allocate(plisttr(3,self%ntpair,nproc),nnls(self%nat,self%nat))
   allocate(nntmp(self%nat),npid(nproc))
   npid = 0

   lrcut2 = self%lrcut*self%lrcut
   srcut2 = self%srcut*self%srcut

   self%nnsas=0
   self%nnlists=0
!$omp parallel default(none) &
!$omp&         shared ( self,xyz,lrcut2,srcut2 ) &
!$omp&         private( i1,i2,x,y,z,dr2,ip,ip2,thrid,nntmp,nnls ) &
!$omp&         shared ( plisttr, npid )
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   thrid=1
!$ thrid=omp_get_thread_num() + 1
!$omp do
   do kk=1,self%ntpair
      i1=self%ppind(1,kk)
      i2=self%ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      self%ddpair(2,kk)=x
      self%ddpair(3,kk)=y
      self%ddpair(4,kk)=z
      self%ddpair(1,kk)=sqrt(dr2)
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
   do i1=1,self%nat
      do i2=1,nntmp(i1)
         self%nnlists(self%nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      self%nnsas(i1)=self%nnsas(i1)+nntmp(i1)
   enddo
!$omp end critical
!$omp end parallel

   self%nnrad=0
   do thrid=1,nproc
      do kk = self%nnrad+1,self%nnrad+npid(thrid)
         self%nnlistr(1:3,kk)=plisttr(1:3,kk-self%nnrad,thrid)
      enddo
      self%nnrad = self%nnrad + npid(thrid)
   enddo

   deallocate(nntmp,nnls)

end subroutine update_nnlist_gbsa_parallel
!> setup a pairlist and compute pair distances of all neighbors
!  within thresholds lrcut and srcut
!  Sequential version.
pure subroutine update_nnlist_gbsa_sequential(self,xyz)
   implicit none
   type(TSolvent),intent(inout) :: self

   real(wp),intent(in) :: xyz(3,self%nat)

   integer kk,i1,i2
   real(wp) rcutn2,lrcut2,srcut2
   real(wp) x,y,z,dr2
   integer ip,ip2
   integer, allocatable :: nntmp(:)
   integer, allocatable :: nnls(:,:)

   allocate(nnls(self%nat,self%nat))
   allocate(nntmp(self%nat))

   lrcut2 = self%lrcut*self%lrcut
   srcut2 = self%srcut*self%srcut

   self%nnsas=0
   self%nnlists=0
   ip=0
   ip2=0
   nntmp=0
   nnls=0
   do kk=1,self%ntpair
      i1=self%ppind(1,kk)
      i2=self%ppind(2,kk)
      x=xyz(1,i1)-xyz(1,i2)
      y=xyz(2,i1)-xyz(2,i2)
      z=xyz(3,i1)-xyz(3,i2)
      dr2=x**2+y**2+z**2
      self%ddpair(2,kk)=x
      self%ddpair(3,kk)=y
      self%ddpair(4,kk)=z
      self%ddpair(1,kk)=sqrt(dr2)
      if(dr2.lt.lrcut2) then
         ip = ip + 1
         self%nnlistr(1,ip)=i1
         self%nnlistr(2,ip)=i2
         self%nnlistr(3,ip)=kk
         if(dr2.lt.srcut2) then
            nntmp(i1) = nntmp(i1) + 1
            nntmp(i2) = nntmp(i2) + 1
            nnls(nntmp(i1),i1)=i2
            nnls(nntmp(i2),i2)=i1
         endif
      endif
   enddo
   self%nnrad = ip
   do i1=1,self%nat
      do i2=1,nntmp(i1)
         self%nnlists(self%nnsas(i1)+i2,i1)=nnls(i2,i1)
      enddo
      self%nnsas(i1)=self%nnsas(i1)+nntmp(i1)
   enddo

   deallocate(nntmp,nnls)

end subroutine update_nnlist_gbsa_sequential

pure subroutine update_dist_gbsa(self,xyz)
   implicit none
   type(TSolvent),intent(inout) :: self

   real(wp),intent(in) :: xyz(3,self%nat)

   integer i1,i2,kk

   do kk = 1, self%ntpair
      i1=self%ppind(1,kk)
      i2=self%ppind(2,kk)
      self%ddpair(2:4,kk)=xyz(1:3,i1)-xyz(1:3,i2)
      self%ddpair(1,kk)=sqrt(self%ddpair(2,kk)**2+ &
      &                      self%ddpair(3,kk)**2+ &
      &                      self%ddpair(4,kk)**2)
   enddo

end subroutine update_dist_gbsa

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine getADet(nAtom, xyz, rad, aDet)
   use xtb_mctc_math, only : matDet3x3

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   !> Shape descriptor of the structure
   real(wp), intent(out) :: aDet

   integer :: iat
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3)
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do

   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))

end subroutine getADet


subroutine getADetDeriv(nAtom, xyz, rad, kEps, qvec, gradient)
   use xtb_mctc_math, only : matDet3x3

   !> Number of atoms
   integer, intent(in) :: nAtom

   !> Cartesian coordinates
   real(wp), intent(in) :: xyz(:, :)

   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   real(wp), intent(in) :: kEps
   real(wp), intent(in) :: qvec(:)

   !> Molecular gradient
   real(wp), intent(inout) :: gradient(:, :)

   integer :: iat
   real(wp) :: r2, rad2, rad3, totRad3, vec(3), center(3), inertia(3, 3), aDet
   real(wp) :: aDeriv(3, 3), qtotal
   real(wp), parameter :: tof = 2.0_wp/5.0_wp, unity(3, 3) = reshape(&
      & [1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 1.0_wp], &
      & [3, 3])

   qtotal = 0.0_wp
   totRad3 = 0.0_wp
   center(:) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      totRad3 = totRad3 + rad3
      center(:) = center + xyz(:, iat) * rad3
      qtotal = qtotal + qvec(iat)
   end do
   center = center / totRad3

   inertia(:, :) = 0.0_wp
   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      r2 = sum(vec**2)
      inertia(:, :) = inertia + rad3 * ((r2 + tof*rad2) * unity &
         & - spread(vec, 1, 3) * spread(vec, 2, 3))
   end do
   aDet = sqrt(matDet3x3(inertia)**(1.0_wp/3.0_wp)/(tof*totRad3))

   aDeriv(:, :) = reshape([&
      & inertia(1,1)*(inertia(2,2)+inertia(3,3))-inertia(1,2)**2-inertia(1,3)**2, &
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,2)*inertia(3,3)-inertia(1,3)*inertia(2,3), & ! xy
      & inertia(2,2)*(inertia(1,1)+inertia(3,3))-inertia(1,2)**2-inertia(2,3)**2, &
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(1,3)*inertia(2,2)-inertia(1,2)*inertia(3,2), & ! xz
      & inertia(1,1)*inertia(2,3)-inertia(1,2)*inertia(1,3), & ! yz
      & inertia(3,3)*(inertia(1,1)+inertia(2,2))-inertia(1,3)**2-inertia(2,3)**2],&
      & shape=[3, 3]) * (250.0_wp / (48.0_wp * totRad3**3 * aDet**5)) &
      & * (-0.5_wp * kEps * qtotal**2 / aDet**2)

   do iat = 1, nAtom
      rad2 = rad(iat) * rad(iat)
      rad3 = rad2 * rad(iat)
      vec(:) = xyz(:, iat) - center
      gradient(:, iat) = gradient(:, iat) + rad3 * matmul(aderiv, vec)
   end do

end subroutine getADetDeriv


end module xtb_solv_gbobc
