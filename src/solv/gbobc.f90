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
   use xtb_mctc_blas, only : mctc_dot, mctc_gemv
   use xtb_solv_born, only : compute_bornr
   use xtb_solv_gbsa, only : TBorn, addGradientHBond, getADet, addADetDeriv, &
      & addHBondDeriv, compute_fhb, getDebyeHueckel, update_nnlist_gbsa
   use xtb_solv_kernel, only : gbKernel, addBornMatSaltStill, addBornMatStill, &
      & addBornMatP16, addGradientSaltStill, addGradientStill, addGradientP16, &
      & addBornDerivSaltStill, addBornDerivStill
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_model, only : TSolvModel, newBornModel
   use xtb_solv_sasa, only : compute_numsa
   use xtb_solv_state, only : solutionState
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvent, only : allocate_gbsa, deallocate_gbsa
   implicit none

   public :: lsalt
   public :: initGBSA,new_gbsa
   public :: gshift
   public :: ionst,ion_rad
   public :: TBorn
   public :: allocate_gbsa,deallocate_gbsa
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

subroutine initGBSA(env,gfn_method,verbose,input)
   use xtb_mctc_strings
   use xtb_readin
   implicit none
   character(len=*), parameter :: source = 'solv_gbobc_initGBSA'
   type(TEnvironment), intent(inout) :: env
   integer, intent(in) :: gfn_method
   logical, intent(in) :: verbose
   type(TSolvInput), intent(in) :: input

   integer :: i,fix,inum,ich
   real(wp) :: rad
   real(wp) :: gamma_in, tmp(94), gstate, dum
   character(len=:),allocatable :: fname
   logical ex
   character(len=80) a80
   type(gbsa_parameter) :: gfn_solvent

   if (gfn_method.gt.1) then
      fname = '.param_gbsa2_'//trim(input%solvent)
   else if (gfn_method.eq.0) then
      fname = '.param_gbsa0_'//trim(input%solvent)
   else
      fname = '.param_gbsa_'//trim(input%solvent)
   endif
   fname = xfind(fname)
   if (verbose) then
      write(env%unit,*) 'Solvent             : ', trim(input%solvent)
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
         select case(lowercase(trim(input%solvent)))
         case default
            call env%error('solvent : '//trim(input%solvent)//&
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
         end select
      !else if (gfn_method.eq.0) then
            !call env%error('solvent : '//trim(input%solvent)//&
               !' not parametrized for GFN0-xTB Hamiltonian',source)
      else
         select case(lowercase(trim(input%solvent)))
         case default
            call env%error('solvent : '//trim(input%solvent)//&
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
         end select
      endif
      if (verbose) then
         write(env%unit,'(1x,"Using internal parameter file:",1x,a)') fname
      end if
   endif

   if (input%alpb) gfn_solvent%alpha = 0.571412_wp

   call new_gbsa_model(gbm,gfn_solvent,input%state,input%temperature,input%nAng)

   gbm%kernel = input%kernel

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

   if(mode.eq.solutionState%reference) then ! gsolv=reference option in COSMOTHERM
      !               RT*(ln(ideal gas mol volume)+ln(rho/M))
      gstate=(gbm%temp*8.31451/1000./4.184)* &
      &      (log(24.79_wp*gbm%temp/298.15)+ &
      &       log(1000.0_wp*gbm%rhos/gbm%smass))
      gbm%gshift=(gbm%gshift+gstate)*kcaltoau
   elseif(mode.eq.solutionState%gsolv)then !gsolv option in COSMOTHERM to which it was fitted
      gbm%gshift=gbm%gshift*kcaltoau
   elseif(mode.eq.solutionState%mol1bar)then ! 1 bar gas/ 1 M solution is not implemented in COSMOTHERM although its the canonical choice
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

subroutine new_gbsa(self,env,n,at)
   use xtb_solv_lebedev
   implicit none
   type(TBorn), intent(inout) :: self
   type(TEnvironment), intent(inout) :: env

   integer, intent(in) :: n
   integer, intent(in) :: at(n)

   integer i,j,k,iAng
   integer ierr
   real(wp) minvdwr, maxvdwr
   real(wp) maxrasasa
   real(wp) r
   real(wp), allocatable :: xang(:),yang(:),zang(:),wang(:)

   ! get some space
   call allocate_gbsa(self,n,gbm%nangsa)
   allocate(self%bornMat(n, n))

   self%bornScale = gbm%c1
   self%keps = gbm%keps
   self%kappa = gbm%kappa
   self%lsalt = gbm%lsalt
   self%lhb = gbm%lhb
   self%kernel = gbm%kernel
   self%dielectricConst = gbm%epsv
   self%ionRad = gbm%ion_rad
   self%alpbet = gbm%alpbet

   ! initialize the vdw radii array
   self%at = at
   maxvdwr=0.0_wp
   minvdwr=1000.0_wp
   do i=1,self%nat
      self%vdwr(i)=gbm%rvdw(self%at(i))*aatoau
      self%rho(i)=self%vdwr(i)*gbm%sx(self%at(i))
      self%svdw(i)=self%vdwr(i)-gbm%soset
      maxvdwr=max(maxvdwr,self%vdwr(i))
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

   do i = 1, self%nat
      self%hbmag(i) = gbm%hb_mag(self%at(i))
   end do

end subroutine new_gbsa


pure subroutine update_dist_gbsa(nat,ntpair,ppind,xyz,ddpair)

   integer, intent(in) :: nat
   integer, intent(in) :: ntpair
   integer, intent(in) :: ppind(:, :)
   real(wp),intent(in) :: xyz(:, :)
   real(wp), intent(out) :: ddpair(:, :)

   integer i1,i2,kk

   do kk = 1, ntpair
      i1=ppind(1,kk)
      i2=ppind(2,kk)
      ddpair(2:4,kk)=xyz(1:3,i1)-xyz(1:3,i2)
      ddpair(1,kk)=sqrt(ddpair(2,kk)**2+ &
         &              ddpair(3,kk)**2+ &
         &              ddpair(4,kk)**2)
   enddo

end subroutine update_dist_gbsa


end module xtb_solv_gbobc
