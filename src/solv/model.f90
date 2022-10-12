
! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Generic solvation model used to create actual solvation calculators
module xtb_solv_model
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_constants, only : fourpi
   use xtb_mctc_convert, only : aatoau, autoaa, kcaltoau, autokcal
   use xtb_mctc_strings, only : lowercase
   use xtb_mctc_systools, only : rdpath
   use xtb_param_vdwradd3, only : vanDerWaalsRadD3
   use xtb_param_vdwradcosmo, only : vanDerWaalsRadCosmo
   use xtb_solv_gbsa, only : TBorn, init_ => init
   use xtb_solv_cosmo, only : TCosmo, init_ => init
   use xtb_solv_input, only : TSolvInput
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_state, only : solutionState, getStateShift
   use xtb_type_environment, only : TEnvironment
   use xtb_type_solvation, only : TSolvation
   use ieee_arithmetic, only : ieee_value,ieee_positive_inf
   implicit none
   private

   public :: TSolvModel, init, info, newSolvationModel, newBornModel


   !> Payload to create actual solvation models
   type, extends(TSolvInput) :: TSolvModel

      !> Meta data: parameter file name
      character(len=:), allocatable :: paramFile

      !> Dielectric constant of the solvent
      real(wp) :: dielectricConst

      !> Molar mass of a solvent molecule (g/mol)
      real(wp) :: molarMass

      !> Density of the solvent
      real(wp) :: density

      !> Free energy shift added to the solvation state specific shift
      real(wp) :: freeEnergyShift

      !> Scaling factor for the Born radii
      real(wp) :: bornScale

      !> Offset parameter for the Born radii integration
      real(wp) :: bornOffset

      !> Probe radius for the integration of the solvent accessible surface area
      real(wp) :: probeRad

      !> Analytical Linearized Poisson Boltzmann alpha constant (zero for GBSA)
      real(wp) :: alpha

      !> Van-der-Waals radii used for various purposes
      real(wp), allocatable :: vdwRad(:)

      !> Surface tension scaling
      real(wp), allocatable :: surfaceTension(:)

      !> Dielectric descreening parameter for GBOBC integrator
      real(wp), allocatable :: descreening(:)

      !> Hydrogen bond strength
      real(wp), allocatable :: hBondStrength(:)

   end type TSolvModel


   !> Initialize solvation model from input
   interface init
      module procedure :: initSolvModel
   end interface init


   !> Define the local data structure for the parameters
   type :: gbsa_parameter
      real(wp) :: epsv = 0.0_wp
      real(wp) :: smass = 0.0_wp
      real(wp) :: rhos = 0.0_wp
      real(wp) :: c1 = 0.0_wp
      real(wp) :: rprobe = 0.0_wp
      real(wp) :: gshift = 0.0_wp
      real(wp) :: soset = 0.0_wp
      real(wp) :: alpha = 0.0_wp
      real(wp) :: gamscale(94) = 0.0_wp
      real(wp) :: sx(94) = 0.0_wp
      real(wp) :: tmp(94) = 0.0_wp
   end type gbsa_parameter

   include 'param_gbsa_acetone.fh'
   include 'param_gbsa_acetonitrile.fh'
   include 'param_gbsa_benzene.fh'
   include 'param_gbsa_ch2cl2.fh'
   include 'param_gbsa_chcl3.fh'
   include 'param_gbsa_cs2.fh'
   include 'param_gbsa_dmso.fh'
   include 'param_gbsa_ether.fh'
   include 'param_gbsa_h2o.fh'
   include 'param_gbsa_methanol.fh'
   include 'param_gbsa_thf.fh'
   include 'param_gbsa_toluene.fh'
   include 'param_gbsa_dmf.fh'
   include 'param_gbsa_nhexan.fh'

   include 'param_alpb_acetone.fh'
   include 'param_alpb_acetonitrile.fh'
   include 'param_alpb_aniline.fh'
   include 'param_alpb_benzaldehyde.fh'
   include 'param_alpb_benzene.fh'
   include 'param_alpb_ch2cl2.fh'
   include 'param_alpb_chcl3.fh'
   include 'param_alpb_cs2.fh'
   include 'param_alpb_dioxane.fh'
   include 'param_alpb_dmf.fh'
   include 'param_alpb_dmso.fh'
   include 'param_alpb_ether.fh'
   include 'param_alpb_ethylacetate.fh'
   include 'param_alpb_furane.fh'
   include 'param_alpb_hexadecane.fh'
   include 'param_alpb_hexane.fh'
   include 'param_alpb_nitromethane.fh'
   include 'param_alpb_octanol.fh'
   include 'param_alpb_phenol.fh'
   include 'param_alpb_thf.fh'
   include 'param_alpb_toluene.fh'
   include 'param_alpb_water.fh'
   include 'param_alpb_woctanol.fh'
   include 'param_alpb_methanol.fh'
   include 'param_alpb_ethanol.fh'

   include 'param_cosmo.fh'
   include 'param_cosmo_inf.fh'

   !> Solvent density (g/cm^3) and molar mass (g/mol)
   real(wp), parameter :: molcm3toau = 8.92388e-2_wp

   !> Surface tension (mN/m=dyn/cm)
   real(wp), parameter :: surfaceTension = 1.0e-5_wp
   real(wp), parameter :: mNmtokcal = 4.0305201015221386e-4_wp
   real(wp), parameter :: kcaltomNm = 1.0_wp/mNmtokcal
   real(wp), parameter :: automNm = autokcal * kcaltomNm

   real(wp),parameter :: lrcut = 35.0_wp * aatoau
   real(wp),parameter :: srcut = 2.0_wp * aatoau


contains


!> Initialize solvation model from input
subroutine initSolvModel(self, env, input, level)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_initSolvModel'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Data structure
   type(TSolvModel), intent(out) :: self

   !> Input data
   type(TSolvInput), intent(in) :: input

   !> Method level for solvation model
   integer, intent(in) :: level

   character(len=:), allocatable :: solvent
   logical :: exist

   self%TSolvInput = input
   call normalizeSolventName(solvent, input%solvent)

   if (input%alpb) then
      call getParamFile(env, solvent, 'alpb', level, self%paramFile)
   end if

   if (input%cosmo) then
      call getParamFile(env, solvent, 'cosmo', level, self%paramFile)
   end if

   if (.not.allocated(self%paramFile) .and. (.not.input%cosmo)) then
      call getParamFile(env, solvent, 'gbsa', level, self%paramFile)
   end if

   if (allocated(self%paramFile)) then
      call readParamFile(self, env)
   else
      call loadInternalParam(self, env, solvent, level)
   end if

   self%freeEnergyShift = self%freeEnergyShift + getStateShift(self%state, &
      & self%temperature, self%density, self%molarMass)
end subroutine initSolvModel


!> Try to find a normalized version of an input solvent name
subroutine normalizeSolventName(solvent, input)

   !> Name of the solvent (normalized)
   character(len=:), allocatable, intent(out) :: solvent

   !> Name of the solvent (unnormalized)
   character(len=*), intent(in) :: input

   solvent = lowercase(input)

end subroutine normalizeSolventName


!> Load internal parameters, if no parameter file is present
subroutine loadInternalParam(self, env, solvent, level)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_loadInternalParam'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Data structure
   type(TSolvModel), intent(inout) :: self

   !> Name of the solvent (normalized)
   character(len=*), intent(in) :: solvent

   !> Method level for solvation model
   integer, intent(in) :: level

   !> I/O handling

   integer :: iostat

   type(gbsa_parameter), allocatable :: param, stub


   select case(self%kernel)
   case (gbKernel%still)
      if (any(level == [2, 0])) then
         self%paramFile = 'internal GFN2-xTB/GBSA'
         select case(solvent)
         case('acetone');      param = gfn2_acetone
         case('acetonitrile'); param = gfn2_acetonitrile
         case('benzene');      param = gfn2_benzene
         case('ch2cl2','dichlormethane'); param = gfn2_ch2cl2
         case('chcl3','chloroform');      param = gfn2_chcl3
         case('cs2');          param = gfn2_cs2
         case('dmso');         param = gfn2_dmso
         case('ether');        param = gfn2_ether
         case('h2o','water');  param = gfn2_h2o
         case('methanol');     param = gfn2_methanol
         case('thf');          param = gfn2_thf
         case('toluene');      param = gfn2_toluene
         case('dmf');          param = gfn2_dmf
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn2_nhexan
         end select
      else
         self%paramFile = 'internal GFN1-xTB/GBSA'
         select case(solvent)
         case('acetone');      param = gfn1_acetone
         case('acetonitrile'); param = gfn1_acetonitrile
         case('benzene');      param = gfn1_benzene
         case('ch2cl2','dichlormethane'); param = gfn1_ch2cl2
         case('chcl3','chloroform');      param = gfn1_chcl3
         case('cs2');          param = gfn1_cs2
         case('dmso');         param = gfn1_dmso
         case('ether');        param = gfn1_ether
         case('h2o','water');  param = gfn1_h2o
         case('methanol');     param = gfn1_methanol
         case('thf');          param = gfn1_thf
         case('toluene');      param = gfn1_toluene
         end select
      end if

      if (.not.allocated(param)) then
         call env%error("solvent: '"//solvent//"' is not parametrized", source)
         return
      end if

   case (gbKernel%p16)
      if (level == 2) then
         self%paramFile = 'internal GFN2-xTB/ALPB'
         select case(solvent)
         case('acetone');      param = gfn2_alpb_acetone
         case('acetonitrile'); param = gfn2_alpb_acetonitrile
         case('aniline');      param = gfn2_alpb_aniline
         case('benzaldehyde');      param = gfn2_alpb_benzaldehyde
         case('benzene');      param = gfn2_alpb_benzene
         case('dioxane');      param = gfn2_alpb_dioxane
         case('ethylacetate');      param = gfn2_alpb_ethylacetate
         case('furane');      param = gfn2_alpb_furane
         case('hexadecane');      param = gfn2_alpb_hexadecane
         case('nitromethane');      param = gfn2_alpb_nitromethane
         case('octanol');      param = gfn2_alpb_octanol
         case('woctanol');      param = gfn2_alpb_woctanol
         case('phenol');      param = gfn2_alpb_phenol 
         case('ch2cl2','dichlormethane'); param = gfn2_alpb_ch2cl2
         case('chcl3','chloroform');      param = gfn2_alpb_chcl3
         case('cs2');          param = gfn2_alpb_cs2
         case('dmso');         param = gfn2_alpb_dmso
         case('ether');        param = gfn2_alpb_ether
         case('h2o','water');  param = gfn2_alpb_water
         case('methanol');     param = gfn2_alpb_methanol 
         case('thf');          param = gfn2_alpb_thf
         case('toluene');      param = gfn2_alpb_toluene
         case('dmf');          param = gfn2_alpb_dmf
         case('ethanol');      param = gfn2_alpb_ethanol
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn2_alpb_hexane
         end select
      elseif (level == 0) then
         self%paramFile = 'internal GFN-FF/ALPB'
         select case(solvent)
         case('acetone');      param = gfnff_alpb_acetone
         case('acetonitrile'); param = gfnff_alpb_acetonitrile
         case('aniline');      param = gfnff_alpb_aniline
         case('benzaldehyde');      param = gfnff_alpb_benzaldehyde
         case('benzene');      param = gfnff_alpb_benzene
         case('dioxane');      param = gfnff_alpb_dioxane
         case('ethylacetate');      param = gfnff_alpb_ethylacetate
         case('furane');      param = gfnff_alpb_furane
         case('hexadecane');      param = gfnff_alpb_hexadecane
         case('nitromethane');      param = gfnff_alpb_nitromethane
         case('octanol');      param = gfnff_alpb_octanol
         case('woctanol');      param = gfnff_alpb_woctanol
         case('phenol');      param = gfnff_alpb_phenol 
         case('ch2cl2','dichlormethane'); param = gfnff_alpb_ch2cl2
         case('chcl3','chloroform');      param = gfnff_alpb_chcl3
         case('cs2');          param = gfnff_alpb_cs2
         case('dmso');         param = gfnff_alpb_dmso
         case('ether');        param = gfnff_alpb_ether
         case('h2o','water');  param = gfnff_alpb_water
         case('methanol');     param = gfnff_alpb_methanol
         case('ethanol');      param = gfnff_alpb_ethanol
         case('thf');          param = gfnff_alpb_thf
         case('toluene');      param = gfnff_alpb_toluene
         case('dmf');          param = gfnff_alpb_dmf
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfnff_alpb_hexane
      end select
      else
         self%paramFile = 'internal GFN1-xTB/ALPB'
         select case(solvent)
         case('acetone');      param = gfn1_alpb_acetone
         case('acetonitrile'); param = gfn1_alpb_acetonitrile
         case('aniline');      param = gfn1_alpb_aniline
         case('benzaldehyde');      param = gfn1_alpb_benzaldehyde
         case('benzene');      param = gfn1_alpb_benzene
         case('dioxane');      param = gfn1_alpb_dioxane
         case('ethylacetate');      param = gfn1_alpb_ethylacetate
         case('furane');      param = gfn1_alpb_furane
         case('hexadecane');      param = gfn1_alpb_hexadecane
         case('nitromethane');      param = gfn1_alpb_nitromethane
         case('octanol');      param = gfn1_alpb_octanol
         case('woctanol');      param = gfn1_alpb_woctanol
         case('phenol');      param = gfn1_alpb_phenol 
         case('ch2cl2','dichlormethane'); param = gfn1_alpb_ch2cl2
         case('chcl3','chloroform');      param = gfn1_alpb_chcl3
         case('cs2');          param = gfn1_alpb_cs2
         case('dmso');         param = gfn1_alpb_dmso
         case('ether');        param = gfn1_alpb_ether
         case('h2o','water');  param = gfn1_alpb_water
         case('methanol');     param = gfn1_alpb_methanol
         case('ethanol');      param = gfn1_alpb_ethanol
         case('thf');          param = gfn1_alpb_thf
         case('toluene');      param = gfn1_alpb_toluene
         case('dmf');          param = gfn1_alpb_dmf
         case('nhexan','n-hexan','nhexane','n-hexane','hexane');
            param = gfn1_alpb_hexane
         end select
      end if
      
      if (self%cosmo) then
         if (allocated(param)) then
            self%paramFile = self%paramFile//"+COSMO"
            stub = param
            param = gfn_cosmo
            param%epsv = stub%epsv
            param%smass = stub%smass
            param%rhos = stub%rhos
            param%gamscale = stub%gamscale
         else if ((solvent .eq. "inf") .or. (solvent .eq. "infinity")) then
            self%paramFile = "internal GFN-xTB/COSMO"
            param = gfn_cosmo_inf
            param%epsv = ieee_value(param%epsv,ieee_positive_inf)
         else
            self%paramFile = "internal GFN-xTB/COSMO"
            param = gfn_cosmo
            read(solvent,*,iostat=iostat) param%epsv
            if (iostat .ne. 0) then
               Call env%error(solvent// " is neither a solvent name nor a dielectric constant.",source)
               return
            end if
            if (param%epsv .eq. 0) then
               Call env%error("You really chose 0 as your dielectric constant?",source)
               return
            end if
         end if
      end if

      if (.not.allocated(param)) then
         call env%error("solvent: '"//solvent//"' is not parametrized", source)
         return
      end if

   end select

   call paramToModel(self, param)

end subroutine loadInternalParam


subroutine readParamFile(self, env)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_readParamFile'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Data structure
   type(TSolvModel), intent(inout) :: self

   type(gbsa_parameter) :: param
   integer :: unit, stat, statc, ii

   statc = 0
   call open_file(unit, self%paramFile, 'r')
   if (unit == -1) then
      call env%error("Could not open parameter file '"//self%paramFile//"'", &
         & source)
      return
   end if

   read(unit,*,iostat=stat) param%epsv
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%smass
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%rhos
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%c1
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%rprobe
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%gshift
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%soset
   if (stat /= 0) statc = statc + 1
   read(unit,*,iostat=stat) param%alpha
   if (stat /= 0) statc = statc + 1
   do ii = 1, 94
      read(unit,*,iostat=stat) param%gamscale(ii), param%sx(ii), param%tmp(ii)
      if (stat /= 0) statc = statc + 1
   enddo

   if (statc /= 0) then
      call env%error("Could not read parameters from '"//self%paramFile//"'", &
         & source)
      return
   end if

   call paramToModel(self, param)

end subroutine readParamFile


subroutine paramToModel(self, param)

   !> Data structure
   type(TSolvModel), intent(inout) :: self

   type(gbsa_parameter), intent(in) :: param

   integer :: ii

   self%dielectricConst = param%epsv
   self%molarMass = param%smass
   self%density = param%rhos/param%smass * molcm3toau
   self%freeEnergyShift = param%gshift * kcaltoau
   self%bornScale = param%c1
   self%bornOffset = param%soset * 0.1_wp*aatoau
   self%probeRad = param%rprobe * aatoau

   if (self%cosmo) then
      self%vdwRad=vanDerWaalsRadCosmo*param%sx
   else
      self%vdwRad=vanDerWaalsRadD3
   end if

   self%surfaceTension = param%gamscale * fourpi*surfaceTension
   self%descreening = param%sx
   allocate(self%hBondStrength(size(param%tmp)))
   if (any(abs(param%tmp).gt.1.e-3_wp)) then
      do ii = 1, size(param%tmp)
         self%hBondStrength(ii) = -(param%tmp(ii)**2) * kcaltoau
      end do
   else
      self%hBondStrength(:) = 0.0_wp
   end if

end subroutine paramToModel


!> Try to find a suitable parameter file for parsing
subroutine getParamFile(env, solvent, model, level, pfile)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_getParamFile'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Name of the solvent (normalized)
   character(len=*), intent(in) :: solvent

   !> Name of the solvation model
   character(len=*), intent(in) :: model

   !> Method level for solvation model
   integer, intent(in) :: level

   !> Parameter file, unallocated if nothing was found
   character(len=:), allocatable, intent(out) :: pfile
   character(len=:), allocatable :: fname

   ! try "param_<model>[level]_<solvent>.txt" first
   select case(level)
   case(2)
      fname = 'param_' // model // '2_' // solvent // '.txt'
   case(1)
      fname = 'param_' // model // '1_' // solvent // '.txt'
   case(0)
      fname = 'param_' // model // '0_' // solvent // '.txt'
   end select

   if (allocated(fname)) then
      call rdpath(env%xtbpath, fname, pfile)
   end if

   if (.not.allocated(pfile)) then
      ! try "param_<model>_<solvent>.txt" next
      fname = 'param_' // model // '_' // solvent // '.txt'
      call rdpath(env%xtbpath, fname, pfile)
   end if

   ! finally try ".param_<model>_<solvent>" which is the old default
   if (.not.allocated(pfile)) then
      select case(level)
      case default
         fname = '.param_' // model // '_' // solvent
      case(2)
         fname = '.param_' // model // '2_' // solvent
      case(0)
         fname = '.param_' // model // '0_' // solvent
      end select

      call rdpath(env%xtbpath, fname, pfile)
   end if

end subroutine getParamFile


!> Print information on the solvation model
subroutine info(self, unit)

   !> Data structure
   class(TSolvModel), intent(in) :: self

   !> Unit for IO
   integer, intent(in) :: unit

   write(unit, '(6x, "*", 1x, a, ":", t40)', advance='no') "Solvation model"
   if (self%cosmo) then
      write(unit, '(a)') "COSMO"
   else
      if (self%alpb) then
         write(unit, '(a)') "ALPB"
      else
         write(unit, '(a)') "GBSA"
      end if
   end if

   write(unit, '(8x, a, t40, a)') "Solvent", self%solvent
   write(unit, '(8x, a, t40, a)') "Parameter file", self%paramFile
   write(unit, '(8x, a, t40, es14.4)') &
      & "Dielectric constant", self%dielectricConst

   write(unit, '(8x, a, t40)', advance='no') "Reference state"
   select case(self%state)
   case default
      write(unit, '(a)') 'gsolv [1 M gas/solution]'
   case(solutionState%reference)
      write(unit, '(a)') 'gsolv=reference [X=1]'
   case(solutionState%mol1bar)
      write(unit, '(a)') 'gsolv [1 bar gas/1 M solution]'
   end select

   write(unit, '(8x, a, t40, es14.4, 1x, a, t60, es14.4, 1x, a)') &
      & "Free energy shift", self%freeEnergyShift, "Eh", &
      & self%freeEnergyShift * autokcal, "kcal/mol"
   write(unit, '(8x, a, t40, es14.4, 1x, a)') &
      & "Temperature", self%temperature, "K"
   write(unit, '(8x, a, t40, es14.4, 1x, a)') &
      & "Density", self%density / (molcm3toau/self%molarMass), "kg/L"
   write(unit, '(8x, a, t40, es14.4, 1x, a)') &
      & "Solvent mass", self%molarMass, "g/mol"

   if (self%cosmo) then
      write(unit, '(8x, a, t40, es14.4)') &
         & "vdW Radii scaling", self%bornScale
   else
      write(unit, '(8x, a, t40)', advance='no') "Interaction kernel"
      select case(self%kernel)
      case default
         write(unit, '(i0, 1x, a)') self%kernel, '(internal error)'
      case(gbKernel%still)
         write(unit, '(a)') 'Still'
      case(gbKernel%p16)
         write(unit, '(a)') 'P16'
      end select

      write(unit, '(8x, a, t40, es14.4)') &
         "Born radius scaling (c1)", self%bornScale
      write(unit, '(8x, a, t40, a)') "Born radii integrator", "GBOBC"
      write(unit, '(8x, a, t40, es14.4, 1x, a, t60, es14.4, 1x, a)') &
         "Born offset", self%bornOffset, "a0", self%bornOffset/autoaa, "AA"

      write(unit, '(8x, a, t40)', advance='no') "H-bond correction"
      if (any(self%hBondStrength < 0.0_wp)) then
         write(unit, '(a)') "true"
      else
         write(unit, '(a)') "false"
      end if

      write(unit, '(8x, a, t40)', advance='no') "Ion screening"
      if (self%ionStrength > 0.0_wp) then
         write(unit, '(a)') "true"
      else
         write(unit, '(a)') "false"
      end if
   end if

   if (allocated(self%surfaceTension)) then
      write(unit, '(8x, a, t40, es14.4, 1x, a, t60, es14.4, 1x, a)') &
         "Surface tension", surfaceTension, "Eh", surfaceTension*automNm, "dyn/cm"
   end if

   write(unit, '(8x, a, t40, i14, 1x, a)') &
      "Grid points", self%nAng, "per atom"

end subroutine info


subroutine newSolvationModel(self, env, model, num)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_newSolvationModel'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Data structure
   class(TSolvModel), intent(in) :: self

   !> Solvation model to create
   class(TSolvation), allocatable, intent(out) :: model

   !> Atomic numbers
   integer, intent(in) :: num(:)

   type(TCosmo), allocatable :: cosmo
   type(TBorn), allocatable :: born

   if (self%cosmo) then
      allocate(cosmo)
      call init_(cosmo, env, num, self%dielectricConst, self%nAng, self%bornScale, &
         & self%vdwRad, self%surfaceTension, self%probeRad, srcut)
      call move_alloc(cosmo, model)
   else
      allocate(born)
      call init_(born, env, num, self%vdwRad, self%dielectricConst, &
         & self%freeEnergyShift, self%descreening, self%bornScale, &
         & self%bornOffset, self%surfaceTension, self%probeRad, lrcut, srcut, &
         & self%nAng, self%hBondStrength, self%temperature, self%kernel, &
         & self%alpb, self%ionStrength, self%ionRad)
      call move_alloc(born, model)
   end if

end subroutine newSolvationModel


subroutine newBornModel(self, env, model, num)

   !> Error source for traceback
   character(len=*), parameter :: source = 'solv_model_newSolvationModel'

   !> Computation environment
   type(TEnvironment), intent(inout) :: env

   !> Data structure
   class(TSolvModel), intent(in) :: self

   !> Solvation model to create
   type(TBorn), intent(out) :: model

   !> Atomic numbers
   integer, intent(in) :: num(:)

   call init_(model, env, num, self%vdwRad, self%dielectricConst, &
      & self%freeEnergyShift, self%descreening, self%bornScale, &
      & self%bornOffset, self%surfaceTension, self%probeRad, lrcut, srcut, &
      & self%nAng, self%hBondStrength, self%temperature, self%kernel, &
      & self%alpb, self%ionStrength, self%ionRad)

end subroutine newBornModel


end module xtb_solv_model
