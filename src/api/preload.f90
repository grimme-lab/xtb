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

!> Preloading module of xtb to change the global state of the library.
module xtb_api_preload
   use iso_c_binding
   use xtb_mctc_accuracy, only : wp
   use xtb_api_utils
   use xtb_type_param, only : scc_parameter, TxTBParameter
   use xtb_type_environment, only : TEnvironment, init
   use xtb_paramset
   implicit none

   type(scc_parameter) :: global_parameter

contains

!> Load GFN-xTB parametrisation and create a lock, the file name is optional.
!  The lock (gfn_method) is maintained by the shared libary.
integer(c_int) function load_xtb_parameters_api(gfn, filename) &
      &  result(status) bind(C, name='load_xTB_parameters')
   use xtb_setparam, only : gfn_method
   use xtb_readparam, only : readParam
   !> GFN level to be loaded
   integer(c_int), intent(in) :: gfn
   !> File name for the parametrisation
   character(kind=c_char), intent(in), optional :: filename(*)

   type(TEnvironment) :: env
   character(len=:), allocatable :: fnv
   type(TxTBParameter) :: globpar
   integer :: ipar
   logical :: exist

   call init(env)

   !$omp critical(xtb_load_api)
   status = 1
   gfn_method = -1

   call mctc_init('peeq',10,.true.)
   call mctc_mute

   if (present(filename)) then
      call c_string_convert(fnv, filename)
      inquire(file=fnv, exist=exist)
      if (exist) then
         open(file=fnv, newunit=ipar)
         call readParam(env, ipar, globpar, .true.)
         close(ipar)
         status = 0
      end if
   else
      select case(gfn)
      case(0_c_int)
         call load_xtb_parameters(env, '.param_gfn0.xtb', globpar, status)
      case(1_c_int)
         call load_xtb_parameters(env, '.param_gfn.xtb', globpar, status)
      case(2_c_int)
         call load_xtb_parameters(env, '.param_gfn2.xtb', globpar, status)
      case default
         status = 2
      end select
   end if

   if (status == 0) then
      select case(gfn)
      case(0_c_int)
         call set_gfn0_parameter(global_parameter, globpar)
      case(1_c_int)
         call set_gfn1_parameter(global_parameter, globpar)
      case(2_c_int)
         call set_gfn2_parameter(global_parameter, globpar)
      case default
         status = 2
      end select
   end if

   if (status == 0) then
      gfn_method = gfn
   end if
   !$omp end critical(xtb_load_api)

end function load_xtb_parameters_api

subroutine load_xtb_parameters(env, p_fnv, globpar, status)
   use xtb_mctc_systools
   use xtb_aoparam, only : use_parameterset
   use xtb_readparam, only : readParam
   type(TEnvironment), intent(inout) :: env
   character(len=*), intent(in) :: p_fnv
   integer, intent(out) :: status
   character(len=:), allocatable :: xtbpath
   character(len=:), allocatable :: fnv
   type(TxTBParameter), intent(out) :: globpar
   integer :: ipar
   logical :: exist
   status = 1

   ! we will try an internal parameter file first to avoid IO
   call use_parameterset(p_fnv, globpar, exist)
   if (exist) then
      status = 0
   else ! no luck, we have to fire up some IO to get our parameters
      call rdvar('XTBPATH', xtbpath)
      if (len(xtbpath) > 0) then
         ! let's check if we can find the parameter file
         call rdpath(xtbpath, p_fnv, fnv, exist)
      end if
      ! maybe the user provides a local parameter file, this was always
      ! an option in `xtb', so we will give it a try
      if (.not.exist) fnv = p_fnv
      inquire(file=fnv, exist=exist)
      if (exist) then
         open(file=fnv, newunit=ipar)
         call readParam(env, ipar, globpar, .true.)
         close(ipar)
         status = 0
      end if
   end if

end subroutine load_xtb_parameters

!> allows loading custom solvent parameters in the shared library
function gbsa_model_preload_api &
      &  (epsv_in,smass_in,rhos_in,c1_in,rprobe_in,gshift_in,soset_in,dum_in, &
      &   gamscale_in,sx_in,tmp_in) &
      &  result(status) bind(C,name='GBSA_model_preload')

   use xtb_solv_gbobc, only : load_custom_parameters

   !> Dielectric data
   real(c_double), intent(in) :: epsv_in
   real(wp), allocatable :: epsv
   !> Molar mass (g/mol)
   real(c_double), intent(in) :: smass_in
   real(wp), allocatable :: smass
   !> Solvent density (g/cm^3)
   real(c_double), intent(in) :: rhos_in
   real(wp), allocatable :: rhos
   !> Born radii
   real(c_double), intent(in) :: c1_in
   real(wp), allocatable :: c1
   !> Atomic surfaces
   real(c_double), intent(in) :: rprobe_in
   real(wp), allocatable :: rprobe
   !> Gshift (gsolv=reference vs. gsolv)
   real(c_double), intent(in) :: gshift_in
   real(wp), allocatable :: gshift
   !> offset parameter (fitted)
   real(c_double), intent(in) :: soset_in
   real(wp), allocatable :: soset
   real(c_double), intent(in) :: dum_in
   real(wp), allocatable :: dum
   !> Surface tension (mN/m=dyn/cm)
   real(c_double), intent(in) :: gamscale_in(94)
   real(wp), allocatable :: gamscale(:)
   !> dielectric descreening parameters
   real(c_double), intent(in) :: sx_in(94)
   real(wp), allocatable :: sx(:)
   real(c_double), intent(in) :: tmp_in(94)
   real(wp), allocatable :: tmp(:)
   integer(c_int) :: status

   call c_get(epsv,epsv_in)
   call c_get(smass,smass_in)
   call c_get(rhos,rhos_in)
   call c_get(c1,c1_in)
   call c_get(rprobe,rprobe_in)
   call c_get(gshift,gshift_in)
   call c_get(soset,soset_in)
   call c_get(dum,dum_in)
   call c_get(gamscale,gamscale_in)
   call c_get(sx,sx_in)
   call c_get(tmp,tmp_in)

   call load_custom_parameters(epsv=epsv)
   call load_custom_parameters(smass=smass)
   call load_custom_parameters(rhos=rhos)
   call load_custom_parameters(c1=c1)
   call load_custom_parameters(rprobe=rprobe)
   call load_custom_parameters(gshift=gshift)
   call load_custom_parameters(soset=soset)
   call load_custom_parameters(dum=dum)
   call load_custom_parameters(gamscale=gamscale)
   call load_custom_parameters(sx=sx)
   call load_custom_parameters(tmp=tmp)

   status = 0

end function gbsa_model_preload_api


end module xtb_api_preload
