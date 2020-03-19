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

!> Enduser bindings for the xtb API. This module tries to provide the most
!  convenient interface to commonly used xtb methods.
module xtb_api_interface
   use iso_c_binding
   use xtb_api_structs
   use xtb_api_utils
   use xtb_api_preload
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_type_environment, only : TEnvironment, init
   use xtb_xtb_data
   use xtb_xtb_gfn0
   use xtb_xtb_gfn1
   use xtb_xtb_gfn2

   implicit none

contains

!> Unsafe xTB calculation, this calculation mode requires to setup at least
!  the parametrisation and the molecular structure data.
!
!  In case null-pointers are passed instead of the basisset and wavefunction,
!  those will be constructed on-the-fly and deleted afterwards.
!  Properties like energy, gradient and stress tensor are returned directly,
!  while other properties can be obtained from querying the wavefunction.
integer(c_int) function xtb_calculation_api &
      & (c_mol, c_param, c_basis, c_wfn, c_pcem, c_opt, c_output, &
      &  c_energy, c_gradient, c_stress) &
      & result(status) bind(C, name="xTB_calculation")
   use xtb_setparam, only : gfn_method
   use xtb_type_environment, only : TEnvironment, init
   use xtb_type_molecule
   use xtb_type_basisset
   use xtb_type_wavefunction
   use xtb_type_options
   use xtb_type_data
   use xtb_type_pcem
   use xtb_scf
   use xtb_peeq
   !> Opaque pointer to the molecular structure data
   type(c_ptr), value, intent(in) :: c_mol
   !> Actual molecular structure data
   type(TMolecule), pointer :: mol
   !> Opaque pointer to the paramatrisation, currently not used
   type(c_ptr), value, intent(in) :: c_param
   !> Opaque pointer to the tight binding basisset
   type(c_ptr), value, intent(in) :: c_basis
   !> Actual basisset, constructed on-the-fly if non is given above
   type(TBasisset), pointer :: basis
   !> Opaque pointer to the tight binding wavefunction
   type(c_ptr), value, intent(in) :: c_wfn
   !> Actual wavefunction, constructed on-the-fly if non is given above
   type(TWavefunction), pointer :: wfn
   !> Opaque pointer to the external potential, currently not used
   type(c_ptr), value, intent(in) :: c_pcem
   !> Options for this calculation run as transparent structure
   type(c_scc_options), intent(in) :: c_opt
   !> Converted representation of the calculation options
   type(scc_options) :: opt
   !> File name for output, can also be "-" for standard output
   character(kind=c_char), intent(in) :: c_output(*)
   !> Converted representation of the file name
   character(len=:), allocatable :: output
   !> Pointer to store the total energy (in Hartree)
   real(c_double), intent(out), optional :: c_energy
   !> Pointer to store the molecular gradient (in Hartree/Bohr)
   real(c_double), intent(out), optional :: c_gradient(3, *)
   !> Pointer to store the stress tensor (in Hartree/Bohr³)
   real(c_double), intent(out), optional :: c_stress(3, 3)

   type(TEnvironment) :: env
   type(TxTBData) :: xtbData
   real(wp) :: energy, sigma(3, 3), egap
   real(wp), allocatable :: gradient(:, :)
   type(scc_results) :: res
   type(tb_pcem) :: pcem
   integer(c_int) :: stat_basis
   integer :: iunit
   logical :: exist, sane, exitRun

   call init(env)
   status = 1

   ! perform some sanity checks first, if they fail, quickly return
   if (.not.(check_xTB_init() .and. c_associated(c_mol))) return
   call c_f_pointer(c_mol, mol)
   if (mol%n <= 0) return

   select case(gfn_method)
   case(0)
      call initGFN0(xtbData)
   case(1)
      call initGFN1(xtbData)
   case(2)
      call initGFN2(xtbData)
   end select

   energy = 0.0_wp
   sigma = 0.0_wp
   allocate(gradient(3, mol%n), source=0.0_wp)

   ! start converting transparent C-types
   opt = c_opt

   ! open a unit for IO
   call c_open_file(iunit, c_output, opt%prlevel)

   ! print the xtb banner with version number and compilation date
   if (opt%prlevel > 2) then
      call xtb_header(iunit)
      call disclamer(iunit)
      call citation(iunit)
   endif

   ! handle basisset
   if (c_associated(c_basis)) then
      call c_f_pointer(c_basis, basis)
   else
      allocate(basis)
      call new_xtb_basisset(mol, basis, stat_basis)
      if (stat_basis /= 0) then
         ! don't call finalize here
         deallocate(basis)
         status = stat_basis
         return
      end if
   end if

   ! handle wavefunction
   if (c_associated(c_wfn)) then
      call c_f_pointer(c_wfn, wfn)
      if (.not.verify_xtb_wavefunction(basis, wfn)) then
         call finalize
         return
      end if
   else
      allocate(wfn)
      call wfn%allocate(basis%n, basis%nshell, basis%nao)
      wfn%nel = nint(sum(mol%z) - mol%chrg)
      wfn%nopen = mol%uhf
      if (mod(wfn%nopen, 2) == 0 .and. mod(wfn%nel, 2) /= 0) wfn%nopen = 1
      if (mod(wfn%nopen, 2) /= 0 .and. mod(wfn%nel, 2) == 0) wfn%nopen = 0
   end if

   ! perform actual calculation
   select case(gfn_method)
   case(0)
      call peeq &
         & (env, mol, wfn, basis, global_parameter, xtbData, &
         &  egap, opt%etemp, opt%prlevel, .false., opt%ccm, opt%acc, &
         &  energy, gradient, sigma, res)
   case(1)
      call scf &
         & (env, mol, wfn, basis, global_parameter, pcem, xtbData, &
         &  egap, opt%etemp, opt%maxiter, opt%prlevel, .false., .false., opt%acc, &
         &  energy, gradient, res)
   case(2)
      call scf &
         & (env, mol, wfn, basis, global_parameter, pcem, xtbData, &
         &  egap, opt%etemp, opt%maxiter, opt%prlevel, .false., .false., opt%acc, &
         &  energy, gradient, res)
   case default
      status = 3
      call finalize
      return
   end select

   call env%check(exitRun)
   if (exitRun) then
      call env%show("Single point calculator terminated")
      call finalize
      status = 4
      return
   end if

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then
      status = 4
      call finalize
      return
   endif

   ! return values to C
   if (present(c_energy)) c_energy = energy
   if (present(c_gradient)) c_gradient(1:3, 1:mol%n) = gradient
   if (present(c_stress)) c_stress = sigma / mol%volume

   ! cleanup and return
   call finalize
   status = 0

contains
subroutine finalize
   if (iunit /= stdout) close(iunit)
   if (.not.c_associated(c_basis)) then
      call basis%deallocate
      deallocate(basis)
      nullify(basis)
   end if
   if (.not.c_associated(c_wfn)) then
      call wfn%deallocate
      deallocate(wfn)
      nullify(wfn)
   end if
end subroutine finalize
end function xtb_calculation_api

function peeq_api &
      &   (natoms,attyp,charge,uhf,coord,lattice,pbc,opt_in,file_in, &
      &    etot,grad,stress,glat) &
      &    result(status) bind(C,name="GFN0_PBC_calculation")

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_options

   use xtb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   real(c_double), intent(in) :: lattice(3,3)
   logical(c_bool),intent(in) :: pbc(3)
   type(c_peeq_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: stress(3,3)
   real(c_double),intent(out) :: glat(3,3)

   type(TMolecule)    :: mol
   type(peeq_options)   :: opt
   type(TEnvironment) :: env

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane, exitRun
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: lattice_gradient(3,3)
   real(wp) :: stress_tensor(3,3)
   real(wp),allocatable :: gradient(:,:)

   status = load_xtb_parameters_api(0_c_int)
   if (status /= 0) return

   call init(env)

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   ! open a unit for IO
   call c_open_file(iunit, file_in, opt%prlevel)

   if (opt%prlevel > 2) then
      ! print the xtb banner with version number and compilation date
      call xtb_header(iunit)
      ! make sure you cannot blame us for destroying your computer
      call disclamer(iunit)
      ! how to cite this program
      call citation(iunit)
   endif

   ! ====================================================================
   !  STEP 3: aquire the molecular structure and fill with data from C
   ! ====================================================================
   mol = TMolecule(natoms, attyp, coord, charge, uhf, lattice, pbc)

   ! ====================================================================
   !  STEP 4: reserve some memory
   ! ====================================================================
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   ! ====================================================================
   !  STEP 5: call the actual Fortran API to perform the calculation
   ! ====================================================================
   ! shut down fatal errors from the MCTC library, so it will not kill the caller
   call mctc_mute

   call gfn0_calculation &
      (iunit,env,opt,mol,hl_gap,energy,gradient,stress_tensor,lattice_gradient)

   call env%check(exitRun)
   if (exitRun) then
      call env%show("Single point calculator terminated")
      call finalize
      status = 1
      return
   end if

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   call c_return(etot, energy)
   call c_return(grad, gradient)
   call c_return(glat, lattice_gradient)
   call c_return(stress, stress_tensor)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      deallocate(gradient)
      if (iunit.ne.stdout) close(iunit)
   end subroutine finalize
end function peeq_api

function gfn2_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in, &
      &    etot,grad,dipole,q,dipm,qp,wbo) &
      &    result(status) bind(C,name="GFN2_calculation")

   use xtb_type_options

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: q(natoms)
   real(c_double),intent(out) :: wbo(natoms,natoms)
   real(c_double),intent(out) :: dipm(3,natoms)
   real(c_double),intent(out) :: qp(6,natoms)
   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: dipole(3)
   type(TxTBData) :: xtbData

   status = load_xtb_parameters_api(2_c_int)
   if (status /= 0) return

   call mctc_init('peeq',10,.true.)

   call initGFN2(xtbData)
   status = gfn12_calc_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData,etot,grad,dipole,q,wbo,dipm,qp)

end function gfn2_api

function gfn1_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,etot,grad,dipole,q,wbo) &
      &    result(status) bind(C,name="GFN1_calculation")

   use xtb_type_options

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: q(natoms)
   real(c_double),intent(out) :: wbo(natoms,natoms)
   real(c_double),intent(out) :: dipole(3)
   type(TxTBData) :: xtbData

   status = load_xtb_parameters_api(1_c_int)
   if (status /= 0) return

   call mctc_init('peeq',10,.true.)

   call initGFN1(xtbData)
   status = gfn12_calc_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData,etot,grad,dipole,q,wbo)

end function gfn1_api

function gfn12_calc_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData,etot,grad,dipole,q,wbo,dipm,qp) &
      &    result(status)

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_options
   use xtb_type_pcem
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_data

   use xtb_setparam, only : ngrida, gfn_method
   use xtb_scf
   use xtb_solv_gbobc

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)
   type(TxTBData), intent(in) :: xtbData

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: q(natoms)
   real(c_double),intent(out) :: wbo(natoms,natoms)
   real(c_double),intent(out) :: dipole(3)
   real(c_double),intent(out),optional :: dipm(3,natoms)
   real(c_double),intent(out),optional :: qp(6,natoms)

   type(TMolecule) :: mol
   type(scc_options) :: opt
   type(TEnvironment) :: env
   type(tb_pcem) :: pcem
   type(TWavefunction) :: wfn
   type(TBasisset) :: basis
   type(scc_results) :: res

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane, exitRun
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)
   integer(c_int) :: stat_basis

   call init(env)

   ! convert the options from C struct to actual Fortran type
   opt = opt_in

   ! open a unit for IO
   call c_open_file(iunit, file_in, opt%prlevel)

   if (opt%prlevel > 2) then
      ! print the xtb banner with version number and compilation date
      call xtb_header(iunit)
      ! make sure you cannot blame us for destroying your computer
      call disclamer(iunit)
      ! how to cite this program
      call citation(iunit)
   endif

   ! aquire the molecular structure and fill with data from C
   mol = TMolecule(natoms, attyp, coord, charge, uhf)

   ! reserve some memory
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   ! setup solvent model
   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none"
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! call the actual Fortran API to perform the calculation
   call new_xtb_basisset(mol, basis, stat_basis)
   if (stat_basis /= 0) then
      status = stat_basis
      call finalize
      return
   end if

   call wfn%allocate(basis%n, basis%nshell, basis%nao)
   wfn%nel = nint(sum(mol%z) - mol%chrg)
   wfn%nopen = mol%uhf
   if (mod(wfn%nopen, 2) == 0 .and. mod(wfn%nel, 2) /= 0) wfn%nopen = 1
   if (mod(wfn%nopen, 2) /= 0 .and. mod(wfn%nel, 2) == 0) wfn%nopen = 0

   call eeq_guess_wavefunction(env, mol, wfn, xtbData)

   call scf &
      & (env, mol, wfn, basis, global_parameter, pcem, xtbData, &
      &  hl_gap, opt%etemp, opt%maxiter, opt%prlevel, .false., .false., opt%acc, &
      &  energy, gradient, res)

   call env%check(exitRun)
   if (exitRun) then
      call env%show("Single point calculator terminated")
      call finalize
      status = 1
      return
   end if

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   ! finally return values to C
   call c_return(etot, energy)
   call c_return(grad, gradient)
   call c_return(q, wfn%q)
   call c_return(wbo, wfn%wbo)
   call c_return(dipole, sum(wfn%dipm,dim=2) + matmul(mol%xyz,wfn%q))
   if (present(dipm)) call c_return(dipm, wfn%dipm)
   if (present(qp)) call c_return(qp, wfn%qp)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call wfn%deallocate
      call basis%deallocate
      deallocate(gradient)
      if (iunit.ne.stdout) close(iunit)
   end subroutine finalize
end function gfn12_calc_impl

function gfn0_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN0_calculation")

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_options

   use xtb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_peeq_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)

   type(TMolecule)    :: mol
   type(peeq_options)    :: opt
   type(TEnvironment) :: env

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane, exitRun
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: dum(3,3)
   real(wp),allocatable :: gradient(:,:)

   status = load_xtb_parameters_api(0_c_int)
   if (status /= 0) return

   call mctc_init('peeq',10,.true.)

   call init(env)

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   ! open a unit for IO
   call c_open_file(iunit, file_in, opt%prlevel)

   if (opt%prlevel > 2) then
      ! print the xtb banner with version number and compilation date
      call xtb_header(iunit)
      ! make sure you cannot blame us for destroying your computer
      call disclamer(iunit)
      ! how to cite this program
      call citation(iunit)
   endif

   ! ====================================================================
   !  STEP 3: aquire the molecular structure and fill with data from C
   ! ====================================================================
   mol = TMolecule(natoms, attyp, coord, charge, uhf)

   ! ====================================================================
   !  STEP 4: reserve some memory
   ! ====================================================================
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   ! ====================================================================
   !  STEP 5: call the actual Fortran API to perform the calculation
   ! ====================================================================
   ! shut down fatal errors from the MCTC library, so it will not kill the caller
   call mctc_mute

   call gfn0_calculation &
      (iunit,env,opt,mol,hl_gap,energy,gradient,dum,dum)

   call env%check(exitRun)
   if (exitRun) then
      call env%show("Single point calculator terminated")
      call finalize
      status = 1
      return
   end if

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   call c_return(etot, energy)
   call c_return(grad, gradient)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      deallocate(gradient)
      if (iunit.ne.stdout) close(iunit)
   end subroutine finalize
end function gfn0_api

function gfn2_pcem_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN2_QMMM_calculation")

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_options
   use xtb_type_pcem
   use xtb_type_wavefunction

   use xtb_aoparam

   use xtb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int), intent(in) :: npc
   real(c_double), intent(in) :: pc_q(npc)
   integer(c_int), intent(in) :: pc_at(npc)
   real(c_double), intent(in) :: pc_gam(npc)
   real(c_double), intent(in) :: pc_coord(3,npc)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: pc_grad(3,npc)
   type(TxTBData) :: xtbData

   status = load_xtb_parameters_api(2_c_int)
   if (status /= 0) return

   call mctc_init('peeq',10,.true.)

   call initGFN2(xtbData)
   status = gfn12_pcem_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad)

end function gfn2_pcem_api

function gfn1_pcem_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN1_QMMM_calculation")

   use xtb_type_options

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int), intent(in) :: npc
   real(c_double), intent(in) :: pc_q(npc)
   integer(c_int), intent(in) :: pc_at(npc)
   real(c_double), intent(in) :: pc_gam(npc)
   real(c_double), intent(in) :: pc_coord(3,npc)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: pc_grad(3,npc)
   type(TxTBData) :: xtbData

   status = load_xtb_parameters_api(1_c_int)
   if (status /= 0) return

   call mctc_init('peeq',10,.true.)

   call initGFN1(xtbData)
   status = gfn12_pcem_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad)

end function gfn1_pcem_api

function gfn12_pcem_impl &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,xtbData, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status)

   use xtb_type_molecule
   use xtb_type_param
   use xtb_type_options
   use xtb_type_pcem
   use xtb_type_wavefunction
   use xtb_type_basisset
   use xtb_type_data

   use xtb_aoparam
   use xtb_setparam, only : gfn_method, ngrida
   use xtb_scf
   use xtb_solv_gbobc

   use xtb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   integer(c_int), intent(in) :: uhf
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)
   type(TxTBData), intent(in) :: xtbData

   integer(c_int), intent(in) :: npc
   real(c_double), intent(in) :: pc_q(npc)
   integer(c_int), intent(in) :: pc_at(npc)
   real(c_double), intent(in) :: pc_gam(npc)
   real(c_double), intent(in) :: pc_coord(3,npc)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: pc_grad(3,npc)

   type(TMolecule) :: mol
   type(scc_options) :: opt
   type(TWavefunction) :: wfn
   type(TBasisset) :: basis
   type(TEnvironment) :: env
   type(tb_pcem) :: pcem

   character(len=:),allocatable :: outfile

   integer(c_int) :: stat_basis
   type(scc_results) :: res
   integer  :: iunit
   logical  :: sane, exitRun
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   call init(env)

   ! convert the options from C struct to actual Fortran type
   opt = opt_in

   ! open a unit for IO
   call c_open_file(iunit, file_in, opt%prlevel)

   if (opt%prlevel > 2) then
      ! print the xtb banner with version number and compilation date
      call xtb_header(iunit)
      ! make sure you cannot blame us for destroying your computer
      call disclamer(iunit)
      ! how to cite this program
      call citation(iunit)
   endif

   ! aquire the molecular structure and fill with data from C
   mol = TMolecule(natoms, attyp, coord, charge, uhf)

   call pcem%allocate(npc)
   pcem%q   = pc_q
   pcem%xyz = pc_coord
   where(pc_at > 0)
      pcem%gam = gam(pc_at)
   elsewhere
      pcem%gam = pc_gam
   endwhere

   ! reserve some memory
   allocate(gradient(3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp

   ! setup solvent model
   lgbsa = len_trim(opt%solvent).gt.0 .and. opt%solvent.ne."none"
   if (lgbsa) then
      call init_gbsa(iunit,trim(opt%solvent),0,opt%etemp,gfn_method,ngrida)
   endif

   ! call the actual Fortran API to perform the calculation
   call new_xtb_basisset(mol, basis, stat_basis)
   if (stat_basis /= 0) then
      status = stat_basis
      call finalize
      return
   end if

   call wfn%allocate(basis%n, basis%nshell, basis%nao)
   wfn%nel = nint(sum(mol%z) - mol%chrg)
   wfn%nopen = mol%uhf
   if (mod(wfn%nopen, 2) == 0 .and. mod(wfn%nel, 2) /= 0) wfn%nopen = 1
   if (mod(wfn%nopen, 2) /= 0 .and. mod(wfn%nel, 2) == 0) wfn%nopen = 0

   call eeq_guess_wavefunction(env, mol, wfn, xtbData)

   call scf &
      & (env, mol, wfn, basis, global_parameter, pcem, xtbData, &
      &  hl_gap, opt%etemp, opt%maxiter, opt%prlevel, .false., .false., opt%acc, &
      &  energy, gradient, res)

   call env%check(exitRun)
   if (exitRun) then
      call env%show("Single point calculator terminated")
      call finalize
      status = 1
      return
   end if

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   ! finally return values to C
   call c_return(etot, energy)
   call c_return(grad, gradient)
   call c_return(pc_grad, pcem%grd)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call pcem%deallocate
      call wfn%deallocate
      deallocate(gradient)
      if (iunit.ne.stdout) close(iunit)
   end subroutine finalize
end function gfn12_pcem_impl

function gbsa_calculation_api &
      &  (natoms,attyp,coord, &
      &   solvent_in,reference,temperature,method,grid_size,file_in, &
      &   brad,sasa) &
      &  result(status) bind(C,name='GBSA_calculation')

   use xtb_mctc_constants

   use xtb_type_molecule
   use xtb_type_options

   use xtb_solv_gbobc

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: coord(3,natoms)

   character(kind=c_char), intent(in) :: file_in(*)
   character(kind=c_char), intent(in) :: solvent_in(*)
   integer(c_int), intent(in) :: reference
   real(c_double), intent(in) :: temperature
   integer(c_int), intent(in) :: method
   integer(c_int), intent(in) :: grid_size
   integer(c_int) :: status

   real(c_double), intent(out) :: brad(natoms)
   real(c_double), intent(out) :: sasa(natoms)

   integer :: iunit
   logical :: sane, exitRun
   character(len=:), allocatable :: outfile
   character(len=:), allocatable :: solvent

   type(TMolecule)    :: mol
   type(TSolvent)     :: gbsa
   type(TEnvironment) :: env

   call mctc_init('gbobc',10,.true.)
   call init(env)

   call c_string_convert(solvent,solvent_in)
   call c_string_convert(outfile,file_in)

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = stdout
      endif
   else
      iunit = stdout
   endif

   ! shut down fatal errors from the MCTC library, so it will not kill the caller
   call mctc_mute

   call init_gbsa(iunit,solvent,reference,temperature,method,grid_size)

   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call mol%update

   call new_gbsa(gbsa,mol%n,mol%at)
   ! initialize the neighbor list
   call update_nnlist_gbsa(gbsa,mol%xyz,.false.)
   ! compute Born radii
   call compute_brad_sasa(gbsa,mol%xyz)

   call mctc_sanity(sane)
   if (.not.sane) then
      call finalize
      status = 1
      return
   endif

   call c_return(brad, gbsa%brad)
   call c_return(sasa, gbsa%sasa*fourpi/gbsa%gamsasa)

   call finalize
   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call deallocate_gbsa(gbsa)
      if (iunit.ne.stdout) call close_file(iunit)
   end subroutine finalize
end function gbsa_calculation_api

end module xtb_api_interface
