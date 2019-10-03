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

module c_api
   use iso_c_binding
   use iso_fortran_env, only: wp => real64, istdout => output_unit

   implicit none

   !> some overloading for convenience
   interface c_return
      module procedure :: c_return_double_0d
      module procedure :: c_return_double_1d
      module procedure :: c_return_double_2d
   end interface c_return

   interface c_return_alloc
      module procedure :: c_return_double_0d_alloc
      module procedure :: c_return_double_1d_alloc
      module procedure :: c_return_double_2d_alloc
   end interface c_return_alloc

   interface c_get
      module procedure :: c_get_double_0d_alloc
      module procedure :: c_get_double_1d_alloc
      module procedure :: c_get_double_2d_alloc
      module procedure :: c_get_double_0d
      module procedure :: c_get_double_1d
      module procedure :: c_get_double_2d
      module procedure :: c_get_int_0d
      module procedure :: c_get_int_1d
      module procedure :: c_get_int_2d
   end interface c_get

contains

function peeq_api &
      &   (natoms,attyp,charge,uhf,coord,lattice,pbc,opt_in,file_in, &
      &    etot,grad,stress,glat) &
      &    result(status) bind(C,name="GFN0_PBC_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(gfn_parameter)  :: gfn
   type(peeq_options)   :: opt
   type(tb_environment) :: env

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: lattice_gradient(3,3)
   real(wp) :: stress_tensor(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! set periodicity of system
   mol%lattice = lattice
   mol%pbc = pbc
   mol%npbc = count(pbc)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

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
      (iunit,env,opt,mol,gfn,hl_gap,energy,gradient,stress_tensor,lattice_gradient)

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
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function peeq_api

function gfn2_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in, &
      &    etot,grad,dipole,q,dipm,qp,wbo) &
      &    result(status) bind(C,name="GFN2_calculation")

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(tb_wavefunction):: wfn
   type(gfn_parameter)  :: gfn
   type(scc_options)    :: opt
   type(tb_environment) :: env
   type(tb_pcem)        :: pcem

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

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

   call gfn2_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   call c_return(q, wfn%q)
   call c_return(dipm, wfn%dipm)
   call c_return(qp, wfn%qp)
   call c_return(wbo, wfn%wbo)
   call c_return(dipole, sum(wfn%dipm,dim=2) + matmul(mol%xyz,wfn%q))

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call wfn%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn2_api

function gfn1_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,etot,grad,dipole,q,wbo) &
      &    result(status) bind(C,name="GFN1_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem
   use tbdef_wavefunction

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(gfn_parameter)  :: gfn
   type(scc_options)    :: opt
   type(tb_environment) :: env
   type(tb_pcem)        :: pcem
   type(tb_wavefunction):: wfn

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

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

   call gfn1_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   call c_return(q, wfn%q)
   call c_return(wbo, wfn%wbo)
   call c_return(dipole, sum(wfn%dipm,dim=2) + matmul(mol%xyz,wfn%q))

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call wfn%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn1_api

function gfn0_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN0_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(gfn_parameter)  :: gfn
   type(peeq_options)    :: opt
   type(tb_environment) :: env

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: dum(3,3)
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

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
      (iunit,env,opt,mol,gfn,hl_gap,energy,gradient,dum,dum)

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
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn0_api

function gfn2_pcem_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN2_QMMM_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem
   use tbdef_wavefunction

   use aoparam

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(gfn_parameter)  :: gfn
   type(scc_options)    :: opt
   type(tb_wavefunction):: wfn
   type(tb_environment) :: env
   type(tb_pcem)        :: pcem

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

   call pcem%allocate(npc)
   pcem%q   = pc_q
   pcem%xyz = pc_coord
   where(pc_at > 0)
      pcem%gam = gam(pc_at)
   elsewhere
      pcem%gam = pc_gam
   endwhere

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

   call gfn2_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   call c_return(pc_grad, pcem%grd)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call pcem%deallocate
      call wfn%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn2_pcem_api

function gfn1_pcem_api &
      &   (natoms,attyp,charge,uhf,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN1_QMMM_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem
   use tbdef_wavefunction

   use aoparam

   use tb_calculators

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

   type(tb_molecule)    :: mol
   type(gfn_parameter)  :: gfn
   type(scc_options)    :: opt
   type(tb_environment) :: env
   type(tb_pcem)        :: pcem
   type(tb_wavefunction):: wfn

   character(len=:),allocatable :: outfile

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp),allocatable :: gradient(:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   call c_string_convert(outfile, file_in)

   if (outfile.ne.'-'.and.opt%prlevel > 0) then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

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
   call mol%allocate(natoms)
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   call c_get(mol%chrg, charge, 0.0_wp)
   call c_get(mol%uhf, uhf, 0)

   ! update everything from xyz and lattice
   call mol%set_nuclear_charge
   call mol%set_atomic_masses
   call mol%update

   call pcem%allocate(npc)
   pcem%q   = pc_q
   pcem%xyz = pc_coord
   where(pc_at > 0)
      pcem%gam = gam(pc_at)
   elsewhere
      pcem%gam = pc_gam
   endwhere

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

   call gfn1_calculation &
      (iunit,env,opt,mol,gfn,pcem,wfn,hl_gap,energy,gradient)

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
   call c_return(pc_grad, pcem%grd)

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call pcem%deallocate
      call wfn%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn1_pcem_api

!> allows loading custom solvent parameters in the shared library
function gbsa_model_preload_api &
      &  (epsv_in,smass_in,rhos_in,c1_in,rprobe_in,gshift_in,soset_in,dum_in, &
      &   gamscale_in,sx_in,tmp_in) &
      &  result(status) bind(C,name='GBSA_model_preload')

   use gbobc, only: load_custom_parameters

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

function gbsa_calculation_api &
      &  (natoms,attyp,coord, &
      &   solvent_in,reference,temperature,method,grid_size,file_in, &
      &   brad,sasa) &
      &  result(status) bind(C,name='GBSA_calculation')

   use mctc_constants

   use tbdef_molecule
   use tbdef_options

   use gbobc

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
   logical :: sane
   character(len=:), allocatable :: outfile
   character(len=:), allocatable :: solvent

   type(tb_molecule)    :: mol
   type(tb_solvent)     :: gbsa
   type(tb_environment) :: env

   call mctc_init('gbobc',10,.true.)
   call env%setup

   call c_string_convert(solvent,solvent_in)
   call c_string_convert(outfile,file_in)

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
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
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gbsa_calculation_api

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
subroutine c_return_double_0d(c_array, f_array)
   real(c_double), intent(out), target :: c_array
   real(wp), intent(in) :: f_array
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_0d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran value has been calculated
subroutine c_return_double_0d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array
   real(wp), allocatable, intent(in) :: f_array
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_0d_alloc

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated
subroutine c_return_double_1d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:)
   real(wp), intent(in) :: f_array(:)
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_1d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated
subroutine c_return_double_1d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:)
   real(wp), allocatable, intent(in) :: f_array(:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_1d_alloc

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated (2D version)
subroutine c_return_double_2d(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:,:)
   real(wp), intent(in) :: f_array(:,:)
   if (c_associated(c_loc(c_array))) then
      c_array = f_array
   endif
end subroutine c_return_double_2d

!> optional return to a c_ptr in case it is not a null pointer and the
!  Fortran array has been populated (2D version)
subroutine c_return_double_2d_alloc(c_array, f_array)
   real(c_double), intent(out), target :: c_array(:,:)
   real(wp), allocatable, intent(in) :: f_array(:,:)
   if (c_associated(c_loc(c_array)) .and. allocated(f_array)) then
      c_array = f_array
   endif
end subroutine c_return_double_2d_alloc

!> convert vector of chars to deferred size character
subroutine c_string_convert(f_string, c_string)
   character(c_char), dimension(*), intent(in) :: c_string
   character(len=:), allocatable, intent(out) :: f_string
   integer :: i
   i = 0
   f_string = ''
   do
      i = i+1
      if (c_string(i).eq.c_null_char) exit
      f_string = f_string//c_string(i)
   enddo
end subroutine c_string_convert

subroutine c_get_double_0d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array
   real(wp), allocatable, intent(out) :: f_array
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_0d_alloc

subroutine c_get_double_1d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:)
   real(wp), allocatable, intent(out) :: f_array(:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_1d_alloc

subroutine c_get_double_2d_alloc(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:,:)
   real(wp), allocatable, intent(out) :: f_array(:,:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_2d_alloc

subroutine c_get_double_0d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array
   real(wp), intent(out) :: f_array
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_0d

subroutine c_get_double_1d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array(:)
   real(wp), intent(out) :: f_array(:)
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_1d

subroutine c_get_double_2d(f_array, c_array, default)
   real(c_double), intent(in), target :: c_array(:,:)
   real(wp), intent(out) :: f_array(:,:)
   real(wp), intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_double_2d

subroutine c_get_int_0d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array
   integer, intent(out) :: f_array
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_0d

subroutine c_get_int_1d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array(:)
   integer, intent(out) :: f_array(:)
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_1d

subroutine c_get_int_2d(f_array, c_array, default)
   integer(c_int), intent(in), target :: c_array(:,:)
   integer, intent(out) :: f_array(:,:)
   integer, intent(in) :: default
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   else
      f_array = default
   endif
end subroutine c_get_int_2d

end module c_api
