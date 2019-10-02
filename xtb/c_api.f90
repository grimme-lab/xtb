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

   !> some overloading for convience
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
      module procedure :: c_get_double_0d
      module procedure :: c_get_double_1d
      module procedure :: c_get_double_2d
   end interface c_get

contains

function peeq_api &
      &   (natoms,attyp,charge,coord,lattice,pbc,opt_in,file_in,etot,grad,glat) &
      &    result(status) bind(C,name="GFN0_PBC_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

   use tb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   real(c_double), intent(in) :: coord(3,natoms)
   real(c_double), intent(in) :: lattice(3,3)
   logical(c_bool),intent(in) :: pbc(3)
   type(c_peeq_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
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
   mol%chrg = charge

   ! update everything from xyz and lattice
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

   call peeq_calculation &
      (iunit,env,opt,mol,gfn,hl_gap,energy,gradient,lattice_gradient)

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
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
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
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
   mol%chrg = charge
   call mol%calculate_distances

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
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
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
      &   (natoms,attyp,charge,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN1_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

   use tb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   real(c_double), intent(in) :: coord(3,natoms)
   type(c_scc_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)

   type(tb_molecule)    :: mol
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
   mol%chrg = charge
   call mol%calculate_distances

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
      (iunit,env,opt,mol,gfn,pcem,hl_gap,energy,gradient)

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
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
end function gfn1_api

function gfn0_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN0_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

   use tb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
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
   mol%chrg = charge
   call mol%calculate_distances

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
      (iunit,env,opt,mol,gfn,hl_gap,energy,gradient)

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
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
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
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
   mol%chrg = charge
   call mol%calculate_distances

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
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
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
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn2_pcem_api

function gfn1_pcem_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN1_QMMM_calculation")

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

   use aoparam

   use tb_calculators

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
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
   mol%chrg = charge
   call mol%calculate_distances

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
      (iunit,env,opt,mol,gfn,pcem,hl_gap,energy,gradient)

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
      call finalize

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   pc_grad = pcem%grd

   call finalize

   status = 0

contains
   subroutine finalize
      call mol%deallocate
      call pcem%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)
   end subroutine finalize
end function gfn1_pcem_api

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

subroutine c_get_double_0d(f_array, c_array)
   real(c_double), intent(in), target :: c_array
   real(wp), allocatable, intent(out) :: f_array
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_0d

subroutine c_get_double_1d(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:)
   real(wp), allocatable, intent(out) :: f_array(:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_1d

subroutine c_get_double_2d(f_array, c_array)
   real(c_double), intent(in), target :: c_array(:,:)
   real(wp), allocatable, intent(out) :: f_array(:,:)
   if (c_associated(c_loc(c_array))) then
      f_array = c_array
   endif
end subroutine c_get_double_2d

end module c_api
