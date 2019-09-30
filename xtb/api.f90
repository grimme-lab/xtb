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

function peeq_api &
      &   (natoms,attyp,charge,coord,lattice,pbc,opt_in,file_in,etot,grad,glat) &
      &    result(status) bind(C,name="GFN0_PBC_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

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
   
   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      call mol%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   glat = lattice_gradient

   call mol%deallocate
   deallocate(gradient)
   if (iunit.ne.istdout) call close_file(iunit)

   status = 0

end function peeq_api

function gfn2_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
      &    etot,grad,dipole,q,dipm,qp,wbo) &
      &    result(status) bind(C,name="GFN2_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

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

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      call mol%deallocate
      deallocate(gradient)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   q = wfn%q
   dipm = wfn%dipm
   qp = wfn%qp
   wbo = wfn%wbo
   dipole = sum(wfn%dipm,dim=2) + matmul(mol%xyz,wfn%q)

   call mol%deallocate
   deallocate(gradient)

   status = 0

end function gfn2_api

function gfn1_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN1_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

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

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      call mol%deallocate
      deallocate(gradient)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient

   call mol%deallocate
   deallocate(gradient)

   status = 0

end function gfn1_api

function gfn0_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in,etot,grad) &
      &    result(status) bind(C,name="GFN0_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

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

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      call mol%deallocate
      deallocate(gradient)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient

   call mol%deallocate
   deallocate(gradient)

   status = 0

end function gfn0_api

function gfn2_pcem_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN2_QMMM_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

   use aoparam

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

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      (iunit,env,opt,mol,gfn,pcem,hl_gap,energy,gradient)

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
      call mol%deallocate
      deallocate(gradient)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   pc_grad = pcem%grd

   call mol%deallocate
   deallocate(gradient)

   status = 0

end function gfn2_pcem_api

function gfn1_pcem_api &
      &   (natoms,attyp,charge,coord,opt_in,file_in, &
      &    npc,pc_q,pc_at,pc_gam,pc_coord,etot,grad,pc_grad) &
      &    result(status) bind(C,name="GFN1_QMMM_calculation")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use tbdef_molecule
   use tbdef_param
   use tbdef_options
   use tbdef_pcem

   use aoparam

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

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

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
      call mol%deallocate
      deallocate(gradient)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   pc_grad = pcem%grd

   call mol%deallocate
   deallocate(gradient)

   status = 0

end function gfn1_pcem_api

function eeq_api &
      &   (natoms,attyp,charge,coord,lattice,pbc,opt_in,file_in, &
      &    qout,dipole,etot,grad,glat) &
      &    result(status) bind(C,name="EEQ_charges")

   use iso_fortran_env, wp => real64, istdout => output_unit
   use iso_c_binding

   use mctc_econv

   use tbdef_molecule
   use tbdef_param
   use tbdef_options

   use ncoord
   use eeq_model
   use pbc_tools

   implicit none

   integer(c_int), intent(in) :: natoms
   integer(c_int), intent(in) :: attyp(natoms)
   real(c_double), intent(in) :: charge
   real(c_double), intent(in) :: coord(3,natoms)
   real(c_double), intent(in) :: lattice(3,3)
   logical(c_bool),intent(in) :: pbc(3)
   type(c_eeq_options), intent(in) :: opt_in
   character(kind=c_char),intent(in) :: file_in(*)

   integer(c_int) :: status

   real(c_double),intent(out) :: etot
   real(c_double),intent(out) :: grad(3,natoms)
   real(c_double),intent(out) :: glat(3,3)
   real(c_double),intent(out) :: qout(natoms)
   real(c_double),intent(out) :: dipole(3)

   type(tb_molecule)    :: mol
   type(chrg_parameter) :: chrgeq
   type(eeq_options)    :: opt
   type(tb_environment) :: env

   character(len=:),allocatable :: outfile

   integer, parameter :: wsc_rep(3) = [1,1,1] ! FIXME

   integer  :: iunit
   logical  :: sane
   integer  :: i
   real(wp) :: energy
   real(wp) :: hl_gap
   real(wp) :: sigma(3,3),inv_lat(3,3)
   real(wp) :: lattice_gradient(3,3)
   real(wp),allocatable :: gradient(:,:)
   real(wp),allocatable :: cn(:)
   real(wp),allocatable :: dcndr(:,:,:)
   real(wp),allocatable :: dcndL(:,:,:)
   real(wp),allocatable :: q(:)
   real(wp),allocatable :: dqdr(:,:,:)
   real(wp),allocatable :: dqdL(:,:,:)

   ! ====================================================================
   !  STEP 1: setup the environment variables
   ! ====================================================================
   call mctc_init('peeq',10,.true.)
   call env%setup

   i = 0
   outfile = ''
   do
      i = i+1
      if (file_in(i).eq.c_null_char) exit
      outfile = outfile//file_in(i)
   enddo

   if (outfile.ne.'-') then
      call open_file(iunit,outfile,'w')
      if (iunit.eq.-1) then
         iunit = istdout
      endif
   else
      iunit = istdout
   endif

   ! ====================================================================
   !  STEP 2: convert the options from C struct to actual Fortran type
   ! ====================================================================
   opt = opt_in

   ! ====================================================================
   !  STEP 3: aquire the molecular structure and fill with data from C
   ! ====================================================================
   call mol%allocate(natoms)
   ! set periodicity of system
   mol%lattice = lattice
   mol%pbc = pbc
   mol%npbc = 0
   do i = 1,3
      if (mol%pbc(i)) mol%npbc = mol%npbc + 1
   enddo
   ! get atomtypes, coordinates and total charge
   mol%at = attyp
   mol%xyz = coord
   mol%chrg = charge
   if (mol%npbc > 0) then
      call generate_wsc(mol,mol%wsc,wsc_rep)
   endif

   ! update everything from xyz and lattice
   call mol%update

   ! ====================================================================
   !  STEP 4: reserve some memory
   ! ====================================================================
   allocate(gradient(3,mol%n),q(mol%n),dqdr(3,mol%n,mol%n+1),dqdL(3,3,mol%n+1), &
      &     cn(mol%n),dcndr(3,mol%n,mol%n),dcndL(3,3,mol%n))
   energy = 0.0_wp
   gradient = 0.0_wp
   sigma = 0.0_wp

   ! ====================================================================
   !  STEP 5: call the actual Fortran API to perform the calculation
   ! ====================================================================
   ! shut down fatal errors from the MCTC library, so it will not kill the caller
   call mctc_mute

   if (mol%npbc > 0) then
      call get_erf_cn(mol,cn,dcndr,dcndL,thr=900.0_wp)
      call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)
   else
      call get_erf_cn(mol,cn,dcndr,thr=900.0_wp)
      call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)
   endif

   call new_charge_model_2019(chrgeq,mol%n,mol%at)

   call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma, &
      &            .false.,.true.,.true.)

   if (mol%npbc > 0) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,lattice_gradient)
   endif

   ! check if the MCTC environment is still sane, if not tell the caller
   call mctc_sanity(sane)
   if (.not.sane) then !! it's stark raving mad and on fire !!

      ! at least try to destroy the molecule structure data
      call mol%deallocate
      deallocate(gradient)
      if (iunit.ne.istdout) call close_file(iunit)

      status = 1
      return
   endif

   ! ====================================================================
   !  STEP 6: finally return values to C
   ! ====================================================================
   etot = energy
   grad = gradient
   glat = lattice_gradient
   qout = q
   dipole = matmul(mol%xyz,q)

   call mol%deallocate
   deallocate(gradient)
   if (iunit.ne.istdout) call close_file(iunit)

   status = 0

end function eeq_api

!function eeeq_api &
!      &   (natoms,attyp,charge,coord,lattice,pbc,opt_in,file_in, &
!      &    chi,gam,kappa,alpha,beta,dpol, &
!      &    qout,dout,dipole,etot,grad,glat) &
!      &    result(status) bind(C,name="EEQ_multipoles")
!
!   use iso_fortran_env, wp => real64, istdout => output_unit
!   use iso_c_binding
!
!   use mctc_econv
!
!   use tbdef_molecule
!   use tbdef_param
!   use tbdef_options
!
!   use ncoord
!   use eeq_model
!   use pbc_tools
!
!   implicit none
!
!   integer(c_int), intent(in) :: natoms
!   integer(c_int), intent(in) :: attyp(natoms)
!   real(c_double), intent(in) :: charge
!   real(c_double), intent(in) :: coord(3,natoms)
!   real(c_double), intent(in) :: lattice(3,3)
!   logical(c_bool),intent(in) :: pbc(3)
!   type(c_eeq_options), intent(in) :: opt_in
!   character(kind=c_char),intent(in) :: file_in(*)
!
!   real(c_double), intent(in) :: chi(*)
!   real(c_double), intent(in) :: gam(*)
!   real(c_double), intent(in) :: kappa(*)
!   real(c_double), intent(in) :: alpha(*)
!   real(c_double), intent(in) :: beta(*)
!   real(c_double), intent(in) :: dpol(*)
!
!   integer(c_int) :: status
!
!   real(c_double),intent(out) :: etot
!   real(c_double),intent(out) :: grad(3,natoms)
!   real(c_double),intent(out) :: glat(3,3)
!   real(c_double),intent(out) :: qout(natoms)
!   real(c_double),intent(out) :: dout(3,natoms)
!   real(c_double),intent(out) :: dipole(3)
!
!   type(tb_molecule)    :: mol
!   type(chrg_parameter) :: chrgeq
!   type(eeq_options)    :: opt
!   type(tb_environment) :: env
!
!   character(len=:),allocatable :: outfile
!
!   integer, parameter :: wsc_rep(3) = [1,1,1] ! FIXME
!
!   integer  :: iunit
!   logical  :: sane
!   integer  :: i
!   real(wp) :: energy
!   real(wp) :: hl_gap
!   real(wp) :: sigma(3,3),inv_lat(3,3)
!   real(wp) :: lattice_gradient(3,3)
!   real(wp),allocatable :: gradient(:,:)
!   real(wp),allocatable :: cn(:)
!   real(wp),allocatable :: dcndr(:,:,:)
!   real(wp),allocatable :: dcndL(:,:,:)
!   real(wp),allocatable :: dipm(:,:)
!   real(wp),allocatable :: q(:)
!   real(wp),allocatable :: dqdr(:,:,:)
!   real(wp),allocatable :: dqdL(:,:,:)
!
!   ! ====================================================================
!   !  STEP 1: setup the environment variables
!   ! ====================================================================
!   call mctc_init('peeq',10,.true.)
!   call env%setup
!
!   i = 0
!   outfile = ''
!   do
!      i = i+1
!      if (file_in(i).eq.c_null_char) exit
!      outfile = outfile//file_in(i)
!   enddo
!
!   if (outfile.ne.'-') then
!      call open_file(iunit,outfile,'w')
!      if (iunit.eq.-1) then
!         iunit = istdout
!      endif
!   else
!      iunit = istdout
!   endif
!
!   ! ====================================================================
!   !  STEP 2: convert the options from C struct to actual Fortran type
!   ! ====================================================================
!   opt = opt_in
!
!   ! ====================================================================
!   !  STEP 3: aquire the molecular structure and fill with data from C
!   ! ====================================================================
!   call mol%allocate(natoms)
!   ! set periodicity of system
!   mol%lattice = lattice
!   mol%pbc = pbc
!   mol%npbc = 0
!   do i = 1,3
!      if (mol%pbc(i)) mol%npbc = mol%npbc + 1
!   enddo
!   ! get atomtypes, coordinates and total charge
!   mol%at = attyp
!   mol%xyz = coord
!   mol%chrg = charge
!   if (mol%npbc > 0) then
!      call generate_wsc(mol,mol%wsc,wsc_rep)
!   endif
!
!   ! update everything from xyz and lattice
!   call mol%update
!
!   ! ====================================================================
!   !  STEP 4: reserve some memory
!   ! ====================================================================
!   allocate(gradient(3,mol%n),q(mol%n),dqdr(3,mol%n,mol%n+1),dqdL(3,3,mol%n+1), &
!      &     dipm(3,mol%n),cn(mol%n),dcndr(3,mol%n,mol%n),dcndL(3,3,mol%n))
!   energy = 0.0_wp
!   gradient = 0.0_wp
!   sigma = 0.0_wp
!
!   ! ====================================================================
!   !  STEP 5: call the actual Fortran API to perform the calculation
!   ! ====================================================================
!   ! shut down fatal errors from the MCTC library, so it will not kill the caller
!   call mctc_mute
!
!   if (mol%npbc > 0) then
!      call get_erf_cn(mol,cn,dcndr,dcndL,thr=900.0_wp)
!      call dncoord_logcn(mol%n,cn,dcndr,dcndL,cn_max=8.0_wp)
!   else
!      call get_erf_cn(mol,cn,dcndr,thr=900.0_wp)
!      call dncoord_logcn(mol%n,cn,dcndr,cn_max=8.0_wp)
!   endif
!
!   call chrgeq%allocate(mol%n,extended=.true.)
!   do i = 1, mol%n
!      chrgeq%en   (i) = chi  (mol%at(i))
!      chrgeq%gam  (i) = gam  (mol%at(i))
!      chrgeq%kappa(i) = kappa(mol%at(i))
!      chrgeq%alpha(i) = alpha(mol%at(i))
!      chrgeq%beta (i) = beta (mol%at(i))
!      chrgeq%dpol (i) = dpol (mol%at(i))
!   enddo
!
!   call eeq_multieq(mol,chrgeq,cn,q,dipm,energy)
!   !call eeq_chrgeq(mol,chrgeq,cn,dcndr,dcndL,q,dqdr,dqdL,energy,gradient,sigma, &
!   !   &            .false.,.true.,.true.)
!
!   if (mol%npbc > 0) then
!      inv_lat = mat_inv_3x3(mol%lattice)
!      call sigma_to_latgrad(sigma,inv_lat,lattice_gradient)
!   endif
!
!   ! check if the MCTC environment is still sane, if not tell the caller
!   call mctc_sanity(sane)
!   if (.not.sane) then !! it's stark raving mad and on fire !!
!
!      ! at least try to destroy the molecule structure data
!      call mol%deallocate
!      deallocate(gradient)
!      if (iunit.ne.istdout) call close_file(iunit)
!
!      status = 1
!      return
!   endif
!
!   ! ====================================================================
!   !  STEP 6: finally return values to C
!   ! ====================================================================
!   etot = energy
!   grad = gradient
!   glat = lattice_gradient
!   qout = q
!   dout = dipm
!   dipole = matmul(mol%xyz,q) + sum(dipm,dim=2)
!
!   call mol%deallocate
!   deallocate(gradient,q,dqdr,dqdL,dipm,cn,dcndr,dcndL)
!   if (iunit.ne.istdout) call close_file(iunit)
!
!   status = 0
!
!end function eeeq_api
