!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
! 2D => 3D structure converter, uses GFN-FF MD + OPT to generate a 3D structure
! of an 2 dimensional input structure.
!------------------------------------------------------------------------------
module xtb_gfnff_convert
contains

subroutine struc_convert( &
         & env,restart,mol,chk,egap,et,maxiter,maxcycle,&
         & etot,g,sigma)
  use xtb_mctc_accuracy, only : wp
  use xtb_gfnff_param
  use xtb_gfnff_setup
  use xtb_disp_dftd3param
  use xtb_type_environment
  use xtb_type_molecule
  use xtb_type_restart
  use xtb_gfnff_calculator, only : TGFFCalculator, newGFFCalculator
  use xtb_type_data
  use xtb_restart
  use xtb_setmod
  use xtb_setparam
  use xtb_dynamic
  use xtb_geoopt
  use xtb_readin, only : xfind
  implicit none
! Dummy -----------------------------------------------------------------------
  type(TEnvironment),intent(inout)            :: env
  type(TMolecule),intent(inout)               :: mol
  type(TRestart),intent(inout)                :: chk
  integer,intent(in)                          :: maxiter
  integer,intent(in)                          :: maxcycle
  real(wp),intent(inout)                      :: etot
  real(wp),intent(in)                         :: et
  real(wp),intent(inout)                      :: egap
  real(wp),intent(inout)                      :: g(3,mol%n)
  real(wp),intent(inout)                      :: sigma(3,3)
  logical,intent(in)                          :: restart
  character(len=:),allocatable                :: fnv
! Stack -----------------------------------------------------------------------
  type(TGFFCalculator)                        :: calc
  integer                                     :: ich
  integer                                     :: idum
  integer                                     :: mode_input
  real(wp)                                    :: time_in
  real(wp)                                    :: temp_in
  real(wp)                                    :: step_in
  real(wp)                                    :: dump_in
  real(wp)                                    :: hmass_in
  logical                                     :: exist
  logical                                     :: fail
  integer, allocatable :: opt_in
  character(len=*),parameter                  :: p_fname_param_gfnff = '.param_gfnff.xtb'
! loop geoopt -----------------------------------------------------------------
  integer                                     :: i, j
  integer, parameter                          :: num_shift_runs = 3
  real(wp), allocatable                       :: shift(:)
  real(wp), allocatable                       :: mol_xyz_arr(:,:,:)
  real(wp)                                    :: etot_arr(num_shift_runs)
  real(wp)                                    :: sign_threshold
  type(TMolecule)                             :: mol_shifted
  allocate(shift(mol%n), source=0.0_wp)
  allocate(mol_xyz_arr(3, mol%n, num_shift_runs), source=0.0_wp)
  etot_arr = 0.0_wp
!------------------------------------------------------------------------------
! set up force field
  call struc_convert_header(env%unit)
  if (allocated(set%opt_engine)) then
    opt_in = set%opt_engine
  end if
  set%opt_engine = p_engine_rf
  mode_input = set%mode_extrun
  set%mode_extrun = p_ext_gfnff
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  call newGFFCalculator(env, mol, calc, fnv, restart, gffVersion%harmonic2020)

!===============================
! Set Block
  time_in  = set%time_md
  set%time_md  = 5.0_wp              ! short 5 ps MD to move it in 3D
  temp_in  = set%temp_md
  set%temp_md  = 298.0_wp            ! md temperature 298 K
  step_in  = set%tstep_md
  set%tstep_md = 2.5_wp              ! md time step 2.5 fs
  dump_in  = set%dump_md2
  set%dump_md2 = 100.0_wp            ! md dump 100 fs
  hmass_in = set%md_hmass
  set%md_hmass = 4.0_wp              ! md hydrogen mass
  set%nvt_md = .true.                ! md thermostat
!===============================
  if (allocated(set%opt_logfile)) then
    fnv = set%opt_logfile
  else
    deallocate(fnv)
  endif
  set%opt_logfile = 'convert.log'
!------------------------------------------------------------------------------
! force field geometry optimization
  ! loop runs 3 geoopt with different shifts in the new 3rd coordinate
  ! and then keeps the mol%xyz with lowest etot for md
  mol_shifted = mol
  do i=1, num_shift_runs
    if (i.ne.42) then  ! if in case you want an opt with the 2D init -> e.g. use if (i.ne.1)
      mol_shifted = mol
      ! create array with random shifts: 0<= shift <1
      call RANDOM_NUMBER(shift)
      ! change signs of shifts randomly
      do j=1, size(shift)
        call RANDOM_NUMBER(sign_threshold)
        if (sign_threshold.lt.0.5_wp) shift(j) = -1.0_wp*shift(j)
      enddo
      mol_shifted%xyz(3,:) = mol_shifted%xyz(3,:) + shift  ! apply shifts
    endif
    
    call geometry_optimization &
        &     (env,mol_shifted,chk,calc,   &
        &      egap,set%etemp,maxiter,maxcycle,etot,g,sigma,p_olev_crude,.false.,.true.,fail)
    mol_xyz_arr(:,:,i) = mol_shifted%xyz  ! store optimized xyz
    etot_arr(i) = etot                    ! store energy etot
  enddo

  mol%xyz(:,:) = mol_xyz_arr(:,:,minloc(etot_arr, DIM=1))  ! keep xyz with lowest etot

    if (allocated(fnv)) then
      set%opt_logfile = fnv
    else
      deallocate(set%opt_logfile)
    endif
    write(*,*)
!------------------------------------------------------------------------------
! force field md simulation
  idum = 0
  call md                &
      &   (env,mol,chk,calc, &
      &    egap,set%etemp,maxiter,etot,g,sigma,0,set%temp_md,idum)
!------------------------------------------------------------------------------
! set all back to input
  set%time_md  = time_in
  set%temp_md  = temp_in
  set%tstep_md = step_in
  set%dump_md2 = dump_in
  set%md_hmass = hmass_in
  if (allocated(opt_in)) then
    set%opt_engine = opt_in
  else
    deallocate(set%opt_engine)
  end if
!------------------------------------------------------------------------------
  write(*,*)
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,'(10x,"|           2D => 3D conversion done!             |")')
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,*)
  set%mode_extrun = mode_input
  mol%info%two_dimensional=.false.
  call gfnff_param_dealloc(calc%topo)

end subroutine struc_convert
end module xtb_gfnff_convert
