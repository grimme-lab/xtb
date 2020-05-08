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
subroutine struc_convert( &
         & env,restart,mol,wfn,calc,egap,et,maxiter,maxcycle,&
         & etot,g,sigma)
  use xtb_mctc_accuracy, only : wp
  use gff_param
  use xtb_disp_dftd3param
  use xtb_type_environment
  use xtb_type_molecule
  use xtb_type_wavefunction
  use xtb_type_calculator
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
  type(TWavefunction),intent(inout)           :: wfn
  type(TCalculator),intent(in)                :: calc
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
  character(len=*),parameter                  :: p_fname_param_gfnff = '.param_gfnff.xtb'
!------------------------------------------------------------------------------
! set up force field
  call struc_convert_header
  mode_input = mode_extrun
  mode_extrun = p_ext_gfnff
  if (.not.allocated(fnv)) fnv=xfind(p_fname_param_gfnff)
  if (.not.allocated(reference_c6)) call copy_c6(reference_c6)
  call open_file(ich,fnv,'r')
  exist = ich .ne. -1
  if (exist) then
     write(*,'(10x,"2D input recognized, conversion to 3D using GFN-FF:")')
     call gfnff_read_param(ich)
     call close_file(ich)
  end if
  call gfnff_setup(verbose,restart,mol,p_ext_gfnff)
!===============================
! Set Block
  ffmode   = -1                  ! just simple rep+harm. on bonds FF
  time_in  = time_md
  time_md  = 5.0_wp              ! short 5 ps MD to move it in 3D
  temp_in  = temp_md
  temp_md  = 298.0_wp            ! md temperature 298 K
  step_in  = tstep_md
  tstep_md = 2.5_wp              ! md time step 2.5 fs
  dump_in  = dump_md2
  dump_md2 = 100.0_wp            ! md dump 100 fs
  hmass_in = md_hmass
  md_hmass = 4.0_wp              ! md hydrogen mass
  nvt_md = .true.                ! md thermostat
!===============================
  if (allocated(opt_logfile)) then
    fnv = opt_logfile
  else
    deallocate(fnv)
  endif
  opt_logfile = 'convert.log'
!------------------------------------------------------------------------------
! force field geometry optimization
  call geometry_optimization &
      &     (env,mol,wfn,calc,   &
      &      egap,etemp,maxiter,maxcycle,etot,g,sigma,p_olev_crude,.false.,.true.,fail)
  if (allocated(fnv)) then
    opt_logfile = fnv
  else
    deallocate(opt_logfile)
  endif
  write(*,*)
!------------------------------------------------------------------------------
! force field md simulation
  idum = 0
  call md                &
      &   (env,mol,wfn,calc, &
      &    egap,etemp,maxiter,etot,g,sigma,0,temp_md,idum)
!------------------------------------------------------------------------------
! set all back to input
  time_md  = time_in
  temp_md  = temp_in
  tstep_md = step_in
  dump_md2 = dump_in
  md_hmass = hmass_in
!------------------------------------------------------------------------------
  write(*,*)
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,'(10x,"|           2D => 3D conversion done!             |")')
  write(*,'(10x," ------------------------------------------------- ")')
  write(*,*)  
  mode_extrun = mode_input
  mol%struc%two_dimensional=.false.
  call gfnff_param_dealloc()

end subroutine struc_convert
