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

!! ========================================================================
!  *    WELCOME TO THE   C O N T R O L   &   I N P U T   MODULE IN XTB    *
!! ------------------------------------------------------------------------
!  I use a singleton pattern for our input parser to ensure that you cannot
!  override your input data, therefore it is strictly necessary to follow a
!  strict hierachy on input sources. We decided that command arguments have
!  the highest priority, than follows the various input files and finally a
!  global configuration file is read in. The program structure is set up in
!  exactly that way, so this module has to behave like a singleton. Fortran
!  provides the save keyword, which is used in each set subroutine for this
!  purpose.
!! ------------------------------------------------------------------------
!  The syntax is:
!> # logical instructions require a single statement
!> $fit
!> $samerand
!> # special logicals are chrg and spin
!> $chrg <int>
!> $spin <int>
!> # data instructons are opened by a keyword, which by default has no effect
!> $gfn
!>    version=<int>
!>    d4=<bool>
!> $scc
!>    temp=<real>
!>    broydamp=<real>
!> ...
!> $end # this is optional
!  unknown instructions are ignored and can be used to store additional info
!  that does not affect the outcome of the calculation
!
!  generally we distinguish between, ints, bools and reals. List's of arbitrary
!  combinations those can be archived by the `parse' subroutine from the MCTC lib.
!  Be carefully when reading arrays of unknown size, because you cannot backspace
!  on the input file. Multiline input is also not supported, but we have
!  arbitrary long lines to compensate for this.
!
!  To implement a new feature start by adding the variables to the header of
!  setparam.f90 (here), think about what you are doing, this module is integrated
!  in almost every important unit of xtb, so be careful when chosing variable
!  names (dum, nat, xyz and i usually bad ideas). You should use verbose names.
!
!  Next step is to define a default and add the printout in `write_set', add
!  the necessary information about the feature in the man page of `xcontrol'
!  and a short description in the manual. SERIOUSLY, do this first.
!  For all features that are incompatible with the singleton nature of `setparam'
!  use `constrain_param' and read the introduction there carefully.
!
!  Now you can add the necessary line in the `readflags' loop of `rdcontrol'
!  and define your own set subroutine here. If you add an option to an existing
!  instruction just add the case in the appropriate and set subroutine and
!  use a saved logical to guard it (I called them set<int>).
!
!  If you have any questions on this module, feel free to ask me or another dev
!  SAW: ehlert@thch.uni-bonn.de
!
module set_module
   use iso_fortran_env, only : wp => real64

   use readin, only : mirror_line,get_value

   use setparam

   implicit none

   private :: wp,mirror_line,get_value

   character,private,parameter :: flag = '$'
   character,private,parameter :: colon = ':'
   character,private,parameter :: space = ' '
   character,private,parameter :: equal = '='
   character,private,parameter :: hash = '#'
   character,private,parameter :: dot = '.'
   character,private,parameter :: comma = ','
   character(len=*),private,parameter :: flag_end = flag//'end'

!  Using allocatable arrays of dynamic length strings is only possible
!  with a lot of hacks, so we use good'ol fixed size stack arrays.
!  Let's choose something different from 42 that is not dividable by 10... ;)
!  Happy debugging!
   integer,private,parameter :: p_str_length = 48
   integer,private,parameter :: p_arg_length = 24

   public

contains

subroutine write_set(ictrl)

   implicit none

   integer,intent(in) :: ictrl

!  write the coord file, charge and spin information first
   write(ictrl,'(a,"chrg",1x,i0)') flag,ichrg
   write(ictrl,'(a,"spin",1x,i0)') flag,nalphabeta

!  was the fit-flag set?
   if (fit) write(ictrl,'(a,"fit")') flag
   if (samerand) write(ictrl,'(a,"samerand")') flag

   call write_set_gfn(ictrl)
   call write_set_scc(ictrl)
   call write_set_opt(ictrl)
   call write_set_thermo(ictrl)
   call write_set_md(ictrl)
   call write_set_siman(ictrl)
   call write_set_hess(ictrl)
   call write_set_gbsa(ictrl)
   call write_set_modef(ictrl)
   call write_set_cube(ictrl)
   call write_set_symmetry(ictrl)
   call write_set_embedding(ictrl)
   call write_set_write(ictrl)
   call write_set_external(ictrl)
   call write_set_stm(ictrl)
   call write_set_path(ictrl)
   call write_set_split(ictrl)
   call write_set_wall(ictrl)
   call write_set_fix(ictrl)
   call write_set_constrain(ictrl)
   call write_set_scan(ictrl)
   call write_set_metadyn(ictrl)
   call write_set_reactor(ictrl)

end subroutine write_set

!  write *all* informations contained in the setcommon to xcontrol
subroutine write_set_gfn(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"gfn")') flag
   write(ictrl,'(3x,"method=",i0)') gfn_method
!   select case(gfn_method)
!   case(p_method_gfn0xtb); write(ictrl,'(3x,"method=gfn0-xtb")')
!   case(p_method_gfn1xtb); write(ictrl,'(3x,"method=gfn1-xtb")')
!   case(p_method_gfn2xtb); write(ictrl,'(3x,"method=gfn2-xtb")')
!   end select
   if (gfn_method.gt.2) &
      write(ictrl,'(3x,"d4=",a)') bool2string(newdisp)
   write(ictrl,'(3x,"scc=",a)') bool2string(solve_scc)
   write(ictrl,'(3x,"periodic=",a)') bool2string(periodic)
end subroutine write_set_gfn

subroutine write_set_scc(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"scc")') flag
   write(ictrl,'(3x,"maxiterations=",i0)') maxscciter
   write(ictrl,'(3x,"temp=",g0)') eTemp
   write(ictrl,'(3x,"broydamp=",g0)') broydamp
end subroutine write_set_scc

subroutine write_set_opt(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"opt")') flag
   write(ictrl,'(3x,"engine=")',advance='no')
   select case(opt_engine)
   case default;            write(ictrl,'("unknown")')
   case(p_engine_rf);       write(ictrl,'("rf")')
   case(p_engine_lbfgs);    write(ictrl,'("lbfgs")')
   case(p_engine_inertial); write(ictrl,'("inertial")')
   end select
   if (allocated(opt_outfile)) &
      write(ictrl,'(3x,"output=",a)')  opt_outfile
   if (allocated(opt_logfile)) &
      write(ictrl,'(3x,"logfile=",a)') opt_logfile
   write(ictrl,'(3x,"optlevel=",a)') int2optlevel(optset%optlev)
   write(ictrl,'(3x,"microcycle=",i0)') optset%micro_opt
   write(ictrl,'(3x,"maxcycle=",i0)') optset%maxoptcycle
   write(ictrl,'(3x,"maxdispl=",g0)') optset%maxdispl_opt
   write(ictrl,'(3x,"hlow=",g0)') optset%hlow_opt
   write(ictrl,'(3x,"hessian=")',advance='no')
   select case(mhset%model)
   case default;          write(ictrl,'("lindh-d2")')
   case(p_modh_unit);     write(ictrl,'("unit")')
   case(p_modh_old);      write(ictrl,'("old")')
   case(p_modh_lindh);    write(ictrl,'("lindh")')
   case(p_modh_lindh_d2); write(ictrl,'("lindh-d2")')
   case(p_modh_swart);    write(ictrl,'("swart")')
   end select
   write(ictrl,'(3x,"s6=",g0)') mhset%s6
   write(ictrl,'(3x,"kstretch=",g0)') mhset%kr
   write(ictrl,'(3x,"kbend   =",g0)') mhset%kf
   write(ictrl,'(3x,"ktorsion=",g0)') mhset%kt
   write(ictrl,'(3x,"koutofp =",g0)') mhset%ko
   write(ictrl,'(3x,"kvdw    =",g0)') mhset%kd
   write(ictrl,'(3x,"kes     =",g0)') mhset%kq
   write(ictrl,'(3x,"rcut    =",g0)') sqrt(mhset%rcut)
   write(ictrl,'(3x,"ts=",i0)') bool2int(tsopt)
   write(ictrl,'(3x,"tsroot=",i0)') tsroot
   write(ictrl,'(3x,"exact rf=",g0)') bool2string(optset%exact_rf)
end subroutine write_set_opt

subroutine write_set_thermo(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   integer :: i,idum
   write(ictrl,'(a,"thermo")') flag
   if (thermotemp(idum+1).ne.0.0_wp) idum = size(thermotemp,1)
   write(ictrl,'(3x,"temp=")',advance='no')
   ! now we print all but the last argument to avoid a trailing comma
   do i = 1, nthermo-1
      write(ictrl,'(g0,",")',advance='no') thermotemp(i)
   enddo
   write(ictrl,'(g0)') thermotemp(nthermo)
   write(ictrl,'(3x,"sthr=",g0)') thermo_sthr
end subroutine write_set_thermo

subroutine write_set_md(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"md")') flag
   write(ictrl,'(3x,"temp=",g0)') temp_md
   write(ictrl,'(3x,"time=",g0)') time_md
   write(ictrl,'(3x,"dump=",g0)') dump_md2
   write(ictrl,'(3x,"velo=",i0)') bool2int(velodump)
   write(ictrl,'(3x,"nvt=",i0)') bool2int(nvt_md)
   write(ictrl,'(3x,"skip=",g0)') skip_md
   write(ictrl,'(3x,"step=",g0)') tstep_md
   write(ictrl,'(3x,"hmass=",i0)') md_hmass
   write(ictrl,'(3x,"shake=",i0)') shake_mode
   write(ictrl,'(3x,"sccacc=",g0)') accu_md
end subroutine write_set_md

subroutine write_set_siman(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"siman")') flag
   write(ictrl,'(3x,"dump=",g0)') dump_md
   write(ictrl,'(3x,"n=",i0)') ntemp_siman
   write(ictrl,'(3x,"ewin=",g0)') ewin_conf
   write(ictrl,'(3x,"temp=",g0)') Tend_siman
   write(ictrl,'(3x,"enan=",i0)') bool2int(enan_siman)
   write(ictrl,'(3x,"check=",i0)') bool2int(.not.check_rmsd)
end subroutine write_set_siman

subroutine write_set_hess(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"hess")') flag
   write(ictrl,'(3x,"sccacc=",g0)') accu_hess
   write(ictrl,'(3x,"step=",g0)') step_hess
end subroutine write_set_hess

subroutine write_set_gbsa(ictrl)
   use gbobc, only: lgbsa,ionst,ion_rad
   implicit none
   integer,intent(in) :: ictrl
   if (lgbsa) then
      write(ictrl,'(a,"gbsa")') flag
      if (allocated(solvent)) write(ictrl,'(3x,"solvent=",a)') solvent
      write(ictrl,'(3x,"ion_st=",g0)') ionst
      write(ictrl,'(3x,"ion_rad=",g0)') ion_rad
      write(ictrl,'(3x,"gbsagrid=")',advance='no')
      select case(ngrida)
      case(p_angsa_normal);   write(ictrl,'(a)') "normal"
      case(p_angsa_tight);    write(ictrl,'(a)') "tight"
      case(p_angsa_verytight);write(ictrl,'(a)') "verytight"
      case(p_angsa_extreme);  write(ictrl,'(a)') "extreme"
      case default;           write(ictrl,'(i0)') ngrida
      end select
   endif
end subroutine write_set_gbsa

subroutine write_set_modef(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"modef")') flag
   write(ictrl,'(3x,"n=",i0)') mode_nscan
   write(ictrl,'(3x,"step=",g0)') mode_step
   write(ictrl,'(3x,"updat=",g0)') mode_updat
   write(ictrl,'(3x,"local=",i0)') mode_local
   write(ictrl,'(3x,"vthr=",g0)') mode_vthr
   write(ictrl,'(3x,"prj=",g0)') mode_prj
   write(ictrl,'(3x,"mode=",g0)') mode_follow
end subroutine write_set_modef

subroutine write_set_cube(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"cube")') flag
   write(ictrl,'(3x,"step=",g0)') cube_step
   write(ictrl,'(3x,"pthr=",g0)') cube_pthr
end subroutine write_set_cube

subroutine write_set_symmetry(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"symmetry")') flag
   write(ictrl,'(3x,"desy=",g0)') desy
   write(ictrl,'(3x,"maxat=",i0)') maxatdesy
end subroutine write_set_symmetry

subroutine write_set_embedding(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"embedding")') flag
   write(ictrl,'(3x,"at=",i0)')  pcem_dummyatom
   write(ictrl,'(3x,"es=",a)')   bool2string(pcem_l_es)
   !write(ictrl,'(3x,"aes=",a)')  bool2string(pcem_l_aes)
   !write(ictrl,'(3x,"disp=",a)') bool2string(pcem_l_disp)
   !write(ictrl,'(3x,"dipm=",a)') bool2string(pcem_l_dipm)
   !write(ictrl,'(3x,"qp=",a)')   bool2string(pcem_l_qp)
   !write(ictrl,'(3x,"cn=",a)')   bool2string(pcem_l_cn)
   !write(ictrl,'(3x,"atm=",a)')  bool2string(pcem_l_atm)
end subroutine write_set_embedding

subroutine write_set_write(ictrl)
   use readin, only : bool2int,bool2string
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"write")') flag
   if (allocated(property_file)) &
   write(ictrl,'(3x,"output file=",a)')      property_file
   write(ictrl,'(3x,"esp=",a)')              bool2string(pr_esp)
   if (allocated(esp_gridfile)) write(ictrl,'(3x,"gridfile=",a)')         esp_gridfile
   write(ictrl,'(3x,"mos=",a)')              bool2string(pr_molden_input)
   write(ictrl,'(3x,"gbw=",a)')              bool2string(pr_gbw)
   write(ictrl,'(3x,"tm mos=",a)')           bool2string(pr_tmmos)
   write(ictrl,'(3x,"tm basis=",a)')         bool2string(pr_tmbas)
   write(ictrl,'(3x,"lmo=",a)')              bool2string(pr_lmo)
   write(ictrl,'(3x,"density=",a)')          bool2string(pr_density)
   write(ictrl,'(3x,"spin population=",a)')  bool2string(pr_spin_population)
   write(ictrl,'(3x,"spin density=",a)')     bool2string(pr_spin_density)
   write(ictrl,'(3x,"fod=",a)')              bool2string(pr_fod)
   write(ictrl,'(3x,"fod population=",a)')   bool2string(pr_fod_pop)
   write(ictrl,'(3x,"wiberg=",a)')           bool2string(pr_wiberg)
   write(ictrl,'(3x,"wbo fragments=",a)')    bool2string(pr_wbofrag)
   write(ictrl,'(3x,"dipole=",a)')           bool2string(pr_dipole)
   write(ictrl,'(3x,"charges=",a)')          bool2string(pr_charges)
   write(ictrl,'(3x,"mulliken=",a)')         bool2string(pr_mulliken)
   write(ictrl,'(3x,"orbital energies=",a)') bool2string(pr_eig)
   write(ictrl,'(3x,"inertia=",a)')          bool2string(pr_moments)
   write(ictrl,'(3x,"distances=",a)')        bool2string(pr_distances)
   write(ictrl,'(3x,"angles=",a)')           bool2string(pr_angles)
   write(ictrl,'(3x,"torsions=",a)')         bool2string(pr_torsions)
   write(ictrl,'(3x,"final struct=",a)')     bool2string(pr_finalstruct)
   write(ictrl,'(3x,"geosum=",a)')           bool2string(pr_geosum)
   write(ictrl,'(3x,"stm=",a)')              bool2string(pr_stm)
   write(ictrl,'(3x,"modef=",a)')            bool2string(pr_modef)
   write(ictrl,'(3x,"gbsa=",a)')             bool2string(pr_gbsa)
end subroutine write_set_write

subroutine write_set_external(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   write(ictrl,'(a,"external")') flag
   if (allocated(ext_orca%path))          &
      write(ictrl,'(3x,"orca bin=",a)')        ext_orca%executable
   if (allocated(ext_orca%input_string))  &
      write(ictrl,'(3x,"orca input line=",a)') ext_orca%input_string
   if (allocated(ext_orca%input_file))    &
      write(ictrl,'(3x,"orca input file=",a)') ext_orca%input_file
   if (allocated(ext_mopac%path))         &
      write(ictrl,'(3x,"mopac bin=",a)')       ext_mopac%executable
   if (allocated(ext_mopac%input_string)) &
      write(ictrl,'(3x,"mopac input=",a)')     ext_mopac%input_string
   if (allocated(ext_mopac%input_file))   &
      write(ictrl,'(3x,"mopac file=",a)')      ext_mopac%input_file
   if (allocated(ext_turbo%path))         &
      write(ictrl,'(3x,"turbodir=",a)')        ext_turbo%path
end subroutine write_set_external

subroutine write_set_stm(ictrl)
   implicit none
   integer,intent(in) :: ictrl
   if (.not.pr_stm) return
   write(ictrl,'(a,"stm")') flag
   write(ictrl,'(3x,"broadening=",g0,1x,"#",1x,a)') stm_alp,"in eV"
   write(ictrl,'(3x,"current=",g0,1x,"#")')         stm_targ
   write(ictrl,'(3x,"grid=",g0,1x,"#",1x,a)')       stm_grid,"in au"
   write(ictrl,'(3x,"thr=",g0)')                    stm_thr
   write(ictrl,'(3x,"potential=",g0,1x,"#",1x,a)')  stm_pot,"in V"
end subroutine write_set_stm

subroutine write_set_path(ictrl)
   use tbdef_atomlist
   implicit none
   integer,intent(in) :: ictrl
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer :: i
   write(ictrl,'(a,"path")') flag
   write(ictrl,'(3x,"nrun=",i0)')  pathset%nrun
   write(ictrl,'(3x,"nopt=",i0)')  pathset%nopt
   write(ictrl,'(3x,"anopt=",i0)') pathset%anopt
   write(ictrl,'(3x,"kpush=",g0)') pathset%kpush
   write(ictrl,'(3x,"kpull=",g0)') pathset%kpull
   write(ictrl,'(3x,"alp=",g0)')   pathset%alp
   if (allocated(pathset%fname)) &
      write(ictrl,'(3x,"product=",a)') pathset%fname
   if (pathset%nat > 0) then
      call atl%new(pathset%atoms(:pathset%nat))
      call atl%to_string(string)
      write(ictrl,'(3x,"atoms:",1x,a)') string
   endif

end subroutine write_set_path

subroutine write_set_reactor(ictrl)
   use tbdef_atomlist
   implicit none
   integer,intent(in) :: ictrl
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer :: i
   write(ictrl,'(a,"reactor")') flag
   write(ictrl,'(3x,"max=",i0)')     reactset%nmax
   write(ictrl,'(3x,"density=",g0,1x,"# in kg/L")') &
                                     reactset%dens
   write(ictrl,'(3x,"kpush=",g0)')   reactset%kpush
   write(ictrl,'(3x,"alp=",g0)')     reactset%alp
   if (reactset%nat > 0) then
      call atl%new(reactset%atoms(:reactset%nat))
      call atl%to_string(string)
      write(ictrl,'(3x,"atoms:",1x,a)') string
   endif
end subroutine write_set_reactor


subroutine write_set_constrain(ictrl)
   use scanparam
   use splitparam
   use mctc_econv, only : autoaa
   use mctc_constants, only : pi
   implicit none
   integer, intent(in) :: ictrl
   integer :: i,idum
   if (nconstr.eq.0) return

   write(ictrl,'(a,"constrain")') flag
   write(ictrl,'(3x,"force constant=",g0)') fcconstr
   do i = 1, nconstr
      ! instead of simply saving the kind of constraint while reading,
      ! we used some really obvious integer code system depending on
      ! the number of elements in the atconstr-Array.
      ! Expected to break at some point, add a FIXME and wait for complaints...
      idum = atconstr(1,i)
      if (atconstr(2,i).gt.0) idum = 1
      if (atconstr(3,i).gt.0) idum = 2
      if (atconstr(4,i).gt.0) idum = 3
      select case(idum)
      case default ! this is bad, we don't want this to happen
         call raise('E','This is an internal error, please report this!',1)
      case(-2)
         write(ictrl,'(3x,"center:",1x,g0,1x,i0,1x,"# force constant")') &
            valconstr(i),iatf1
      case(-1)
         write(ictrl,'(3x,"z:",1x,g0,1x,"# in Ångström")') &
            valconstr(i)*autoaa
      case(0)
         write(ictrl,'(3x,"cma:",1x,g0,1x,"# in Ångström")') &
            valconstr(i)*autoaa
      case(1)
         write(ictrl,'(3x,"distance:",1x,2(i0,",",1x),g0,1x,"# in Ångström")') &
            atconstr(1,i),atconstr(2,i),valconstr(i)*autoaa
      case(2)
         write(ictrl,'(3x,"angle:",1x,3(i0,",",1x),g0,1x,"# in Degree")') &
            atconstr(1,i),atconstr(2,i),atconstr(3,i),valconstr(i)*180.0_wp/pi
      case(3)
         write(ictrl,'(3x,"dihedral:",1x,4(i0,",",1x),g0,1x,"# in Degree")') &
            atconstr(1,i),atconstr(2,i),atconstr(3,i),atconstr(4,i), &
            valconstr(i)*180.0_wp/pi
      end select
   enddo

end subroutine write_set_constrain

subroutine write_set_scan(ictrl)
   use scanparam
   use mctc_econv, only : autoaa
   use mctc_constants, only : pi
   implicit none
   integer, intent(in) :: ictrl
   integer :: i,idum
   if (nscan.eq.0) return

   write(ictrl,'(a,"scan")') flag
   select case(scan_mode)
   case default ! this should never happen...
      call raise('E','This is an internal error, please report this!',1)
   case(p_scan_sequential)
      write(ictrl,'(3x,"mode=sequential")')
   case(p_scan_concerted)
      write(ictrl,'(3x,"mode=concerted")')
   end select
   do i = 1, nscan
      ! instead of simply saving the kind of constraint while reading,
      ! we used some really obvious integer code system depending on
      ! the number of elements in the atconstr-Array.
      ! Expected to break at some point, add a FIXME and wait for complaints...
      idum = atconstr(1,i)
      if (atconstr(2,i).gt.0) idum = 1
      if (atconstr(3,i).gt.0) idum = 2
      if (atconstr(4,i).gt.0) idum = 3
      if (idum.le.1) then
         write(ictrl,'(3x,i0,":",1x,2(g0,",",1x),i0,1x,"# in Ångström")') &
            scan_list(i)%iconstr,scan_list(i)%valscan(1)*autoaa, &
            scan_list(i)%valscan(scan_list(i)%nscan)*autoaa,scan_list(i)%nscan
      else
         write(ictrl,'(3x,i0,":",1x,2(g0,",",1x),i0,1x,"# in Degree")') &
            scan_list(i)%iconstr,scan_list(i)%valscan(1)*180.0_wp/pi, &
            scan_list(i)%valscan(scan_list(i)%nscan)*180.0_wp/pi,scan_list(i)%nscan
      endif
   enddo

end subroutine write_set_scan

subroutine write_set_fix(ictrl)
   use tbdef_atomlist
   use fixparam
   implicit none
   integer, intent(in) :: ictrl
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer :: i
   if (fixset%n.eq.0 .and. freezeset%n.eq.0) return

   write(ictrl,'(a,"fix")') flag
   if (fixset%n > 0) then
      call atl%new(fixset%atoms(:fixset%n))
      call atl%to_string(string)
      write(ictrl,'(3x,"atoms:",1x,a)') string
   endif

   write(ictrl,'(3x,"freeze frequency=",g0,1x,"# in atomic units")') freezeset%fc
   if (freezeset%n > 0) then
      call atl%new(freezeset%atoms(:freezeset%n))
      call atl%to_string(string)
      write(ictrl,'(3x,"freeze:",1x,a)') string
   endif

end subroutine write_set_fix

subroutine write_set_split(ictrl)
   use tbdef_atomlist
   use splitparam
   implicit none
   integer, intent(in) :: ictrl
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer :: i
   if ((iatf1.eq.0).and.(iatf2.eq.0)) return

   write(ictrl,'(a,"split")') flag
   do i = 1, maxval(splitlist)
      call atl%new
      call atl%add(splitlist.eq.i)
      call atl%to_string(string)
      write(ictrl,'(3x,"fragment:",1x,i0,",",a)') i, string
   enddo

end subroutine write_set_split

subroutine write_set_wall(ictrl)
   use sphereparam
   implicit none
   integer, intent(in) :: ictrl
   integer :: i,j,nlist

   write(ictrl,'(a,"wall")') flag
   write(ictrl,'(3x,"potential=")',advance='no')
   if (spherepot_type.eq.p_type_polynomial) then
      write(ictrl,'(a)') 'polynomial'
   else
      write(ictrl,'(a)') 'logfermi'
   endif
   write(ictrl,'(3x,"alpha=",i0)') sphere_alpha
   write(ictrl,'(3x,"beta=",g0)') sphere_beta
   write(ictrl,'(3x,"temp=",g0)') sphere_temp
   write(ictrl,'(3x,"autoscale=",g0)') sphere_autoscale
   write(ictrl,'(3x,"axisshift=",g0)') sphere_shift

   if (number_walls.eq.0) return

   do i = 1, number_walls
      if (all(abs(wpot(i)%radius-wpot(i)%radius(1)).lt.1e-6_wp)) then
         write(ictrl,'(3x,"sphere:",1x,g0,",",1x)',advance='no') &
            wpot(i)%radius(1)
      else
         write(ictrl,'(3x,"ellipsoid:",1x,3(g0,",",1x))',advance='no') &
            wpot(i)%radius
      endif
      if (allocated(wpot(i)%list)) then
         nlist = size(wpot(i)%list,1)
         do j = 1, nlist-1
            write(ictrl,'(i0,",",1x)',advance='no') wpot(i)%list(j)
         enddo
         write(ictrl,'(i0,1x,"# radius in Bohr")') wpot(i)%list(nlist)
      else
         write(ictrl,'("all # radius in Bohr")')
      endif
   enddo

end subroutine write_set_wall

subroutine write_set_metadyn(ictrl)
   use tbdef_atomlist
   use fixparam
   implicit none
   integer, intent(in) :: ictrl
   type(tb_atomlist) :: atl
   character(len=:), allocatable :: string
   integer :: i

   if (metaset%maxsave.eq.0) return

   write(ictrl,'(a,"metadyn")') flag
   write(ictrl,'(3x,"save=",i0)')   metaset%maxsave
   write(ictrl,'(3x,"factor=",1x,g0)') metaset%global_factor
   if (any((metaset%factor - metaset%global_factor).ne.0.0_wp)) then
      do i = 1, metaset%maxsave
         if ((metaset%factor(i) - metaset%global_factor).ne.0.0_wp) then
            write(ictrl,'(3x,"modify factor:",1x,i0,",",1x,g0)') &
               i,metaset%factor(i)
         endif
      enddo
   endif
   if (allocated(metaset%fname)) write(ictrl,'(3x,"coord=",a)') metaset%fname
   if (metaset%nat > 0) then
      call atl%new(metaset%atoms(:metaset%nat))
      call atl%to_string(string)
      write(ictrl,'(3x,"atoms:",1x,a)') string
   endif

end subroutine write_set_metadyn

subroutine open_set(ictrl,fname)
   use mctc_timings, only : prtimestring
   implicit none
   integer,intent(out) :: ictrl
   character(len=*),intent(in) :: fname
   character(len=:),allocatable :: buffer
   logical :: exist
   integer :: i

   call open_file(ictrl,fname,'w')

!  write the command line which has called this program
   call get_command(length=i)
   allocate(character(len=i) :: buffer)
   call get_command(buffer)
   write(ictrl,'(a,"cmd",x,a)') flag,buffer

!  write the date when this program was started, use timings for this job
   write(ictrl,'(a,"date",x,a)') flag,prtimestring('S')

end subroutine open_set

subroutine close_set(ictrl)
   integer,intent(in) :: ictrl
   write(ictrl,'("'//flag_end//'")')
   call close_file(ictrl)
end subroutine close_set

subroutine rdcontrol(fname,copy_file)
   use iso_fortran_env, only : output_unit,iostat_end
   use readin, only : find_new_name
   use splitparam,  only : maxfrag
   use scanparam,   only : maxconstr,maxscan
   use sphereparam, only : maxwalls
   implicit none
   character(len=*),intent(in)  :: fname
   character(len=:),allocatable :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   character(len=:),allocatable :: newname
   logical,intent(in),optional  :: copy_file
   integer :: i
   integer :: id
   integer :: ic
   integer :: ie
   integer :: ncount
   integer :: copy
   integer :: err
   logical :: exist
   logical :: do_copy

   if (present(copy_file)) then
      do_copy=copy_file
   else
      do_copy=.false.
   endif

   call open_file(id,fname,'r')
   if (id.eq.-1) then
      call raise('S',"could not find '"//fname//"'",1)
      return
   endif

   if (do_copy) then
      newname = find_new_name(fname)
!     newunit would be nice to have here, but it would always open to a
!     negative number, currently I'm checking for a negative number in
!     in mirror_line to avoid copying to another unit, so I go with unit 41
      call open_file(copy,newname,'w')
   else
      copy = -1 ! deactivate copy in mirror_line
   endif


!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on id (dammit Turbomole format)
   call mirror_line(id,copy,line,err)
   readflags: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         select case(line(2:))
         ! logical
         case('fit'     ); call set_fit;      call mirror_line(id,copy,line,err)
         case('samerand'); call set_samerand; call mirror_line(id,copy,line,err)
         case('cma'     ); call set_cma;      call mirror_line(id,copy,line,err)
         ! data
         case('cube'     ); call rdblock(set_cube,    line,id,copy,err,ncount)
         case('write'    ); call rdblock(set_write,   line,id,copy,err,ncount)
         case('gfn'      ); call rdblock(set_gfn,     line,id,copy,err,ncount)
         case('scc'      ); call rdblock(set_scc,     line,id,copy,err,ncount)
         case('opt'      ); call rdblock(set_opt,     line,id,copy,err,ncount)
         case('hess'     ); call rdblock(set_hess,    line,id,copy,err,ncount)
         case('md'       ); call rdblock(set_md,      line,id,copy,err,ncount)
         case('siman'    ); call rdblock(set_siman,   line,id,copy,err,ncount)
         case('modef'    ); call rdblock(set_modef,   line,id,copy,err,ncount)
         case('gbsa'     ); call rdblock(set_gbsa,    line,id,copy,err,ncount)
         case('embedding'); call rdblock(set_pcem,    line,id,copy,err,ncount)
         case('thermo'   ); call rdblock(set_thermo,  line,id,copy,err,ncount)
         case('external' ); call rdblock(set_external,line,id,copy,err,ncount)
         case('symmetry' ); call rdblock(set_symmetry,line,id,copy,err,ncount)
         case('metadyn'  ); call rdblock(set_metadyn, line,id,copy,err,ncount)
         case('path'     ); call rdblock(set_path,    line,id,copy,err,ncount)
         case('reactor'  ); call rdblock(set_reactor, line,id,copy,err,ncount)
         case('stm'      ); call rdblock(set_stm,     line,id,copy,err,ncount)
         ! data + user data which is read later, but we start counting here
         case('fix'      ); call rdblock(set_fix,     line,id,copy,err,ncount)
         case('wall'     ); call rdblock(set_wall,    line,id,copy,err,ncount)
                            maxwalls = maxwalls + ncount
         case('scan'     ); call rdblock(set_scan,    line,id,copy,err,ncount)
                            maxscan = maxscan + ncount
                            maxconstr = maxconstr + ncount ! better save than sorry
         case('constrain'); call rdblock(set_constr,  line,id,copy,err,ncount)
                            maxconstr = maxconstr + ncount
         case('split'    ); call rdblock(set_split,   line,id,copy,err,ncount)
                            maxfrag = maxfrag + ncount
         ! legacy
         case('set'      ); call rdsetbl(set_legacy,line,id,copy,err)
         case default ! unknown keyword -> ignore, we don't raise them
         !  except for chrg and spin which you actually can set here
         !  read them here because select case might not catch them that easy
            if (index(line(2:),'chrg').eq.1) call set_chrg(line(7:))
            if (index(line(2:),'spin').eq.1) call set_spin(line(7:))
!           get a new line
            call mirror_line(id,copy,line,err)
         end select
      else ! not a keyword -> ignore
         call mirror_line(id,copy,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readflags
!     if (index(line,flag_end).ne.0) exit readflags ! compatibility reasons
   enddo readflags

   if (do_copy) call close_file(copy)
   call close_file(id)
end subroutine rdcontrol

subroutine rdsetbl(handler,line,id,copy,err)
   use iso_fortran_env, only : iostat_end
   implicit none
   interface
      subroutine handler(key,val)
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: val
      end subroutine handler
   end interface
   integer,intent(in) :: id
   integer,intent(in) :: copy
   external :: handler
   integer,intent(out) :: err
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   call raise('S','set-blocks will become obsolete in xtb 6.0 and newer',1)
   do
      call mirror_line(id,copy,line,err)
      if (err.eq.iostat_end) exit
      if (index(line,flag).ne.0) exit

      ! find the first space
      ie = index(line,space)
      if ((line.eq.'').or.(ie.eq.0)) cycle
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(key,val)
   enddo

end subroutine rdsetbl

subroutine rdblock(handler,line,id,copy,err,ncount)
   use iso_fortran_env, only : iostat_end
   implicit none
   interface
      subroutine handler(key,val)
      character(len=*),intent(in) :: key
      character(len=*),intent(in) :: val
      end subroutine handler
   end interface
   integer,intent(in) :: id
   integer,intent(in) :: copy
   external :: handler
   integer,intent(out) :: err
   integer,intent(out) :: ncount
   character(len=:),allocatable,intent(out) :: line
   character(len=:),allocatable :: key
   character(len=:),allocatable :: val
   integer :: ie
   ncount = 0
   do
      call mirror_line(id,copy,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return

      ! find the equal sign
      ie = index(line,equal)
      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
      if (ie.eq.0) cycle 
      key = trim(line(:ie-1))
      val = trim(adjustl(line(ie+1:)))
      call handler(key,val)
   enddo

end subroutine rdblock

subroutine set_exttyp(typ)
   implicit none
   character(len=*),intent(in) :: typ
   logical,save :: set = .true.
   if (.not.set) return
   select case(typ)
   case default ! do nothing
      call raise('S',typ//' is no valid exttyp (internal error)',1)

   case('vtb')
      mode_extrun = p_ext_vtb
   case('eht')
      mode_extrun = p_ext_eht
   case('xtb')
      mode_extrun = p_ext_xtb
   case('qmdff')
      mode_extrun = p_ext_qmdff
   case('orca')
      mode_extrun = p_ext_orca
   case('turbomole')
      mode_extrun = p_ext_turbomole
      extcode = 1
      extmode = 1
   case('mopac')
      mode_extrun = p_ext_mopac

   end select
   set = .false.
end subroutine set_exttyp

subroutine set_geopref(typ)
   implicit none
   character(len=*),intent(in) :: typ
   logical,save :: set = .true.
   if (.not.set) return
   select case(typ)
   case default ! do nothing
      call raise('S',typ//' is no valid geometry format (internal error)',1)

   case('sdf')
      geometry_inputfile = p_geo_sdf

   case('xmol','xyz')
      geometry_inputfile = p_geo_xmol

   case('coord','tm','turbomole')
      geometry_inputfile = p_geo_coord

   case('vasp','poscar')
      geometry_inputfile = p_geo_poscar

   end select
   set = .false.
end subroutine set_geopref
 
subroutine set_runtyp(typ)
   implicit none
   character(len=*),intent(in) :: typ
   logical,save :: set = .true.
   if (.not.set) then
      call raise('S','Runtyp already set and locked, please use a composite runtyp instead.',1)
      return
   endif
   select case(typ)
   case default ! do nothing
      call raise('E',typ//' is no valid runtyp (internal error)',1)

   case('nox')
      runtyp = p_run_nox
   case('stda')
      runtyp = p_run_stda

   case('scc')
      runtyp = p_run_scc

   case('grad')
      runtyp = p_run_grad
   case('opt','optts')
      runtyp = p_run_opt

   case('hess')
      runtyp = p_run_hess
   case('ohess')
      runtyp = p_run_ohess

   case('md')
      runtyp = p_run_md
   case('omd')
      runtyp = p_run_omd

   case('siman')
      runtyp = p_run_siman

   case('path')
      runtyp = p_run_path

   case('screen')
      runtyp = p_run_screen

   case('gmd')
      runtyp = p_run_gmd

   case('modef')
      runtyp = p_run_modef

   case('mdopt')
      runtyp = p_run_mdopt

   case('metaopt')
      runtyp = p_run_metaopt
   case('reactor')
      runtyp = p_run_reactor

   case('vip')
      runtyp = p_run_vip
   case('vea')
      runtyp = p_run_vea
   case('vipea')
      runtyp = p_run_vipea
   case('vomega')
      runtyp = p_run_vomega
   case('vfukui')
      runtyp = p_run_vfukui
   end select
   set = .false.
end subroutine set_runtyp

subroutine set_fit
   implicit none
   fit = .true.
end subroutine set_fit

subroutine set_cma
   implicit none
   do_cma_trafo = .true.
end subroutine set_cma

subroutine set_enso_mode
   implicit none
   enso_mode = .true.
end subroutine set_enso_mode

subroutine set_samerand
   implicit none
   samerand = .true.
end subroutine set_samerand

subroutine set_define
   implicit none
   define = .true.
end subroutine set_define

subroutine set_chrg(val)
   implicit none
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   logical,save :: set1 = .true.
   if (set1) then
      if (get_value(val,idum)) then
         ichrg = idum
      else
         call raise('E','Charge could not be read from your argument',1)
      endif
   endif
   set1 = .false.
end subroutine set_chrg

subroutine set_spin(val)
   implicit none
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   logical,save :: set1 = .true.
   if (set1) then
      if (get_value(val,idum)) then
         nalphabeta = idum
      else
         call raise('E','Spin could not be read from your argument',1)
      endif
   endif
   set1 = .false.
end subroutine set_spin

subroutine set_write(key,val)
   !Purpose: set global parameter for writeouts
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1  = .true.
   logical,save :: set2  = .true.
   logical,save :: set3  = .true.
   logical,save :: set4  = .true.
   logical,save :: set5  = .true.
   logical,save :: set6  = .true.
   logical,save :: set7  = .true.
   logical,save :: set8  = .true.
   logical,save :: set9  = .true.
   logical,save :: set10 = .true.
   logical,save :: set11 = .true.
   logical,save :: set12 = .true.
   logical,save :: set13 = .true.
   logical,save :: set14 = .true.
   logical,save :: set15 = .true.
   logical,save :: set16 = .true.
   logical,save :: set17 = .true.
   logical,save :: set18 = .true.
   logical,save :: set19 = .true.
   logical,save :: set20 = .true.
   logical,save :: set21 = .true.
   logical,save :: set22 = .true.
   logical,save :: set23 = .true.
   logical,save :: set24 = .true.
   logical,save :: set25 = .true.
   logical,save :: set26 = .true.
   logical,save :: set27 = .true.
   logical,save :: set28 = .true.
   logical,save :: set29 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by write",1)
   case('esp')
      if (get_value(val,ldum).and.set1)  pr_esp = ldum
      set1 = .false.
   case('mos')
      if (get_value(val,ldum).and.set2)  pr_molden_input = ldum
      set2 = .false.
   case('lmo')
      if (get_value(val,ldum).and.set3)  pr_lmo = ldum
      set3 = .false.
   case('density')
      if (get_value(val,ldum).and.set4)  pr_density = ldum
      set4 = .false.
   case('spin population')
      if (get_value(val,ldum).and.set5)  pr_spin_population = ldum
      set5 = .false.
   case('spin density')
      if (get_value(val,ldum).and.set6)  pr_spin_density = ldum
      set6 = .false.
   case('fod')
      if (get_value(val,ldum).and.set7) then
         pr_fod = ldum
         pr_fod_pop = ldum
      endif
      set7 = .false.
   case('wiberg')
      if (get_value(val,ldum).and.set8)  pr_wiberg = ldum
      set8 = .false.
   case('dipole')
      if (get_value(val,ldum).and.set9)  pr_dipole = ldum
      set9 = .false.
   case('charges')
      if (get_value(val,ldum).and.set10) pr_charges = ldum
      set10 = .false.
   case('mulliken')
      if (get_value(val,ldum).and.set11) pr_mulliken = ldum
      set11 = .false.
   case('orbital energies')
      if (get_value(val,ldum).and.set12) pr_eig = ldum
      set12 = .false.
   case('gridfile')
       if (set13) esp_gridfile = val
      set13 = .false.
   case('stm')
      if (get_value(val,ldum).and.set14) pr_stm = ldum
      set14 = .false.
   case('gbw')
      if (get_value(val,ldum).and.set15) pr_gbw = ldum
      set15 = .false.
   case('tm mos')
      if (get_value(val,ldum).and.set16) pr_tmmos = ldum
      set16 = .false.
   case('tm basis')
      if (get_value(val,ldum).and.set17) pr_tmbas = ldum
      set17 = .false.
   case('json')
      if (get_value(val,ldum).and.set18) pr_json = ldum
      set18 = .false.
   case('distances')
      if (get_value(val,ldum).and.set19) pr_distances = ldum
      set19 = .false.
   case('angles')
      if (get_value(val,ldum).and.set20) pr_angles = ldum
      set20 = .false.
   case('torsions')
      if (get_value(val,ldum).and.set21) pr_torsions = ldum
      set21 = .false.
   case('final struct')
      if (get_value(val,ldum).and.set22) pr_finalstruct = ldum
      set22 = .false.
   case('geosum')
      if (get_value(val,ldum).and.set23) pr_geosum = ldum
      set23 = .false.
   case('moments','inertia')
      if (get_value(val,ldum).and.set24) pr_moments = ldum
      set24 = .false.
   case('modef')
      if (get_value(val,ldum).and.set25) pr_modef = ldum
      set25 = .false.
   case('wbo fragments')
      if (get_value(val,ldum).and.set26) pr_wbofrag = ldum
      set26 = .false.
   case('output file')
      if (set27) property_file = val
      set27 = .false.
   case('fod population')
      if (get_value(val,ldum).and.set28) pr_fod_pop = ldum
      set28 = .false.
   case('gbsa')
      if (get_value(val,ldum).and.set29) pr_gbsa = ldum
      set29 = .false.
   end select
end subroutine set_write

subroutine set_pcem(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1  = .true.
   logical,save :: set2  = .true.
   logical,save :: set3  = .true.
   logical,save :: set4  = .true.
   logical,save :: set5  = .true.
   logical,save :: set6  = .true.
   logical,save :: set7  = .true.
   logical,save :: set8  = .true.
   logical,save :: set9  = .true.
   logical,save :: set10 = .true.
   logical,save :: set11 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by embedding",1)
   case('at')
      if (get_value(val,idum).and.set1) pcem_dummyatom = idum
      set1 = .false.
   case('es')
      if (get_value(val,ldum).and.set2) pcem_l_es = ldum
      set2 = .false.
   case('aes')
      if (get_value(val,ldum).and.set3) pcem_l_aes = ldum
      set3 = .false.
   case('disp')
      if (get_value(val,ldum).and.set4) pcem_l_disp = ldum
      set4 = .false.
   case('dipm')
      if (get_value(val,ldum).and.set5) pcem_l_dipm = ldum
      set5 = .false.
   case('qp')
      if (get_value(val,ldum).and.set6) pcem_l_qp = ldum
      set6 = .false.
   case('cn')
      if (get_value(val,ldum).and.set7) pcem_l_cn = ldum
      set7 = .false.
   case('atm')
      if (get_value(val,ldum).and.set8) pcem_l_atm = ldum
      set8 = .false.
   case('interface')
      if (set9) then
      select case(val)
      case default
         call raise('S',"Unknown interface value '"//val//"' is ignored",1)
      case('legacy')
         pcem_interface = p_pcem_legacy
      case('orca')
         pcem_interface = p_pcem_orca
         pcem_orca = .true.
      end select
      endif
      set9 = .false.
   case('input')
      if (set10) pcem_file = val
      set10 = .false.
   case('gradient')
      if (set11) pcem_grad = val
      set11 = .false.
   end select
end subroutine set_pcem

subroutine set_gfn(key,val)
   use mctc_strings, only : lowercase
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by gfn",1)
   case('version','method')
      if (key.eq.'version') &
         call raise('S',"Don't use the 'version' key, since it is confusing",1)
      if (get_value(val,idum).and.set1) then
         if ((idum.ge.0).and.(idum.le.2)) then ! actually, this looks stupid...
            gfn_method = idum
         elseif (idum.eq.3) then
            gfn_method = 2
            call raise('S','Please, request GFN2-xTB with method=2!',1)
         else
            call raise('S','We have not yet made a GFN'//val//'-xTB method',1)
         endif
      endif
      set1 = .false.
   case('d4')
      if (get_value(val,ldum).and.set2) newdisp = ldum
      set2 = .false.
   case('scc')
      if (get_value(val,ldum).and.set3) solve_scc = ldum
      set3 = .false.
   case('periodic')
      if (get_value(val,ldum).and.set4) periodic = ldum
      set4 = .false.
   end select
end subroutine set_gfn

subroutine set_scc(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by scc",1)
   case('temp')
      if (get_value(val,ddum).and.set1) eTemp = ddum
      set1 = .false.
   case('broydamp')
      if (get_value(val,ddum).and.set2) broydamp = ddum
      set2 = .false.
   case('guess')
      if (.not.set3) return
      if (val.eq.'gasteiger') then
         guess_charges = p_guess_gasteiger
      else if (val.eq.'goedecker') then
         guess_charges = p_guess_goedecker
      else if (val.eq.'sad') then
         guess_charges = p_guess_sad
      else if (val.eq.'multieq') then
         guess_charges = p_guess_multieq
      endif
      set3 = .false.
   case('maxiterations')
      if (get_value(val,idum).and.set4) then
         if (idum.le.0) then
            call raise('S','negative SCC-Iterations make no sense',1)
         else
            maxscciter = idum
         endif
      endif
      set4 = .false.
   end select
end subroutine set_scc

subroutine set_opt(key,val)
   use mctc_strings
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1  = .true.
   logical,save :: set2  = .true.
   logical,save :: set3  = .true.
   logical,save :: set4  = .true.
   logical,save :: set5  = .true.
   logical,save :: set6  = .true.
   logical,save :: set7  = .true.
   logical,save :: set8  = .true.
   logical,save :: set9  = .true.
   logical,save :: set10 = .true.
   logical,save :: set11 = .true.
   logical,save :: set12 = .true.
   logical,save :: set13 = .true.
   logical,save :: set14 = .true.
   logical,save :: set15 = .true.
   logical,save :: set16 = .true.
   logical,save :: set17 = .true.
   logical,save :: set18 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by opt",1)
   case('optlevel')
      if (set1) then
         optset%optlev = optlevel2int(val)
      endif
      set1 = .false.
   case('microcycle')
      if (get_value(val,idum).and.set2) optset%micro_opt = idum
      set2 = .false.
   case('maxcycle')
      if (get_value(val,idum).and.set3) optset%maxoptcycle = idum
      set3 = .false.
   case('maxdispl')
      if (get_value(val,ddum).and.set4) optset%maxdispl_opt = ddum
      set4 = .false.
   case('hlow')
      if (get_value(val,ddum).and.set5) optset%hlow_opt = ddum
      set5 = .false.
   case('s6','s6opt')
      if (get_value(val,ddum).and.set6) mhset%s6 = ddum
      set6 = .false.
   case('ts')
      if (get_value(val,ldum).and.set7) tsopt = ldum
      set7 = .false.
   case('tsroot')
      if (get_value(val,idum).and.set8) tsroot = idum
      set8 = .false.
   case('kstretch','kr')
      if (get_value(val,ddum).and.set9)  mhset%kr = ddum
      set9 = .false.
   case('kbend',   'kf')
      if (get_value(val,ddum).and.set10) mhset%kf = ddum
      set10 = .false.
   case('ktorsion','kt')
      if (get_value(val,ddum).and.set11) mhset%kt = ddum
      set11 = .false.
   case('koutofp','ko')
      if (get_value(val,ddum).and.set12) mhset%ko = ddum
      set12 = .false.
   case('kvdw','kd')
      if (get_value(val,ddum).and.set13) mhset%kd = ddum
      set13 = .false.
   case('hessian')
      if (set14) then
         select case(val)
         case("lindh");    mhset%model = p_modh_lindh
         case("lindh-d2"); mhset%model = p_modh_lindh_d2
         case("swart");    mhset%model = p_modh_swart
         case("old");      mhset%model = p_modh_old
         case("unit");     mhset%model = p_modh_unit
         end select
      endif
      set14 = .false.
   case('kes','kq')
      if (get_value(val,ddum).and.set15) mhset%kq = ddum
      set15 = .false.
   case('rcut')
      if (get_value(val,ddum).and.set16) mhset%rcut = ddum*ddum
      set16 = .false.
   case('exact rf')
      if (get_value(val,ldum).and.set17) optset%exact_rf = ldum
      set17 = .false.
   case('engine')
      if (set18) then
         select case(lowercase(val))
         case default; call raise('S',"engine '"//val//"' is not implemented",1)
         case('rf','ancopt');      opt_engine = p_engine_rf
         case('lbfgs','l-ancopt'); opt_engine = p_engine_lbfgs
         case('inertial','fire');  opt_engine = p_engine_inertial
         end select
      endif
      set18 = .false.
   case('output')
      if (.not.allocated(opt_outfile)) opt_outfile = val
   case('logfile')
      if (.not.allocated(opt_logfile)) opt_logfile = val
   end select
end subroutine set_opt

function optlevel2int(level) result(ilvl)
   implicit none
   character(len=*),intent(in)  :: level
   integer :: ilvl

   select case(level)
   case default
      ilvl = p_olev_normal
   case('crude','-3')
      ilvl = p_olev_crude
   case('sloppy','-2')
      ilvl = p_olev_sloppy
   case('loose','-1')
      ilvl = p_olev_loose
   case('lax')
      ilvl = p_olev_lax
   case('normal','0')
      ilvl = p_olev_normal
   case('tight','1')
      ilvl = p_olev_tight
   case('verytight','vtight','2')
      ilvl = p_olev_vtight
   case('extreme','3')
      ilvl = p_olev_extreme
   end select

end function optlevel2int

function int2optlevel(ilvl) result(level)
   implicit none
   character(len=:),allocatable :: level
   integer,         intent(in)  :: ilvl

   select case(ilvl)
   case default
      level = 'normal'
   case(-3)
      level = 'crude'
   case(-2)
      level = 'sloppy'
   case(-1)
      level = 'loose'
   case(-4) ! FIXME
      level = 'lax'
   case(0)
      level = 'normal'
   case(1)
      level = 'tight'
   case(2)
      level = 'verytight'
   case(3)
      level = 'extreme'
   end select

end function int2optlevel

subroutine set_thermo(key,val)
   use mctc_strings, only : parse
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: i
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by thermo",1)
   case('temp')
      if (.not.set1) return ! we could read this twice... but we don't do
      thermotemp = 0.0
      nthermo = 0
      call parse(val,comma,argv,narg)
      idum = 0
      do i = 1, narg
         if (get_value(trim(argv(i)),ddum)) then
            if (ddum.le.0.0_wp) then ! doesn't make sense, skip garbage input
               call raise('S',"A temperature of "//trim(argv(i))//" K sounds strange to me",1)
               cycle
            endif
            idum = idum + 1 ! use only readable arguments
            if (idum.gt.size(thermotemp,1)) exit ! don't overflow
            thermotemp(nthermo+idum) = ddum
         endif
      enddo
      nthermo = nthermo+idum
      set1 = .false.
   case('sthr')
      if (get_value(val,ddum).and.set2) thermo_sthr = ddum
      set2 = .false.
   end select
end subroutine set_thermo

subroutine set_md(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   logical,save :: set7 = .true.
   logical,save :: set8 = .true.
   logical,save :: set9 = .true.
   logical,save :: set10= .true.
   logical,save :: set11= .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by md",1)
   case('temp')
      if (get_value(val,ddum).and.set1) temp_md = ddum
      set1 = .false.
   case('time')
      if (get_value(val,ddum).and.set2) time_md = ddum
      set2 = .false.
   case('dump')
      if (get_value(val,ddum).and.set3) dump_md2 = ddum
      set3 = .false.
   case('velo')
!      if (get_value(val,idum).and.set4) then
!         if (idum.eq.1) then
!            velodump = .true.
!         else if (idum.eq.0) then
!            velodump = .false.
!         endif
!      endif
      if (get_value(val,ldum).and.set4) velodump = ldum
      set4 = .false.
   case('nvt')
!      if (get_value(val,idum).and.set5) then
!         if (idum.eq.1) then
!            nvt_md = .true.
!         else if (idum.eq.0) then
!            nvt_md = .false.
!         endif
!      endif
      if (get_value(val,ldum).and.set5) nvt_md = ldum
      set5 = .false.
   case('skip')
      if (get_value(val,idum).and.set6) skip_md = idum
      set6 = .false.
   case('step')
      if (get_value(val,ddum).and.set7) tstep_md = ddum
      set7 = .false.
   case('hmass')
      if (get_value(val,idum).and.set8) md_hmass = idum
      set8 = .false.
   case('shake')
      if (get_value(val,idum).and.set9) then
         if (idum.eq.2) then
            shake_md = .true.
            xhonly = .false.
            shake_mode = 2
         else if (idum.eq.1) then
            shake_md = .true.
            xhonly = .true.
            shake_mode = 1
         else if (idum.eq.0) then
            shake_md = .false.
            xhonly = .false.
            shake_mode = 0
         else if(idum.eq.3) then
            shake_md = .true.
            shake_mode = 3
         endif
      endif
      set9 = .false.
   case('sccacc')
      if (get_value(val,ddum).and.set10) accu_md = ddum
      set10 = .false.
   case('restart')
      if (get_value(val,ldum).and.set11) restart_md = ldum
      set11 = .false.
   end select
end subroutine set_md

subroutine set_siman(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by siman",1)
   case('dump')
      if (get_value(val,ddum).and.set1) dump_md = ddum
      set1 = .false.
   case('n')
      if (get_value(val,idum).and.set2) ntemp_siman = idum
      set2 = .false.
   case('ewin')
      if (get_value(val,ddum).and.set3) ewin_conf = ddum
      set3 = .false.
   case('temp')
      if (get_value(val,ddum).and.set4) Tend_siman = ddum
      set4 = .false.
   case('enan')
!      if (get_value(val,idum).and.set5) then
!         if (idum.eq.1) then
!            enan_siman = .true.
!         else if (idum.eq.0) then
!            enan_siman = .false.
!         endif
!      endif
      if (get_value(val,ldum).and.set5) enan_siman = ldum
      set5 = .false.
   case('check')
      if (get_value(val,idum).and.set6) then
         if (idum.eq.1) then
            check_rmsd = .false.
         else if (idum.eq.0) then
            check_rmsd = .true.
         endif
      endif
      set6 = .false.
   end select
end subroutine set_siman

subroutine set_hess(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by hess",1)
   case('sccacc')
      if (get_value(val,ddum).and.set1) accu_hess = ddum
      set1 = .false.
   case('step')
      if (get_value(val,ddum).and.set2) step_hess = ddum
      set2 = .false.
   end select
end subroutine set_hess

subroutine set_reactor(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by reactor",1)
   case('kpush')
      if (get_value(val,ddum).and.set1) reactset%kpush = ddum
      set1 = .false.
   case('alp')
      if (get_value(val,ddum).and.set2) reactset%alp   = ddum
      set2 = .false.
   case('max')
      if (get_value(val,idum).and.set3) reactset%nmax  = idum
      set3 = .false.
   case('density')
      if (get_value(val,ddum).and.set4) reactset%dens  = ddum
      set4 = .false.
   end select
end subroutine set_reactor


subroutine set_gbsa(key,val)
   use gbobc, only: lsalt,ionst,ion_rad,lgbsa
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
      select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by gbsa",1)
   case('solvent')
      if (set1 .and. val.ne.'none') then
         solvent = val
         lgbsa=.true.
      endif
      set1 = .false.
   case('ion_st')
      if (get_value(val,ddum).and.set2) then
         ionst = ddum
         if (ionst.gt.0.0_wp) lsalt = .true.
      endif
      set2 = .false.
   case('ion_rad')
      if (get_value(val,ddum).and.set3) ion_rad = ddum
      set3 = .false.
   case('gbsagrid')
      if (set4) then
         select case(val)
         case default
            if (get_value(val,idum)) then
               if (any(idum.eq.ldgrids)) then
                  if (idum < p_angsa_normal) &
                     call raise('S',"Small SASA grids can lead to numerical instabilities",1)
                  ngrida = idum
               else
                  call raise('S',"There is no "//val//" Lebedev grid",1)
               endif
            endif
         case('normal');    ngrida = p_angsa_normal
         case('tight');     ngrida = p_angsa_tight
         case('verytight'); ngrida = p_angsa_verytight
         case('extreme');   ngrida = p_angsa_extreme
         endselect
      endif
      set4 = .false.
   end select
end subroutine set_gbsa

subroutine set_modef(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   logical,save :: set7 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by modef",1)
   case('n')
      if (get_value(val,idum).and.set1) mode_nscan = idum
      set1 = .false.
   case('step')
      if (get_value(val,ddum).and.set2) mode_step = ddum
      set2 = .false.
   case('updat')
      if (get_value(val,ddum).and.set3) mode_updat = ddum
      set3 = .false.
   case('local')
      if (get_value(val,idum).and.set4) mode_local = idum
      set4 = .false.
   case('vthr')
      if (get_value(val,ddum).and.set5) mode_vthr = ddum
      set5 = .false.
   case('prj')
      if (get_value(val,idum).and.set6) mode_prj = idum
      set6 = .false.
   case('mode')
      if (get_value(val,idum).and.set7) mode_follow = idum
      set7 = .false.
   end select
end subroutine set_modef

subroutine set_cube(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by cube",1)
   case('step')
      if (get_value(val,ddum).and.set1) cube_step = ddum
      set1 = .false.
   case('pthr')
      if (get_value(val,ddum).and.set2) cube_pthr = ddum
      set2 = .false.
   case('cal')
      call raise('S',"the key 'cal' has been removed",1)
!      if (get_value(val,idum).and.set3) cube_cal = idum
!      set3 = .false.
   end select
end subroutine set_cube

subroutine set_stm(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by stm",1)
   case('broadening')
      if (get_value(val,ddum).and.set1) stm_alp = ddum
      set1 = .false.
   case('current')
      if (get_value(val,ddum).and.set2) stm_targ = ddum
      set2 = .false.
   case('grid')
      if (get_value(val,ddum).and.set3) stm_grid = ddum
      set3 = .false.
   case('thr')
      if (get_value(val,ddum).and.set4) stm_thr = ddum
      set4 = .false.
   case('potential')
      if (get_value(val,ddum).and.set5) stm_pot = ddum
      set5 = .false.
   end select
end subroutine set_stm

subroutine set_symmetry(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by symmetry",1)
   case('desy')
      if (get_value(val,ddum).and.set1) desy = ddum
      set1 = .false.
   case('maxat')
      if (get_value(val,idum).and.set2) maxatdesy = idum
      set2 = .false.
   end select
end subroutine set_symmetry

subroutine set_external(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   logical,save :: set7 = .true.
   logical,save :: set8 = .true.
   select case(key)
   case default
      call raise('S',"the key '"//key//"' is not recognized by external",1)
   case('orca bin')
      if (set1) ext_orca%executable = val
      set1 = .false.
   case('orca input line')
      if (set2) ext_orca%input_string = val
      set2 = .false.
   case('orca input file')
      if (set3) ext_orca%input_file = val
      set3 = .false.
   case('mopac bin')
      if (set4) ext_mopac%executable = val
      set4 = .false.
   case('mopac input')
      if (set5) ext_mopac%input_string = val
      set5 = .false.
   case('mopac file')
      if (set6) ext_mopac%input_file = val
      set6 = .false.
   case('turbodir')
      if (set7) ext_turbo%path = val
      set7 = .false.
   end select
end subroutine set_external

subroutine set_fix(key,val)
   use fixparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by fix",1)
   case('freeze frequency')
      if (get_value(val,ddum).and.set1) freezeset%fc = ddum
      set1 = .false.
   end select

end subroutine set_fix

subroutine set_constr(key,val)
   use scanparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by constrain",1)
   case('force constant')
      if (get_value(val,ddum).and.set1) fcconstr = ddum
      set1 = .false.
   !  this are global constrains, you should only chose one, I think, everything
   !  else doesn't really make sense, but as user it is your own responsibility
   !  so you can set all logicals for whatever reason
   case('all bonds')
      if (get_value(val,ldum).and.set2) lconstr_all_bonds = ldum
      set2 = .false.
   case('all angles')
      if (get_value(val,ldum).and.set3) lconstr_all_angles = ldum
      set3 = .false.
   case('all torsions')
      if (get_value(val,ldum).and.set4) lconstr_all_torsions = ldum
      set4 = .false.
   case('reference')
      if (set5) potset%fname = val
      set5 = .false.
   end select

end subroutine set_constr

subroutine set_metadyn(key,val)
   use fixparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum

   integer  :: i
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv

   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.

   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by metadyn",1)
   case('save')
      if (get_value(val,idum).and.set1) metaset%maxsave = idum
      set1 = .false.
   case('width','alp')
      if (get_value(val,ddum).and.set2) metaset%width = ddum
      set2 = .false.
   case('factor','kpush')
      if (get_value(val,ddum).and.set3) metaset%global_factor = ddum
      set3 = .false.
   case('coord')
      if (set4) metaset%fname = val
      set4 = .false.
   end select

end subroutine set_metadyn

subroutine set_path(key,val)
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum

   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   logical,save :: set7 = .true.

   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by path",1)
   case('nrun')
      if (get_value(val,idum).and.set1) pathset%nrun = idum
      set1 = .false.
   case('nopt')
      if (get_value(val,idum).and.set2) pathset%nopt = idum
      set2 = .false.
   case('anopt')
      if (get_value(val,idum).and.set3) pathset%anopt = idum
      set3 = .false.
   case('kpush')
      if (get_value(val,ddum).and.set4) pathset%kpush = ddum
      set4 = .false.
   case('kpull')
      if (get_value(val,ddum).and.set5) pathset%kpull = ddum
      set5 = .false.
   case('alp')
      if (get_value(val,ddum).and.set6) pathset%alp = ddum
      set6 = .false.
   case('product')
      if (set7) then
         inquire(file=val,exist=ldum)
         if (.not.ldum) then
            call raise('E',"Could not find: '"//val//"' in $path/product",1)
         endif
         pathset%fname = val
      endif
      set7 = .false.
   end select

end subroutine set_path


! this is a dummy routine
subroutine set_scan(key,val)
   use scanparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by scan",1)
   case('mode')
      if (val.eq.'sequential') then
         scan_mode = p_scan_sequential
      else if (val.eq.'concerted') then
         scan_mode = p_scan_concerted
      endif
      set1 = .false.
   end select

end subroutine set_scan

subroutine set_wall(key,val)
   use sphereparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   logical,save :: set1 = .true.
   logical,save :: set2 = .true.
   logical,save :: set3 = .true.
   logical,save :: set4 = .true.
   logical,save :: set5 = .true.
   logical,save :: set6 = .true.
   logical,save :: set7 = .true.
   logical,save :: set8 = .true.
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by fix",1)
   case('potential')
      if (.not.set1) return
      if (val.eq.'polynomial') then
         spherepot_type = p_type_polynomial
      else if (val.eq.'logfermi') then
         spherepot_type = p_type_logfermi
      endif
      set1 = .false.
   case('alpha')
      if (get_value(val,idum).and.set2) sphere_alpha = idum
      set2 = .false.
   case('beta')
      if (get_value(val,ddum).and.set3) sphere_beta = ddum
      set3 = .false.
   case('temp')
      if (get_value(val,ddum).and.set4) sphere_temp = ddum
      set4 = .false.
   case('autoscale')
      if (get_value(val,ddum).and.set5) sphere_autoscale = ddum
      set5 = .false.
   case('axisshift')
      if (get_value(val,ddum).and.set6) sphere_shift = ddum
      set6 = .false.
   case('auto')
      if (.not.set7) return
      select case(val)
      case('distance')
         sphereauto = p_auto_dist
      case('density')
         sphereauto = p_auto_dens
      case('force')
         sphereauto = p_auto_force
      end select
      set7 = .false.
   case('center')
      if (.not.set8) return
      select case(val)
      case('zero')
         spherecent = p_cent_zero
      case('mass')
         spherecent = p_cent_com
      case('geo')
         spherecent = p_cent_cog
      end select
      set8 = .false.
   end select

end subroutine set_wall

! this is a dummy routine
subroutine set_split(key,val)
   use splitparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
end subroutine set_split

subroutine set_legacy(key,val)
   use mctc_strings, only : parse
   use sphereparam
   implicit none
   character(len=*),intent(in) :: key
   character(len=*),intent(in) :: val
   integer  :: err
   integer  :: idum
   real(wp) :: ddum
   logical  :: ldum
   integer  :: narg
   character(len=p_str_length),dimension(p_arg_length) :: argv
   select case(key)
   case default ! do nothing
      call raise('S',"the key '"//key//"' is not recognized by set",1)
   case('runtyp');      call set_runtyp(val)
   case('chrg','charge'); call set_chrg(val)
   case('uhf');           call set_spin(val)
   case('restartmd','mdrestart'); call set_md('restart','1')
   case('samerand');    call set_samerand
   case('hessf');       continue ! used later in read_userdata
!   case('atomlist-');
   case('fragment1');   continue ! used later in read_userdata
   case('fragment2');   continue ! used later in read_userdata
   case('constrxyz');   continue ! used later in read_userdata
!   case('constrainel');
   case('constrainalltors','constralltors'); call set_constr('all torsions','true')
   case('constrainallbo','constralltbo'); call set_constr('all bonds','true')
!   case('constrain');
!   case('scan');
   case('fit');         call set_fit
   case('optlev');      call set_opt('optlevel',val)
   case('gfnver');      call set_gfn('method',val)
   case('ellips');      continue ! used later in read_userdata
                        maxwalls = maxwalls + 1
!   case('ellipsoi');
   case('sphere');      continue ! used later in read_userdata
                        maxwalls = maxwalls + 1
   case('fix');         continue ! used later in read_userdata
   case('fixfc');       call set_fix('force constant',val)
   case('springexp');   call set_fix('spring exponent',val)
   case('constrfc');    call set_constr('force constant',val)
   case('hlowopt');     call set_opt('hlow',val)
   case('s6opt');       call set_opt('s6',val)
   case('microopt');    call set_opt('microcycle',val)
   case('maxopt');      call set_opt('maxcycle',val)
   case('maxdispl');    call set_opt('maxdispl',val)
   case('mdtime');      call set_md('time',val)
   case('mdtemp');      call set_md('temp',val)
   case('etemp');       call set_scc('temp',val)
   case('broydamp');    call set_scc('broydamp',val)
   case('nsiman');      call set_siman('n',val)
   case('mddump');      call set_siman('dump',val)
   case('mddumpxyz');   call set_md('dump',val)
   case('mdskip');      call set_md('skip',val)
   case('shake');       call set_md('shake',val)
   case('md_hmass');    call set_md('hmass',val)
   case('tend');        call set_siman('temp',val)
   case('mdstep');      call set_md('step',val)
   case('velodump');    call set_md('velo',val)
   case('sccmd');       call set_md('sccacc',val)
   case('scchess');     call set_hess('sccacc',val)
   case('stephess');    call set_hess('step',val)
   case('nvt');         call set_md('nvt',val)
   case('enan');        call set_siman('enan',val)
   case('mode_n');      call set_modef('n',val)
   case('mode_step');   call set_modef('step',val)
   case('mode_vthr');   call set_modef('vthr',val)
   case('mode_updat');  call set_modef('updat',val)
   case('mode_local');  call set_modef('local',val)
   case('mode_prj');    call set_modef('prj',val)
   case('path_kpush');  call set_path('kpush',val)
   case('path_kpull');  call set_path('kpull',val)
   case('path_alp');    call set_path('alp',val)
   case('path_nopt');   call set_path('nopt',val)
   case('path_anopt');  call set_path('anopt',val)
   case('path_nrun');   call set_path('nrun',val)
   case('metadyn') ! just for the records, I was against implementing it that way
      call parse(val,space,argv,narg) ! need to parse for spaces...
      if (narg <3) then
         call raise('S','deprecated $set/metadyn keyword broken by user input',1)
         return ! okay, you screwed up, let's get outta here
      endif
      call set_metadyn('factor',trim(argv(1)))
      call set_metadyn('width', trim(argv(2)))
      call set_metadyn('save',  trim(argv(3)))

   case('atomlist+'); continue ! use later in constrain_param
   case('atomlist-'); call raise('S',"$set/atomlist- is not implemented",1)
   case('cube_pthr');   call set_cube('pthr',val)
   case('cube_step');   call set_cube('step',val)
   case('cube_cal'); call raise('S',"The key 'cube_cal' has been removed from $set",1)
   case('gbsa');        call set_gbsa('solvent',val)
   case('ion_st');      call set_gbsa('ion_st',val)
   case('ion_rad');     call set_gbsa('ion_rad',val)
   case('gbsagrid');    call set_gbsa('gbsagrid',val)
   case('ewin_conf');   call set_siman('ewin',val)
   case('desy');        call set_symmetry('desy',val)
   case('thermo');      call set_thermo('temp',val)
   case('thermo_sthr'); call set_thermo('sthr',val)
   case('desymaxat');   call set_symmetry('maxat',val)
!   case('ex_open_HS');
!   case('ex_open_LS');
!   case('orca_mpi');
!   case('orca_exe');
!   case('orca_line');
   case('check_equal'); call set_siman('check',val)
   case('stm_alp');     call set_stm('broadening',val)
   case('stm_targ');    call set_stm('current',val)
   case('stm_grid');    call set_stm('grid',val)
   case('stm_thr');     call set_stm('thr',val)
   case('stm_pot');     call set_stm('potential',val)
   end select

end subroutine set_legacy

end module set_module
