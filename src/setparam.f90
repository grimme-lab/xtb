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

!! ------------------------------------------------------------------------
!  The set documentation has been moved to set_module.f90
!  here you can define all global variables but no more I/O-stuff
!! ------------------------------------------------------------------------
module xtb_setparam
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_kernel, only : gbKernel
   use xtb_solv_input, only : TSolvInput
   use xtb_type_setvar

   implicit none

   private :: wp

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

!! ------------------------------------------------------------------------
!  GFN Hamiltonian: 1=JCTC 2017 version, monopole ES, 2=multipole ES/D4(atm),
   integer,parameter :: p_method_gfn0xtb = 0
   integer,parameter :: p_method_gfn1xtb = 1
   integer,parameter :: p_method_gfn2xtb = 2

!  Geometry input type
   integer,parameter :: p_geo_coord  = 1
   integer,parameter :: p_geo_xmol   = 2
   integer,parameter :: p_geo_sdf    = 3
   integer,parameter :: p_geo_poscar = 4

   ! Initial guess charges
   integer, parameter :: p_guess_sad = 0
   integer, parameter :: p_guess_gasteiger = 1
   integer, parameter :: p_guess_goedecker = 2
   integer, parameter :: p_guess_multieq = 3

   integer, parameter :: p_olev_crude   = -3
   integer, parameter :: p_olev_sloppy  = -2
   integer, parameter :: p_olev_loose   = -1
   integer, parameter :: p_olev_lax     = -4
   integer, parameter :: p_olev_normal  =  0
   integer, parameter :: p_olev_tight   =  1
   integer, parameter :: p_olev_vtight  =  2
   integer, parameter :: p_olev_extreme =  3
   ! geometry optimization backend
   integer, parameter :: p_engine_rf        = 1
   integer, parameter :: p_engine_lbfgs     = 2
   integer, parameter :: p_engine_inertial  = 3
   integer, parameter :: p_engine_pbc_lbfgs = 4

   integer, parameter :: p_modh_read     = -2
   integer, parameter :: p_modh_unit     = -1
   integer, parameter :: p_modh_lindh    =  1
   integer, parameter :: p_modh_lindh_d2 =  2
   integer, parameter :: p_modh_swart    =  3
   integer, parameter :: p_modh_old      =  4
   integer, parameter :: p_modh_gff      =  5

!  Choose the grid for the GBSA
   integer, parameter :: ldgrids(32) = &
      &[  6,  14,  26,  38,  50,  74,  86, 110, 146, 170, &
      & 194, 230, 266, 302, 350, 434, 590, 770, 974,1202, &
      &1454,1730,2030,2354,2702,3074,3470,3890,4334,4802, &
      &5294,5810]
   integer, parameter :: p_angsa_normal    = ldgrids(12)
   integer, parameter :: p_angsa_tight     = ldgrids(19)
   integer, parameter :: p_angsa_verytight = ldgrids(23)
   integer, parameter :: p_angsa_extreme   = ldgrids(32)

   ! interface mode
   integer,parameter :: p_pcem_legacy = 1
   integer,parameter :: p_pcem_orca = 2
   
   type oniom_settings
      
      !> inner region charge
      integer  :: innerchrg
      
      !> cut high order(>1) covalent bonds
      logical :: ignore_topo = .false.
      
      !> derived mode
      logical :: derived = .false.
      
      !> dummy execution to check inner region geo and chrg
      logical :: cut_inner = .false.

      !> explicite charges (inner:outer)
      logical :: fixed_chrgs= .false.
      
      !> mute external output (ORCA,TURBOMOLE)
      logical :: silent = .false.
      
      !> print optimization logs for inner region calculations
      logical :: logs = .false.

      !> if saturate outer region
      logical :: outer = .false.
      
      !> log units
      integer:: ilog1, ilog2
  
   end type oniom_settings

   type qm_external
      
      character(len=:), allocatable :: path
      
      !> absolute path to executable
      character(len=:), allocatable :: executable
      
      !> external input
      character(len=:), allocatable :: input_file
      
      !> alternative for input_file
      character(len=:), allocatable :: input_string
      
      !> molecular structure file
      character(len=:), allocatable :: str_file
      
      !> if input_file exist
      logical :: exist
      
      !> special case of the oniom embedding 
      logical :: oniom=.false.


   end type qm_external

   type TPTBSetup
      !> Do PTB additionally to the normal run in hessian
      logical :: ptb_in_hessian = .false.
      !> Electronic structure method for the energetic hessian part
      character(len=:), allocatable :: hessmethod
      !> temperature for Raman (in K)
      real(wp):: raman_temp = 298.15_wp
   
      !> incident laser wavelength for Raman (in nm)
      real(wp):: raman_lambda = 19435.0_wp
   end type TPTBSetup

   integer, parameter :: p_elprop_beta = 2
   integer, parameter :: p_elprop_alpha = 1
   integer, parameter :: p_elprop_dipole = 0

   integer, parameter :: p_ext_vtb       = -1
   integer, parameter :: p_ext_eht       =  0
   integer, parameter :: p_ext_xtb       =  1
   integer, parameter :: p_ext_driver    =  3
   integer, parameter :: p_ext_turbomole =  4
   integer, parameter :: p_ext_orca      =  5
   integer, parameter :: p_ext_mopac     = 12
   integer, parameter :: p_ext_gfnff     = 13
   integer, parameter :: p_ext_oniom     = 14
   integer, parameter :: p_ext_iff       = 15
   integer, parameter :: p_ext_tblite    = 16
   integer, parameter :: p_ext_ptb       = 17
   integer, parameter :: p_ext_mcgfnff   = 18

   integer, parameter :: p_run_prescc  =   1
   integer, parameter :: p_run_scc    =   2
   integer, parameter :: p_run_grad   =   3
   integer, parameter :: p_run_opt    =   4
   integer, parameter :: p_run_hess   =   5
   integer, parameter :: p_run_ohess  =   7
   integer, parameter :: p_run_bhess  =  71
   integer, parameter :: p_run_md     =   6
   integer, parameter :: p_run_omd    =   8
   integer, parameter :: p_run_path   =  10
   integer, parameter :: p_run_screen =  11
   integer, parameter :: p_run_modef  =  13
   integer, parameter :: p_run_mdopt  =  14
   integer, parameter :: p_run_metaopt=  15
   integer, parameter :: p_run_vip    = 100
   integer, parameter :: p_run_vea    = 101
   integer, parameter :: p_run_vipea  = 102
   integer, parameter :: p_run_vfukui = 103
   integer, parameter :: p_run_vomega = 104

   integer,private :: idum

   type :: TSet
   integer  :: gfn_method = -1
   integer  :: maxscciter = 250
   real(wp) :: acc = 1.0_wp
   logical  :: newdisp = .true.
   logical  :: solve_scc = .true.
   logical  :: periodic = .false.
   logical  :: optcell = .true.

!  Geometry input type
   integer  :: geometry_inputfile = p_geo_coord

!! ------------------------------------------------------------------------
!  electronic SCC temperature
   real(wp) :: etemp = 300.0_wp
!  damping in Broyden SCC procedure (0.05 for critical cases, autoadjusted)
   real(wp) :: broydamp = 0.40_wp
   real(wp) :: dispscale = 1.0_wp

   integer  :: guess_charges = p_guess_goedecker

!  shift molecule to center of mass
   logical  :: do_cma_trafo = .false.

!  static homogenous external electric field in a.u.
   real(wp) :: efield(3) = [0.0_wp, 0.0_wp, 0.0_wp]

! linear dependencies overlap cut-off stuff
   real(wp) :: lidethr = 0.00001_wp   ! cut-off threshold for small overlap eigenvalues

!! ------------------------------------------------------------------------
!  default optimization level
!  crude = -3,     sloppy = -2,      loose = -1,      normal = 0,
!  tight = 1,      verytight = 2,    extreme = 3
   character(len=:), allocatable :: opt_outfile
   character(len=:), allocatable :: opt_logfile
   integer, allocatable :: opt_engine

   !> ANCopt settings
   type(ancopt_setvar) :: optset = ancopt_setvar (&
      optlev = p_olev_normal, &
      micro_opt = 20, &  ! increased during opt.
      maxoptcycle = 0, & ! det. in ancopt routine if not read in
      maxdispl_opt = 1.000_wp, &
      hlow_opt = 0.010_wp, & ! 0.002 is too small
      average_conv = .false.)
   
   type(modhess_setvar) :: mhset = modhess_setvar (&
      model = p_modh_old, &
!  force constants for stretch, bend and torsion
      kr = 0.4000_wp, &
      kf = 0.1300_wp, &
      kt = 0.0075_wp, &
      ko = 0.0000_wp, &
      kd = 0.0000_wp, &
      kq = 0.0000_wp, &
!  cutoff for constructing Hessian
      rcut = 70.0_wp, &
!  dispersion scaling in ANC generation
      s6 = 20.0_wp)

!! ------------------------------------------------------------------------
!  temp. for thermostatistical calc. (could be more than just one T)
   integer  :: nthermo = 1
   real(wp) :: thermotemp(50) = [298.15_wp, spread(0.0_wp, 1, 49)]
!  rotor cut-off (cm-1) in thermo  (was 100 cm-1 previously)
   real(wp) :: thermo_sthr = 50.0_wp
!  threshold (cm-1) for inverting imaginary modes
   real(wp) :: thermo_ithr = -20.0_wp
   ! frequency scaling for therostatistical calculation
   real(wp) :: thermo_fscal = 1.0_wp

!! ------------------------------------------------------------------------
!  MD thermostat/initial siman/GBSA temperature
   real(wp) :: temp_md = 298.15_wp
!  MD run time in ps
   real(wp) :: time_md = 50.0_wp
!  dump (=optimize) structure in siman every mddump fs
   real(wp) :: dump_md = 1000.0_wp ! scoord
!  MD dump step in fs for xyz output file, MUST BE .eq. mdstep for power
!  IR spectra
   real(wp) :: dump_md2 = 50.0_wp ! molden,xyzfile
!  set to 1 if dumps (trj file) should contain velocities
   logical  :: velodump = .false.
!  use thermostat (=1, =0 for NVE)
   logical  :: nvt_md = .true.
!  skip interval in -mdav, -mdopt
   integer  :: skip_md = 500 ! mdopt, mdav
!  MD time step in fs (automatically determined if < 0),
!  could be 4-5 fs with shake =2, md_hmass=4
   real(wp) :: tstep_md = 4.0_wp
!  increase hydrogen mass to this value in amu (at const. tot. mass)
!  allowing large time steps (=0 off)
   integer  :: md_hmass = 4
!  shake on (=0: off which is default) for X-H bonds only (=1),
!  or all bonds (=2) or user defined bonds (=3)
   integer  :: shake_mode = 2
   logical  :: shake_md = .true.
   logical  :: xhonly = .true.
   logical  :: honly = .false.
   logical :: forcewrrestart = .false.

!! ------------------------------------------------------------------------
!  target rmsd value for bhess run in Ångström
   real(wp) :: target_rmsd = 0.1_wp

!! ------------------------------------------------------------------------
!  number of siman annealing blocks
   integer  :: ntemp_siman = 3
!  energy window (kcal) for considering conformers
   real(wp) :: ewin_conf = 20.0_wp
!  highest siman annealing temperature (very system specific)
   real(wp) :: Tend_siman = 1000.0_wp
!  include enantiomers in siman (=1)
   logical  :: enan_siman = .false.
!  SCC accuracy level in MD. Every 10th step the SCC is properly converged
!  at sccconv=1.0. sccmd should be < 5 in critical cases, effects may show
!  up as bad thermostating
   real(wp) :: accu_md = 2.0_wp

!! ------------------------------------------------------------------------
!  SCC accuracy level in Hessian runs
   real(wp) :: accu_hess = 0.3_wp
!  Cartesian displacement increment for numerical Hessian
   real(wp) :: step_hess = 0.005_wp
   ! Scaling factor for the hessian elements
   real(wp) :: scale_hess = 1.0_wp

!  switch on gbsa for solvent if second argument is a valid solvent name
   type(TSolvInput) :: solvInput

!! ------------------------------------------------------------------------
!  of points along normal mode path scan
   integer  :: mode_nscan = 31
!  step lengths for scan (should be around 1 because its adjusted
!  internally to mode mass and FC)
   real(wp) :: mode_step = 1.0_wp
!  update search mode with a fraction of the displacement at every step
!  (0.0 means no update, 0.1-0.2 is a good choice)
   real(wp) :: mode_updat = 0.2_wp
!  use canonical normal modes (=0) or Pipek-Mezey localized ones (=1)
   integer  :: mode_local = 0
!  threshold up to which frequency modes are used for mode based conformer
!  search (def. is 300)
   real(wp) :: mode_vthr = 0.0_wp
!  number of second mode which should be projected out in mode following
!  (normally = 7 ie the TS mode which is fixed then)
   integer  :: mode_prj = 0
!  set by -modef via cmdline
   integer  :: mode_follow = 7

!! ------------------------------------------------------------------------
!  biased path finder based on RMSD criteria
!! ------------------------------------------------------------------------
   type(path_setvar) :: pathset

!! ------------------------------------------------------------------------
!  nano reactor based on RMSD biasing potential
!! ------------------------------------------------------------------------
   type(reactor_setvar) :: reactset

!! ------------------------------------------------------------------------
!  grid spacing for cube file
   real(wp) :: cube_step = 0.4_wp
!  density matrix neglect threshold
   real(wp) :: cube_pthr = 0.05_wp
!  cube boundary offset
   real(wp) :: cube_boff = 3.0_wp
!! ------------------------------------------------------------------------
!  PRINTOUT
!! ------------------------------------------------------------------------
   character(len=:),allocatable  :: property_file
   logical  :: pr_esp = .false.
   character(len=:),allocatable  :: esp_gridfile
   character(len=10) :: lmoinfo_fname='xtblmoinfo'
   logical  :: pr_molden_input = .false.
   logical  :: pr_lmo = .false.
   logical  :: pr_local = .true.
   logical  :: pr_density = .false.
   logical  :: pr_spin_population = .false.
   logical  :: pr_spin_density = .false.
   logical  :: pr_fod = .false.
   logical  :: pr_fod_pop = .false.
   logical  :: pr_wiberg = .true.
   logical  :: pr_wbofrag = .false.
   logical  :: pr_charges = .true.
   logical  :: pr_dipole = .true.
   logical  :: pr_mulliken = .true.
   logical  :: pr_eig = .true.
   logical  :: pr_gbw = .false.
   logical  :: pr_tmmos = .false.
   logical  :: pr_tmbas = .false.
   logical  :: pr_json = .false.
   logical  :: pr_distances = .true.
   logical  :: pr_angles = .false.
   logical  :: pr_torsions = .false.
   logical  :: pr_geosum = .true.
   logical  :: pr_finalstruct = .true.
   logical  :: pr_moments = .true.
   logical  :: pr_modef = .false.
   logical  :: pr_gbsa = .false.
   logical  :: pr_nmtm = .false.
   logical  :: pr_dftbp_hessian_out = .false.

!! ------------------------------------------------------------------------
!  point group symmetrization threshold
   real(wp) :: desy = 0.1_wp
!  point group determination skipped if # atoms > this value
!  (i.e. desymaxat 0 switches it off)
   integer  :: maxatdesy = 200

!  compare molecules in ensemble for removing doubles by RMSD (=0)
!  or rot.const.(=1)
   logical  :: check_rmsd = .true.

!! ------------------------------------------------------------------------
!  (point) charge embedding stuff
!! ------------------------------------------------------------------------
   integer  :: pcem_dummyatom = 7 ! nitrogen
   integer  :: pcem_interface = p_pcem_legacy
   ! pcharge input file
   character(len=:),allocatable :: pcem_file
   character(len=:),allocatable :: pcem_grad
   logical  :: pcem_orca   = .false.
!  controls which interactions included in the Fockian depend on the
!  external point charges
   logical  :: pcem_l_es   = .true.
!  external point charges are included in the anisotropic electrostatics
   logical  :: pcem_l_aes  = .false.
!  external point charges are included as type A non-additive dispersion
   logical  :: pcem_l_disp = .false.
!  external point charges can be anisotropic and will be included
!  in the AES (GFN2 only, works only with pcem_l_aes = .true.)
   logical  :: pcem_l_dipm = .false.
   logical  :: pcem_l_qp   = .false.
!  external point charges can have a coordination number
   logical  :: pcem_l_cn   = .false.
!  external point charges are included in the ATM calculation for
!  type B non-additive dispersion
   logical  :: pcem_l_atm  = .false.

!! ------------------------------------------------------------------------
!  STM images
!! ------------------------------------------------------------------------
   logical  :: pr_stm     = .false.
   real(wp) :: stm_alp    =1.5_wp
   real(wp) :: stm_targ   =1.e-4_wp
   real(wp) :: stm_grid   =0.5_wp
   real(wp) :: stm_pot    =0.0_wp
   real(wp) :: stm_thr    =1.0_wp

!  exchange correction scaling factor for HS case 0.3, for LS case -1.4
   real(wp) :: ex_open ! set to 0.5/-0.5 in .xtbrc, respectively

!! ------------------------------------------------------------------------
!  ONIOM
!! ------------------------------------------------------------------------
   type(oniom_settings) :: oniom_settings

!! ------------------------------------------------------------------------
!  External settings
!! ------------------------------------------------------------------------
   type(qm_external) :: ext_driver
   type(qm_external) :: ext_orca
   type(qm_external) :: ext_turbo
   type(qm_external) :: ext_mopac

!! ------------------------------------------------------------------------
!  information about molecule
!! ------------------------------------------------------------------------
   integer  :: ichrg = 0
   logical  :: clichrg = .false.
   integer  :: nalphabeta = 0

!  cannot be set by .xtbrc/setblock
   integer  :: modflag(50) = 0
   integer  :: tsroot = 0
   integer  :: extcode = 0
   integer  :: extmode = 0
   integer  :: mode_extrun = 1 ! xtb is default
!  integer  :: dummyint ! not used
   integer  :: runtyp = 2 ! SCC by default
   integer  :: elprop = 0 ! dipole by default
   logical  :: rdset = .false.

   ! ENSO (ENergic SOrting something algorithm) compatibility mode
   logical  :: enso_mode = .false.

   logical  :: restart_md = .false.
   logical  :: fit = .false. ! write fit data in scf.f
   logical  :: tsopt = .false.
   logical  :: mdrtrconstr = .false. ! not used
!  initialize at each start the RNG if .false.
   logical  :: samerand = .false.
!  just check the input, don't do calculations
   logical  :: define = .false.
!  printlevel for the main program
   logical  :: silent = .false.
   logical  :: verbose = .false.
   logical  :: veryverbose = .false.
   logical  :: ceasefiles  = .false.

!  character(len=80) :: inputname = ''
   character(len= 4) :: pgroup = 'C1  '
!! ------------------------------------------------------------------------
   !> PTB settings
   type(TPTBSetup) :: ptbsetup
   !> GFN-FF manual setup of nb list via xcontrol
   !  ffnb(42,i) stores the number of neighbors of atom i
   integer, allocatable :: ffnb(:,:)
   end type TSet

   type(TSet) :: set

   type(env_setvar) :: xenv

   character(len=:),allocatable :: molnameline
   character(len=:),allocatable :: commentline

contains


subroutine initrand
   implicit none
   integer :: i,j
   integer,allocatable :: iseed(:)
   integer :: imagic = 41
   if (set%samerand) then
      call random_seed(size=j)
      allocate(iseed(j), source = imagic)
      do i = 1, j
         iseed(i) = iseed(i)+j
      enddo
      call random_seed(put=iseed)
      deallocate(iseed)
   else
      call random_seed()
   endif
end subroutine initrand

function get_namespace(string) result(name)
   use xtb_mctc_global, only : persistentEnv
   implicit none
   character(len=*),intent(in)  :: string
   character(len=:),allocatable :: name
   if (string(1:1).eq.'/') then
      name = string
      return
   endif
   if (allocated(persistentEnv%io%namespace)) then
      if (string(1:1).eq.dot) then
         name = dot//persistentEnv%io%namespace//string
      else
         name = persistentEnv%io%namespace//dot//string
      endif
   else
      name = string
   endif
end function get_namespace

end module xtb_setparam
