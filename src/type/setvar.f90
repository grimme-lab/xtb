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

module xtb_type_setvar
   use xtb_mctc_accuracy, only : wp
   implicit none

!  NOTE: some set variables here are nothing but placeholders and are
!        currently NOT IN USE! Only variables declared as public will
!        be accessable in setparam.f90 and elsewhere.

   public :: env_setvar
   public :: ancopt_setvar
   public :: modhess_setvar
   public :: fix_setvar
   public :: constr_setvar
   public :: md_setvar
   public :: metadyn_setvar
   public :: reactor_setvar
   public :: path_setvar

   private
   integer,private :: idum

   type :: env_setvar
      character(len=:),allocatable :: home
      character(len=:),allocatable :: path
      character(len=:),allocatable :: namespace
      character(len=:),allocatable :: directory
   end type env_setvar

!! ------------------------------------------------------------------------
!  runtyp
!! ------------------------------------------------------------------------
   type :: run_setvar
   end type run_setvar

!! ------------------------------------------------------------------------
!  GFN Hamiltonian
!! ------------------------------------------------------------------------
   type :: gfn_setvar
      integer  :: gfn_method = -1
      integer  :: maxscciter = 250
      logical  :: newdisp = .true.
      logical  :: solve_scc = .true.
      logical  :: periodic = .false.
   end type gfn_setvar

!! ------------------------------------------------------------------------
!  SCC calculation
!! ------------------------------------------------------------------------
   type :: scc_setvar
!     electronic SCC temperature
      real(wp) :: etemp = 300.0_wp
!     damping in Broyden SCC procedure (0.05 for critical cases, autoadjusted)
      real(wp) :: broydamp = 0.40_wp
   end type scc_setvar

   !>  approximate normal coordinate rational function optimizer 
   type :: ancopt_setvar
      
      !> default optimization level
      !> crude = -3,     sloppy = -2,      loose = -1,      normal = 0,
      !> tight = 1,      verytight = 2,    extreme = 3
      integer  :: optlev = 0
      
      !> number of opt. cycles before new ANC are made by model Hessian
      integer  :: micro_opt = 0

      !> total number of opt. cycles, 0 means automatically determined
      integer  :: maxoptcycle = 0 
      
      !> maximum coordinate displacement in ancopt
      real(wp) :: maxdispl_opt = 0.0_wp
      
      !> lowest force constant in ANC generation (should be > 0.005)
      real(wp) :: hlow_opt = 0.0_wp
      
      logical :: exact_rf = .false.

      !> average energy and gradient before checking for convergence 
      !> to accelerate numerically noisy potential energy surfaces
      logical :: average_conv = .false.

   end type ancopt_setvar

   type modhess_setvar
      integer  :: model = 0
      
      !> cutoff for constructing internal coordinates
      real(wp) :: rcut = 0.0_wp
      
      !> dispersion scaling in ANC generation
      real(wp) :: s6 = 0.0_wp
   
      !> force constants for stretch, bend and torsion
      real(wp) :: kr = 0.0_wp
      real(wp) :: kf = 0.0_wp
      real(wp) :: kt = 0.0_wp
      real(wp) :: ko = 0.0_wp
      real(wp) :: kd = 0.0_wp
      real(wp) :: kq = 0.0_wp
   
   end type modhess_setvar

!! ------------------------------------------------------------------------
!  Constraining potential
!! ------------------------------------------------------------------------
   type :: fix_setvar
      integer  :: n    = 0
      integer  :: nval = 0
      real(wp) :: fc   = 0.0_wp
      real(wp),allocatable :: expo(:)
      integer, allocatable :: atoms(:)
      real(wp),allocatable :: val(:)
   contains
      procedure :: allocate => allocate_fix
      procedure :: deallocate => deallocate_fix
   end type fix_setvar

   type :: constr_setvar
      integer  :: n  = 0
      real(wp) :: fc = 0.0_wp
      character(len=:),allocatable :: fname
      real(wp),allocatable :: xyz(:,:)
      integer, allocatable :: lookup(:)
      integer, allocatable :: typeid(:)
      type(fix_setvar) :: pos
      type(fix_setvar) :: dist
      type(fix_setvar) :: angle
      type(fix_setvar) :: dihedral
   contains
      procedure :: allocate => allocate_constr
      procedure :: deallocate => deallocate_constr
   end type constr_setvar

!! ------------------------------------------------------------------------
!  Hessian
!! ------------------------------------------------------------------------
   type :: hessian_setvar
!     SCC accuracy level in Hessian runs
      real(wp) :: accu_hess = 0.3_wp
!     Cartesian displacement increment for numerical Hessian
      real(wp) :: step_hess = 0.005_wp
   end type hessian_setvar

!! ------------------------------------------------------------------------
!  thermostatistical corrections
!! ------------------------------------------------------------------------
   type :: thermo_setvar  
!     temp. for thermostatistical calc. (could be more than just one T)
      real(wp) :: thermotemp(50) = (/298.15_wp,(0.0_wp,idum=2,50)/)
!     rotor cut-off (cm-1) in thermo  (was 100 cm-1 previously)
      real(wp) :: thermo_sthr = 50.0_wp
   end type thermo_setvar


!! ------------------------------------------------------------------------
!  Molecular dynamics (MD)
!! ------------------------------------------------------------------------
   type :: md_setvar
!     MD thermostat/initial siman/GBSA temperature
      real(wp) :: temp = 298.15_wp
!     MD run time in ps
      real(wp) :: time = 50.0_wp
!     dump (=optimize) structure in siman every mddump fs
      real(wp) :: dump = 1000.0_wp ! scoord
!     MD dump step in fs for xyz output file, MUST BE .eq. mdstep for power
!     IR spectra
      real(wp) :: dump2 = 50.0_wp ! molden,xyzfile
!     set to 1 if dumps (trj file) should contain velocities
      logical  :: velodump = .false.
!     use thermostat (=1, =0 for NVE)
      logical  :: nvt = .true.
!     skip interval in -mdav, -mdopt
      integer  :: skip = 500 ! mdopt, mdav
!     MD time step in fs (automatically determined if < 0),
!     could be 4-5 fs with shake =2, md_hmass=4
      real(wp) :: tstep = 4.0_wp
!     increase hydrogen mass to this value in amu (at const. tot. mass)
!     allowing large time steps (=0 off)
      integer  :: hmass = 4
!     shake on (=0: off which is default) for X-H bonds only (=1),
!     or all bonds (=2) or user defined bonds (=3)
      integer  :: shake_mode = 2
      logical  :: shake_md = .true.
      logical  :: xhonly = .true.
      logical  :: honly = .false.
!     SCC accuracy level in MD. Every 10th step the SCC is properly converged
!     at sccconv=1.0. sccmd should be < 5 in critical cases, effects may show
!     up as bad thermostating
      real(wp) :: accu_md = 2.0_wp
   end type md_setvar

!! ------------------------------------------------------------------------
!  simulated annealing (SIMAN)
!! ------------------------------------------------------------------------
   type :: siman_setvar
!     number of siman annealing blocks
      integer  :: ntemp_siman = 3
!     energy window (kcal) for considering conformers
      real(wp) :: ewin_conf = 20.0_wp
!     highest siman annealing temperature (very system specific)
      real(wp) :: Tend_siman = 1000.0_wp
!     include enantiomers in siman (=1)
      logical  :: enan_siman = .false.
   end type siman_setvar

!! ------------------------------------------------------------------------
!  Mode following
!! ------------------------------------------------------------------------
   type :: modef_setvar
!     of points along normal mode path scan
      integer  :: nscan = 31
!     step lengths for scan (should be around 1 because its adjusted
!     internally to mode mass and FC)
      real(wp) :: step = 1.0_wp
!     update search mode with a fraction of the displacement at every step
!     (0.0 means no update, 0.1-0.2 is a good choice)
      real(wp) :: updat = 0.2_wp
!     use canonical normal modes (=0) or Pipek-Mezey localized ones (=1)
      integer  :: local = 0
!     threshold up to which frequency modes are used for mode based conformer
!     search (def. is 300)
      real(wp) :: vthr = 0.0_wp
!     number of second mode which should be projected out in mode following
!     (normally = 7 ie the TS mode which is fixed then)
      integer  :: prj = 0
!     set by -modef via cmdline
      integer  :: follow = 7
   end type modef_setvar

!! ------------------------------------------------------------------------
!  RMSD based meta dynamic feature
!! ------------------------------------------------------------------------
   type :: metadyn_setvar
      integer  :: maxsave = 0
      integer  :: nstruc = 0
      real(wp) :: global_factor = 0.0_wp
      real(wp),allocatable :: factor(:)
      real(wp) :: global_width = 1.0_wp
      real(wp),allocatable :: width(:)
      real(wp) :: ramp = 0.03_wp
      integer  :: nat = 0
      logical  :: static = .true.
      integer, allocatable :: atoms(:)
      real(wp),allocatable :: xyz(:,:,:)
      character(len=:),allocatable :: fname
   contains
      procedure :: allocate => allocate_metadyn
      procedure :: deallocate => deallocate_metadyn
   end type metadyn_setvar

!! ------------------------------------------------------------------------
!  biased path finder based on RMSD criteria
!! ------------------------------------------------------------------------
   type :: path_setvar
!     number of runs with different settings (<1 means auto setup)
      integer  :: nrun =3  
!     approx. no. of opt. (=points) on path 
      integer  :: nopt =50 
!     max # of "an" opt steps (2-4 is reasonable)
      integer  :: anopt=3  
!     push (educt)   Gaussian RSMD prefactor  
      real(wp) :: kpush= 0.05_wp
!     pull (product) Gaussian RSMD prefactor  
      real(wp) :: kpull=-0.04_wp
!     Gaussian RSMD width
      real(wp) :: alp  = 0.70_wp
!     file name for product structure
      character(len=:),allocatable :: fname
!     default pull strength on path point
      real(wp) :: ppull = 0.05_wp
!     atom list
      integer  :: nat = 0
      integer, allocatable :: atoms(:)
   contains
      procedure :: allocate => allocate_path
      procedure :: deallocate => deallocate_path
   end type path_setvar

!! ------------------------------------------------------------------------
!  nano reactor based on RMSD biasing potential
!! ------------------------------------------------------------------------
   type :: reactor_setvar
!     number of saved images for RMSD potential
      integer  :: nmax =50 
!     push (educt) Gaussian RSMD prefactor  
      real(wp) :: kpush= 0.05_wp
!     Gaussian RSMD width
      real(wp) :: alp  = 0.70_wp
!     density of system
      real(wp) :: dens  = 1.0_wp
!     atom list
      integer  :: nat = 0
      integer, allocatable :: atoms(:)
   contains
      procedure :: allocate => allocate_reactor
      procedure :: deallocate => deallocate_reactor
   end type reactor_setvar

contains

subroutine allocate_metadyn(self,nat,nstruc)
   implicit none
   class(metadyn_setvar) :: self
   integer, intent(in) :: nat
   integer, intent(in) :: nstruc
   self%maxsave = nstruc
   call self%deallocate
   allocate( self%atoms(nat),          source = 0 )
   allocate( self%factor(nstruc),      source = 0.0_wp )
   allocate( self%width(nstruc),       source = 0.0_wp )
   allocate( self%xyz  (3,nat,nstruc), source = 0.0_wp )
end subroutine allocate_metadyn

subroutine deallocate_metadyn(self)
   class(metadyn_setvar) :: self
   if(allocated(self%atoms))  deallocate( self%atoms )
   if(allocated(self%factor)) deallocate( self%factor )
   if(allocated(self%width))  deallocate( self%width )
   if(allocated(self%xyz  ))  deallocate( self%xyz   )
end subroutine deallocate_metadyn

subroutine allocate_path(self,nat)
   implicit none
   class(path_setvar)  :: self
   integer, intent(in) :: nat
   call self%deallocate
   allocate( self%atoms(nat),          source = 0 )
end subroutine allocate_path

subroutine deallocate_path(self)
   class(path_setvar) :: self
   if(allocated(self%atoms))  deallocate( self%atoms )
end subroutine deallocate_path

subroutine allocate_reactor(self,nat)
   implicit none
   class(reactor_setvar) :: self
   integer, intent(in) :: nat
   call self%deallocate
   allocate( self%atoms(nat),          source = 0 )
end subroutine allocate_reactor

subroutine deallocate_reactor(self)
   class(reactor_setvar) :: self
   if(allocated(self%atoms))  deallocate( self%atoms )
end subroutine deallocate_reactor

subroutine allocate_fix(self,nat,nval,fc,expo)
   implicit none
   class(fix_setvar)   :: self
   integer, intent(in) :: nat
   integer, intent(in),optional :: nval
   real(wp),intent(in),optional :: fc
   real(wp),intent(in),optional :: expo
   call self%deallocate
   if (present(nval)) self%nval = nval
   if (present(fc))   self%fc = fc
   if (present(nval) .and. present(expo)) then
      allocate(self%expo(nval))
      self%expo(:) = expo
   end if
   allocate( self%atoms(nat), source = 0 )
   if (present(nval)) allocate( self%val(nval),  source = 0.0_wp )
end subroutine allocate_fix

subroutine deallocate_fix(self)
   class(fix_setvar) :: self
   self%n = 0
   self%fc = 0.0_wp
   if(allocated(self%expo)) deallocate( self%expo )
   if(allocated(self%atoms)) deallocate( self%atoms )
   if(allocated(self%val))   deallocate( self%val )
end subroutine deallocate_fix

subroutine allocate_constr(self,nat,nval,fc,expo)
   implicit none
   class(constr_setvar) :: self
   integer, intent(in)  :: nat
   integer, intent(in)  :: nval
   real(wp),intent(in),optional :: fc
   real(wp),intent(in),optional :: expo
   call self%deallocate
   self%n = nval
   if (present(fc)) self%fc = fc
   allocate( self%lookup(nval), source = 0 )
   allocate( self%typeid(nval), source = 0 )
   call self%pos     %allocate(nat,nat*(nat+1)/2)
   call self%dist    %allocate(nval*2,nval,expo=expo)
   call self%angle   %allocate(nval*3,nval)
   call self%dihedral%allocate(nval*4,nval)
   if (present(fc)) then
      self%pos     %fc = fc
      self%dist    %fc = fc
      self%angle   %fc = fc
      self%dihedral%fc = fc
   endif
end subroutine allocate_constr

subroutine deallocate_constr(self)
   class(constr_setvar) :: self
   call self%pos     %deallocate
   call self%dist    %deallocate
   call self%angle   %deallocate
   call self%dihedral%deallocate
   self%n = 0
   self%fc = 0.0_wp
   if(allocated(self%lookup)) deallocate( self%lookup )
   if(allocated(self%typeid)) deallocate( self%typeid )
end subroutine deallocate_constr

end module xtb_type_setvar
