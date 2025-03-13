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

!> implementation of different relaxation procedures, every relax algorithm
!  comes with its own parameter class (most are general parameters like thresholds)
module xtb_relaxation_engine
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_accuracy, only : wp, sp
   use xtb_mctc_convert, only : fstoau, amutoau
   use xtb_mctc_fileTypes, only : fileType
   use xtb_type_environment, only : TEnvironment
   use xtb_optimizer, only : convergence_log, load_turbomole_log
   use xtb_bfgs
   use xtb_david2
   implicit none
   private :: wp
   !> precision of the rational function step (usually single instead of double)

   !> generic optimization options bundle, contains mostly convergence thresholds
   !  and IO informations
   type relax_options
      !> printlevel
      integer  :: printlevel = 2
      !> single point accuracy
      real(wp) :: acc = 1.0_wp
      !> maximum displacement for MD step
      real(wp) :: max_displacement = 0.2_wp
      !> write optimization log for every nth step
      integer  :: logstep = 1
      !> unit for the optimization log
      integer  :: ilog = -1
      !> file type for optimization log
      integer  :: ftype = 1
      !> threshold for the converence of the energy change of the system
      real(wp) :: e_thr = 5.0e-6_wp
      !> threshold for the converence of the gradient norm of the system
      real(wp) :: g_thr = 1.0e-3_wp
      !> threshold for the converence of the displacement norm of the system
      real(wp) :: d_thr = 1.0e-3_wp
      !> threshold for the converence of the maximum compounent of the gradient
      real(wp) :: g_max = 1.0e-3_wp
      !> threshold for the converence of the maximum compounent of the displacement
      real(wp) :: d_max = 1.0e-3_wp
   end type relax_options

   !> settings for the fast inertial relaxation engine
   type, extends(relax_options) :: fire_options
      !> number of steps for resetting ANC's
      integer :: micro_cycle = 20
      !> time step for the MD propagation
      real(wp) :: time_step = 0.25_wp * fstoau
      !> time step for the MD propagation
      real(wp) :: lat_time_step = 0.25_wp * fstoau
      !> parameter for adjusting velocities
      real(wp) :: astart = 0.1_wp
      !> scaling parameter for astart
      real(wp) :: fa = 0.99_wp
      !> step to start incrementing propagation time length
      integer  :: nmin = 5
      !> increment of time step
      real(wp) :: finc = 1.1_wp
      !> decrement of time step at internal reset
      real(wp) :: fdec = 0.5_wp
      !> maximium number of timestep
      real(wp) :: dtmax = 1.25_wp * fstoau
      !> power threshold for internal reset
      real(wp) :: pcut = 0.0_wp
      !> mass of particles
      real(wp) :: mass = 1.0_wp * amutoau
      !> mass of lattice parameters
      real(wp) :: lat_mass = 1.0_wp * amutoau
      !> use preconditioning with (model) hessian
      logical :: precon = .false.
      !> use BFGS update of hessian
      logical :: update = .false.
   end type fire_options

   !> settings for the low memory BFGS optimizer
   type, extends(relax_options) :: lbfgs_options
      !> lifetime of the approximate normal coordinate system
      integer :: micro_cycle = 20
      !> convergence threshold for the iterative RF solver
      real(sp) :: rfthr = 1.0e-8_sp
      !> number of images to consider in the updating procedure
      integer :: memory = 20
      !> lower eigenvalue cutoff for (model) hessian
      real(wp) :: hlow = 0.02_wp
      !> highest eigenvalue cutoff for (model) hessian
      real(wp) :: hmax = 5.0_wp
      !> increment for microcycles
      real(wp) :: cycle_inc = 1.1_wp
   end type lbfgs_options

   !> use profiler while performing relaxation
   logical, private, parameter :: profile = .true.

   !> format for printing reals in setup section using scientific notation
   character(len=*), private, parameter :: scifmt = &
      '(10x,":",3x,a,e21.7,1x,a,1x,":")'
   !> format for printing reals in setup section
   character(len=*), private, parameter :: dblfmt = &
      '(10x,":",3x,a,f17.7,5x,a,1x,":")'
   !> format for printing ints in setup section
   character(len=*), private, parameter :: intfmt = &
      '(10x,":",3x,a,i17,  5x,a,1x,":")'
   !> format for printing strings in setup section
   character(len=*), private, parameter :: chrfmt = &
      '(10x,":",3x,a,a17,      11x,":")'

   !> format for printlevel==2 cycle header
   character(len=*), private, parameter :: cyclefmt = &
      '(/,72("."),/,30(".")," CYCLE",i5,1x,30("."),/,72("."))'

contains

!> frontend implementation of the fast inertial relaxation engine
subroutine fire &
      &   (env,ilog,mol,chk,calc, &
      &    optlevel,maxstep,energy,egap,gradient,sigma,printlevel,fail)

   use xtb_mctc_convert

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_timer

   use xtb_setparam
   use xtb_fixparam

   use xtb_optimizer
   use xtb_axis
   use xtb_hessian
   use xtb_lsrmsd
   use xtb_pbc_tools

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "relaxation_engine_fire"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   integer, intent(in) :: ilog

   type(TMolecule), intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
   !> optimization level
   integer, intent(in) :: optlevel
   !> maximum number of steps
   integer, intent(in) :: maxstep
   !> energy of the system, contains energy of initial configuration on start
   !  and energy of final configuration at exit
   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: egap
   !> gradient of the system, contains gradient of initial configuration on start
   !  and gradient of final configuration at exit
   real(wp), intent(inout) :: gradient(:,:)
   real(wp), intent(inout) :: sigma(3,3)
   integer, intent(in) :: printlevel
   logical, intent(out) :: fail

   logical :: minpr
   logical :: pr
   logical :: debug
   logical :: converged
   logical :: linear
   integer :: iter
   integer :: nvar
   integer :: nat3
   integer :: thisstep
   integer :: maxcycle

   type(TMolecule) :: molopt

   real(wp) :: time
   real(wp) :: a,b,c
   real(wp) :: U(3,3), x_center(3), y_center(3), rmsdval
   real(wp) :: estart,esave
   real(wp), allocatable :: xyz0(:,:)
   real(wp), allocatable :: xyzopt(:,:)
   real(wp), allocatable :: grmsd(:,:)
   real(wp), allocatable :: velocities(:,:)
   !  packed (model) hessian
   real(wp), allocatable :: hessp(:)
   !  dummy for projection
   real(wp), allocatable :: pmode(:,:)

   !  settings for the fast inertial relaxation engine
   type(fire_options) :: opt

   !  optimize cell parameters
   logical :: optcell
   real(wp) :: lat_gradient(3,3)
   real(wp) :: lat_velocities(3,3)
   real(wp) :: inv_lat(3,3)

   !  timer for profiling (only for printlevel > 1)
   type(tb_timer) :: timer

   ! its an atom, we are done here
   if (mol%n.eq.1) return
   ! setup the timer for a profiling run
   if (profile) call timer%new(8,.false.)
   if (profile) call timer%measure(1,'setup')

   ! settings for optimizer, defaults from ase/optimize/precon/fire.py
   opt = fire_options(logstep = 1, &
      &               dtmax = 1.25_wp * fstoau, time_step = 0.25_wp * fstoau, &
      &               micro_cycle = set%optset%micro_opt, &
      &               max_displacement = set%optset%maxdispl_opt, &
      &               printlevel = printlevel, &
      &               precon = .false., update = .false.)

   ! obtain the thresholds from the optlevel, possible optlevels are currently
   call get_optthr(mol%n,optlevel,opt%e_thr,opt%g_thr,maxcycle,opt%acc)
   ! check for user input regarding the maximum number of cycles
   if (maxstep > 0) maxcycle = maxstep
   ! deactivate microcyles (seems to be expensive with only small gains)

   ! open the logfile, the log is bound to unit 942, so we cannot use newunit
   ! and have to hope that nobody else is currently occupying this identifier
   opt%ilog = ilog
   if (mol%npbc > 0) then
      opt%ftype = fileType%vasp
   else
      opt%ftype = fileType%xyz
   endif
   !call open_file(opt%ilog,'xtbopt.log','w')
   ! write starting structure to log
   if (opt%ilog.ne.-1) then
      call writeMolecule(mol, opt%ilog,format=opt%ftype, energy=energy, &
         & gnorm=norm2(gradient))
   endif

   ! get memory
   nat3 = 3*mol%n
   nvar = nat3
   if(fixset%n.gt.0) then ! exact fixing
      nvar=nat3-3*fixset%n-3
      if(nvar.le.0) nvar=1
   endif

   allocate( velocities(3,mol%n), pmode(nvar,1), hessp(nat3*(nat3+1)/2), &
      &      xyzopt(3,mol%n), source = 0.0_wp )
   ! set defaults
   iter = 0
   converged = .false.
   optcell = mol%npbc > 0
   minpr = opt%printlevel > 0
   pr    = opt%printlevel > 1
   debug = opt%printlevel > 2

   ! initial velocities
   velocities = -opt%time_step * gradient/opt%mass
   if (optcell) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,lat_gradient)
      lat_velocities = -opt%lat_time_step * lat_gradient/opt%lat_mass
   endif

   call axis(mol%n,mol%at,mol%xyz,a,b,c)
   linear = c.lt.1.0e-10_wp

   ! print a nice summary with all settings and thresholds of FIRE
   if(pr)then
      write(env%unit,'(a)') &
         "      ----------------------------------------------------------- ",&
         "     |              Fast Inertial Relaxation Engine              |",&
         "      ----------------------------------------------------------- "
      write(env%unit,'(/,10x,51("."))')
      write(env%unit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(env%unit,'(10x,":",49("."),":")')
      write(env%unit,scifmt) "energy convergence",opt%e_thr,           "Eh   "
      write(env%unit,scifmt) "grad. convergence ",opt%g_thr,           "Eh/a0"
      write(env%unit,chrfmt) "optimization log  ",bool2string(opt%ilog.ne.-1)
      write(env%unit,intfmt) "writing optlog all",opt%logstep,         "steps"
      write(env%unit,intfmt) "Hessian reset     ",opt%micro_cycle,     "steps"
      write(env%unit,intfmt) "time step scaling ",opt%nmin,            "steps"
      write(env%unit,dblfmt) "maximium step len.",opt%max_displacement,"a0   "
      write(env%unit,dblfmt) "initial time step ",opt%time_step*autofs,"fs   "
      write(env%unit,dblfmt) "maximum time step ",opt%dtmax*autofs,    "fs   "
      write(env%unit,dblfmt) "power reset thr   ",opt%pcut,            "au   "
      write(env%unit,dblfmt) "alpha parameter   ",opt%astart,          "     "
      write(env%unit,chrfmt) "BFGS update       ",bool2string(opt%update)
      write(env%unit,chrfmt) "preconditioning   ",bool2string(opt%precon)
      write(env%unit,'(10x,51("."))')
   endif

   if (profile) call timer%measure(1)

   estart = energy
   thisstep = opt%micro_cycle
   molopt = mol

   if (.not.pr.and.minpr) write(env%unit,'(a6,a14,a16,a16,a15,a6)') &
      &          "cycle", "energy", "change", "gnorm", "step", "conv?"
! ======================================================================
   Precon_microiter: do while (.not.converged .and. iter.lt.maxcycle)
! ======================================================================
   if (profile) call timer%measure(7,"model hessian")
   if (opt%precon) then
      if (minpr) write(env%unit,'(" * calculating model hessian...")')
      call modhes(env,calc,set%mhset,molopt%n,molopt%xyz,molopt%at,hessp,pr)
      if(fixset%n.gt.0)then
         ! exact fixing
         call trproj(molopt%n,molopt%n*3,molopt%xyz,hessp,.false.,-1,pmode,1)
      else
         if (.not.linear) &
         ! normal
         call trproj(molopt%n,molopt%n*3,molopt%xyz,hessp,.false.,0,pmode,1) 
      endif

   endif
   if (profile) call timer%measure(7)
   esave = energy
   xyz0 = molopt%xyz

   call inertial_relax &
      &   (env,iter,thisstep,opt,molopt, &
      &    chk,calc,energy,egap,gradient,sigma,hessp,velocities, &
      &    lat_velocities,optcell,converged,fail,timer)

   thisstep = min(ceiling(thisstep*opt%finc),2*opt%micro_cycle)

   call rmsd(molopt%n,xyz0,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)

   if (.not.converged.and.pr) then
      write(env%unit,'(" * RMSD in coord.:",f14.7,1x,"a0")',advance='no') rmsdval
      write(env%unit,'(5x,"energy gain",e16.7,1x,"Eh")') energy-esave
   endif
! ======================================================================
   enddo Precon_microiter
! ======================================================================

   if (converged) then
      if(pr) then
         call rmsd(mol%n,mol%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"CYCLES"
         write(env%unit,'(72("-"))')
         write(env%unit,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
            energy-estart, (energy-estart)*autokcal
         write(env%unit,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
            rmsdval, rmsdval*autoaa
         if (profile) then
            write(env%unit,'(1x,"total power (kW/mol):",F18.7,1x,"(step)",F10.4,1x,"(real)")') &
               & (energy-estart)*autokJ/iter, (energy-estart)*autokJ/timer%get()
         endif
         write(env%unit,'(72("-"))')
      endif
   else
      ! not converging in the given cycles is a FAILURE, we should make this clearer
      ! This is still no ERROR, since we want the geometry written afterwards
      if(pr) then
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN",iter,"CYCLES"
      endif
   endif

   ! save optimized geometry
   mol = molopt

   ! we cannot be sure that the geometry was written in the last optimization
   ! step due to the logstep > 1, so we append the last structure to the optlog
   if (mod(iter-1,opt%logstep).ne.0 .and. opt%ilog.ne.-1) then
      call writeMolecule(mol, opt%ilog, format=opt%ftype, energy=energy, &
         & gnorm=norm2(gradient))
   endif
   !call close_file(opt%ilog)

   ! finally flush the timer to the output unit
   if (minpr.and.profile) call timer%write(env%unit,'FIRE')

end subroutine fire

!> frontend implementation of the low memory/linear scaling (by taste)
!  approximate normal coordinate rational function optimizer (L-ANCopt)
subroutine l_ancopt &
      &   (env,ilog,mol,chk,calc, &
      &    optlevel,maxcycle_in,energy,egap,gradient,sigma,printlevel,fail)

   use xtb_mctc_convert
   use xtb_mctc_lapack, only : lapack_syev

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_xtb_calculator
   use xtb_gfnff_calculator
   use xtb_extern_turbomole, only : TTMCalculator
   use xtb_type_data
   use xtb_type_timer

   use xtb_setparam
   use xtb_fixparam

   use xtb_optimizer
   use xtb_axis
   use xtb_hessian
   use xtb_lsrmsd
   use xtb_detrotra, only : detrotra4

   use xtb_gfnff_fraghess

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "relaxation_engine_l_ancopt"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   integer, intent(in) :: ilog

   type(TMolecule), intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc
   !> optimization level
   integer, intent(in) :: optlevel
   !> maximum number of optimization cycles
   integer, intent(in) :: maxcycle_in
   !> energy of the system, contains energy of initial configuration on start
   !  and energy of final configuration at exit
   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: egap
   !> gradient of the system, contains gradient of initial configuration on start
   !  and gradient of final configuration at exit
   real(wp), intent(inout) :: gradient(:,:)
   real(wp), intent(inout) :: sigma(3,3)
   integer, intent(in) :: printlevel
   logical, intent(out) :: fail

   logical :: minpr
   logical :: pr
   logical :: debug
   logical :: converged
   logical :: linear
   integer :: iter
   integer :: nvar
   integer :: maxcycle
   integer :: thiscycle

   type(TMolecule) :: molopt

   real(wp) :: a,b,c
   real(wp) :: U(3,3), x_center(3), y_center(3), rmsdval
   real(wp) :: estart,esave
   real(wp), allocatable :: xyzopt(:,:)
   real(wp), allocatable :: grmsd(:,:)
   !  packed (model) hessian
   real(wp), allocatable :: hessp(:)
   !  dummy for projection
   real(wp), allocatable :: pmode(:,:)

   !  ANC generation
   integer :: i,j,ij,k
   integer :: nat3
   integer :: lwork
   integer :: info
   integer :: itry
   real(wp) :: thr
   real(wp) :: damp
   real(wp) :: edum
   real(sp), allocatable :: hess(:,:)
   real(sp), allocatable :: aux(:)
   real(sp), allocatable :: eig(:)
   real(wp), allocatable :: trafo(:,:)
   real(wp), allocatable :: hdiag(:)
   real(wp), allocatable :: xyz0(:,:)
   type(convergence_log), allocatable :: avconv

   !  settings for the LBFGS optimizer
   type(lbfgs_options) :: opt

   !  determine wheter hessian is fragmented or not
   logical :: fragmented_hessian

   !  timer for profiling (only for printlevel > 1)
   type(tb_timer) :: timer

   ! its an atom, we are done here
   if (mol%n.eq.1) return
   ! setup the timer for a profiling run
   if (profile) call timer%new(8,.false.)
   if (profile) call timer%measure(1,'setup')

   ! settings for optimizer, defaults from opt.f
   opt = lbfgs_options(logstep = 1, &
      &                printlevel = printlevel, &
      &                micro_cycle = set%optset%micro_opt, &
      &                memory = set%optset%micro_opt, &
      &                max_displacement = set%optset%maxdispl_opt, &
      &                hlow = set%optset%hlow_opt )

   ! obtain the thresholds from the optlevel, possible optlevels are currently
   call get_optthr(mol%n,optlevel,opt%e_thr,opt%g_thr,maxcycle,opt%acc)
   ! check for user input regarding the maximum number of cycles
   if (maxcycle_in > 0) maxcycle = maxcycle_in

   ! Activate averaged convergence criterium
   if (set%optset%average_conv) then
      select type(calc)
      class is(TTMCalculator)
         avconv = load_turbomole_log(maxcycle)
         if (avconv%nlog > 0 .and. pr) then
            write(env%unit, '(a, 1x, i0, 1x, a)') &
               "Convergence averaging initialized with", avconv%nlog, "entries"
         end if
      class default
         avconv = convergence_log(maxcycle)
      end select
   end if

   ! open the logfile, the log is bound to unit 942, so we cannot use newunit
   ! and have to hope that nobody else is currently occupying this identifier
   opt%ilog = ilog
   if (mol%npbc > 0) then
      opt%ftype = fileType%vasp
   else
      opt%ftype = fileType%xyz
   endif
   !call open_file(opt%ilog,'xtbopt.log','w')
   if (opt%ilog.ne.-1) then
      call writeMolecule(mol, opt%ilog, format=opt%ftype, energy=energy, &
         & gnorm=norm2(gradient))
   endif

   ! get memory, allocate single and double precision arrays separately
   nat3 = 3*mol%n
   allocate( pmode(nat3,1), hessp(nat3*(nat3+1)/2), trafo(nat3,nat3), &
      &      xyz0(3,mol%n), xyzopt(3,mol%n), &
      &      source = 0.0_wp )
   allocate( hess(nat3,nat3), eig(nat3), source = 0.0_sp )
   ! set defaults
   iter = 0
   converged = .false.
   minpr = opt%printlevel > 0
   pr    = opt%printlevel > 1
   debug = opt%printlevel > 2
   fragmented_hessian = set%mode_extrun.eq.p_ext_gfnff.and.mol%n.gt.500

   if (set%mode_extrun.eq.p_ext_gfnff) opt%hlow = 0.02
   if (fragmented_hessian) then
      if (mol%n .gt. 2000) opt%hlow = 2.0e-2_wp + (mol%n - 2000.0_wp) * 5.0e-6_wp
      opt%hlow=min(opt%hlow,0.05_wp)
   end if   

   call axis(mol%n,mol%at,mol%xyz,a,b,c)
   linear = c.lt.1.0e-10_wp

   ! different DOF in case of frag hess
   nat3 = 3 * mol%n
   if (fragmented_hessian) then
      nvar = nat3
   else
      nvar = nat3 - 6
      if(linear) nvar = nat3 - 5

      if(fixset%n.gt.0) then ! exact fixing
         nvar=nat3-3*fixset%n-3
         if(nvar.le.0) nvar=1
      endif
   end if

   allocate( hdiag(nvar), source = 0.0_wp )

   ! print a nice summary with all settings and thresholds of ANCopt
   if(pr)then
      write(env%unit,'(a)') &
         "      ----------------------------------------------------------- ",&
         "     |                       L-ANC optimizer                     |",&
         "      ----------------------------------------------------------- "
      write(env%unit,'(10x,a)') &
         "low memory version by Sebastian Ehlert"
      write(env%unit,'(/,10x,51("."))')
      write(env%unit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(env%unit,'(10x,":",49("."),":")')
      write(env%unit,scifmt) "energy convergence",opt%e_thr,           "Eh   "
      write(env%unit,scifmt) "grad. convergence ",opt%g_thr,           "Eh/a0"
      write(env%unit,chrfmt) "optimization log  ",bool2string(opt%ilog.ne.-1)
      write(env%unit,intfmt) "writing optlog all",opt%logstep,         "steps"
      write(env%unit,intfmt) "Hessian reset     ",opt%micro_cycle,     "steps"
      write(env%unit,dblfmt) "maximium step len.",opt%max_displacement,"a0   "
      write(env%unit,intfmt) "LBFGS images      ",opt%memory,          "     "
      write(env%unit,'(10x,51("."))')
   endif

   if (profile) call timer%measure(1)

   estart = energy
   thiscycle = opt%micro_cycle
   molopt = mol

   if (.not.pr.and.minpr) write(env%unit,'(a6,a14,a16,a16,a15,a6)') &
      &          "cycle", "energy", "change", "gnorm", "step", "conv?"
! ======================================================================
   ANC_microiter: do while (.not.converged .and. iter.lt.maxcycle)
! ======================================================================
   if (profile) call timer%measure(2,"model hessian")
   if (minpr) write(env%unit,'(" * calculating model hessian...")')
   call modhes(env,calc,set%mhset,molopt%n,molopt%xyz,molopt%at,hessp,pr)

   ! Project translation, rotation and fixed atoms
    if(fixset%n.gt.0)then
        call trproj(molopt%n,nat3,molopt%xyz,hessp,.false., -1  ,pmode,1)     ! exact fixing
    else
        if (.not.linear) &
        call trproj(molopt%n,nat3,molopt%xyz,hessp,.false.,0,pmode,1)     ! normal
    endif


   if (profile) call timer%measure(2)
   if (profile) call timer%measure(3,"ANC generation")
   ! blowup hessian
   do i = 1, 3*molopt%n
      do j = 1, i
         ij = i*(i-1)/2 + j
         hess(i,j) = hessp(ij)
         hess(j,i) = hessp(ij)
      enddo
   enddo
   ! diagonalize the hessian
   select type(calc)
   type is(TGFFCalculator)
   if (fragmented_hessian) then
      if (calc%topo%nsystem.gt.1) then
         write(env%unit,'(" * fragmented diagonalization...",1x,i0,1x,"fragments")') calc%topo%nsystem
         call frag_hess_diag(mol%n,hess,eig,calc%topo%ispinsyst, &
            & calc%topo%nspinsyst,calc%topo%nsystem)
       else if (calc%topo%nsystem.eq.1) then
          lwork  = 1 + 6*nat3 + 2*nat3**2
          allocate(aux(lwork))
          call lapack_syev ('V','U',nat3,hess,nat3,eig,aux,lwork,info)
          deallocate(aux)
       end if
   else
      lwork  = 1 + 6*nat3 + 2*nat3**2
      allocate(aux(lwork))
      call lapack_syev ('V','U',nat3,hess,nat3,eig,aux,lwork,info)
      deallocate(aux)
   end if
   class default
      lwork  = 1 + 6*nat3 + 2*nat3**2
      allocate(aux(lwork))
      call lapack_syev ('V','U',nat3,hess,nat3,eig,aux,lwork,info)
      deallocate(aux)
   end select

   if (.not. fragmented_hessian) then
      call detrotra4(linear,mol,hess,eig)
   end if

   select type(calc)
   class default
      thr = 1.0e-11_wp
      ! shift all non-zero eigenvalues by
      edum = minval(eig)
      damp = max(opt%hlow - edum,0.0_wp)
      where(abs(eig).gt.thr) eig = eig+damp
   type is(TGFFCalculator)
      thr = 1.0e-11_wp
      ! shift all non-zero eigenvalues by
      edum = minval(eig)
      damp = max(opt%hlow - edum,0.0_wp)
   end select

   if(pr)then
      write(env%unit,*) 'Shifting diagonal of input Hessian by ', damp
      write(env%unit,*) 'Lowest  eigenvalues of input Hessian'
      write(env%unit,'(6F12.6)')(eig(i),i=1,min(18,nat3))
      write(env%unit,*) 'Highest eigenvalues'
      write(env%unit,'(6F12.6)')(eig(i),i=nat3-5,nat3)
      write(env%unit,*)
   endif

! initialize hessian for opt.
   do itry = 1, 4
      trafo = 0.0_wp
      j = 0
      k = 0
      hessp = 0.0_wp
! take largest (positive) first
      do i=nat3,1,-1
         if(abs(eig(i)).gt.thr .and. k.lt.nvar)then
            k=k+1
            trafo(1:nat3,k)=hess(1:nat3,i)
            hessp(k+k*(k-1)/2)=min(max(eig(i),opt%hlow),opt%hmax)
         endif
      enddo

      if (k.ne.nvar) thr=thr*0.1_wp
   enddo
   if(k.ne.nvar) then
      write(env%unit,*)'k=',k,'  nvar=',nvar
      write(env%unit,*) 'ANC generation failed'
      fail=.true.
      return
   endif

   ! sort
   call sort(nat3,nvar,hessp,trafo)
   ! and invert diagonal for LBFGS
   do i = 1, nvar
      j = i*(i-1)/2 + i
      hdiag(i) = 1.0_wp / hessp(j)
   enddo

   ! reset approximate normal coordinate system
   xyz0 = molopt%xyz
   esave = energy
   if (profile) call timer%measure(3)

   call lbfgs_relax &
      &   (env,iter,thiscycle,opt,molopt, &
      &    chk,calc,energy,egap,gradient,sigma,nvar,hdiag,trafo,xyz0, &
      &    converged,fail,timer,avconv)

   thiscycle = min(ceiling(thiscycle*opt%cycle_inc),2*opt%micro_cycle)

   call rmsd(molopt%n,xyz0,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)

   if (.not.converged.and.pr) then
      write(env%unit,'(" * RMSD in coord.:",f14.7,1x,"a0")',advance='no') rmsdval
      write(env%unit,'(5x,"energy gain",e16.7,1x,"Eh")') energy-esave
   endif
! ======================================================================
   enddo ANC_microiter
! ======================================================================

   if (converged) then
      if(minpr) then
         call rmsd(mol%n,mol%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"CYCLES"
         write(env%unit,'(72("-"))')
         write(env%unit,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
            energy-estart, (energy-estart)*autokcal
         write(env%unit,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
            rmsdval, rmsdval*autoaa
         if (profile) then
            write(env%unit,'(1x,"total power (kW/mol):",F18.7,1x,"(step)",F10.4,1x,"(real)")') &
               & (energy-estart)*autokJ/iter, (energy-estart)*autokJ/timer%get()
         endif
         write(env%unit,'(72("-"))')
      endif
   else
      ! not converging in the given cycles is a FAILURE, we should make this clearer
      ! This is still no ERROR, since we want the geometry written afterwards
      if(minpr) then
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN",iter,"CYCLES"
      endif
   endif

   ! save optimized geometry
   mol = molopt

   if (profile) call timer%measure(7,"optimization log")
   ! we cannot be sure that the geometry was written in the last optimization
   ! step due to the logstep > 1, so we append the last structure to the optlog
   if (mod(iter,opt%logstep).ne.1.and.opt%ilog.ne.-1) then
      call writeMolecule(mol, opt%ilog, format=opt%ftype, energy=energy, &
         & gnorm=norm2(gradient))
   end if
   !call close_file(opt%ilog)
   if (profile) call timer%measure(7)

   ! finally flush the timer to the output unit
   if (minpr.and.profile) call timer%write(env%unit,'L-ANCopt')

end subroutine l_ancopt

!> updates displacement using the formula given by Nocedal, generally known
!  as limited memory BFGS algorithm
subroutine lbfgs_step(iter,memory,nvar,gradient,glast,displacement, &
      &               hdiag,s,y,rho)
   implicit none
   !> current iteration step
   integer, intent(in) :: iter
   !> memory limit for the LBFGS update
   integer, intent(in) :: memory
   !> 3*natoms
   integer, intent(in) :: nvar
   !> current gradient
   real(wp), intent(in) :: gradient(:)
   !> gradient on the last point
   real(wp), intent(in) :: glast(:)
   !> on input displacement from the last step, on exit new displacement vector
   real(wp), intent(inout) :: displacement(:)
   !> *inverse* Hessian in diagonal form
   real(wp), intent(in) :: hdiag(:)

   !> LBFGS scratch array of former displacements
   real(wp), intent(inout) :: s(:,:)
   !> LBFGS scratch array of former gradient changes
   real(wp), intent(inout) :: y(:,:)
   !> LBFGS scratch array of dot products between s and y
   real(wp), intent(inout) :: rho(:)

   real(wp), allocatable :: d(:)
   real(wp), allocatable :: q(:)
   real(wp), allocatable :: a(:)
   real(wp) :: b

   real(wp), external :: ddot

   integer :: thisiter
   integer :: lastiter
   integer :: mem
   integer :: i

   allocate( q(nvar), d(nvar), a(memory), source = 0.0_wp )

   thisiter = mod(iter-1, memory) + 1
   s(:,thisiter) = displacement
   y(:,thisiter) = gradient - glast

   b = ddot(nvar,s(:,thisiter),1,y(:,thisiter),1)

   rho(thisiter) = 1.0_wp / max(abs(b), epsilon(1.0_wp))

   q = gradient
   !do mem = 1, min(memory-1, iter-1)
      !i = mod(thisiter - mem - 1, memory) + 1
   do mem = iter, max(1,iter-memory), -1
      i = mod(mem-1,memory)+1
      a(i) = rho(i) * ddot(nvar,s(:,i),1,q,1)
      q = q - a(i) * y(:,i)
   enddo

   d = Hdiag * q

   !lastiter = mod(thisiter, memory) + 1
   !do mem = 0, min(memory-1,iter-1)
      !i = mod(lastiter + mem - 1, memory) + 1
   do mem = max(1,iter-memory), iter
      i = mod(mem-1,memory)+1
      b = rho(i) * ddot(nvar,y(:,i),1,d,1)
      d = d + s(:,i) * (a(i) - b)
   enddo

   displacement = -d
end subroutine lbfgs_step

!> backend implementation of the low memory BFGS algorithm, this implementation
!  is augmented with a coordinate transformation in approximate normal coordinates
subroutine lbfgs_relax &
      &   (env,iter,maxcycle,opt,mol, &
      &    chk,calc,energy,egap,g_xyz,sigma,nvar,hdiag,trafo,xyz0, &
      &    converged,fail,timer,avconv)

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_timer

   use xtb_setparam

   use xtb_optimizer

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "relaxation_engine_lbfgs_relax"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol

   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc

   !> settings for the low memory BFGS
   type(lbfgs_options), intent(in) :: opt
   !> current iteration
   integer, intent(inout) :: iter
   !> maximum number of steps to be performed by this optimizer
   integer, intent(in) :: maxcycle
   !> convergence information
   logical, intent(out) :: converged
   !> internal error
   logical, intent(out) :: fail
   !> energy of the system, contains energy of initial configuration on start
   !  and energy of final configuration at exit
   real(wp), intent(inout) :: energy
   !> HOMO-LUMO gap
   real(wp), intent(inout) :: egap
   !> gradient of the system, contains gradient of initial configuration on start
   !  and gradient of final configuration at exit
   real(wp), intent(inout) :: g_xyz(:,:)
   real(wp), intent(inout) :: sigma(3,3)
   !> dimension of approximate normal coordinate system
   integer, intent(in) :: nvar
   !> diagonal hessian
   real(wp), intent(in) :: hdiag(:)
   !> transformation matrix for cartesian coordinates to ANC's
   real(wp), intent(in) :: trafo(:,:)
   !> start geometry used to generate ANC's
   real(wp), intent(in) :: xyz0(:,:)
   !> timer for profiling the relaxation procedure
   type(tb_timer), intent(inout) :: timer
   !> Average convergence memory
   type(convergence_log), intent(inout), optional :: avconv

   type(scc_results) :: res

   integer :: icycle
   logical :: minpr
   logical :: pr
   logical :: debug
   integer :: memory
   integer :: i,j,ij

   real(wp) :: alpha
   real(wp) :: gnorm
   real(wp) :: elast
   real(wp) :: ediff
   real(wp) :: depred
   real(wp) :: max_displacement
   real(wp) :: step_length
   real(wp), allocatable :: displacement(:)
   real(wp), allocatable :: g_anc(:)
   real(wp), allocatable :: glast(:)

   !  LBFGS variables
   real(wp), allocatable :: lbfgs_s(:,:)
   real(wp), allocatable :: lbfgs_y(:,:)
   real(wp), allocatable :: lbfgs_rho(:)

   !> approximate normal coordinate system
   real(wp), allocatable :: anc(:)
   !  RF variables
   integer :: nvar1,npvar,npvar1
   real(sp) :: dsnrm
   real(sp), external :: sdot
   real(sp), allocatable :: eaug(:)
   real(sp), allocatable :: Uaug(:,:)
   real(sp), allocatable :: Aaug(:)

   !  BLAS interface
   real(wp), external :: ddot

   !  displacement in coordinates is lower than given threshold
   logical  :: dconverged
   logical  :: econverged
   logical  :: gconverged

   minpr = opt%printlevel > 0
   pr    = opt%printlevel > 1
   debug = opt%printlevel > 2

   fail = .false.
   converged  = .false.
   dconverged = .false.
   econverged = .false.
   gconverged = .false.

   memory = min(opt%memory,maxcycle)

   allocate( displacement(nvar), g_anc(nvar), glast(nvar), anc(nvar), &
      &      source = 0.0_wp )

   if (profile) call timer%measure(4,'coordinate transformation')
   ! transform cartesian gradient in ANC's
   call dgemv('t',3*mol%n,nvar,1.0_wp,trafo,3*mol%n,g_xyz,1,0.0_wp,g_anc,1)
   ! get current gradient norm
   gnorm = sqrt(ddot(nvar,g_anc,1,g_anc,1))
   if (profile) call timer%measure(4)
   if (profile) call timer%measure(5,"Rational function")

   ! we need an initial guess for the displacement, to get one, we solve
   !   ⎛ H  g ⎞ ⎛ dx ⎞     ⎛ dx ⎞
   !   ⎝ g  0 ⎠ ⎝  1 ⎠ = λ ⎝  1 ⎠
   ! 'cause it has to a damn good one!
   nvar1  = nvar+1 ! dimension of RF calculation
   npvar = nvar*nvar1/2
   npvar1 = nvar1*(1+nvar1)/2
   allocate(Uaug(nvar1,1),eaug(nvar1),Aaug(npvar1), source = 0.0_sp)
   ! generate the augmented Hessian first
   Aaug = 0.0_sp
   ! at the very beginning, the Hessian is diagonal in the ANC coordinates
   do i = 1, nvar
      Aaug(i*(1+i)/2) = 1.0_wp / hdiag(i)
   enddo
   Aaug(npvar+1:npvar1-1) = g_anc
   ! do steepest decent guess
   Uaug(:,1)=[-real(g_anc,sp),1.0_sp]
   dsnrm = sqrt(sdot(nvar1,Uaug,1,Uaug,1))
   Uaug  = Uaug/dsnrm
   ! refine guess with davidson
   !call sdavid(debug,.false.,nvar1,1,opt%rfthr,Aaug,Uaug,eaug)
   call solver_sdavidson(nvar1,opt%rfthr,Aaug,Uaug,eaug,fail,debug)
   if (abs(Uaug(nvar1,1)) .lt. 1.0e-10_sp) error stop "Internal RF error"
   displacement = Uaug(:nvar,1) / Uaug(nvar1,1)
   ! discard memory allocation in favor of LBFGS
   deallocate( Aaug, Uaug, eaug )
   if (profile) call timer%measure(5)

   ! now get the memory for the LBFGS
   allocate( lbfgs_s(nvar,memory), lbfgs_y(nvar,memory), lbfgs_rho(memory), &
      &      source = 0.0_wp )

   opt_cycle: do icycle = 1, maxcycle
      iter = iter+1

      if (pr) write(env%unit,cyclefmt) iter

      ! rescale the steplength based on the gradient norm
      if (gnorm.lt.0.002_wp) then
         alpha = 1.5_wp
      elseif (gnorm.lt.0.0006_wp) then
         alpha = 2.0_wp
      elseif (gnorm.lt.0.0003_wp) then
         alpha = 3.0_wp
      else
         alpha = 1.0_wp
      endif

      ! check if step is too large
      step_length = alpha * sqrt(ddot(nvar,displacement,1,displacement,1))
      ! the old ANCopt used to cut of the displacement, which will immediately
      ! destroy the LBFGS procedure, therefore we *have* to rescale
      if (step_length > opt%max_displacement) then
         if (debug) write(env%unit,'(" * rescaling step by",f14.7)') &
            &             opt%max_displacement / step_length
         displacement = opt%max_displacement * displacement / step_length
      endif

      ! check for convergence in step length (only printout)
      dconverged = step_length < opt%d_thr

      ! update ANC's
      anc = anc + displacement * alpha
      ! transform from ANC to cartesian coordinates
      if (profile) call timer%measure(4)
      call trfp2xyz(nvar,3*mol%n,anc,xyz0,trafo,mol%xyz)
      if (profile) call timer%measure(4)

      ! save energy and gradient from previous calculation
      glast = g_anc
      elast = energy

      ! get singlepoint energy
      if (profile) call timer%measure(6,"singlepoint calculation")
      call calc%singlepoint(env,mol,chk,opt%printlevel-1,.true.,energy,g_xyz,sigma,egap,res)
      !call gfnff_eg(.false.,mol%n,charge,attyp,xyz,q,.true.,g_xyz,energy)
      if (profile) call timer%measure(6)
      call env%check(fail)
      if (fail) then
         call env%error('SCF not converged, aborting...', source)
         return
      endif

      ! transform cartesian gradient in ANC coordinate system
      if (profile) call timer%measure(4)
      call dgemv('t',3*mol%n,nvar,1.0_wp,trafo,3*mol%n,g_xyz,1,0.0_wp,g_anc,1)
      if (profile) call timer%measure(4)

      if (profile) call timer%measure(7,"optimization log")
      if (opt%ilog.ne.-1) then
         if (mod(icycle,opt%logstep).eq.1.or.opt%logstep.eq.1) then
            call writeMolecule(mol, opt%ilog,format=opt%ftype,energy=res%e_total, &
               & gnorm=res%gnorm)
         end if
      endif
      !call wrlog(mol%n,xyz,attyp,energy,gnorm,.false.)
      if (profile) call timer%measure(7)

      if (present(avconv)) then
         call avconv%set_eg_log(energy, gnorm)
         energy = avconv%get_averaged_energy()
         gnorm = avconv%get_averaged_gradient()
         if (pr) then
            write(env%unit,'("av. E:",1x,f14.7,1x,"->",1x,f14.7)') &
               avconv%elog(avconv%nlog), energy
            write(env%unit,'("av. G:",1x,f14.7,1x,"->",1x,f14.7)') &
               avconv%glog(avconv%nlog), gnorm
         end if
      end if

      ! check for convergence of the energy change
      ediff = energy - elast
      econverged = ediff < 0.0_wp .and. abs(ediff) < opt%e_thr

      ! check for convergence of the gradient norm
      gnorm = sqrt(ddot(nvar,g_anc,1,g_anc,1))
      gconverged = gnorm < opt%g_thr

      if (pr) then
         write(env%unit,'(" * total energy  :",f14.7,1x,"Eh")',advance='no') energy
         write(env%unit,'(5x,"change   ",e18.7,1x,"Eh")')                    ediff
         write(env%unit,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no') gnorm
         write(env%unit,'(2x,"converged?",3x,3(1x,a,"=",l1))') &
            "E",econverged,"G",gconverged,"D",dconverged
         write(env%unit,'(3x,"step length   :",f14.7,1x,"a0")') step_length
      else if (minpr) then
         write(env%unit,'(i6,f14.7,1x,"(",e14.7,")",2(1x,f14.7),1x,"(",3l1,")")') &
            & iter, energy, ediff, gnorm, step_length, &
            & econverged, gconverged, dconverged
      endif

      ! check for convergence of the minimization
      converged = econverged.and.gconverged
      if (converged) then
         exit opt_cycle
      endif

      ! update the inverse Hessian and obtain next step
      if (profile) call timer%measure(8,"hessian update")
      call lbfgs_step(icycle,memory,nvar,g_anc,glast,displacement, &
         &            hdiag,lbfgs_s,lbfgs_y,lbfgs_rho)
      if (profile) call timer%measure(8)

   enddo opt_cycle

end subroutine lbfgs_relax

!> backend implementation of the fast inertial relaxation engine
subroutine inertial_relax &
      &   (env,iter,maxstep,opt,mol, &
      &    chk,calc,energy,egap,gradient,sigma,hessp,velocities, &
      &    lat_velocities,optcell,converged,fail,timer)

   use xtb_mctc_convert

   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_timer

   use xtb_setparam

   use xtb_pbc_tools

   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "relaxation_engine_inertial_relax"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol
   type(TRestart),intent(inout) :: chk
   class(TCalculator), intent(inout) :: calc

   !> settings for the fast inertial relaxation engine
   type(fire_options), intent(in) :: opt
   !> current iteration
   integer, intent(inout) :: iter
   !> maximum number of steps to be performed by this relaxation engine
   integer, intent(in) :: maxstep
   !> convergence information
   logical, intent(out) :: converged
   !> internal error
   logical, intent(out) :: fail
   !> energy of the system, contains energy of initial configuration on start
   !  and energy of final configuration at exit
   real(wp), intent(inout) :: energy
   real(wp), intent(inout) :: egap
   !> gradient of the system, contains gradient of initial configuration on start
   !  and gradient of final configuration at exit
   real(wp), intent(inout) :: gradient(:,:)
   real(wp), intent(inout) :: sigma(3,3)
   !> packed (model) hessian for preconditioner
   real(wp), intent(inout) :: hessp(:)
   !> velocities of the atoms
   real(wp), intent(inout) :: velocities(:,:)
   logical,  intent(in)    :: optcell
   real(wp), intent(inout) :: lat_velocities(3,3)
   !> timer for profiling the relaxation procedure
   type(tb_timer), intent(inout) :: timer

   type(scc_results) :: res

   integer :: istep
   logical :: minpr
   logical :: pr
   logical :: debug
   integer :: jstep

   real(wp) :: power
   real(wp) :: elast
   real(wp) :: ediff
   real(wp) :: gnorm
   real(wp) :: speed
   real(wp) :: step_length
   real(wp) :: time_step
   real(wp) :: max_tstep
   real(wp) :: apar
   real(wp) :: gpinvg
   real(wp), allocatable :: displacement(:,:)
   real(wp), allocatable :: pinvg(:,:)
   real(wp), allocatable :: scratch(:,:)

   !  cell optimization
   integer  :: lstep
   real(wp) :: lat_step_length
   real(wp) :: lat_displacement(3,3)
   real(wp) :: lat_gradient(3,3)
   real(wp) :: inv_lat(3,3)
   real(wp) :: lat_gnorm
   real(wp) :: lat_power
   real(wp) :: lat_speed
   real(wp) :: lat_time_step
   real(wp) :: lat_apar

   !  BLAS interface
   real(wp), external :: ddot

   !  displacement in coordinates is lower than given threshold
   logical  :: dconverged
   logical  :: econverged
   logical  :: gconverged
   logical  :: logthis

   allocate( displacement(3,mol%n), pinvg(3,mol%n), scratch(3,mol%n), &
      &      source = 0.0_wp )

   minpr = opt%printlevel > 0
   pr    = opt%printlevel > 1
   debug = opt%printlevel > 2

   fail = .false.
   converged  = .false.
   dconverged = .false.
   econverged = .false.
   gconverged = .false.

   time_step = opt%time_step
   max_tstep = opt%dtmax
   apar = opt%astart
   if (optcell) then
      lat_apar = opt%astart
      lat_time_step = opt%lat_time_step
   endif

   ! get current gradient norm
   gnorm = sqrt(ddot(3*mol%n,gradient,1,gradient,1))
   ! get current speed
   speed = sqrt(ddot(3*mol%n,velocities,1,velocities,1))
   ! get power based on v·F = -v·g
   power = -ddot(3*mol%n,velocities,1,gradient,1)/opt%mass
   if (optcell) then
      inv_lat = mat_inv_3x3(mol%lattice)
      call sigma_to_latgrad(sigma,inv_lat,lat_gradient)
      lat_gnorm = sqrt(ddot(9,lat_gradient,1,lat_gradient,1))
      lat_speed = sqrt(ddot(9,lat_velocities,1,lat_velocities,1))
      lat_power = -ddot(9,lat_velocities,1,lat_gradient,1)/opt%lat_mass
   endif

   jstep = 0
   lstep = 0
   md_step: do istep = 1, maxstep
      iter = iter+1
      jstep = jstep+1
      lstep = lstep+1

      if (pr) write(env%unit,cyclefmt) iter

      if (profile) call timer%measure(6,"preconditioner")
      if (opt%precon) &
         call get_preconditioner_metric(mol%n*3,hessp,gradient,pinvg,fail)
      if (profile) call timer%measure(6)

      if (profile) call timer%measure(2,"velocities")
      ! update velocities based on current forces
      if (power > opt%pcut) then
         if (.not.fail .and. opt%precon) then
            gpinvg = -ddot(3*mol%n,gradient,1,pinvg,1)
            call dspmv('u',3*mol%n,1.0_wp,hessp,velocities,1,0.0_wp,scratch,1)
            velocities = (1.0_wp - apar) * velocities &
               &       + apar * pinvg/gpinvg * ddot(3*mol%n,scratch,1,velocities,1)
         else
            velocities = (1.0_wp - apar) * velocities &
               &       - apar * gradient/gnorm * speed
         endif
         if (jstep > opt%nmin) then
            time_step = min(time_step * opt%finc, max_tstep)
            apar = apar*opt%fa
         endif
      else
         write(env%unit,'(" * reset velocities")')
         ! internal reset of current velocities (TODO: return to driver?)
         velocities = 0.0_wp
         apar = opt%astart
         time_step = time_step*opt%fdec
         jstep = 0
      endif
      if (opt%precon) then
         velocities = velocities + time_step * pinvg/opt%mass
      else
         velocities = velocities - time_step * gradient/opt%mass
      endif

      if (optcell) then
         if (lat_power > opt%pcut) then
            lat_velocities = (1.0_wp - lat_apar) * lat_velocities &
               &       - lat_apar * lat_gradient/lat_gnorm * lat_speed
            if (lstep > opt%nmin) then
               lat_time_step = min(lat_time_step * opt%finc, max_tstep)
               lat_apar = lat_apar*opt%fa
            endif
         else
            write(env%unit,'(" * reset lattice velocities")')
            ! internal reset of current velocities (TODO: return to driver?)
            lat_velocities = 0.0_wp
            lat_apar = opt%astart
            lat_time_step = lat_time_step*opt%fdec
            lstep = 0
         endif

         lat_velocities = lat_velocities - lat_time_step * lat_gradient/opt%lat_mass
      endif

      ! get new speed
      speed = sqrt(ddot(3*mol%n,velocities,1,velocities,1))
      if (optcell) then
         lat_speed = sqrt(ddot(9,lat_velocities,1,lat_velocities,1))
      endif

      if (profile) call timer%measure(2)
      if (profile) call timer%measure(3,"propagation")

      ! get propagation from velocities
      displacement = time_step * velocities
      if (optcell) then
         lat_displacement = lat_time_step * lat_velocities
      endif

      ! renormalize if step taken is to large
      step_length = sqrt(ddot(3*mol%n,displacement,1,displacement,1))
      if (step_length > opt%max_displacement) then
         write(env%unit,'(" * rescaling displacement step")')
         displacement = opt%max_displacement * displacement / step_length
         !velocities = opt%max_displacement * velocities / step_length
      endif
      !print*,velocities
      !print*,gradient
      !print*,displacement
      if (optcell) then
         lat_step_length = sqrt(ddot(9,lat_displacement,1,lat_displacement,1))
         if (lat_step_length > opt%max_displacement) then
            write(env%unit,'(" * rescaling lattice displacement")')
            lat_displacement = opt%max_displacement*lat_displacement/lat_step_length
         endif
      endif

      dconverged = step_length < opt%d_thr .and. &
         & (.not.optcell .or. lat_step_length < opt%d_thr)

      ! update positions
      if (optcell) then
         mol%lattice = mol%lattice + lat_displacement
         call mol%wrap_back
      endif
      mol%xyz = mol%xyz + displacement
      ! save last energy
      elast = energy
      ! save last gradient
      scratch = gradient

      if (profile) call timer%measure(3)
      if (profile) call timer%measure(4,"singlepoint calculation")
      ! get singlepoint energy
      call calc%singlepoint(env,mol,chk,opt%printlevel-1,.true.,energy,gradient,sigma,egap,res)
      if (profile) call timer%measure(4)
      call env%check(fail)
      if (fail) then
         call env%error('SCF not converged, aborting...', source)
         return
      endif
      if (profile) call timer%measure(5,"log and printout")

      ! check if we write to the optimization log this step
      logthis = opt%ilog.ne.-1 .and. mod(istep-1,opt%logstep).eq.0
      if (logthis) then
         call writeMolecule(mol, opt%ilog,format=opt%ftype,energy=res%e_total, &
            & gnorm=res%gnorm)
      endif

      ! check for convergence of the energy change
      ediff = energy - elast
      econverged = ediff < 0.0_wp .and. abs(ediff) < opt%e_thr

      ! get new power
      power = -ddot(3*mol%n,velocities,1,gradient,1)/opt%mass
      ! check for convergence of the gradient norm
      gnorm = sqrt(ddot(3*mol%n,gradient,1,gradient,1))
      if (optcell) then
         inv_lat = mat_inv_3x3(mol%lattice)
         call sigma_to_latgrad(sigma,inv_lat,lat_gradient)
         lat_power = -ddot(9,lat_velocities,1,lat_gradient,1)/opt%lat_mass
         lat_gnorm = sqrt(ddot(9,lat_gradient,1,lat_gradient,1))
      endif

      gconverged = gnorm < opt%g_thr .and. &
         & (.not.optcell .or. lat_gnorm < opt%g_thr)

      if (pr) then
         write(env%unit,'(" * total energy  :",f14.7,1x,"Eh")',advance='no')   energy
         write(env%unit,'(5x,"change    :",e16.7,1x,"Eh")')                    ediff
         write(env%unit,'(3x,"gradient norm :",f14.7,1x,"Eh/a0")',advance='no')gnorm
         write(env%unit,'(2x,"converged?",3x,3(1x,a,"=",l1))') &
            "E",econverged,"G",gconverged,"D",dconverged
         write(env%unit,'(3x,"time step     :",f14.7,1x,"fs")',advance='no') &
            time_step * autofs
         write(env%unit,'(5x,"power     :",e16.7)') power
         write(env%unit,'(3x,"step length   :",f14.7,1x,"a0")',advance='no') &
            step_length
         write(env%unit,'(5x,"speed     :",e16.7)') speed
         if (optcell) then
            write(env%unit,'(3x,"time step     :",f14.7,1x,"fs")',advance='no') &
               lat_time_step * autofs
            write(env%unit,'(5x,"cell power:",e16.7)') lat_power
            write(env%unit,'(3x,"step length   :",f14.7,1x,"a0")',advance='no') &
               lat_step_length
            write(env%unit,'(5x,"cell speed:",e16.7)') lat_speed
         endif
      else if (minpr) then
         write(env%unit,'(i6,f14.7,1x,"(",e14.7,")",2(1x,f14.7),1x,"(",3l1,")")') &
            & iter, energy, ediff, gnorm, step_length, &
            & econverged, gconverged, dconverged
      endif
      if (profile) call timer%measure(5)

      ! check for convergence of the minimization
      converged = econverged.and.gconverged
      if (converged) then
         exit md_step
      endif

      ! update the Hessian
      ! we perform the Hessian update after the convergence check, this avoids
      ! doing an update on a converged geometry, but still updates at the end
      ! of the microcycle (if you drop the Hessian afterwards, this is wasted)
      if (profile) call timer%measure(8,"hessian update")
      if (opt%update) &
         call bfgs(3*mol%n,gnorm,gradient,scratch,displacement,hessp)
      if (profile) call timer%measure(8)

   enddo md_step

end subroutine inertial_relax

!> Calculates predicted energy change according to the second order model.
function predict_energy_change(nat3,grad,displ,hess) result(depred)

   implicit none

   !> 3*natoms
   integer, intent(in) :: nat3
   !> Gradient vector
   real(wp), dimension(nat3), intent(in) :: grad
   !> Displacement vector
   real(wp), dimension(nat3), intent(in) :: displ
   !> Hessian matrix stored as lower triangle
   real(wp), dimension(nat3*(nat3+1)/2), intent(in) :: hess

   !> Predicted energy change
   real(wp) :: depred

   real(wp), dimension(nat3) :: hdx
   real(wp), external :: ddot
   real(wp) :: gtmp, htmp

   call dspmv('u',nat3,0.5_wp,hess,displ,1,0.0_wp,hdx,1)

   gtmp   = ddot(nat3,displ,1,grad,1)

   htmp   = ddot(nat3,displ,1,hdx,1)

   depred = htmp + gtmp

end function predict_energy_change

!> obtain preconditioner by solving P·x = y.
!  P·g⁻¹ is saved in pinvg on exit
subroutine get_preconditioner_metric(nvar,hess,grad,pinvg,fail)
   implicit none
   !> calculation dimension (usually 3*natoms)
   integer, intent(in) :: nvar
   real(wp), intent(in) :: hess(:)
   real(wp), intent(in) :: grad(:,:)
   real(wp), intent(out) :: pinvg(:,:)
   logical, intent(out) :: fail

   integer :: i,ii,j,ij
   integer :: info
   integer :: lwork
   integer, allocatable :: ipiv(:)
   real(wp) :: test(1)
   real(wp), allocatable :: work(:)
   real(wp), allocatable :: htmp(:,:)

   real(wp), parameter :: c_stab = 0.1_wp

   allocate( ipiv(nvar), source = 0 )
   allocate( htmp(nvar,nvar), source = 0.0_wp )

   ! blowup packed hessian (avoids destroying Hessian in the following step)
   do i = 1, nvar
      ii = i*(i-1)/2
      do j = 1, i-1
         ij = ii + j
         htmp(i,j) = hess(ij)
         htmp(j,i) = hess(ij)
      enddo
      ! add stabilization term to diagonal (not tracked by BFGS!)
      htmp(i,i) = hess(ii+i) + c_stab
   enddo

   ! initialize right hand side
   pinvg = -grad

   call dsysv('u',nvar,1,htmp,nvar,ipiv,pinvg,nvar,test,-1,info)
   if (info.eq.0) then
      lwork = int(test(1))
      allocate( work(lwork), source = 0.0_wp )

      call dsysv('u',nvar,1,htmp,nvar,ipiv,pinvg,nvar,work,lwork,info)

      deallocate(work,ipiv)
   endif

   fail = info.ne.0
end subroutine get_preconditioner_metric

pure function bool2string(bool) result(string)
   logical,intent(in) :: bool
   character(len=:),allocatable :: string
   if (bool) then
      string = 'true'
   else
      string = 'false'
   endif
end function bool2string

end module xtb_relaxation_engine
