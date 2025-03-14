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

module xtb_optimizer
   use xtb_io_writer, only : writeMolecule
   use xtb_mctc_accuracy, only : wp, sp
   use xtb_mctc_fileTypes, only : fileType
   use xtb_type_environment, only : TEnvironment
   use xtb_extern_turbomole, only : TTMCalculator
   use xtb_bfgs
   use xtb_hessian, only : numhess
   use xtb_david2
   implicit none

   !> time profiling
   logical,private,parameter :: profile = .true.

   type :: convergence_log
      integer :: nlog
      real(wp), allocatable :: elog(:)
      real(wp), allocatable :: glog(:)
   contains
      procedure :: set_eg_log
      procedure :: get_averaged_energy
      procedure :: get_averaged_gradient
   end type convergence_log
   interface convergence_log
      module procedure new_convergence_log
   end interface convergence_log

contains

!> construct optimizer settings from optimization level
subroutine get_optthr(n,olev,ethr,gthr,maxcycle,acc)
   use xtb_setparam
   implicit none
   integer, intent(in)  :: n
   integer, intent(in)  :: olev
   real(wp),intent(out) :: ethr
   real(wp),intent(out) :: gthr
   integer, intent(out) :: maxcycle
   real(wp),intent(out) :: acc
   
   select case(olev)
   ! very approximate = crude !
   case(p_olev_crude)
      ethr   = 5.d-4
      gthr   = 1.d-2
      maxcycle=n
      acc=3.00d0
   
   ! approximate = sloopy !
   case(p_olev_sloppy)
      ethr   = 1.d-4
      gthr   = 6.d-3
      maxcycle=n
      acc=3.00d0
   
   ! loose !
   case(p_olev_loose)
      ethr   = 5.d-5
      gthr   = 4.d-3
      maxcycle=n*2
      acc=2.00d0

   ! for DCOSMO-RS opts with TM i.e. between loose and normal !
   case(p_olev_lax)
      ethr   = 2.d-5
      gthr   = 2.5d-3
      maxcycle=n*2
      acc=2.00d0
   
   ! normal !
   case default
      ethr  = 5.d-6
      gthr  = 1.d-3
      maxcycle=n*3
      acc=1.0d0
   
   ! tight !
   case(p_olev_tight)
      ethr   = 1.d-6
      gthr   = 8.d-4
      maxcycle=n*5
      acc=0.20d0

   ! very tight !
   case(p_olev_vtight)
      ethr   = 1.d-7
      gthr   = 2.d-4
      maxcycle=n*20
      acc=0.05d0

   ! extreme !
   case(p_olev_extreme)
      ethr   = 5.d-8
      gthr   = 5.d-5
      maxcycle=n*20
      acc=0.01d0
   end select
   
   maxcycle=min(maxcycle,10000)
   maxcycle=max(maxcycle,200)

end subroutine get_optthr

!> Approximate Normal Coordinate rational function optimizer 
subroutine ancopt(env,ilog,mol,chk,calc, &
      &           egap,et,maxiter,maxcycle_in,etot,g,sigma,tight,pr,fail)
   
   use xtb_mctc_convert
   use xtb_mctc_la

   use xtb_type_molecule
   use xtb_type_anc
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_timer

   use xtb_setparam
   use xtb_fixparam

   use xtb_setmod, only : int2optlevel

   use xtb_axis, only : axis2
   use xtb_hessian, only : trproj,rdhess
   use xtb_readin
   use xtb_lsrmsd

   implicit none
   
   !> traceback for error handling 
   character(len=*), parameter :: source = "optimizer_ancopt"
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular structure data
   type(TMolecule), intent(inout) :: mol
   
   !> optimization level
   integer, intent(in)    :: tight
   
   !> max number of SCC cycles
   integer, intent(in)    :: maxiter
   
   !> max number of optimization cycles
   integer, intent(in)    :: maxcycle_in
   
   !> wrapper for changing info during SCC
   type(TRestart),intent(inout) :: chk

   !> polymorphic calculator instance
   class(TCalculator), intent(inout) :: calc
   
   !> electronic energy
   real(wp) :: eel

   !> total energy
   real(wp),intent(inout) :: etot
   
   !> electronic temperature
   real(wp),intent(in)    :: et
   
   !> HOMO-LUMO gap
   real(wp),intent(inout) :: egap
   
   !> gradients
   real(wp),intent(inout) :: g(3,mol%n)
   
   !> strain derivatives
   real(wp),intent(inout) :: sigma(3,3)
   
   !> printlevel
   logical, intent(in)    :: pr
   
   !> optimization failure
   logical, intent(out)   :: fail

!-----------------!   
! local variables !
!-----------------!   
   
   type(TMolecule) :: molopt
   type(scc_results) :: res
   type(tb_anc) :: anc
   type(tb_timer) :: timer

   !> maximum coordinate displacement
   real(wp) :: maxdispl

   !> lowest value for force constant
   real(wp) :: hlow

   !> highest value for force constant
   real(wp) :: hmax 

   !> dispersion scaling
   real(wp) :: s6

   !> total energy from initial single point callculation
   real(wp) :: estart

   !> energy & gradient convergence thresholds
   real(wp) :: ethr, gthr
   
   
   !> number of modes followed
   integer :: modef
   
   !> max number of opt cycles
   integer :: maxopt

   !> max number of opt cycles before generation of new ANC
   integer :: maxmicro
   
   !> algorithm for updating Hessian  
   integer :: iupdat ! 0 = BFGS, 1 = Powell
   
   !> current iteration 
   integer :: iter

   !> number of atomic coordinates
   integer :: nat3

   !> deegrees of freedom
   integer :: nvar

   !> work array size
   integer :: lwork, liwork

   !> step size for numerical Hessian
   real(wp) :: step_hess

   real(wp) :: step,amu2au,au2cm,dumi,dumj,damp,edum,thr,aaa,bbb
   real(wp) :: energy,acc,rij(3),t1,t0,w1,w0,ccc

   integer  :: n3,i,j,k,l,jjj,ic,jc,ia,ja,ii,jj,info
   integer  :: nread,itry,iii
   integer  :: id,ihess,error
   integer, intent(in)  :: ilog
   integer, external    :: lin

   !> Hessian matrix
   real(wp),allocatable :: h (:,:), hess_tmp(:,:)

   real(wp),allocatable :: b (:,:)
   real(wp),allocatable :: pmode(:,:)
   real(wp),allocatable :: grmsd(:,:)

   !> force constants
   real(wp),allocatable :: fc(:)

   !> eigenvalues
   real(wp),allocatable :: eig(:)

   real(wp),allocatable :: aux(:)
   real(wp),allocatable :: hess(:)
   integer, allocatable :: iwork(:)
   integer, allocatable :: totsym(:)
   type(convergence_log), allocatable :: avconv
   
   real(wp) :: U(3,3), x_center(3), y_center(3), rmsdval
   real(wp) :: esave
   logical :: restart,ex,converged
   
   !> if linear molecule
   logical :: linear

   !> dipole
   real(wp),allocatable :: dipgrad(:,:)
   integer, allocatable :: list(:)
   
   !> formatting strings for output 
   character(len=*),parameter :: scifmt = &
      '(10x,":",3x,a,e22.7,1x,a,1x,":")'
   character(len=*),parameter :: dblfmt = &
      '(10x,":",3x,a,f18.7,5x,a,1x,":")'
   character(len=*),parameter :: intfmt = &
      '(10x,":",3x,a,i18,      10x,":")'
   character(len=*),parameter :: chrfmt = &
      '(10x,":",3x,a,a18,      10x,":")'
  
   !> error string
   character(len=128) :: errStr
   
   !> debug mode
   logical, parameter :: debug(2) = [.false.,.false.]
   character(len=9):: hessfmt

   ! print ANCopt header !
   call ancopt_header(env%unit,set%veryverbose)
   
   if(mol%n.eq.1) return ! skip optimization for 1 atom 

   ! performance profiling timer !
   if (profile) call timer%new(8,.false.)
   if (profile) call timer%measure(1,'optimizer setup')
   
   ! defaults !
   iter = 0
   fail  = .false.
   modef = 0
   iupdat= 0 
   hmax  = 5.0_wp
   nat3 = 3 * mol%n
   maxdispl = set%optset%maxdispl_opt  ! def: 1.0
   hlow = set%optset%hlow_opt          ! def: 0.01 
   s6   = set%mhset%s6                 ! def: 20.0
   maxmicro=set%optset%micro_opt       ! def: 20 
   estart = etot
   call get_optthr(mol%n,tight,ethr,gthr,maxopt,acc)
   
   ! if provided by user !
   if(maxcycle_in.gt.0)then
      maxopt=maxcycle_in
   endif
   if(maxopt.lt.maxmicro) maxmicro=maxopt
   
   ! use Powell update if TS optimization ! 
   if(set%tsopt)then
      hlow=max(hlow,0.250d0)
      iupdat=1
   endif
   
   ! energy and gradient averaging !
   if (set%optset%average_conv) then ! default: .false.
      select type(calc)
      class is(TTMCalculator)
         avconv = load_turbomole_log(maxopt)
         if (avconv%nlog > 0 .and. pr) then
            write(env%unit, '(a, 1x, i0, 1x, a)') &
               "Convergence averaging initialized with", avconv%nlog, "entries"
         end if
      class default
         avconv = convergence_log(maxopt)
      end select
   end if

   ! determine if linear molecule via rotational constants !
   call axis2(mol%n,mol%xyz,aaa,bbb,ccc,dumi,dumj)
   if (ccc.lt.1.d-10) then
      linear = .true.
      nvar = nat3 - 5
   else
      linear = .false.
      nvar = nat3 - 6
   endif

   ! exact fixing case ! 
   if(fixset%n.gt.0) then 
      nvar = nat3 - 3*fixset%n - 3
      if(nvar.le.0) nvar=1
   endif

   ! adjust if restrated !
   call open_binary(id,'.xtbtmpmode','r')
   if(id.ne.-1)then
      read (id) modef
      allocate(pmode(nat3,modef))
      read (id) pmode
      call close_file(id)
      nvar=nvar-modef
   else
      allocate(pmode(nat3,1)) ! dummy allocated
   endif

   ! print ANCopt settings !
   if(pr)then
      write(env%unit,'(/,10x,51("."))')
      write(env%unit,'(10x,":",22x,a,22x,":")') "SETUP"
      write(env%unit,'(10x,":",49("."),":")')
      write(env%unit,chrfmt) "optimization level", int2optlevel(tight)
      write(env%unit,intfmt) "max. optcycles    ", maxopt
      write(env%unit,intfmt) "ANC micro-cycles  ", maxmicro
      write(env%unit,intfmt) "degrees of freedom", nvar
      
      if (modef>0) then
         write(env%unit,intfmt) "# mode follow     ", modef
      endif
      
      write(env%unit,'(10x,":",49("."),":")')
      
      if (set%optset%exact_rf) then
         write(env%unit,chrfmt) "RF solver         ", "spevx"
      else
         write(env%unit,chrfmt) "RF solver         ", "davidson"
      endif
      
      write(env%unit,chrfmt) "write xtbopt.log  ", bool2string(ilog.ne.-1)
      
      if (linear) then
         write(env%unit,chrfmt) "linear (good luck)", bool2string(linear)
      else
         write(env%unit,chrfmt) "linear?           ", bool2string(linear)
      endif
      
      write(env%unit,scifmt) "energy convergence", ethr,    "Eh  "
      write(env%unit,scifmt) "grad. convergence ", gthr,    "Eh/α"
      write(env%unit,dblfmt) "maximum RF displ. ", maxdispl,"    "
      write(env%unit,scifmt) "Hlow (freq-cutoff)", hlow,    "    "
      write(env%unit,dblfmt) "Hmax (freq-cutoff)", hmax,    "    "
      write(env%unit,dblfmt) "S6 in model hess. ", s6,      "    "
      write(env%unit,'(10x,51("."))')
   endif

   ! work arrays !
   lwork  = 1 + 6*nat3 + 2*nat3**2
   liwork = 8 * nat3

   allocate(h(nat3,nat3),fc(nat3*(nat3+1)/2),eig(nat3))

   ! read in Hessian from file !
   if (set%mhset%model == p_modh_read) then
      call open_file(ihess, 'hessian', 'r')
      if (ihess == -1) then
         call env%error("Could not read in hessian as requested.", source)
         return
      endif
      call read_hessian(ihess, nat3, h, error)
      if (error /= 0) then
         call env%error("Could not read hessian from file.", source)
         return
      endif

      ! do not reset the hessian
      maxmicro = maxopt
      ex = .true.
   else
      ex = .false.
   endif

   call anc%allocate(mol%n,nvar,hlow,hmax) ! allocate ANC
   molopt = mol ! copy molecular information
   if (profile) call timer%measure(1) ! start opt timer

! ======================================================================
   ANC_microiter: do
! ======================================================================

!----------------------------------------------------------------!
!--------------------- Hessian generation -----------------------!
!----------------------------------------------------------------!

   if (profile) call timer%measure(2,'model hessian') ! start timer for Hessian
   
   ! initial guess for Hessian !
   if (.not.ex)then ! normal case
      
      if (pr) &
         &  write(env%unit,'(/,''generating ANC from exact Hessian ...'')')
      
      call modhes(env, calc, set%mhset, molopt%n, molopt%xyz, molopt%at, fc, pr)   ! def: Lindh model (1995)
      call env%check(fail)
      if (fail) then
         call env%error("Calculation of model hessian failed", source)
         return
      end if

      thr=1.d-11
   
   ! read in Hessian from file !
   else 
      
      if(pr) &
         &  write(env%unit,'(/,''generating ANC from read Hessian ...'')')
      
      ! pack Hessian !
      k=0
      do i=1,nat3
         do j=1,i
            k=k+1
            fc(k)=h(j,i)
         enddo
      enddo

      thr=1.d-10

   endif

   ! project out translational and rotational modes !
   if(modef.eq.0)then
      if(fixset%n.gt.0)then
         call trproj(molopt%n,nat3,molopt%xyz,fc,.false., -1  ,pmode,1)     ! exact fixing
      else
         if (.not.linear) &
         call trproj(molopt%n,nat3,molopt%xyz,fc,.false.,modef,pmode,1)     ! normal
      endif
   else
      call trproj(molopt%n,nat3,molopt%xyz,fc,.false.,modef,pmode,modef) ! NMF
   endif

   ! blow up Hessian !
   k=0
   do i=1,nat3
      do j=1,i
         k=k+1
         h(i,j)=fc(k)
         h(j,i)=fc(k)
      enddo
   enddo

   if (debug(2) .and. nat3 <= 30) then !######## DEBUG ########
      write(env%unit,'(/,''Hessian matrix'')')
      write(hessfmt,'(a,i0,a)') '(', nat3, 'F10.6)'
      write(env%unit,hessfmt) (h(:,i), i=1,nat3)
   endif
   
   if (profile) call timer%measure(2) ! stop timer for model Hessian

!----------------------------------------------------------------!
!---------------------- ANC generation --------------------------!
!----------------------------------------------------------------!

   if (profile) call timer%measure(3,'ANC generation') ! start timer for ANC generation   

   ! diagonalize Hessian and sort eigenvalues !
   call anc%new(env%unit,molopt%xyz,h,pr,linear)  
   if (profile) call timer%measure(3) ! stop timer for ANC generation
   esave = etot ! save energy 

!----------------------------------------------------------------!
!---------------------- Optimization ----------------------------!
!----------------------------------------------------------------!
   call relax(env,iter,molopt,anc,restart,maxmicro,maxdispl,ethr,gthr, &
      & iii,chk,calc,egap,acc,et,maxiter,iupdat,etot,g,sigma,ilog,pr,fail, &
      & converged,timer,set%optset%exact_rf,avconv)

   ! check for errors !
   call env%check(fail)
   if (fail) then
      call env%error("Could not relax/optimize structure", source)
      return
   endif

   ! dynamically adjust number of micro iterations !
   maxmicro=min(int(maxmicro*1.1),2*set%optset%micro_opt)

   ! assess the optimization by RMSD change !
   call rmsd(molopt%n,anc%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)

   ! this comes close to a goto, but it's not a goto ... it's even worse !
   if (restart.and.iter.lt.maxopt) then
      if (pr) then
         write(env%unit,'(" * RMSD in coord.:",f14.7,1x,"α")',advance='no') rmsdval
         write(env%unit,'(6x,"energy gain",e16.7,1x,"Eh")') etot-esave
      end if
      cycle ANC_microiter
   endif
   exit  ANC_microiter
! ======================================================================
   enddo ANC_microiter
! ======================================================================

   ! convergence report !
   if (converged) then
      if(pr) then
         call rmsd(mol%n,mol%xyz,molopt%xyz,1,U,x_center,y_center,rmsdval,.false.,grmsd)
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "GEOMETRY OPTIMIZATION CONVERGED AFTER",iter,"ITERATIONS"
         write(env%unit,'(72("-"))')
         write(env%unit,'(1x,"total energy gain   :",F18.7,1x,"Eh",F14.4,1x,"kcal/mol")') &
            etot-estart, (etot-estart)*autokcal
         write(env%unit,'(1x,"total RMSD          :",F18.7,1x,"a0",F14.4,1x,"Å")') &
            rmsdval, rmsdval*autoaa
         if (profile) then
            write(env%unit,'(1x,"total power (kW/mol):",F18.7,1x,"(step)",F10.4,1x,"(real)")') &
               & (etot-estart)*autokJ/iter, (etot-estart)*autokJ/timer%get()
         endif
         write(env%unit,'(72("-"))')
      endif
   else
      ! not converging in the given cycles is a FAILURE, we should make this clearer !
      ! This is still no ERROR, since we want the geometry written afterwards !
      if(pr) then
         write(env%unit,'(/,3x,"***",1x,a,1x,i0,1x,a,1x,"***",/)') &
            "FAILED TO CONVERGE GEOMETRY OPTIMIZATION IN",iter,"ITERATIONS"
      endif
   endif

   mol = molopt ! copy optimized geometry back to molecule

   if (pr.and.profile) call timer%write(env%unit,'ANCopt') ! write timer report

   ! cleanup !
   if (profile) call timer%deallocate
   if (allocated(pmode))  deallocate(pmode)
   if (allocated(h))      deallocate(h)
   if (allocated(fc))     deallocate(fc)
   call anc%deallocate

end subroutine ancopt


!> @brief actual optimization procedure
!!
!! note that this subroutines implements the main loop of the optimization
!! procedure INSIDE the mirco iterator which is implemented by the calling
!! routine, the only way to leave this subroutine is by an EARLY return,
!! without going over the RESTART logical at the end of the subroutine.
!* I have warned you, be careful not to break anything.
subroutine relax(env,iter,mol,anc,restart,maxcycle,maxdispl,ethr,gthr, &
      &          ii,chk,calc, &
      &          egap,acc_in,et,maxiter,iupdat,etot,g,sigma,ilog,pr,fail,converged, &
      &          timer,exact,avconv)

   use xtb_mctc_blas, only : blas_gemv
   use xtb_type_molecule
   use xtb_type_anc
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_type_timer

   implicit none

   !> source of errors in the main program unit
   character(len=*), parameter :: source = "optimizer_relax"

   !> calculation environment
   type(TEnvironment), intent(inout) :: env

   !> molecular structure data
   type(TMolecule),    intent(inout) :: mol
   
   !> timer instance
   type(tb_timer),       intent(inout) :: timer
   
   !> ANC instance
   type(tb_anc),         intent(inout) :: anc
   
   !> WFN 
   type(TRestart),intent(inout) :: chk
   
   !> polymorphic calculator instance
   class(TCalculator), intent(inout) :: calc
   
   !> number of micro iterations 
   integer, intent(in)    :: maxiter

   !> Hessian update algorithm
   integer, intent(in)    :: iupdat

   !> iunit for log file
   integer, intent(in)    :: ilog

   !> temperature
   real(wp),intent(in)    :: et
   
   !> total energy
   real(wp),intent(inout) :: etot
   
   !> SCC accuracy
   real(wp),intent(in)    :: acc_in

   !> gradients
   real(wp),intent(inout) :: g(3,mol%n)

   !> strain derivatives
   real(wp),intent(inout) :: sigma(3,3)

   !> HOMO-LUMO gap
   real(wp),intent(inout) :: egap

   !> conv check
   logical, intent(out)   :: converged

   !> solver type
   logical, intent(in)    :: exact
   
   !> averaging
   type(convergence_log), intent(inout), optional :: avconv

   !> SCC results storage
   type(scc_results) :: res

   !> dimensions for RFO
   integer  :: nvar1,npvar,npvar1

   !> local booleans to control optimization
   logical  :: restart, first, pr, fail
   logical  :: econverged
   logical  :: gconverged
   logical  :: lowered

   integer  :: maxcycle, iter, prlevel
   integer  :: i, j, ii, jj, k, lwork, info, m, idum, imax(3)
   real(wp) :: energy,ethr,gthr,dsnrm,maxdispl,t0,w0,t1,w1
   real(wp) :: lambda,gnorm,dnorm,ddot,eold,xdum,estart,acc,e_in
   real(wp) :: depred,echng,dummy,maxd,alp,gchng,smallreal,gnold
   real(wp),allocatable :: gold(:)
   real(wp),allocatable :: displ(:), gint(:)

   !> RFO eigenvalues
   real(sp),allocatable :: eaug(:)

   !> RFO eigenvectors
   real(sp),allocatable :: Uaug(:,:)

   !> Augmented Hessian
   real(sp),allocatable :: Aaug(:)

   real(sp) :: r4dum,sdot
   parameter (r4dum=1.e-8)
   parameter (smallreal=1.d-14)

!----------------------------------------------------------------!
!--------------------- Initialization ---------------------------!
!----------------------------------------------------------------!
   

   ! set printlevel !
   if (pr) then 
      prlevel = 1
   else 
      prlevel=0
   endif
   
   ! initialize variables !
   gnorm  = 0.0_wp
   depred = 0.0_wp
   echng  = 0.0_wp
   alp    = 1.0_wp
   maxd   = maxdispl
   acc    = acc_in
   energy = etot
   e_in   = etot
   first  = .true.
   converged = .false.

   nvar1  = anc%nvar+1         ! dimension of RF calculation
   npvar  = anc%nvar*(nvar1)/2 ! packed size of Hessian (note the abuse of nvar1!)
   npvar1 = nvar1*(nvar1+1)/2  ! packed size of augmented Hessian
   
   allocate( gold(anc%nvar), displ(anc%nvar), gint(anc%nvar), source = 0.0_wp )
   allocate( Uaug(nvar1,1),eaug(nvar1),Aaug(npvar1) )

!! ========================================================================
   main_loop: do ii=1,maxcycle
!! ========================================================================
   
   iter=iter+1 ! iteration counter
   if(pr) &
      write(env%unit,'(/,72("."),/,30(".")," CYCLE",i5,1x,30("."),/,72("."))')iter

   ! save values from previous iteration !
   gold = gint
   gnold= gnorm
   eold = energy

   ! dE via 2-nd order Taylor expansion !
   ! E = E0 + delta * G + delta^2 * H !
   if (ii > 1) &
      call prdechng(anc%nvar,gold,displ,anc%hess,depred)
  
   ! displace cartesian coordinates !
   if (profile) call timer%measure(4,'coordinate transformation')
   call anc%get_cartesian(mol%xyz)
   if (profile) call timer%measure(4)

   ! single point + analytical gradients !
   if (profile) call timer%measure(5,'single point calculation')
   g = 0.0_wp
   call calc%singlepoint(env,mol,chk,prlevel,iter.eq.1,energy,g,sigma,egap,res)
   if (profile) call timer%measure(5)
   
   ! something went wrong in SCC or diag !
   call env%check(fail)
   if (fail) then
      call env%error('SCF not converged, aborting...', source)
      return
   endif
   if (energy.gt.1.d42) then
      call env%error('energy is bogus! aborting...', source)
      fail=.true.
      return
   endif
   
   ! write geometry to log file !
   if (profile) call timer%measure(6,'optimization log')
   call writeMolecule(mol, ilog, format=fileType%xyz, energy=res%e_total, &
      & gnorm=res%gnorm)
   if (profile) call timer%measure(6)


   ! transform xyz (g) to internal gradient (gint) !
   if (profile) call timer%measure(4)
   call dgemv('t',anc%n3,anc%nvar,1.0_wp,anc%B,anc%n3,g,1,0.0_wp,gint,1)
   if (profile) call timer%measure(4)

   ! check Euclidean norm of the normal gradient !
   gnorm = norm2(gint)
   if(gnorm.gt.500.) then   
      call env%error('|grad| > 500, something is totally wrong!', source)
      fail=.true.
      return
   endif

   ! average energy and gradient !
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

   ! adjust SCC accuracy dynamically !
   if(gnorm    .lt.0.004)then
      acc=acc_in
   elseif(gnorm.lt.0.02)then
      acc=3.0d0*acc_in
   elseif(gnorm.gt.0.02)then
      acc=6.0d0*acc_in
   endif

   first =.false. ! distinguish between first and subsequent iterations

   ! calculate the change in energy and gradient norm !
   gchng = gnorm -gnold
   echng = energy-eold

   ! check for convergence !
   econverged = abs(echng).lt.ethr
   gconverged = gnorm.lt.gthr
   lowered    = echng.lt.0.0_wp

   ! print energy and gradient norm of a current cycle !
   if(pr) then
      write(env%unit,'(" * total energy  :",f14.7,1x,"Eh")',advance='no')   energy
      write(env%unit,'(5x,"change   ",e18.7,1x,"Eh")')                      echng
      write(env%unit,'(3x,"gradient norm :",f14.7,1x,"Eh/α")',advance='no') gnorm
      write(env%unit,'(3x,"predicted",e18.7)',advance='no')                 depred
      write(env%unit,'(1x,"("f7.2"%)")')         (depred-echng)/(echng+1e-34_wp)*100
   endif
   
   ! check 0 energy case !
   if ( energy .eq. 0 ) then
      call env%error('external program error', source)
      return
   end if

   if(ii.eq.1) estart=energy ! save energy of first iteration

   ! adjust step size dynamically !
   if(gnorm.lt.0.002)then 
      alp = 1.5d0       ! 1.5
   elseif(gnorm.lt.0.0006)then
      alp = 2.0d0       ! 2
   elseif(gnorm.lt.0.0003)then
      alp = 3.0d0       ! 3
   else
      alp = 1.0d0
   endif

!! ------------------------------------------------------------------------
!  Update the hessian
!! ------------------------------------------------------------------------

   if (profile) call timer%measure(7,'hessian update')
   if(ii.gt.1)then
      if(iupdat.eq.0)then
         call bfgs  (anc%nvar,gnorm,gint,gold,displ,anc%hess)
      else
         call powell(anc%nvar,gnorm,gint,gold,displ,anc%hess)
      endif
   endif
   if (profile) call timer%measure(7)

! ------------------------------------------------------------------------
!  the following lines define the rational function (RF) method
! ------------------------------------------------------------------------
   if (profile) call timer%measure(8,"rational function")
!  funny that implementing a RF algorithm takes 5 lines when done right... SAW
!  As it is now condensed down to a few lines, I feel the urge to explain
!  how the calculation is actually done. We solve this:
!   ⎛ H  g ⎞ ⎛ dx ⎞     ⎛ dx ⎞
!   ⎝ g  0 ⎠ ⎝  1 ⎠ = λ ⎝  1 ⎠
!  first augment hessian by gradient, we keep everything nicely packed, no blowup
   Aaug(1:npvar)          = anc%hess ! pack hessian
   Aaug(npvar+1:npvar1-1) = gint ! augment with gradient
   Aaug(npvar1)           = 0.0_sp ! add zero
   
   ! chose your favourite solver !
   if (exact .or. nvar1.lt.50) then
      call solver_sspevx(nvar1,r4dum,Aaug,Uaug,eaug,fail)
   else
      ! steepest decent guess for displacement !
      if (ii.eq.1) then
         Uaug(:,1)=[-real(gint(1:anc%nvar),sp),1.0_sp]
         dsnrm = sqrt(sdot(nvar1,Uaug,1,Uaug,1))
         Uaug  = Uaug/dsnrm
      endif
      call solver_sdavidson(nvar1,r4dum,Aaug,Uaug,eaug,fail,.false.)
      if (fail) & ! retry with better solver
      call solver_sspevx(nvar1,r4dum,Aaug,Uaug,eaug,fail)
   endif
   
   ! error handling after RF calculation !
   if (fail .or. abs(Uaug(nvar1,1)).lt.1.e-10) then
      call env%error("internal rational function error", source)
      return
   end if

   ! divide by last element to get the displacement vector !
   displ(1:anc%nvar) = Uaug(1:anc%nvar,1)/Uaug(nvar1,1)
   
   ! check if step is too large, cut off everything that is too large !
   do j=1,anc%nvar
      if(abs(displ(j)).gt.maxd) then
         if(displ(j) < 0) displ(j)=-maxd
         if(displ(j) > 0) displ(j)= maxd
      endif
   enddo
   
   ! now some output !
   dsnrm=sqrt(ddot(anc%nvar,displ,1,displ,1)) ! displacement norm
   if(pr)then
      
      ! find the largest displacement !
      gold = abs(displ)
      imax(1) = maxloc(gold,1); gold(imax(1)) = 0.0_wp
      imax(2) = maxloc(gold,1); gold(imax(2)) = 0.0_wp
      imax(3) = maxloc(gold,1)

      write(env%unit,'(3x,"displ. norm   :",f14.7,1x,"α")',advance='no') &
         dsnrm*alp
      write(env%unit,'(6x,"lambda   ",e18.7)') eaug(1) ! eigenvalue of the RF 
      write(env%unit,'(3x,"maximum displ.:",f14.7,1x,"α")',advance='no') &
         abs(displ(imax(1)))*alp
      write(env%unit,'(6x,"in ANC''s ",3("#",i0,", "),"...")') imax
   endif
   
   if (profile) call timer%measure(8)

   ! 2nd: exit and redo hessian (internal restart) !
   if(ii.gt.2.and.dsnrm.gt.2.0) then
      if (pr) write(*,*) 'exit because of too large step'
      exit main_loop
   endif

   ! new coordinates !
   anc%coord = anc%coord + displ * alp

   ! conv check !
   if(abs(echng).lt.ethr.and.gnorm.lt.gthr.and.echng.lt.1.0e-10_wp) then
      restart=.false.
      converged = .true.
      etot=energy
      return
   endif
!! ========================================================================
   enddo main_loop
!! ========================================================================
   
   ! cleanup !
   if (allocated(Uaug))  deallocate(Uaug)
   if (allocated(eaug))  deallocate(eaug)
   if (allocated(Aaug))  deallocate(Aaug)

   ! post micro cycle processing !
   restart=.true.
   etot=energy
   call anc%get_cartesian(mol%xyz)

end subroutine relax

pure subroutine solver_ssyevx(n,thr,A,U,e,fail)
   use xtb_mctc_lapack, only : lapack_syevx
   implicit none
   integer, intent(in)    :: n
   real(sp),intent(in)    :: thr
   real(sp),intent(inout) :: A(:,:)
   real(sp),intent(inout) :: U(:,:)
   real(sp),intent(inout) :: e(:)
   logical, intent(out)   :: fail

   integer :: i,j,k
   integer :: lwork,info
   real(sp),allocatable :: work(:)
   integer, allocatable :: iwork(:)
   integer, allocatable :: ifail(:)
   real(sp) :: dum

   lwork = 1+8*n+n**2
   allocate(iwork(5*n),work(lwork),ifail(n))

   j=1
   call lapack_syevx('V','I','U',n,A,n,dum,dum,j,j,thr, &
   &           i,e,U,n,work,lwork,iwork,ifail,info)
   fail = .false.
   if (info.ne.0) fail = .true.

   deallocate(iwork,work,ifail)
end subroutine solver_ssyevx

pure subroutine solver_sspevx(n,thr,A,U,e,fail)
   use xtb_mctc_lapack, only : lapack_spevx
   implicit none
   integer, intent(in)    :: n
   real(sp),intent(in)    :: thr
   real(sp),intent(inout) :: A(:)
   real(sp),intent(inout) :: U(:,:)
   real(sp),intent(inout) :: e(:)
   logical, intent(out)   :: fail

   integer :: i,j,k
   integer :: info
   real(sp),allocatable :: work(:)
   integer, allocatable :: iwork(:)
   integer, allocatable :: ifail(:)
   real(sp) :: dum

   allocate(iwork(5*n),work(8*n),ifail(n))

   j=1
   call lapack_spevx('V','I','U',n,A,dum,dum,j,j,thr,i,e,U,n,work,iwork,ifail,info)
   fail = .false.
   if (info.ne.0) fail = .true.

   deallocate(iwork,work,ifail)
end subroutine solver_sspevx

subroutine findlroot(i0,n,e,el)
   implicit none
   integer, intent(in)  :: i0
   integer, intent(in)  :: n
   real(wp),intent(in)  :: e(n)
   real(wp),intent(out) :: el
   integer :: i

   el=1.d+42
   do i=1,n
      if(i.eq.i0.and.i0.gt.0) cycle
      if(abs(e(i)) .gt. 1.d-10 .and. e(i) .lt. el ) el = e(i)
   enddo

end subroutine findlroot

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine sort(nat3,nvar,hess,b)
   implicit none
   integer, intent(in) :: nat3
   integer, intent(in) :: nvar
   real(wp),intent(inout) :: hess(nvar*(nvar+1)/2)
   real(wp),intent(inout) :: b(nat3,nat3)
   real(wp),allocatable :: edum(:)
   real(wp) :: pp,sc1
   integer  :: ii,k,j,m,i
   allocate( edum(nvar), source = 0.0_wp )

   do k=1,nvar
      edum(k)=hess(k+k*(k-1)/2)
   enddo
! sort
   do ii = 2, nvar
      i = ii - 1
      k = i
      pp= edum(i)
      do j = ii, nvar
         if (edum(j) .gt. pp) cycle
         k = j
         pp= edum(j)
      enddo
      if (k .eq. i) cycle
      edum(k) = edum(i)
      edum(i) = pp
      do m=1,nat3
         sc1=b(m,i)
         b(m,i)=b(m,k)
         b(m,k)=sc1
      enddo
   enddo

   do k=1,nvar
      hess(k+k*(k-1)/2)=edum(k)
   enddo

end subroutine sort

!> calculate predicted energy change
subroutine prdechng(nat3,grad,displ,hess,depred)
!---------------------------------------------------------------------
! Purpose:
! Calculates predicted energy change according to the second order
! model.
!
! Input:
! nat3  - 3*natoms
! hess  - Hessian matrix stored as lower triangle
! grad  - Gradient vector
! displ - Displacement vector
!
! Output:
! depred - Predicted energy change
!---------------------------------------------------------------------
   use xtb_mctc_blas, only : blas_spmv
   implicit none

! Input:
   integer, intent(in) :: nat3
   real(wp),intent(in) :: grad(nat3)
   real(wp),intent(in) :: displ(nat3)
   real(wp),intent(in) :: hess(nat3*(nat3+1)/2)

! Output:
   real(wp),intent(out) :: depred

! Local:
   real(wp), allocatable :: hdx(:)
   real(wp) :: gtmp, htmp

! BLAS functions:
   real(wp), external :: ddot
!---------------------------------------------------------------------
   allocate( hdx(nat3), source = 0.0_wp )

   call blas_spmv('u',nat3,0.5d0,hess,displ,1,0.0d0,hdx,1) ! hdx = 0.5*H*displ

   gtmp   = ddot(nat3,displ,1,grad,1) ! gtmp = grad*displ

   htmp   = ddot(nat3,displ,1,hdx,1) ! htmp = hdx*displ

   depred = htmp + gtmp

end subroutine prdechng



subroutine trfp2xyz(nvar,nat3,p,xyz0,h,dspl)
   implicit none
   integer, intent(in)  :: nat3
   integer, intent(in)  :: nvar
   integer :: nat,icount,i,j,k
   real(wp),intent(in)  :: xyz0(3,nat3/3)
   real(wp),intent(out) :: dspl(3,nat3/3)
   real(wp),intent(in)  :: h(nat3,nat3)
   real(wp),intent(in)  :: p(nvar)
   real(wp) :: dum

   dspl = 0.0d0
   nat=nat3/3

! generate cartesian displacement vector
   do i=1,nvar
      icount=0
      do j=1,nat
         do k=1,3
            icount=icount+1
            dum=h(icount,i)*p(i)
            dspl(k,j)=dspl(k,j)+dum
         enddo
      enddo
   enddo

   dspl=dspl + xyz0

   return
end subroutine trfp2xyz


subroutine prdispl(nvar,displ)
   implicit none
   integer, intent(in)  :: nvar
   real(wp),intent(in)  :: displ(nvar)
   real(wp),allocatable :: er(:)
   integer, allocatable :: merk(:)
   integer  :: i,j,ii,k
   integer  :: ihilf
   real(wp) :: pp
   allocate( er(nvar), source = 0.0_wp )
   allocate( merk(nvar), source = 0 )

   er = abs(displ)

   do i=1,nvar
      merk(i)=i
   enddo
   do ii=2,nvar
      i =ii-1
      k =i
      pp=er(i)
      do j=ii, nvar
         if(er(j).le.pp) cycle
         k=j
         pp=er(j)
      enddo
      if (k.eq.i) cycle
      er(k)=er(i)
      er(i)=pp
      ihilf=merk(i)
      merk(i)=merk(k)
      merk(k)=ihilf
   enddo

   write(*,'(''Largest |displ|/coords:'',5(f8.4,'' ('',i4,'')''))') &
      (er(i),merk(i),i=1,min(3,nvar))

end subroutine prdispl

!> generate model Hessian
subroutine modhes(env, calc, modh, natoms, xyz, chg, Hess, pr)
   use xtb_type_setvar
   use xtb_modelhessian
   use xtb_setparam
   use xtb_type_calculator
   use xtb_gfnff_calculator
!
!       generates a Lindh Model Hessian
!       Chem. Phys. Let. 241(1995) 423-428
!       with personal permission from R.L.
!
   implicit none

   !> Source of errors in the main program unit
   character(len=*), parameter :: source = "optimizer_modhes"

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   !> Calculator
   class(TCalculator), intent(inout) :: calc

   type(modhess_setvar),intent(in) :: modh
   logical, intent(in)  :: pr

   !> other variables
   integer  :: i
   integer  :: nhess
   integer, intent(in)  :: natoms
   real(wp),intent(in)  :: xyz(3,natoms)
   
   !> model hessian
   real(wp),intent(out) :: Hess((natoms*3)*((natoms*3)+1)/2)
   integer, intent(in)  :: chg(natoms)

   ! initialize !
   Hess=0.0_wp

   select type(calc)
   class default ! all calculator cases except GFN-FF
      select case(modh%model) 
      case default
         call env%error("internal error in model hessian!", source)
         return
      case(p_modh_old)
        if (pr) write(env%unit,'(a)') "Using Lindh-Hessian (1995)"
        call ddvopt(xyz, natoms, Hess, chg, modh%s6)
      case(p_modh_lindh_d2)
        if (pr) write(env%unit,'(a)') "Using Lindh-Hessian"
        call mh_lindh_d2(xyz, natoms, Hess, chg, modh)
      case(p_modh_lindh)
        if (pr) write(env%unit,'(a)') "Using Lindh-Hessian (2007)"
        call mh_lindh(xyz, natoms, Hess, chg, modh)
      case(p_modh_swart)
        if (pr) write(env%unit,'(a)') "Using Swart-Hessian"
        call mh_swart(xyz, natoms, Hess, chg, modh)
      end select
   type is(TGFFCalculator) ! GFN-FF case
      select case(modh%model)
      case default
         call env%error("internal error in model hessian!", source)
         return
      case(p_modh_old, p_modh_gff)
         if (pr) write(env%unit,'(a)') "Using GFN-FF Lindh-Hessian"
         call gff_ddvopt(xyz, natoms, Hess, chg, modh%s6, calc%param, calc%topo, calc%neigh)
      case(p_modh_lindh_d2)
        if (pr) write(env%unit,'(a)') "Using Lindh-Hessian"
        call mh_lindh_d2(xyz, natoms, Hess, chg, modh)
      case(p_modh_lindh)
        if (pr) write(env%unit,'(a)') "Using Lindh-Hessian (2007)"
        call mh_lindh(xyz, natoms, Hess, chg, modh)
      case(p_modh_swart)
        if (pr) write(env%unit,'(a)') "Using Swart-Hessian"
        call mh_swart(xyz, natoms, Hess, chg, modh)
      end select
   end select

!  constraints
   call constrhess(natoms,chg,xyz,Hess)

end subroutine modhes

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! H(cart) = Bt * H(int) *B
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine itochess(nvar,nat3,b,hess,fc)
   use xtb_mctc_blas, only : blas_gemm
   implicit none
   integer, intent(in)  :: nvar
   integer, intent(in)  :: nat3
   real(wp),intent(in)  :: b(nat3,nat3)
   real(wp),intent(in)  :: hess(nvar*(nvar+1)/2)
   real(wp),intent(out) :: fc(nat3*(nat3+1)/2)

   integer :: i,j,k
   real(wp),allocatable :: h(:,:)
   real(wp),allocatable :: hb(:,:)
   real(wp),allocatable :: chess(:,:)

   allocate(h(nvar,nvar),hb(nvar,nat3),chess(nat3,nat3))

!  BLAS implementation of the transformation to cartesians
!  unpack BFGS hess
   k=0
   do i=1, nvar
      do j=1,i
         k=k+1
         h(i,j)=hess(k)
         h(j,i)=h(i,j)
      enddo
   enddo

!  calculate H*B (note: B is saved as Bt)
   call blas_gemm ('n','t',nvar,nat3,nvar,1.0d0,h,nvar,b,nat3,0.0d0,hb,nvar)

!  calculate Bt*HB
   call blas_gemm ('n','n',nat3,nat3,nvar,1.0d0,b,nat3,hb,nvar,0.0d0,chess,nat3)

!  pack
   k=0
   do i=1,nat3
      do j=1,i
         k=k+1
         fc(k)=chess(j,i)
      enddo
   enddo

end subroutine itochess

subroutine read_hessian(ihess, ndim, hessian, error)
   use xtb_mctc_systools
   implicit none
   integer, intent(in) :: ihess
   integer, intent(in) :: ndim
   real(wp), intent(out) :: hessian(:, :)
   character(len=:), allocatable :: line
   integer, intent(out) :: error
   integer :: i, j
   error = 0
   do while(error == 0)
      call getline(ihess, line, error)
      if (index(line, '$hessian') == 1) then
         read(ihess, *, iostat=error) &
            & ((hessian(j, i), j = 1, ndim), i = 1, ndim)
         exit
      endif
   enddo
end subroutine read_hessian

subroutine geoconvav(nc, e, g, val, deriv)
   implicit none
   !> total number of E/G points
   integer nc
   !> total energy in Eh
   real*8  e(*)
   !> norm of Cartesian gradient (in TM: |dE/dxyz|)
   real*8  g(*)
   !> av. energy in Eh to be used further
   real*8  val
   !> av. gradient
   real*8  deriv

   integer nav, low
   integer i, j
   parameter (nav=5)    ! average over last nav
   real*8 eav,gav

   ! only apply it if sufficient number of points i.e. a "tail" can exist
   ! with the censo blockl = 8 default, this can first be effective in the second
   if(nc.lt.3*nav)then
      val = e(nc)
      deriv = g(nc)
      return
   endif

   low = max(1,nc-nav+1)
   j = 0
   eav = 0
   do i=nc,low,-1
      j = j + 1
      eav = eav + e(i)
      gav = gav + g(i)
   enddo
   val = eav / float(j)

   low = max(1,nc-nav+1)
   j = 0
   gav = 0
   do i=nc,low,-1
      j = j + 1
      gav = gav + g(i)
   enddo
   ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
   ! DCOSMO-RS gradient for large cases can never (even on average)
   ! become lower than the "-opt normal" thresholds
   deriv=gav / float(j) / 2.d0
end subroutine geoconvav

pure function new_convergence_log(nmax) result(self)
   integer, intent(in) :: nmax
   type(convergence_log) :: self
   self%nlog = 0
   allocate(self%elog(nmax))
   allocate(self%glog(nmax))
end function new_convergence_log

function load_turbomole_log(nmax) result(self)
   use xtb_mctc_resize, only : resize
   use xtb_mctc_systools, only : getline
   integer, intent(in) :: nmax
   type(convergence_log) :: self
   integer :: nlog, ilog
   real(wp) :: etmp, gtmp
   real(wp), allocatable :: elog(:), glog(:)
   character(len=:), allocatable :: line
   logical :: exist
   integer, parameter :: initial_size = 50
   integer :: unit, stat
   nlog = 0
   inquire(file="gradient", exist=exist)
   if (exist) then
      open(file="gradient", newunit=unit)
      allocate(elog(initial_size), glog(initial_size))
      call getline(unit, line, stat)
      do while(stat == 0)
         if (index(line, "cycle") > 0) then
            if (tokenize_cycle(line, etmp, gtmp)) then
               if (nlog >= size(elog)) then
                  call resize(elog, size(elog)+size(elog)/2+1)
                  call resize(glog, size(glog)+size(glog)/2+1)
               end if
               nlog = nlog + 1
               elog(nlog) = etmp
               glog(nlog) = gtmp
            end if
         end if
         call getline(unit, line, stat)
      end do
      close(unit)
   end if
   self = convergence_log(nmax + nlog)
   do ilog = 1, nlog
      call self%set_eg_log(elog(ilog), glog(ilog))
   end do
end function load_turbomole_log

function tokenize_cycle(line, e, g) result(stat)
   character(len=*), intent(in) :: line
   real(wp), intent(out) :: e
   real(wp), intent(out) :: g
   logical :: stat
   integer :: stat1, stat2
   read(line(33:51), *, iostat=stat1) e
   read(line(66:), *, iostat=stat2) g
   stat = stat1 == 0 .and. stat2 == 0
end function tokenize_cycle

pure function get_averaged_energy(self) result(val)
   class(convergence_log), intent(in) :: self
   real(wp) :: eav, val
   integer :: i, j, low
   integer, parameter :: nav = 5

   ! only apply it if sufficient number of points i.e. a "tail" can exist
   ! with the censo blockl = 8 default, this can first be effective in the second
   if(self%nlog.lt.3*nav)then
      val = self%elog(self%nlog)
   else
      low = max(1, self%nlog-nav+1)
      j = 0
      eav = 0
      do i = self%nlog, low, -1
         j = j + 1
         eav = eav + self%elog(i)
      enddo
      val = eav / float(j)
   end if

end function get_averaged_energy

pure function get_averaged_gradient(self) result(deriv)
   class(convergence_log), intent(in) :: self
   real(wp) :: gav, deriv
   integer :: i, j, low
   integer, parameter :: nav = 5

   ! only apply it if sufficient number of points i.e. a "tail" can exist
   ! with the censo blockl = 8 default, this can first be effective in the second
   if(self%nlog.lt.3*nav) then
      deriv = self%glog(self%nlog)
   else
      low = max(1, self%nlog-nav+1)
      j = 0
      gav = 0
      do i = self%nlog, low, -1
         j = j + 1
         gav = gav + self%glog(i)
      enddo
      ! adjust the gradient norm to xtb "conventions" because e.g. a noisy
      ! DCOSMO-RS gradient for large cases can never (even on average)
      ! become lower than the "-opt normal" thresholds
      deriv = gav / float(j) / 2.d0
   end if

end function get_averaged_gradient

pure subroutine set_eg_log(self, e, g)
   use xtb_mctc_resize, only : resize
   class(convergence_log), intent(inout) :: self
   real(wp), intent(in) :: e, g
   if (self%nlog >= size(self%elog)) then
      call resize(self%elog, size(self%elog)+size(self%elog)/2+1)
      call resize(self%glog, size(self%glog)+size(self%glog)/2+1)
   end if
   self%nlog = self%nlog + 1
   self%elog(self%nlog) = e
   self%glog(self%nlog) = g
end subroutine set_eg_log

end module xtb_optimizer
