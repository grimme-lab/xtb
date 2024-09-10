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

module xtb_geoopt
contains

!> wrapper for relaxation engines
!> more info: https://xtb-docs.readthedocs.io/en/latest/optimization.html 
subroutine geometry_optimization &
      &   (env,mol,wfn,calc,egap,et,maxiter,maxcycle_in,etot,g,sigma, &
      &    tight,pr,initial_sp,fail)

   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_io, only : stdout
   use xtb_type_molecule
   use xtb_type_restart
   use xtb_type_calculator
   use xtb_type_data
   use xtb_optimizer
   use xtb_relaxation_engine
   use xtb_single
   ! modules for PBC opt with LBFGS including lattice optimization
   use xtb_pbc_optimizer_driver, only : relax_pbc
   use xtb_pbc_optimizer_lbfgs, only : optimizer_type, lbfgs_optimizer, new_lbfgs_optimizer
   use xtb_pbc_optimizer_lbfgs, only : lbfgs_input
   use xtb_pbc_optimizer_filter_cart, only : cartesian_filter, new_cartesian_filter
   use xtb_setparam

   implicit none
   
   character(len=*), parameter :: source = 'xtb_geoopt'
   
   !> instace of polymorphic calculatior
   class(TCalculator), intent(inout) :: calc
   
   !> calculation environment
   type(TEnvironment), intent(inout) :: env
   
   !> molecular data information
   type(TMolecule), intent(inout) :: mol

   !> wavefunction
   type(TRestart),intent(inout) :: wfn
   
   !> optimization level
   integer, intent(in)    :: tight

   !> max number of SCC cycles
   integer, intent(in)    :: maxiter
   
   !> max number of optimization cycles
   integer, intent(in)    :: maxcycle_in
   
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
   
   !> internal printlevel 
   logical, intent(in)    :: pr
   
   !> optimization convergence
   logical, intent(out)   :: fail
   
   !> perform initial single point 
   logical, intent(in)    :: initial_sp

   type(scc_results) :: res
   logical :: final_sp, exitRun
   integer :: printlevel
   integer :: ilog
   ! for PBC optimization with lbfgs (including lattice optimization)
   logical :: optcell
   type(lbfgs_optimizer) :: lbfgs_opt
   type(lbfgs_input) :: opt_input
   type(cartesian_filter)  :: filter

!----------------!
! Initialization !
!----------------!

   final_sp = pr

   if (pr) then
      call delete_file('NOT_CONVERGED')
      call delete_file('.xtboptok')
      printlevel = 2
   else
      printlevel = 0
   endif
   
   ! create/open optimization log !
   if (.not.allocated(set%opt_logfile)) then
      call open_file(ilog,'xtbopt.log','w')
   else
      if (set%opt_logfile == '-') then
         ilog = stdout
      else
         call open_file(ilog,set%opt_logfile,'w')
         if (pr) write(stdout,'(/,a)') &
            "Optimization log is written to '"//set%opt_logfile//"'"
      endif
   endif
   
   ! create ONIOM optimization logs for low & high methods ! 
   if (set%oniom_settings%logs) then
      call open_file(set%oniom_settings%ilog1,"low.inner_region.log",'w')
      call open_file(set%oniom_settings%ilog2,"high.inner_region.log",'w')
   endif

   call mol%update ! update interatomic and minimum image (periodic) distances 

!-------------!
! Calculation !
!-------------!

   if (initial_sp)  call singlepoint &
         &(env,mol,wfn,calc, &
         & egap,et,maxiter,printlevel-1,.false.,.false.,1.0_wp,etot,g,sigma,res)

   select case(set%opt_engine)
   case(p_engine_rf) ! ANCopt !
      call ancopt &
         &(env,ilog,mol,wfn,calc, &
         & egap,et,maxiter,maxcycle_in,etot,g,sigma,tight,pr,fail)
      final_sp = .true. ! required since ANCopt might perform an untracked displacement 
   
   case(p_engine_lbfgs)
      call l_ancopt & ! L-ANCopt !
         &(env,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   case(p_engine_pbc_lbfgs)
     ! get number of unique species; used in precond_lindh
     !call get_identity(mol%nid, mol%id, mol%at)
     ! create new anc filter
     optcell = mol%npbc > 0 .and. set%optcell
     call new_cartesian_filter(filter, mol, optcell)
     ! create new Limited-memory BFGS optimizer 
     call new_lbfgs_optimizer(lbfgs_opt, env, opt_input, filter)
     ! run optimization
     call relax_pbc(lbfgs_opt, env, mol, wfn, calc, filter, printlevel)
   case(p_engine_inertial)
      call fire & ! FIRE !
         &(env,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   
   end select
  
!-----------------!
! Post-processing !
!-----------------!

   ! check convergence !
   if (pr) then
      if (fail) then
         call touch_file('NOT_CONVERGED')
         call env%warning("Geometry optimization did not converge", source)
      else
         call touch_file('.xtboptok')
      endif
   end if

   call env%check(exitRun) ! error check 

   if (exitRun) then
      call env%rescue("Trying to recover from failed geometry optimization", source)
      fail = .true.
   end if
   
   ! print fine details of optimized geometry !
   if (pr.and.set%pr_finalstruct) then
      write(env%unit,'(''================'')')
      write(env%unit,*) 'final structure:'
      write(env%unit,'(''================'')')
      call writeMolecule(mol, env%unit)
   endif
   if (pr.and.set%pr_geosum) call geosum(mol%n,mol%at,mol%xyz)
   
   ! final SP (usually for all engines) !
   if (final_sp) then
      if (pr) call generic_header(env%unit,'Final Singlepoint',49,10)
      call singlepoint &
         &(env,mol,wfn,calc, &
         & egap,et,maxiter,printlevel,.false.,.false.,1.0_wp,etot,g,sigma,res)
   endif
   
   ! close log file !
   if (ilog.ne.stdout) call close_file(ilog)
   
   ! close ONIOM optlogs !
   if (set%oniom_settings%logs) then
      call close_file(set%oniom_settings%ilog1)
      call close_file(set%oniom_settings%ilog2)
   endif

end subroutine geometry_optimization

end module
