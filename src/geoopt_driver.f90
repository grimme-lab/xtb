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

   use xtb_setparam

   implicit none

   character(len=*), parameter :: source = 'xtb_geoopt'

   !> Calculation environment
   type(TEnvironment), intent(inout) :: env

   type(TMolecule), intent(inout) :: mol
   integer, intent(in)    :: tight
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: maxcycle_in
   type(TRestart),intent(inout) :: wfn
   class(TCalculator), intent(inout) :: calc
   real(wp),intent(inout) :: etot
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: egap
   real(wp),intent(inout) :: g(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)
   logical, intent(in)    :: pr
   logical, intent(out)   :: fail
   logical, intent(in)    :: initial_sp

   type(scc_results) :: res
   logical :: final_sp, exitRun
   integer :: printlevel
   integer :: ilog

   final_sp = pr
   if (pr) then
      call delete_file('NOT_CONVERGED')
      call delete_file('.xtboptok')
      printlevel = 2
   else
      printlevel = 0
   endif

   if (.not.allocated(opt_logfile)) then
      call open_file(ilog,'xtbopt.log','w')
   else
      if (opt_logfile == '-') then
         ilog = stdout
      else
         call open_file(ilog,opt_logfile,'w')
         if (pr) write(stdout,'(/,a)') &
            "Optimization log is written to '"//opt_logfile//"'"
      endif
   endif

   call mol%update

   if (initial_sp) then
      call singlepoint &
         &(env,mol,wfn,calc, &
         & egap,et,maxiter,printlevel-1,.false.,.false.,1.0_wp,etot,g,sigma,res)
   endif

   select case(opt_engine)
   case(p_engine_rf)
      call ancopt &
         &(env,ilog,mol,wfn,calc, &
         & egap,et,maxiter,maxcycle_in,etot,g,sigma,tight,pr,fail)
      ! required since ANCopt might perform an untracked displacement
      final_sp = .true.
   case(p_engine_lbfgs)
      call l_ancopt &
         &(env,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   case(p_engine_inertial)
      call fire &
         &(env,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   end select

   if (pr) then
      if (fail) then
         call touch_file('NOT_CONVERGED')
         call env%warning("Geometry optimization did not converge", source)
      else
         call touch_file('.xtboptok')
      endif
   end if

   call env%check(exitRun)
   if (exitRun) then
      call env%rescue("Trying to recover from failed geometry optimization", source)
      fail = .true.
   end if

   if (pr.and.pr_finalstruct) then
      write(env%unit,'(''================'')')
      write(env%unit,*) 'final structure:'
      write(env%unit,'(''================'')')
      call writeMolecule(mol, env%unit)
   endif
   if (pr.and.pr_geosum) call geosum(mol%n,mol%at,mol%xyz)

   if (final_sp) then
      if (pr) call generic_header(env%unit,'Final Singlepoint',49,10)
      call singlepoint &
         &(env,mol,wfn,calc, &
         & egap,et,maxiter,printlevel,.false.,.false.,1.0_wp,etot,g,sigma,res)
   endif

   if (ilog .ne.stdout) call close_file(ilog)

end subroutine geometry_optimization

end module
