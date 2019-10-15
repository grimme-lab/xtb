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

subroutine geometry_optimization &
      &   (mol,wfn,calc,egap,et,maxiter,maxcycle_in,etot,g,sigma, &
      &    tight,pr,initial_sp,fail)
   use iso_fortran_env, wp => real64, istdout => output_unit

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_calculator
   use tbdef_data

   use optimizer
   use relaxation_engine
   use single
   use write_geometry

   use setparam

   implicit none

   type(tb_molecule), intent(inout) :: mol
   integer, intent(in)    :: tight
   integer, intent(in)    :: maxiter
   integer, intent(in)    :: maxcycle_in
   type(tb_wavefunction),intent(inout) :: wfn
   type(tb_calculator),  intent(in) :: calc
   real(wp),intent(inout) :: etot
   real(wp),intent(in)    :: et
   real(wp),intent(inout) :: egap
   real(wp),intent(inout) :: g(3,mol%n)
   real(wp),intent(inout) :: sigma(3,3)
   logical, intent(in)    :: pr
   logical, intent(out)   :: fail
   logical, intent(in)    :: initial_sp

   type(scc_results) :: res
   logical :: final_sp
   integer :: printlevel
   integer :: iunit, ilog

   final_sp = pr
   if (pr) then
      call delete_file('NOT_CONVERGED')
      call delete_file('.xtboptok')
      printlevel = 2
   else
      printlevel = 0
   endif

   if (.not.allocated(opt_outfile)) then
      iunit = istdout
   else
      if (opt_outfile == '-') then
         iunit = istdout
      else
         call open_file(iunit,opt_outfile,'w')
         if (pr) write(istdout,'(/,a)') &
            "Optimizer printout bound to '"//opt_outfile//"'"
      endif
   endif

   if (.not.allocated(opt_logfile)) then
      call open_file(ilog,'xtbopt.log','w')
   else
      if (opt_logfile == '-') then
         ilog = istdout
      else
         call open_file(ilog,opt_logfile,'w')
         if (pr) write(istdout,'(/,a)') &
            "Optimization log is written to '"//opt_logfile//"'"
      endif
   endif

   call mol%update

   if (initial_sp) then
      call singlepoint &
         &(iunit,mol,wfn,calc, &
         & egap,et,maxiter,printlevel-1,.false.,.false.,1.0_wp,etot,g,sigma,res)
   endif

   select case(opt_engine)
   case(p_engine_rf)
      call ancopt &
         &(iunit,ilog,mol,wfn,calc, &
         & egap,et,maxiter,maxcycle_in,etot,g,sigma,tight,pr,fail)
      ! required since ANCopt might perform an untracked displacement
      final_sp = .true.
   case(p_engine_lbfgs)
      call l_ancopt &
         &(iunit,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   case(p_engine_inertial)
      call fire &
         &(iunit,ilog,mol,wfn,calc, &
         & tight,maxcycle_in,etot,egap,g,sigma,printlevel,fail)
   end select

   if (fail) then
      call touch_file('NOT_CONVERGED')
   else
      call touch_file('.xtboptok')
   endif

   if (pr.and.pr_finalstruct) then
      write(iunit,'(''================'')')
      write(iunit,*) 'final structure:'
      write(iunit,'(''================'')')
      call write_coord(iunit,mol%n,mol%at,mol%xyz)
   endif
   if (pr.and.pr_geosum) call geosum(mol%n,mol%at,mol%xyz)

   if (final_sp) then
      if (pr) call generic_header(iunit,'Final Singlepoint',49,10)
      call singlepoint &
         &(iunit,mol,wfn,calc, &
         & egap,et,maxiter,printlevel,.false.,.false.,1.0_wp,etot,g,sigma,res)
   endif

   if (iunit.ne.istdout) call close_file(iunit)
   if (ilog .ne.istdout) call close_file(ilog)

end subroutine geometry_optimization
