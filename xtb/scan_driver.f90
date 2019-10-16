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

subroutine relaxed_scan(mol,wfx,calc)
   use iso_fortran_env, wp => real64, id => output_unit

   use setparam
   use scanparam

   use tbdef_molecule
   use tbdef_wavefunction
   use tbdef_calculator
   use tbdef_data

   use optimizer, only : wrlog2

   implicit none

   type(tb_molecule), intent(inout) :: mol
   type(tb_calculator),intent(in) :: calc
   type(tb_wavefunction),intent(inout) :: wfx

   integer  :: ilog ! file handle
   real(wp) :: egap
   integer  :: maxiter
   integer  :: maxcycle
   real(wp) :: etot,efix
   real(wp) :: sigma(3,3)
   real(wp),allocatable :: g(:,:)
   real(wp),allocatable :: xyzsave(:,:)
   integer  :: optlevel
   logical  :: pr,fail,reset
   logical  :: exist
   integer  :: i,j,k

   pr = verbose
   reset = .false.
   optlevel = optset%optlev
   maxcycle = optset%maxoptcycle
   maxiter = maxscciter

   allocate ( g(3,mol%n), source = 0.0_wp )
   allocate ( xyzsave(3,mol%n), source = mol%xyz )

   write(id,'(72("="))')
   write(id,'(1x,"RELAXED SCAN")')
   write(id,'(72("-"))')
   write(id,'(1x,"output written to",1x,a)') get_namespace('xtbscan.log')
   call open_file(ilog,'xtbscan.log','w')

!! ========================================================================
!  SEQUENTIAL SCAN
!  -> scan over all constrains one after another
!  -> you can scan a constrain twice, we don't care
   if (scan_mode.eq.p_scan_sequential) then
      k = 0
      do i = 1, nscan
         write(id,'(1x,"scaning constraint",1x,i0,1x,i0)') scan_list(i)%iconstr
         do j = 1, scan_list(i)%nscan
            k = k+1
!           write(id,'(1x,"valconstr:",1x,f12.8)') scan_list%valscan(j)
            valconstr(scan_list(i)%iconstr) = scan_list(i)%valscan(j)
!           write(id,'(i0,1x,i0)') i,j
            if (.not.verbose) &
               write(id,'("... step",1x,i0,1x,"...")') k
            call geometry_optimization &
               &(mol,wfx,calc, &
               & egap,etemp,maxiter,maxcycle,etot,g,sigma,optlevel,pr,.true.,fail)
            efix = 0.0_wp
            call constrpot(mol%n,mol%at,mol%xyz,g,efix)
            if (.not.verbose) then
               write(id,'(" current energy:",1x,f20.8)') etot
               write(id,'("    bias energy:",1x,f20.8)') efix
               write(id,'("unbiased energy:",1x,f20.8)') etot-efix
            endif
            call wrlog2(ilog,mol%n,mol%xyz,mol%at,etot-efix)
         enddo
         if (reset) valconstr(scan_list(i)%iconstr) = scan_list(i)%valconstr
      enddo
!! ========================================================================
!  CONCERTED SCAN
!  -> we scan all constraints at once, therefore we require an equal number
!     of steps for each scan
!  -> you can still scan a constrain twice, which then will lead to
!     abdominations and horrors beyond your imagination, beware.
   else if (scan_mode.eq.p_scan_concerted) then
      if (any(scan_list(1:nscan)%nscan.ne.scan_list(1)%nscan)) then
         ! I'm not happy doing this here, since it should be catched already
         call raise('E','Wrong setup for concerted scan, aborting...',1)
      endif
      do j = 1, scan_list(1)%nscan
         do i = 1, nscan
            valconstr(scan_list(i)%iconstr) = scan_list(i)%valscan(j)
         enddo
         if (.not.verbose) &
            write(id,'("... step",1x,i0,1x,"...")')   j
         call geometry_optimization &
            &(mol,wfx,calc, &
            & egap,etemp,maxiter,maxcycle,etot,g,sigma,optlevel,pr,.true.,fail)
         efix = 0.0_wp
         call constrpot(mol%n,mol%at,mol%xyz,g,efix)
         if (.not.verbose) then
            write(id,'(" current energy:",1x,f20.8)') etot
            write(id,'("    bias energy:",1x,f20.8)') efix
            write(id,'("unbiased energy:",1x,f20.8)') etot-efix
         endif
         call wrlog2(ilog,mol%n,mol%xyz,mol%at,etot-efix)
      enddo
!! ========================================================================
!  BUGGY SCAN
!  if you land here, you used it wrong. Fix it! Now!
   else
      call raise('E','We screwed up in scan_driver.f90, blame us in a bug report',1)
   endif

   call close_file(ilog)
   
end subroutine relaxed_scan
