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

module screening
   use iso_fortran_env, wp => real64

   use dynamic, only : wrc,boltz,xyzsort2

   implicit none

contains

subroutine screen(mol0,wfn,calc,egap,et,maxiter,epot,grd,sigma)

   use mctc_econv, only : autokcal, aatoau

   use tbdef_molecule
   use tbdef_calculator
   use tbdef_wavefunction
   use tbdef_data

   use setparam

   use axis_trafo, only : axis
   use optimizer, only : wrlog2

   implicit none

   type(tb_molecule), intent(inout) :: mol0
   type(tb_wavefunction),intent(inout) :: wfn
   type(tb_calculator),intent(in) :: calc
   integer icall,maxiter
   real(wp) epot,et,egap
   real(wp), intent(inout) :: grd(3,mol0%n)
   real(wp), intent(inout) :: sigma(3,3)
   type(tb_molecule) :: mol

   real(wp),allocatable :: rot (:,:)
   real(wp),allocatable :: de  (:)
   real(wp),allocatable :: ecnf(:)
   real(wp),allocatable :: xyznew(:,:,:)
   real(wp),allocatable :: rmat(:),eread(:)
   real(wp),allocatable :: ecnfnew(:)
   real(wp),allocatable :: er(:),pp(:)
   integer,allocatable :: imass(:)
   integer,allocatable ::double(:)

   real(wp) :: e,gnorm,dens,emin,ethr,ethr2,gg,ss
   real(wp) :: xx(10),r,rmsd,acc,T,beta,A,dum
   integer :: i,j,k,m,i1,i2,iz1,iz2,lin,nst,ncnf,ndum,olev
   integer :: nall,ntemp,ncnf1,icyc,maxoptiter
   character(len=80) :: atmp,line
   character(len=2) :: asym
   logical :: fail,equalrot2,ohbonded,include_enan,enan,ex
   logical :: checkrmsd
   integer :: ich,ilog

   write(*,*)
   write(*,'(7x,''======================================='')')
   write(*,'(7x,''|              S C R E E N            |'')')
   write(*,'(7x,''|          structural screening       |'')')
   write(*,'(7x,''======================================='')')
   write(*,*)

   inquire(file='solvent',exist=ex)
   if(ex)then ! use QMDFF to generate ensemble
      nall=ntemp_siman*time_md*1000/50  ! number of qmdff write outs
      nall=nall/40 ! take every 40th step from qmdff trj
      allocate(xyznew(3,mol0%n,nall),rot(3,0:nall),de(nall), &
         &         ecnf(0:nall),double(0:nall), &
         &         imass(mol0%n),ecnfnew(0:nall))
      call qmdfftoscreen(mol0%n,mol0%at,mol0%xyz,nall,xyznew)
   else  ! the ensemble exists
      atmp='xtbscreen.xyz'
      call cqpath_read_pathfile_parameter(atmp,iz1,iz2,nall)
      if(iz2.ne.mol0%n) call raise('E','read error in screen',1)
      allocate(xyznew(3,mol0%n,nall),rot(3,0:nall),de(nall), &
         &         ecnf(0:nall),double(0:nall),eread(nall), &
         &         imass(mol0%n),ecnfnew(0:nall))
      call cqpath_read_pathfile(atmp,iz1,iz2,nall,xyznew,mol0%at,eread)
   endif

   T     =temp_md
   beta  =1./(T*8.314510/4.184/1000.+1.d-14)
   ! consider two struc. equal if dE(kcal) < this value
   ethr=0.1
   include_enan=.false.
   xyznew = xyznew*aatoau

   ! the input reference
   mol = mol0
   rot = -1
   ecnf=10000.
   i=0
   ecnf(0)=epot
   call axis(mol%n,mol%at,mol%xyz,rot(1,0),rot(2,0),rot(3,0))
   call wrc('scoord.0',mol%n,mol%xyz,mol%at)

   ncnf =nall
   ncnf1=nall
   write(*,*)ncnf,' start structures'

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! opt all
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

   do icyc=1,3

      maxoptiter=0
      if(icyc.eq.1)then
         maxoptiter=10
         ethr2=ewin_conf*2
         olev=-99
      endif
      if(icyc.eq.2)then
         olev=-1
         ethr2=ewin_conf
      endif
      if(icyc.eq.3)then
         olev=0
         ethr2=ewin_conf
      endif

      i=1
      write(*,*)
      write(*,*) 'optimization blocks ',i, ' opt. level ',olev
      do i=1,ncnf

         mol%xyz(1:3,1:mol%n)=xyznew(1:3,1:mol%n,i)

         call geometry_optimization &
            &       (mol,wfn,calc, &
            &        egap,et,maxiter,maxoptiter,ecnf(i),grd,sigma,optset%optlev, &
            &        .false.,.true.,fail)

         if(.not.fail)then
            !call getname1(i,atmp)
            write(atmp,'(''scoord.'',i0)')i
            write(*,'(i4,''  energy : '',f14.7)') i,ecnf(i)
            call axis(mol%n,mol%at,mol%xyz,rot(1,i),rot(2,i),rot(3,i))
            call wrc(atmp,mol%n,mol%xyz,mol%at)
         endif

      enddo

      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! analyse FF ensemble
      !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      double=0
      do i=0,ncnf
         do j=0,i-1
            e=(ecnf(i)-ecnf(j))*autokcal
            if(check_rmsd) then
               call heavyrmsdfile(i,j,dum)
               if(dum.lt.0.2) double(i)=1
            else
               if(equalrot2(i,j,nall,0.05d0,rot).and.abs(e).lt.ethr) then
                  ! check enantiomers
                  if(include_enan)then
                     call heavyrmsdfile(i,j,dum)
                     if(dum.lt.0.2) double(i)=1
                  else
                     double(i)=1
                  endif
               endif
            endif
         enddo
      enddo

      emin=minval(ecnf(1:ncnf))

      m=0
      do i=1,ncnf
         ! opt failed
         if(rot(1,i).lt.0) cycle

         !call getname1(i,atmp)
         write(atmp,'(''scoord.'',i0)')i
         call open_file(ich,atmp,'r')
         read(ich,'(a)')line
         do j=1,mol%n
            read(ich,*) mol%xyz(1:3,j)
         enddo
         call remove_file(ich)

         dum=autokcal*(ecnf(i)-emin)
         if(double(i).eq.0.and.dum.lt.ethr2)then
            m=m+1
            xyznew(1:3,1:mol%n,m)=mol%xyz(1:3,1:mol%n)
            ecnfnew(m) = ecnf(i)
         endif
      enddo
      write(*,'(i3,'' conformers remaining in'',f8.2,'' kcal window'')') &
         &      m,ethr2

      ! sort by energy
      ncnf=m
      call xyzsort2(ncnf,mol%n,ncnf,ecnfnew,xyznew)
      ecnf(1:ncnf)=ecnfnew(1:ncnf)

   enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! output
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
99 write(*,*)
   write(*,*) '==================================================='
   write(*,*) '============= ordered conformer list =============='
   write(*,*) '==================================================='
   write(*,*) 'written to <scoord.nnn> and xtbscreen.log'
   write(*,*) 'dE=0 corresponds to the input structure'

   allocate(er(ncnf+1),pp(ncnf+1))
   call open_file(ilog,'xtbscreen.log','w')
   call wrlog2(ilog,mol0%n,mol0%xyz,mol0%at,ecnf(0))
   j=0
   emin=1.0d+42
   do i=0,ncnf
      er(i+1)=autokcal*(ecnf(i)-ecnf(0))
      !call getname1(i,atmp)
      write(atmp,'(''scoord.'',i0)')i
      if(i.gt.0) then
         call wrc(atmp,mol%n,xyznew(1,1,i),mol%at)
         call wrlog2(ilog,mol%n,xyznew(1,1,i),mol%at,ecnf(i))
      endif
      write(*,'(i4,'' dE (kcal):'',F10.2,F14.6)') i, er(i+1), ecnf(i)
   enddo
   call close_file(ilog)

   call boltz(ncnf+1,T,er,pp)
   A=0
   do i=1,ncnf+1
      A=A+pp(i)*log(pp(i)+1.d-12)
   enddo
   gg=(1./beta)*A
   ss=-1000.*4.184*gg/T
   write(*,'('' S / J/mol K  :'',F10.3)') ss
   write(*,'('' G / kcal/mol :'',F10.3)') gg

   do i=ncnf+1,ncnf1
      !call getname1(i,atmp)
      write(atmp,'(''scoord.'',i0)')i
      call delete_file(trim(atmp))
   enddo

end subroutine screen

! call QMDFF MD to generate conformer ensemble
subroutine qmdfftoscreen(n,at,xyz,nall,xyzall)
   use iso_c_binding, only : c_null_char
   use gbobc, only: lgbsa
   use setparam
   implicit none
   integer  :: n
   integer  :: at(n)
   integer  :: nall
   real(wp) :: xyzall(3,n,nall)
   real(wp) :: xyz   (3,n)
   interface
      function rename(oldpath, newpath) bind(c) result(stat)
         use iso_c_binding
         character(len=*,kind=c_char), intent(in) :: oldpath, newpath
         integer(c_int) :: stat
      end function rename
   end interface

   character(80) :: atmp
   character(2)  :: a2
   integer  :: i,j,k,l,m,ndum,nstep,dumpstep
   real(wp) :: tstart,temp2,temp,xyztmp(3,n)
   integer  :: ilog

   !call wrcoord(16,n,xyz,at,0.0d0,'coord')  ! qmdff needs a coord file
   call wrc('coord',n,xyz,at)
   write(*,*) 'moving hessian to hessian.save'
   !call execute_command_line('mv hessian hessian.save')
   i = rename('hessian'//c_null_char,'hessian.save'//c_null_char)

   temp  =Tend_siman
   tstart=temp_md
   temp2 =tstart
   dumpstep=40

   k=0
   do i=1,ntemp_siman
      write(atmp,'(''qmdff coord -rd -md '',i4,'' -temp '',F6.1, &
         &   '' > tmp'')') idint(time_md),temp2
      if(lgbsa) &
         &   write(atmp,'(''qmdff coord -rd -md '',i4,'' -temp '',F6.1, &
         &   '' -gbsa '',a20,'' > tmp'')') idint(time_md),temp2,solvent
      write(*,*)'QMDFF call:',trim(atmp)
      call execute_command_line(atmp)                                ! run qmdff
      temp2=temp2+(temp-tstart)/(ntemp_siman-1.)
      write(atmp,'(''cp qmdff.trj qmdff.trj.T.'',i1)') i
      call execute_command_line(atmp)                                ! save qmdff trj
      call open_file(ilog,'qmdff.trj','r')
      m=0
10    m=m+1                                            ! read
      read(ilog,*,end=20) ndum
      read(ilog,'(a)') atmp
      do j=1,ndum
         read(ilog,*,end=20) a2,xyztmp(1:3,j)
      enddo
      if(mod(m,dumpstep).eq.0.or.k.eq.0) then
         k=k+1
         if(k.gt.nall) goto 20
         xyzall(1:3,1:ndum,k)=xyztmp(1:3,1:ndum)
      endif
      goto 10
20    call remove_file(ilog)
   enddo

end subroutine qmdfftoscreen

end module screening

