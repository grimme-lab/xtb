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

subroutine qpotini(natoms,iat,xyz)
   use iso_fortran_env, only : wp => real64
   use setparam
   use fixparam
   implicit none
   integer natoms,iat(natoms)
   real(wp) xyz(3,natoms)

   integer i,j,k,nread,at(natoms)
   real(wp) xyzf(3,natoms),rij(3)
   character(80) fname,a80
   logical ex

   if(fixset%n.eq.0)then
      write(*,*)'no old-style atom pair restraining (fixing) potential'
   else
      !       calculate distance matrix for fixed atoms
      k=0
      do i=1,fixset%n
         do j=1,i-1
            k=k+1
            rij=xyz(:,fixset%atoms(j))-xyz(:,fixset%atoms(i))
            fixset%val(k)=sqrt(sum(rij*rij))
         enddo
      enddo
      write(*,'(a,i5,a,F10.4)')'# fixed atoms:',fixset%n,' FC:',fixset%fc
      fixset%fc=fixset%fc/(fixset%n-1)
      write(*,'(a,F10.4)')'Final FC :',fixset%fc
      write(*,'(a,F10.0)')'exponent :',fixset%expo
      write(*,'(''fixed atoms'',15I5)')fixset%atoms(1:fixset%n)
   endif

end subroutine qpotini

subroutine countfix(nat,nfix,fname)
   use iso_fortran_env, only : wp => real64
   use mctc_strings, only : lowercase
   implicit none
   real(wp) xx(10)
   integer n,nn,j,nat,k,nfix
   character(128) line
   character(40) fname
   logical fixed(nat),f
   integer :: ich

   call open_file(ich,fname,'r')
   n=0
   fixed=.false.
   100 read(ich,'(a)',end=300)line
   if(index(line,'$user').ne.0)goto 200
   if(index(line,'$red' ).ne.0)goto 200
   if(index(line,'$end' ).ne.0)goto 200
   call readl(line,xx,nn)
   if(nn.ne.3) goto 100
   n=n+1
   line = lowercase(line)
   call getf(line,f)
   if(f)then
      fixed(n)=.true.
   endif
   goto 100
   200 continue
   nfix=count(fixed)
   call close_file(ich)
   return
   300 continue
   call raise('e','internal error in countfix',1)
end subroutine countfix


subroutine rdfix(nat,n,xyz,iat,fname,fixed)
   use iso_fortran_env, only : wp => real64
   use mctc_strings, only : lowercase
   implicit none
   real(wp) xyz(3,*), xx(10)
   integer nat,iat(nat),n,nn,j,k,fixed(*)
   character(128) line
   character(40) fname
   character(2) asym
   logical f
   integer :: ich

   write(*,'(''Reading '',A,'' ...'')') trim(fname)
   call open_file(ich,fname,'r')
   n=0
   k=0
10 read(ich,'(a)',end=20)line
   if(index(line,'$user').ne.0)goto 30
   if(index(line,'$red' ).ne.0)goto 30
   if(index(line,'$end' ).ne.0)goto 30
   call readl(line,xx,nn)
   if(nn.ne.3) goto 10
   n=n+1
   xyz(1,n)=xx(1)
   xyz(2,n)=xx(2)
   xyz(3,n)=xx(3)
   line = lowercase(line)
   call elem(line,j)
   call getf(line,f)
   if(f)then
      k=k+1
      fixed(k)=n  !save the array position of fixed atoms
   endif
   !        if(j.eq.0) error stop 'error reading element symbol'
   !        iat(n)=j
   goto 10
20 continue

   call close_file(ich)
   call raise('e','internal error in rdfix',1)
30 call close_file(ich)
   return
end subroutine rdfix

subroutine getf(line,fix)
   implicit none
   character(*) line
   CHARACTER(3) e
   logical fix
   integer i,j,k,l,n
   e='   '
   k=1
   fix=.false.
   DO J=1,len(line)
      if (k.gt.2)exit
      N=ICHAR(line(J:J))
      if(len_trim(e).ge.1 .and. n.eq.ichar(' '))exit!break if space after elem. symbol
      if(len_trim(e).ge.1 .and. n.eq.9)exit !break if tab after elem. symbol
      if(n.ge.ichar('a') .and. n.le.ichar('z') )then
         e(k:k)=line(j:j)
         k=k+1
      endif
   enddo
   if(j.eq.len(line))then
      return
   else
      do i=j,len(line)
         n=ichar(line(i:i))
         if (n.eq.ichar('f'))fix=.true.
      enddo
   endif
end subroutine getf

! special routine for the Ln fit for Xrays
subroutine fixmetal(n,iat,xyz)
   use iso_fortran_env, only : wp => real64
   use fixparam
   implicit none
   integer n,iat(n),nb(0:20,n),nn,ifix(n)
   real(wp) xyz(3,n)
   integer i,j,k

   call neighbor(n,xyz,iat,nb)

   ifix=1
   do i=1,n
      if(iat(i).gt.57.and.iat(i).lt.72)then
         nn=nb(20,i)
         ifix(i)=0
         do j=1,nn
            k=nb(j,i)
            ifix(k)=0
         enddo
      endif
   enddo

   fixset%n=0
   do i=1,n
      if(ifix(i).eq.1)then
         fixset%n=fixset%n+1
         fixset%atoms(fixset%n)=i
      endif
   enddo

end subroutine fixmetal


subroutine neighbor(natoms,xyz,iz,nb)
   use iso_fortran_env, only : wp => real64
   implicit none
   integer iz(*),natoms,nb(0:20,*)
   real(wp)  xyz(3,*)

   logical da
   integer iat,i,j,k,nn,ni
   real(wp) dx,dy,dz,r,damp,xn,rr,rco,r2,f,pi,a1
   real(wp) rad(94)
   data rad / &
      & 0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,0.64D0,0.60D0, &
      & 0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,1.00D0,1.01D0, &
      & 2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,1.24D0,1.18D0, &
      & 1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,1.17D0,1.16D0, &
      & 2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,1.36D0,1.34D0, &
      & 1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,1.36D0,1.36D0, &
      & 2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,1.85D0,1.83D0, &
      & 1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,1.74D0,1.64D0, &
      & 1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,1.32D0,1.44D0, &
      & 1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,2.01D0,1.90D0, &
      & 1.84D0,1.83D0,1.80D0,1.80D0 /

   do i=1,natoms
      f=1.3
      k=0
      100    do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r2=dx*dx+dy*dy+dz*dz
            r=sqrt(r2)*0.52917726
            rco=rad(iz(i))+rad(iz(iat))
            if(r.lt.f*rco)then
               k=k+1
               nb(k,i)=iat
            endif
         endif
      enddo
      if(k.lt.1.and.f.lt.1.5)then
         f=f*1.1
         goto 100
      endif
      nb(20,i)=k
   enddo

end subroutine neighbor

