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

      subroutine zmatpr(nat,at,geo,na,nb,nc,molnum)
      use iso_fortran_env, wp => real64, istdout => output_unit
      use mctc_constants
      implicit none
      integer, intent(in) :: nat
      integer, intent(in) :: at(nat)
      integer, intent(in) :: na(nat),nb(nat),nc(nat)
      real(wp),intent(in) :: geo(3,nat)
      character(len=2),external :: asym
      character(len=20) :: filename
      logical  :: ex
      integer  :: i,l,m,n
      integer, intent(in) :: molnum
      integer  :: ich ! file handle
      real(wp) :: bl,ang,dihed

      do i=1,nat
         l=1
         m=1
         n=1
         if(i.eq.1)then
            l=0
            m=0
            n=0
         endif
         if(i.eq.2)then
            m=0
            n=0
         endif
         if(i.eq.3)n=0
         bl=geo(1,i)
         ang=geo(2,i)*180./pi
         dihed=geo(3,i)*180./pi
         if(dihed.gt.180.0d0)dihed=dihed-360.0d0
         write(istdout,'(i4,2x,a2,f12.6,2x,f10.4,2x,f10.4,i6,2i5)')
     .   i,asym(at(i)),bl,ang,dihed,na(i),nb(i),nc(i)
      enddo

      write(filename,'("zmatrix",i0,".zmat")') molnum
      call open_file(ich,trim(filename),'w')

      write(ich,'(a2)') asym(at(1))
      write(ich,'(a2,x,i0,x,f8.3)') asym(at(2)), na(2), geo(1,2)
      write(ich,'(a2,x,i0,x,f8.3,x,i0,x,f8.3)') asym(at(3)), na(3)
     .                      ,geo(1,3),nb(3), geo(2,3)*180./pi

      do i=4,nat
         bl=geo(1,i)
         ang=geo(2,i)*180./pi
         dihed=geo(3,i)*180./pi
         if(dihed.gt.180.0d0)dihed=dihed-360.0d0
         write(ich,'(a2,x,i0,x,f8.3,x,i0,x,f8.3,x,i0,x,f8.3)')
     .                asym(at(i)),na(i),bl,nb(i),ang,nc(i),dihed
      enddo
      write(ich,*)
      close(ich)
      end subroutine zmatpr
