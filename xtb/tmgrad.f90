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

subroutine tmgrad(nat,iat,xyz,g,edum)
   use iso_fortran_env, wp => real64
   use setparam, only : get_namespace
   implicit none
   integer  :: nat
   integer  :: iat(nat)
   real(wp) :: xyz(3,nat)
   real(wp) :: g  (3,nat)
   real(wp) :: xx(10),edum,gsum
   integer  :: i,j,k,nn,stat
   integer  :: ien,jen,igr,jgr
   character(len=:),allocatable :: efile,etemp,gfile,gtemp
   character(len=80) :: a
   character(len=80) :: a1
   character(len=2),external :: esym

   ! write file energy
   a1='$energy'
   j=0
   efile = get_namespace('energy')
   etemp = get_namespace('energy.tmp')
   open(newunit=ien,file=efile)
   open(newunit=jen,file=etemp)
   read(ien,'(a)',end=33)a1
33 continue
   write(jen,'(a)')a1
20 read(ien,'(a)',end=44)a
   if(index(a,'$').eq.0)then
      write(jen,'(a)')a
      call readl(a,xx,nn)
      j=idint(xx(1))
   else
   44 write(jen,'(i6,3F18.11)')j+1,edum,edum,edum
      a='$end'
      write(jen,'(a)')a
      goto 22
   endif
   goto 20
22 continue
   close(ien)
   close(jen)
   call execute_command_line('mv '//etemp//' '//efile)

   gsum=sqrt(sum(g**2))

   ! write file gradient
   a1='$grad'
   gfile = get_namespace('gradient')
   gtemp = get_namespace('gradient.tmp')
   open(newunit=igr,file=gfile)
   open(newunit=jgr,file=gtemp)
   read(igr,'(a)',end=31)a1
31 continue
   write(jgr,'(a)')a1
21 read(igr,'(a)',end=41)a
   if(index(a,'$').eq.0)then
      write(jgr,'(a)')a
      call readl(a,xx,nn)
      if(index(a,'cycle').ne.0)j=idint(xx(1))
   else
   41 write(jgr,'(''  cycle = '',i6,4x,''SCF energy ='',F18.11,3x, &
      &           ''|dE/dxyz| ='',F10.6)')j+1,edum,gsum     
      do i=1,nat
         write(jgr,'(3(F20.14,2x),4x,a2)')xyz(1,i),xyz(2,i),xyz(3,i), &
         &                                esym(iat(i))
      enddo
      do i=1,nat
         write(jgr,'(3D22.13)')g(1,i),g(2,i),g(3,i)
      enddo
      a='$end'
      write(jgr,'(a)')a
      goto 24
   endif
   goto 21
24 continue
   close(igr)
   close(jgr)
   call execute_command_line('mv '//gtemp//' '//gfile)

end subroutine tmgrad

