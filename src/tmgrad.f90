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

subroutine tmgrad(nat,iat,xyz,g,edum)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_symbols, only : toLcSymbol
   use xtb_setparam, only : get_namespace
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

   ! write file energy
   a1='$energy'
   j=0
   efile = get_namespace('energy')
   etemp = get_namespace('energy.tmp')
   open(newunit=ien,file=efile)
   open(newunit=jen,file=etemp)
   read(ien,'(a)')a1
   write(jen,'(a)')a1
   read(ien,'(a)')a
   write(jen,'(i6,3F18.11)')j+1,edum,edum,edum
   a='$end'
   write(jen,'(a)')a
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
   read(igr,'(a)')a1
   write(jgr,'(a)')a1
   read(igr,'(a)')a
   write(jgr,'(''  cycle = '',i6,4x,''SCF energy ='',F18.11,3x, &
   &           ''|dE/dxyz| ='',F10.6)')j+1,edum,gsum
   do i=1,nat
      write(jgr,'(3(F20.14,2x),4x,a2)')xyz(1,i),xyz(2,i),xyz(3,i), &
      &                                toLcSymbol(iat(i))
   enddo
   do i=1,nat
      write(jgr,'(3D22.13)')g(1,i),g(2,i),g(3,i)
   enddo
   a='$end'
   write(jgr,'(a)')a
   close(igr)
   close(jgr)
   call execute_command_line('mv '//gtemp//' '//gfile)

end subroutine tmgrad

