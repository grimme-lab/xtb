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

subroutine printbas(n,at)
  use aoparam
  implicit none

  integer n,at(n),nn(94),i,j,l
  character(30) atmp
  character(3) lnam(0:13)          

  lnam(0)='s'
  lnam(1)='p'
  lnam(2)='d'
  lnam(3)='f'
  lnam(11)='S'
  lnam(12)='sp'
  lnam(13)='spd'

  nn=0
  do i=1,n
     nn(at(i))=nn(at(i))+1
  enddo

  write(*,*)' Z AO/shell   Hii/eV     exponent'
  do i=1,94
     atmp=timestp(i)
     if(nn(i).eq.0      ) cycle
     if(atmp(1:1).eq.'-') cycle
     write(*,'(i3,5x,A30,''  EN:'',F6.3,'' GAM:'',F6.3,''  GM3:'',F7.4)') &
          & i,timestp(i),en(i),gam(i),gam3(i)
     do j=1,ao_n(i)
        l=ao_l(j,i)
        write(*,'(3x,i3,a3,2F12.6)') &
             & ao_pqn(j,i),lnam(l),ao_lev(j,i),ao_exp(j,i)
     enddo
  enddo
end subroutine printbas


subroutine printbas2(n,at)
  use aoparam
  implicit none

  integer n,at(n),nn(94),i,j,l
  character(30) atmp
  character(3) lnam(0:13)          

  lnam(0)='s'
  lnam(1)='p'
  lnam(2)='d'
  lnam(3)='f'
  lnam(11)='S'
  lnam(12)='sp'
  lnam(13)='spd'

  nn=0
  do i=1,n
     nn(at(i))=nn(at(i))+1
  enddo

  write(*,*)' Z AO/shell   Hii/eV     exponent'
  do i=1,94
     atmp=timestp(i)
     if(nn(i).eq.0      ) cycle
     if(atmp(1:1).eq.'-') cycle
     write(*,'(i3,5x,A30,''  EN:'',F6.3,'' GM2:'',F6.3, &
          & ''  GM3:'',F7.4,''  RAES:'',F5.2)') &
          & i,timestp(i),en(i),gam(i),gam3(i),radaes(i)
     do j=1,ao_n(i)
        l=ao_l(j,i)
        write(*,'(3x,i3,a3,2F12.6)') &
             & ao_pqn(j,i),lnam(l),ao_lev(j,i),ao_exp(j,i)
     enddo
  enddo
end subroutine printbas2
