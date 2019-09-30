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

subroutine gfn1_prparam(iunit,n,at,par)
   use iso_fortran_env, wp => real64
   use tbdef_param
   use aoparam
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   type(scc_parameter),intent(in) :: par

   integer nn(94),i,j,l
   character(30) atmp
   character(3) lnam(0:13)          

   write(iunit,'(13x,''KAB for M(3d)-M(3d) :'',f22.4)') 1.1
   write(iunit,'(13x,''KAB for M(4d)-M(4d) :'',f22.4)') 1.2
   write(iunit,'(13x,''KAB for M(5d)-M(5d) :'',f22.4)') 1.2

   write(iunit,'(13x,a,3x,":",10x,F12.4)') &
    "k(s)             ",par%kspd(1),&
    "k(p)             ",par%kspd(2),&
    "k(d)             ",par%kspd(3),&
    "k(f)             ",par%kspd(4),&
    "kEN (H0ij)       ",par%kenscal,&
    "D3 a1            ",par%disp%a1,&
    "D3 a2            ",par%disp%a2,&
    "D3 s6            ",1.0_wp,&
    "D3 s8            ",par%disp%s8,&
    "D3 s9            ",0.0_wp,&
    "alphaj           ",par%alphaj,&
    "XBdamp           ",par%xbdamp,&
    "XBrad            ",par%xbrad
   write(iunit,'(13x,a,3x,":",2x,4F10.4)') &
    "kcnsh            ",par%kcnsh(1:4)/0.01_wp

   write(iunit,'(a)')

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

   write(iunit,*)' Z AO/shell   Hii/eV     exponent'
   do i=1,94
      atmp=timestp(i)
      if(nn(i).eq.0      ) cycle
      if(atmp(1:1).eq.'-') cycle
      write(iunit,'(i3,5x,A30,"  EN:",F6.3," GAM:",F6.3,"  GM3:",F7.4)') &
         & i,timestp(i),en(i),gam(i),gam3(i)
      do j=1,ao_n(i)
         l=ao_l(j,i)
         write(iunit,'(3x,i3,a3,2F12.6)') &
            & ao_pqn(j,i),lnam(l),ao_lev(j,i),ao_exp(j,i)
      enddo
   enddo

end subroutine gfn1_prparam

subroutine gfn2_prparam(iunit,n,at,par)
   use iso_fortran_env, kdp => real64
   use tbdef_param
   use aoparam
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   type(scc_parameter),intent(in) :: par

   integer nn(94),i,j,l
   character(30) atmp
   character(3) lnam(0:13)          

   write(iunit,'(13x,a,":",10x,F12.4)') &
    "k(s)              ",par%kspd(1),&
    "k(p)              ",par%kspd(2),&
    "k(d)              ",par%kspd(3),&
    "k(f)              ",par%kspd(4),&
    "kEN (H0ij)        ",par%kenscal,&
    "D4 a1             ",par%disp%a1,&
    "D4 a2             ",par%disp%a2,&
    "D4 s6             ",par%disp%s6,&
    "D4 s8             ",par%disp%s8,&
    "D4 s9             ",par%disp%s9,&
    "alphaj            ",2.0_kdp
   write(iunit,'(a)')

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

   write(iunit,*)' Z AO/shell   Hii/eV     exponent'
   do i=1,94
      atmp=timestp(i)
      if(nn(i).eq.0      ) cycle
      if(atmp(1:1).eq.'-') cycle
      write(iunit,'(i3,5x,A30,''  EN:'',F6.3,'' GM2:'',F6.3, &
         & ''  GM3:'',F7.4,''  RAES:'',F5.2)') &
         & i,timestp(i),en(i),gam(i),gam3(i),radaes(i)
      do j=1,ao_n(i)
         l=ao_l(j,i)
         write(iunit,'(3x,i3,a3,2F12.6)') &
            & ao_pqn(j,i),lnam(l),ao_lev(j,i),ao_exp(j,i)
      enddo
   enddo

end subroutine gfn2_prparam

subroutine gfn0_prparam(iunit,n,at,par)
   use iso_fortran_env, kdp => real64
   use tbdef_param
   use aoparam
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   type(scc_parameter),intent(in) :: par

   integer nn(94),i,j,l
   character(30) atmp
   character(3) lnam(0:13)          

   write(iunit,'(13x,''KAB for M(3d)-M(3d)    :'',f8.4)') 1.1
   write(iunit,'(13x,''KAB for M(4d)-M(4d)    :'',f8.4)') 1.1
   write(iunit,'(13x,''KAB for M(5d)-M(5d)    :'',f8.4)') 1.1
   write(iunit,'(13x,''KAB for M-M (group XI) :'',f8.4)') 0.9

   write(iunit,'(13x,a,":",10x,F12.4)') &
    "k(s)             ",par%kspd(1),&
    "k(p)             ",par%kspd(2),&
    "k(d)             ",par%kspd(3),&
    "k(f)             ",par%kspd(4),&
    "kEN (H0ij)       ",par%kenscal,&
    "D3 a1            ",par%disp%a1,&
    "D3 a2            ",par%disp%a2,&
    "D3 s6            ",par%disp%s6,&
    "D3 s8            ",par%disp%s8,&
    "D3 s9            ",par%disp%s9,&
    "alphaj           ",par%alphaj,&
    "ken⁴             ",par%xbdamp,&
    "rep dEN          ",par%xbrad,&
    "k(dz)            ",par%kcnsh(1),&
    "SRB shift        ",par%kcnsh(2),&
    "SRB prefact.     ",par%kcnsh(3),&
    "SRB slope        ",par%kcnsh(4)
   write(iunit,'(a)')

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

   write(iunit,*)' Z AO/shell   Hii/eV     exponent'
   do i=1,94
      atmp=timestp(i)
      if(nn(i).eq.0      ) cycle
      if(atmp(1:1).eq.'-') cycle
      write(iunit,'(i3,5x,A30,''  EN:'',F6.3,'' GAM:'',F6.3,''  GM3:'',F7.4)') &
         & i,timestp(i),en(i),gam(i),gam3(i)
      do j=1,ao_n(i)
         l=ao_l(j,i)
         write(iunit,'(3x,i3,a3,2F12.6)') &
            & ao_pqn(j,i),lnam(l),ao_lev(j,i),ao_exp(j,i)
      enddo
   enddo

end subroutine gfn0_prparam
