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

!     *****************************************************************         

pure elemental function asym(i)
   implicit none
   interface
      pure elemental subroutine upper(as)
      character(len=2),intent(inout) :: as
      end subroutine upper
   end interface
   integer,intent(in) :: i
   character(len=2) :: asym
   character(len=2) :: elemnt(118), as
   parameter (elemnt=(/ &
   & 'H ','He', &
   & 'Li','Be','B ','C ','N ','O ','F ','Ne', &
   & 'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
   & 'K ','Ca', &
   & 'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
   &           'Ga','Ge','As','Se','Br','Kr', &
   & 'Rb','Sr', &
   & 'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd', &
   &           'In','Sn','Sb','Te','I ','Xe', &
   & 'Cs','Ba', &
   & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
   & 'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
   &           'Tl','Pb','Bi','Po','At','Rn', &
   & 'Fr','Ra', &
   & 'Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No', &
   & 'Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds','Rg','Cn', &
   &           'Nh','Fl','Mc','Lv','Ts','Og' /))
   if (i.gt.118) then
      asym='XX'
   else
      asym=elemnt(i)
   endif
end function asym

pure elemental subroutine upper(as)
   character(len=2),intent(inout) :: as
   nsp=ichar(' ')
   nd=ichar('A')-ichar('a')
   do i=1,2
      j=ichar(as(i:i))
      if(j.ne.nsp)j=j+nd
      as(i:i)=char(j)
   enddo
end subroutine upper

pure elemental subroutine upper10(as)
   character(len=*),intent(inout) :: as
   nsp=ichar(' ')
   nd=ichar('A')-ichar('a')
   do i=1,10
      j=ichar(as(i:i))
      if(j.ne.nsp)j=j+nd
      as(i:i)=char(j)
   enddo
end subroutine upper10

pure elemental function esym(i)
   integer,intent(in) :: i
   character(len=2) :: esym
   character(len=2),parameter :: elemnt(118) = (/ &
   & 'h ','he', &
   & 'li','be','b ','c ','n ','o ','f ','ne', &
   & 'na','mg','al','si','p ','s ','cl','ar', &
   & 'k ','ca', &
   & 'sc','ti','v ','cr','mn','fe','co','ni','cu','zn', &
   &           'ga','ge','as','se','br','kr', &
   & 'rb','sr', &
   & 'y ','zr','nb','mo','tc','ru','rh','pd','ag','cd', &
   &           'in','sn','sb','te','i ','xe', &
   & 'cs','ba','la', &
   & 'ce','pr','nd','pm','sm','eu','gd','tb','dy','ho','er','tm','yb', &
   & 'lu','hf','ta','w ','re','os','ir','pt','au','hg', &
   &           'tl','pb','bi','po','at','rn', &
   & 'fr','ra','ac', &
   & 'th','pa','u ','np','pu','am','cm','bk','cf','es','fm','md','no', &
   & 'lr','rf','db','sg','bh','hs','mt','ds','rg','cn', &
   &           'nh','fl','mc','lv','ts','og' /)
   esym=elemnt(i)
end function esym
