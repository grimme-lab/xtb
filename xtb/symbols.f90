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

module tbmod_symbols
   implicit none

contains


elemental subroutine symbol_to_number(number, symbol)
   character(len=2), parameter :: pse(118) = [ &
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
      &           'nh','fl','mc','lv','ts','og' ]
   character(len=*), intent(in) :: symbol
   integer, intent(out) :: number
   character(len=2) :: lc_symbol
   integer :: i, j, k, l
   integer, parameter :: offset = iachar('a')-iachar('A')

   number = 0
   lc_symbol = '  '

   k = 0
   do j = 1, len_trim(symbol)
      if (k > 2) exit
      l = iachar(symbol(j:j))
      if (k >= 1 .and. l == iachar(' ')) exit
      if (k >= 1 .and. l == 9) exit
      if (l >= iachar('A') .and. l <= iachar('Z')) l = l + offset
      if (l >= iachar('a') .and. l <= iachar('z')) then
         k = k+1
         lc_symbol(k:k) = achar(l)
      endif
   enddo

   do i = 1, size(pse)
      if (lc_symbol == pse(i)) then
         number = i
         exit
      endif
   enddo

   if (number == 0) then
      select case(lc_symbol)
      case('d ', 't ')
         number = 1
      end select
   end if

end subroutine symbol_to_number

elemental subroutine number_to_symbol(symbol, number)
   integer,intent(in) :: number
   character(len=2), intent(out) :: symbol
   character(len=2), parameter :: pse(118) = [ &
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
   &           'Nh','Fl','Mc','Lv','Ts','Og' ]
   if (number.gt.118) then
      symbol = 'XX'
   else
      symbol = pse(number)
   endif
end subroutine number_to_symbol

elemental function to_number(symbol) result(number)
   character(len=*), intent(in) :: symbol
   integer :: number
   call symbol_to_number(number, symbol)
end function to_number

elemental function to_symbol(number) result(symbol)
   integer,intent(in) :: number
   character(len=2) :: symbol
   call number_to_symbol(symbol, number)
end function to_symbol

end module tbmod_symbols
