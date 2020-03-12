! This file is part of xtb.
!
! Copyright (C) 2019-2020 Sebastian Ehlert
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

!> Handle conversion between element symbols and atomic numbers
module xtb_mctc_symbols
   use xtb_mctc_resize, only : resize
   implicit none
   private

   public :: symbolLength
   public :: symbolToNumber, numberToSymbol, numberToLcSymbol
   public :: toNumber, toSymbol, toLcSymbol, getIdentity


   !> Get chemical identity
   interface getIdentity
      module procedure :: getIdentityNumber
      module procedure :: getIdentitySymbol
   end interface getIdentity


   !> Maximum allowed length of element symbols
   integer, parameter :: symbolLength = 4


   !> Periodic system of elements
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


   !> Lower case version of the periodic system of elements
   character(len=2), parameter :: lcPse(118) = [ &
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


   !> ASCII offset between lowercase and uppercase letters
   integer, parameter :: offset = iachar('a') - iachar('A')


contains


!> Convert element symbol to atomic number
elemental subroutine symbolToNumber(number, symbol)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Atomic number
   integer, intent(out) :: number

   character(len=2) :: lcSymbol
   integer :: i, j, k, l

   number = 0
   lcSymbol = '  '

   k = 0
   do j = 1, len_trim(symbol)
      if (k > 2) exit
      l = iachar(symbol(j:j))
      if (k >= 1 .and. l == iachar(' ')) exit
      if (k >= 1 .and. l == 9) exit
      if (l >= iachar('A') .and. l <= iachar('Z')) l = l + offset
      if (l >= iachar('a') .and. l <= iachar('z')) then
         k = k+1
         lcSymbol(k:k) = achar(l)
      endif
   enddo

   do i = 1, size(lcPse)
      if (lcSymbol == lcPse(i)) then
         number = i
         exit
      endif
   enddo

   if (number == 0) then
      select case(lcSymbol)
      case('d ', 't ')
         number = 1
      end select
   end if

end subroutine symbolToNumber


!> Convert atomic number to element symbol
elemental subroutine numberToSymbol(symbol, number)

   !> Atomic number
   integer, intent(in) :: number

   !> Element symbol
   character(len=2), intent(out) :: symbol

   if (number <= 0 .or. number > size(pse)) then
      symbol = '--'
   else
      symbol = pse(number)
   endif

end subroutine numberToSymbol


!> Convert atomic number to element symbol
elemental subroutine numberToLcSymbol(symbol, number)

   !> Atomic number
   integer, intent(in) :: number

   !> Element symbol
   character(len=2), intent(out) :: symbol

   if (number <= 0 .or. number > size(lcPse)) then
      symbol = '--'
   else
      symbol = lcPse(number)
   endif

end subroutine numberToLcSymbol


!> Convert element symbol to atomic number
elemental function toNumber(symbol) result(number)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Atomic number
   integer :: number

   call symbolToNumber(number, symbol)

end function toNumber


!> Convert atomic number to element symbol
elemental function toSymbol(number) result(symbol)

   !> Atomic number
   integer,intent(in) :: number

   !> Element symbol
   character(len=2) :: symbol

   call numberToSymbol(symbol, number)

end function toSymbol


!> Convert atomic number to element symbol
elemental function toLcSymbol(number) result(symbol)

   !> Atomic number
   integer,intent(in) :: number

   !> Element symbol
   character(len=2) :: symbol

   call numberToLcSymbol(symbol, number)

end function toLcSymbol


!> Get chemical identity from a list of atomic numbers
subroutine getIdentityNumber(nId, identity, number)

   !> Number of unique species
   integer, intent(out) :: nId

   !> Ordinal numbers
   integer, intent(in) :: number(:)

   !> Chemical identity
   integer, intent(out) :: identity(:)

   integer, allocatable :: iTmp(:)
   integer :: nAt, iAt, iId

   nAt = size(identity)
   allocate(iTmp(nAt))
   nId = 0
   do iAt = 1, nAt
      iId = findNumber(iTmp(:nId), number(iAt))
      if (iId == 0) then
         call appendNumber(iTmp, nId, number(iAt))
         iId = nId
      end if
      identity(iAt) = iId
   end do

end subroutine getIdentityNumber


!> Get chemical identity from a list of element symbols
subroutine getIdentitySymbol(nId, identity, symbol)

   !> Number of unique species
   integer, intent(out) :: nId

   !> Element symbols
   character(len=symbolLength), intent(in) :: symbol(:)

   !> Chemical identity
   integer, intent(out) :: identity(:)

   character(len=symbolLength), allocatable :: sTmp(:)
   integer :: nAt, iAt, iId

   nAt = size(identity)
   allocate(sTmp(nAt))
   nId = 0
   do iAt = 1, nAt
      iId = findSymbol(sTmp(:nId), symbol(iAt))
      if (iId == 0) then
         call appendSymbol(sTmp, nId, symbol(iAt))
         iId = nId
      end if
      identity(iAt) = iId
   end do

end subroutine getIdentitySymbol


!> Find element symbol in an unordered list, all entries are required to be unique
pure function findSymbol(list, symbol) result(position)

   !> List of element symbols
   character(len=*), intent(in) :: list(:)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> Position of the symbol in list if found, otherwise zero
   integer :: position
   integer :: iSym

   position = 0
   do iSym = 1, size(list)
      if (symbol == list(iSym)) then
         position = iSym
         exit
      end if
   end do

end function findSymbol


!> Find atomic number in an unordered list, all entries are required to be unique
function findNumber(list, number) result(position)

   !> List of atomic numbers
   integer, intent(in) :: list(:)

   !> Atomic number
   integer, intent(in) :: number

   !> Position of the number in list if found, otherwise zero
   integer :: position
   integer :: iNum

   position = 0
   do iNum = 1, size(list)
      if (number == list(iNum)) then
         position = iNum
         exit
      end if
   end do

end function findNumber


!> Append an element symbol to an unsorted list, to ensure no dublicates search
!> for the element symbol first
pure subroutine appendSymbol(list, nList, symbol)

   !> List of element symbols
   character(len=*), allocatable, intent(inout) :: list(:)

   !> Current occupied size of list
   integer, intent(inout) :: nList

   !> Elements symbol
   character(len=*), intent(in) :: symbol

   if (nList >= size(list)) then
      call resize(list)
   end if

   nList = nList + 1
   list(nList) = symbol

end subroutine appendSymbol


!> Append an atomic number to an unsorted list, to ensure no dublicates search
!> for the atomic number first
pure subroutine appendNumber(list, nList, number)

   !> List of atomic number
   integer, allocatable, intent(inout) :: list(:)

   !> Current occupied size of list
   integer, intent(inout) :: nList

   !> Atomic number
   integer, intent(in) :: number

   if (nList >= size(list)) then
      call resize(list)
   end if

   nList = nList + 1
   list(nList) = number

end subroutine appendNumber


end module xtb_mctc_symbols
