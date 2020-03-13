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

!> Implementation of a map from element symbols/atomic numbers to ids
module xtb_type_identitymap
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_resize, only : resize
   use xtb_mctc_symbols, only : symbolLength
   use xtb_type_molecule, only : TMolecule
   implicit none
   private

   public :: TIdentityMap, init


   !> Maps an identity to all its atoms
   type :: TAtomicMap

      !> Positions in the atomic array
      integer, allocatable :: pos(:)

   end type TAtomicMap


   !> Map from element symbols to identity
   type :: TIdentityMap
      private

      !> Atomic numbers for each id
      integer, allocatable :: num(:)

      !> Element symbols for each id
      character(len=symbolLength), allocatable :: sym(:)

      !> Maps from id to its atoms
      type(TAtomicMap), allocatable :: map(:)

   contains

      !> Write informative printout for the map
      procedure :: writeInfo

      !> Check if symbol or number maps to atoms
      generic :: has => hasSymbol, hasNumber

      !> Check if element symbol maps to atoms
      procedure, private :: hasSymbol

      !> Check if atomic number maps to atoms
      procedure, private :: hasNumber

      !> Get indices
      generic :: get => getIndexSymbol, getIndexNumber

      !> Get indices for an unique element symbol
      procedure, private :: getIndexSymbol

      !> Get indices for a certain atomic number
      procedure, private :: getIndexNumber

      !> Setter functions in atomic arrays
      generic :: set => setRealWithSymbol, setRealWithNumber

      !> Set value in an array using an unique element symbol
      procedure, private :: setRealWithSymbol

      !> Set value in an array using the atomic number
      procedure, private :: setRealWithNumber

   end type TIdentityMap


   !> Initialize identity map
   interface init
      module procedure :: initIdentityMapFromMolecule
      module procedure :: initIdentityMapFromArrays
   end interface init


contains


!> Initialize identity map from molecular structure data
subroutine initIdentityMapFromMolecule(self, mol)

   !> Instance of the map
   type(TIdentityMap), intent(out) :: self

   !> Molecular structure data
   type(TMolecule), intent(in) :: mol

   call init(self, mol%id, mol%sym, mol%at)

end subroutine initIdentityMapFromMolecule


!> Initialize identity map from at atomic numbers and element symbols
subroutine initIdentityMapFromArrays(self, id, sym, num)

   !> Instance of the map
   type(TIdentityMap), intent(out) :: self

   !> Atomic numbers for each id
   integer, intent(in) :: id(:)

   !> Atomic numbers for each atom
   integer, intent(in) :: num(:)

   !> Element symbols for each atom
   character(len=symbolLength), intent(in) :: sym(:)

   integer :: nAt, nId
   integer :: iId, iAt, thisAt, initialSize, nMp
   integer, allocatable :: pos(:)

   nAt = size(id)
   nId = maxval(id)

   allocate(self%map(nId))
   allocate(self%num(nId))
   allocate(self%sym(nId))

   initialSize = nAt / nId + 1
   allocate(pos(initialSize))

   do iId = 1, nId
      pos = 0
      nMp = 0
      thisAt = 0
      do iAt = 1, nAt
         if (id(iAt) == iId) then
            if (thisAt == 0) thisAt = iAt
            if (nMp >= size(pos)) then
               call resize(pos)
            end if
            nMp = nMp + 1
            pos(nMp) = iAt
         end if
      end do
      self%num(iId) = num(thisAt)
      self%sym(iId) = sym(thisAt)
      allocate(self%map(iId)%pos(nMp))
      self%map(iId)%pos(:) = pos(:nMp)
   end do

end subroutine initIdentityMapFromArrays


subroutine writeInfo(self, unit)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Unit for I/O
   integer, intent(in) :: unit

   integer :: iId, iStart, iEnd, iDel
   character(len=:), allocatable :: list

   write(unit, '(a5,1x,a4,1x,a4,3x,a5)') "ID", "Z", "sym.", "atoms"
   do iId = 1, size(self%sym)
      write(unit, '(i5,1x,i4,1x,a4,3x)', advance='no') &
         & iId, self%num(iId), self%sym(iId)
      call serializeAtomicMap(self%map(iId)%pos, list, '-', ', ')
      iStart = 1
      iEnd = min(60, len(list))
      do while(iStart <= len(list))
         if (iEnd >= len(list)) then
            iDel = iEnd
         else
            iDel = index(list(iStart:iEnd), ' ', back=.true.) + iStart - 2
            if (iDel < iStart) then
               iDel = index(list(iEnd:), ' ') + iEnd - 2
               if (iDel == iEnd - 2) iDel = len(list)
            end if
         end if
         write(unit, '(a)') list(iStart:iDel)
         iStart = iDel+1
         iEnd = min(iStart+60, len(list))
         if (iStart < len(list)) then
            write(unit, '(17x)', advance='no')
         end if
      end do
   end do

end subroutine writeInfo


subroutine serializeAtomicMap(list, string, separator, delimiter)

   !> Instance of an atomic map
   integer, intent(in) :: list(:)

   !> Serialized form of the atomic map
   character(len=:), allocatable, intent(out) :: string

   !> Separator for list
   character(len=*), intent(in), optional :: separator

   !> Separator for list
   character(len=*), intent(in), optional :: delimiter

   integer :: ii, jj, kk
   character(len=64) :: buffer
   character(len=:), allocatable :: sep, del

   if (present(separator)) then
      sep = separator
   else
      sep = ', '
   end if

   if (present(delimiter)) then
      del = delimiter
   else
      del = '-'
   end if

   string = ''
   jj = list(1)
   kk = jj
   do ii = 2, size(list)
      if (jj+1 == list(ii)) then
         jj = jj+1
      else
         if (jj-kk == 0) then
            write(buffer, '(i0)') jj
         else if (jj-kk == 1) then
            write(buffer, '(i0,a,i0)') kk, del, jj
         else
            write(buffer, '(i0,a,i0)') kk, sep, jj
         end if
         string = string // trim(buffer) // del
         jj = list(ii)
         kk = jj
      end if
   end do
   if (jj-kk == 0) then
      write(buffer, '(i0)') jj
   else if (jj-kk == 1) then
      write(buffer, '(i0,a,i0)') kk, del, jj
   else
      write(buffer, '(i0,a,i0)') kk, sep, jj
   end if
   string = string // trim(buffer)

end subroutine serializeAtomicMap


!> Check if element symbol maps to atoms
pure function hasSymbol(self, sym) result(exist)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Symbol is present
   logical :: exist

   exist = any(self%sym == sym)

end function hasSymbol


!> Check if atomic number maps to atoms
pure function hasNumber(self, num) result(exist)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Atomic number
   integer, intent(in) :: num

   !> Number is present
   logical :: exist

   exist = any(self%num == num)

end function hasNumber


!> Get indices for an unique element symbol
subroutine getIndexSymbol(self, indx, sym)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Indices for all selected atoms
   integer, allocatable, intent(out) :: indx(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   integer :: iId

   do iId = 1, size(self%sym)
      if (self%sym(iId) == sym) then
         indx = self%map(iId)%pos
         exit
      end if
   end do

end subroutine getIndexSymbol


!> Get indices for a certain atomic number
subroutine getIndexNumber(self, indx, num)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Indices for all selected atoms
   integer, allocatable, intent(out) :: indx(:)

   !> Atomic number
   integer, intent(in) :: num

   integer :: nAt, iId, newSize
   integer, allocatable :: pos(:)

   nAt = 0
   do iId = 1, size(self%num)
      if (self%num(iId) == num) then
         newSize = nAt + size(self%map(iId)%pos)
         if (newSize >= size(pos)) then
            call resize(pos, newSize)
         end if
         pos(nAt+1:newSize) = self%map(iId)%pos
      end if
   end do

   if (nAt > 0) then
      call move_alloc(pos, indx)
   end if

end subroutine getIndexNumber


!> Set value in an array using an unique element symbol
subroutine setRealWithSymbol(self, array, sym, val)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Instance of the array
   real(wp), intent(inout) :: array(:)

   !> Element symbol
   character(len=*), intent(in) :: sym

   !> Value to set
   real(wp), intent(in) :: val

   integer :: iId

   do iId = 1, size(self%sym)
      if (self%sym(iId) == sym) then
         array(self%map(iId)%pos) = val
      end if
   end do

end subroutine setRealWithSymbol


!> Set value in an array using the atomic number
subroutine setRealWithNumber(self, array, num, val)

   !> Instance of the map
   class(TIdentityMap), intent(in) :: self

   !> Instance of the array
   real(wp), intent(inout) :: array(:)

   !> Atomic number
   integer, intent(in) :: num

   !> Value to set
   real(wp), intent(in) :: val

   integer :: iId

   do iId = 1, size(self%num)
      if (self%num(iId) == num) then
         array(self%map(iId)%pos) = val
      end if
   end do

end subroutine setRealWithNumber


end module xtb_type_identitymap
