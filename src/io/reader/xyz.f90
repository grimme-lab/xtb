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

module xtb_io_reader_xyz
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_type_molecule
   use xtb_pbc_tools
   use xtb_readin, only : getline => strip_line
   implicit none
   private

   public :: readMoleculeXYZ


   logical, parameter :: debug = .false.


contains


subroutine readMoleculeXYZ(mol, unit, status, iomsg)
   class(TMolecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer  :: ii, n, iat
   character(len=symbolLength), allocatable :: sym(:)
   real(wp), allocatable :: xyz(:, :)
   real(wp) :: x, y, z
   real(wp) :: conv

   character(len=:),allocatable :: message
   character(len=:),allocatable :: line
   character(len=symbolLength) :: chdum
   integer  :: err

   status = .false.

   conv = aatoau

   read(unit,*,iostat=err) n
   if (err.ne.0) then
      iomsg = "Could not read number of atoms, check format!"
      return
   endif

   if (n.lt.1) then
      iomsg = "Found no atoms, cannot work without atoms!"
      return
   endif

   allocate(sym(n))
   allocate(xyz(3, n))

   ! drop next record
   read(unit,'(a)')

   ii = 0
   do while (ii < n)
      call getline(unit,line,err)
      if (is_iostat_end(err)) exit
      if (err.ne.0) then
         iomsg = "Could not read geometry from Xmol file"
         return
      endif
      if (debug) print'(">",a)',line
      read(line,*,iostat=err) chdum, x, y, z
      if (err.ne.0) then
         iomsg = "Could not parse coordinates from Xmol file"
         return
      endif
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',x, y, z

      iat = toNumber(chdum)
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         ii = ii+1
         sym(ii) = trim(chdum)
         xyz(:,ii) = [x, y, z]*conv
      else
         iomsg = "Unknown element symbol: '"//trim(chdum)//"'"
         return
      endif
   enddo

   if (ii.ne.n) then
      iomsg = "Atom number missmatch in Xmol file"
      return
   endif

   call init(mol, sym, xyz)

   status = .true.

end subroutine readMoleculeXYZ


end module xtb_io_reader_xyz
