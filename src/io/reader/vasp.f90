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

module xtb_io_reader_vasp
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_strings
   use xtb_mctc_systools
   use xtb_type_molecule
   use xtb_pbc_tools
   use xtb_mctc_symbols, only : toNumber, symbolLength
   use xtb_type_molecule
   use xtb_type_vendordata, only : vasp_info
   implicit none
   private

   public :: readMoleculeVasp


   logical, parameter :: debug = .false.


contains


subroutine readMoleculeVasp(mol, unit, status, iomsg)
   type(TMolecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   real(wp) :: lattice(3,3)
   integer, allocatable :: ncount(:)
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   integer :: natoms
   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=2*symbolLength) :: args(256),args2(256)
   real(wp), allocatable :: xyz(:, :)
   character(len=symbolLength), allocatable :: sym(:)

   integer i,j,k,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck

   integer :: iat, idum, err

   status = .false.

   lattice=0

   ncheck=0
   ntype=0
   ! first line contains the symbols of different atom types
   call getline(unit,line,err)
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   if (debug) print'(">",a)',line
   call parse(line,' ',args,ntype)

   ! this line contains the global scaling factor,
   call getline(unit,line,err)
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   if (debug) print'(">",a)',line
   read(line,*,iostat=err) ddum
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   ! the Ang->au conversion is included in the scaling factor
   if (debug) print'("->",g0)',ddum
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(unit,line,err)
      if (err.ne.0) then
         iomsg = "Could not read lattice from POSCAR"
         return
      endif
      if (debug) print'("->",a)',line
      read(line,*,iostat=err) latvec
      if (err.ne.0) then
         iomsg = "Could not read lattice from POSCAR"
         return
      endif
      lattice(:,i) = latvec * scalar
   enddo
   ! Either here are the numbers of each element,
   ! or (>vasp.5.1) here are the element symbols
   call getline(unit,line,err)
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   if (debug) print'(">",a)',line
   ! try to read first element
   read(line,*,iostat=err) idum
   ! CONTCAR files have additional Element line here since vasp.5.1
   if (err.ne.0) then
      call parse(line,' ',args,ntype)
      call getline(unit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) then
         iomsg = "Could not read POSCAR"
         return
      endif
   endif
   call parse(line,' ',args2,nn)
   if (nn.ne.ntype) then
      iomsg = 'Error reading number of atomtypes'
      return
   endif

   allocate(ncount(nn), source = 0)
   ncheck = 0
   do i = 1, nn
      read(args2(i), *, iostat=err) ncount(i)
      iat = toNumber(args(i))
      if (iat < 1 .or. ncount(i) < 1) then
         iomsg = 'unknown element.'
         return
      endif
   enddo

   natoms = sum(ncount)
   allocate(sym(natoms))
   allocate(xyz(3, natoms))

   k = 0
   do i = 1, nn
      do j = 1, ncount(i)
         k = k+1
         sym(k) = trim(args(i))
      enddo
   enddo

   call getline(unit, line, err)
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   if (debug) print'(">",a)', line
   line=adjustl(line)
   if (line(:1).eq.'s' .or. line(:1).eq.'S') then
      selective = .true.
      call getline(unit,line,err)
      if (debug) print'("->",a)', line
      if (err.ne.0) then
         iomsg = "Could not read POSCAR"
         return
      endif
      line = adjustl(line)
   endif

   cartesian=(line(:1).eq.'c' .or. line(:1).eq.'C' .or. &
      &       line(:1).eq.'k' .or. line(:1).eq.'K')
   do i = 1, natoms
      call getline(unit, line, err)
      if (err.ne.0) then
         iomsg = "Could not read geometry from POSCAR"
         return
      endif
      if (debug) print'("-->",a)',line
      read(line, *, iostat=err) coord
      if (err.ne.0) then
         iomsg = "Could not read geometry from POSCAR"
         return
      endif

      if (cartesian) then
         xyz(:,i) = coord*scalar
      else
         xyz(:,i) = matmul(lattice, coord)
      endif

   enddo

   call init(mol, sym, xyz, lattice=lattice)
   ! save information about this POSCAR for later
   mol%vasp = vasp_info(scale=ddum, selective=selective, cartesian=cartesian)

   status = .true.

end subroutine readMoleculeVasp


end module xtb_io_reader_vasp
