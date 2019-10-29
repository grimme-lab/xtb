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

submodule(tbdef_molecule) molecule_reader
   use tbdef_molecule
   implicit none

   interface assignment(=)
      module procedure :: symbol_to_number
   end interface assignment(=)

contains

module subroutine read_molecule_generic(self, unit, format)
   use tbmod_file_utils
   class(tb_molecule), intent(out) :: self
   integer, intent(in) :: unit
   integer, intent(in), optional :: format
   integer :: ftype
   logical :: status
   character(len=:), allocatable :: message
   if (present(format)) then
      ftype = format
   else
      ftype = self%ftype
   endif

   select case(ftype)
   case(p_ftype%xyz)
      call read_molecule_xyz(self, unit, status, iomsg=message)
   case(p_ftype%tmol)
      call read_molecule_tmol(self, unit, status, iomsg=message)
   case(p_ftype%molfile)
      call read_molecule_molfile(self, unit, status, iomsg=message)
   case(p_ftype%sdf)
      call read_molecule_sdf(self, unit, status, iomsg=message)
   case(p_ftype%vasp)
      call read_molecule_vasp(self, unit, status, iomsg=message)
   case(p_ftype%pdb)
      call read_molecule_pdb(self, unit, status, iomsg=message)
   case default
      status = .false.
      message = "coordinate format known"
   end select

   if (.not.status) call raise('E', message, 1)

   self%ftype = ftype

end subroutine read_molecule_generic


subroutine read_molecule_xyz(mol, unit, status, iomsg)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use pbc_tools
   use readin, only : getline => strip_line
   logical, parameter :: debug = .false.
   class(tb_molecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer  :: n, iat
   real(wp) :: xyz(3)
   real(wp) :: conv

   character(len=:),allocatable :: message
   character(len=:),allocatable :: line
   character(len=10) :: chdum
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

   call mol%allocate(n)
   mol%npbc = 0 ! Xmol is always molecular (there are extensions to this...)
   mol%pbc = .false.

   ! drop next record
   read(unit,'(a)')

   n = 0
   do while (n < mol%n)
      call getline(unit,line,err)
      if (is_iostat_end(err)) exit
      if (err.ne.0) then
         iomsg = "Could not read geometry from Xmol file"
         return
      endif
      if (debug) print'(">",a)',line
      read(line,*,iostat=err) chdum, xyz(1), xyz(2), xyz(3)
      if (err.ne.0) then
         iomsg = "Could not parse coordinates from Xmol file"
         return
      endif
      if (debug) print'("->",a)',chdum
      if (debug) print'("->",3g0)',xyz

      iat = chdum
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         n = n+1
         mol%at(n) = iat
         mol%sym(n) = trim(chdum)
         mol%xyz(:,n) = xyz*conv
      else
         iomsg = "Unknown element symbol: '"//trim(chdum)//"'"
         return
      endif
   enddo

   if (n.ne.mol%n) then
      iomsg = "Atom number missmatch in Xmol file"
      return
   endif

   status = .true.
end subroutine read_molecule_xyz


subroutine read_molecule_tmol(mol, unit, status, iomsg)
   use iso_fortran_env, wp => real64
   use readin, getline => strip_line
   use pbc_tools
   implicit none
   type(tb_molecule),intent(inout) :: mol
   integer,intent(in) :: unit !< file handle
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:),allocatable :: line
   integer :: err
   integer :: idum
   logical,parameter :: debug = .false.
   integer :: natoms
   integer :: i
   integer :: npbc
   logical :: exist

   character,parameter :: flag = '$'
   character(len=*),parameter :: flag_end = flag//'end'

   status = .false.

   exist  = .false.
   npbc   = -1
   natoms = -1

!  we need to read this file twice, for a lot of reasons
   rewind(unit) ! FIXME

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on unit (dammit Turbomole format)
   if (debug) print'("first scan")'
   call getline(unit,line,err)
   scanfile: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'coord').eq.2) then
            if (debug) print'("*>",a)',line
            call get_atomnumber(natoms,unit,line,err)
         elseif (index(line,'periodic').eq.2) then
            if (debug) print'("*>",a)',line
            call get_periodic(npbc,unit,line,err)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(unit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(unit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit scanfile
      if (index(line,flag_end).ne.0) exit scanfile
   enddo scanfile

   if (natoms.eq.-1) then
      iomsg = "coord keyword was not found"
      return
   endif
   if (npbc.eq.-1) npbc = 0

   if (natoms.lt.1) then ! this will catch an empty coord data group
      iomsg = 'Found no atoms, cannot work without atoms!'
      return
   endif

   call mol%allocate(natoms)
   mol%npbc = npbc
   do i = 1, npbc
      mol%pbc(i) = .true.
   enddo
   if (mol%npbc == 0) then
      ! this might be wrong
      mol%lattice = reshape([0.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,0.0_wp],[3,3])
      mol%rec_lat = reshape([1.0_wp,0.0_wp,0.0_wp,&
         &                   0.0_wp,1.0_wp,0.0_wp,&
         &                   0.0_wp,0.0_wp,1.0_wp],[3,3])
   endif

   rewind(unit) ! FIXME

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on unit (dammit Turbomole format)
   if (debug) print'("one''n''half read")'
   if (npbc > 0) then
   call getline(unit,line,err)
   readlat: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'lattice').eq.2) then
            if (exist) then
               iomsg = "Multiple definitions of lattice present!"
               return
            endif
            exist = .true.
            if (debug) print'("*>",a)',line
            call get_lattice(unit,line,err,npbc,mol%lattice)
            call dlat_to_cell(mol%lattice,mol%cellpar)
            call dlat_to_rlat(mol%lattice,mol%rec_lat)
            mol%volume = dlat_to_dvol(mol%lattice)
         elseif (index(line,'cell').eq.2) then
            if (exist) then
               iomsg = "Multiple definitions of lattice present!"
               return
            endif
            exist = .true.
            if (debug) print'("*>",a)',line
            call get_cell(unit,line,err,npbc,mol%cellpar)
            call cell_to_dlat(mol%cellpar,mol%lattice)
            call cell_to_rlat(mol%cellpar,mol%rec_lat)
            mol%volume = cell_to_dvol(mol%cellpar)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(unit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(unit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readlat
      if (index(line,flag_end).ne.0) exit readlat
   enddo readlat

   rewind(unit) ! FIXME
   endif

!  read first line before the readloop starts, I have to do this
!  to avoid using backspace on unit (dammit Turbomole format)
   if (debug) print'("second read")'
   call getline(unit,line,err)
   readfile: do
      !  check if there is a $ in the *first* column
      if (index(line,flag).eq.1) then
         if (index(line,'coord').eq.2) then
            if (debug) print'("*>",a)',line
            call get_coord(unit,line,err,mol)
         else
!           get a new line
            if (debug) print'(" >",a)',line
            call getline(unit,line,err)
         endif
      else ! not a keyword -> ignore
         if (debug) print'(" >",a)',line
         call getline(unit,line,err)
      endif
   !  check for end of file, which I will tolerate as alternative to $end
      if (err.eq.iostat_end) exit readfile
      if (index(line,flag_end).ne.0) exit readfile
   enddo readfile

   if (.not.exist .and. npbc > 0) then
      iomsg = "There is no definition of the lattice present!"
      return
   endif

   status = .true.

contains

!> count the number of lines in the coord block -> number of atoms
subroutine get_atomnumber(ncount,unit,line,err)
   implicit none
   integer,intent(in) :: unit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: ncount        !< number of readin lines
   integer,intent(out) :: err           !< status of unit
   ncount = 0
   do
      call getline(unit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
   enddo
end subroutine get_atomnumber

!> get the dimensionality of the system
subroutine get_periodic(npbc,unit,line,err)
   implicit none
   integer,intent(in) :: unit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: npbc          !< number of periodic dimensions
   integer,intent(out) :: err           !< status of unit
   integer :: idum
   if (get_value(line(10:),idum)) then
      if (idum.eq.0 .or. idum.eq.3) then
         npbc = idum
      else
         iomsg = "This kind of periodicity is not implemented"
         return
      endif
   else
      iomsg = "Could not read periodicity from '"//line//"'"
      return
   endif
   call getline(unit,line,err)
end subroutine get_periodic

!> read the lattice vectors from the data group
subroutine get_lattice(unit,line,err,npbc,lat)
   use mctc_econv
   implicit none
   integer,intent(in) :: unit          !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out) :: err           !< status of unit
   integer,intent(in)  :: npbc          !< number of periodic dimensions
   real(wp),intent(out) :: lat(3,3)     !< lattice vectors
   integer :: ncount
   real(wp) :: lvec(3)
   real(wp) :: conv
   if (npbc == 0) then
      iomsg = "lattice data group is present but not periodic?"
      return
   endif
   if (len(line).le.8) then
      conv = 1.0_wp ! default is bohr
   elseif (index(line(9:),'bohr').gt.0) then
      conv = 1.0_wp
   elseif (index(line(9:),'angs').gt.0) then
      conv = aatoau
   else ! fall back to default
      iomsg = "Could not make sense out of unit '"//line(7:)//"'"
      conv = 1.0_wp
   endif
   lat = 0.0_wp
   ncount = 0
   do
      call getline(unit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
      if (ncount.gt.npbc) then
         iomsg = "Input error in lattice data group"
         return
      endif
      read(line,*,iostat=err) lvec(1:npbc)
      if (err.ne.0) then
         iomsg = "Input error in lattice data group"
         return
      endif
      lat(1:npbc,ncount) = lvec(1:npbc)*conv
   enddo
end subroutine get_lattice

!> read the cell parameters and transform to lattice
subroutine get_cell(unit,line,err,npbc,cellpar)
   use mctc_constants
   use mctc_econv
   implicit none
   integer,intent(in)   :: unit        !< file handle
   character(len=:),allocatable :: line !< string buffer
   integer,intent(out)  :: err          !< status of unit
   integer,intent(in)   :: npbc         !< number of periodic dimensions
   real(wp),intent(out) :: cellpar(6)   !< cell parameter
   real(wp) :: conv,vol
   if (npbc == 0) then
      iomsg = "cell data group is present but not periodic?"
      return
   endif

   associate( alen => cellpar(1), blen => cellpar(2), clen => cellpar(3), &
              alp  => cellpar(4), bet  => cellpar(5), gam  => cellpar(6) )

   alen = 1.0_wp; blen = 1.0_wp; clen = 1.0_wp
   alp = 90.0_wp; bet = 90.0_wp; gam = 90.0_wp

   if (len(line).le.5) then
      conv = 1.0_wp ! default is bohr
   elseif (index(line(6:),'bohr').gt.0) then
      conv = 1.0_wp
   elseif (index(line(6:),'angs').gt.0) then
      conv = aatoau
   else ! fall back to default
      iomsg = "Could not make sense out of unit '"//line(7:)//"'"
      conv = 1.0_wp
   endif
   call getline(unit,line,err)
   if (err.eq.iostat_end) return
   if (index(line,flag).ne.0) return
   if (debug) print'(" ->",a)',line
   if (npbc == 3) &
      read(line,*,iostat=err) alen,blen,clen,alp,bet,gam
   if (npbc == 2) &
      read(line,*,iostat=err) alen,blen,             gam
   if (npbc == 1) &
      read(line,*,iostat=err) alen
   if (err.ne.0) then
      iomsg = "Could not read cell from '"//line//"'"
      return
   endif

   if(alen < 0.0d0 .or. blen < 0.0d0 .or. clen < 0.0d0 .or. &
      alp < 0.0d0 .or. bet < 0.0d0 .or. gam < 0.0d0) then
      iomsg = "Negative cell parameters make no sense!"
      return
   endif

   alen = alen*conv
   blen = blen*conv
   clen = clen*conv
   alp  = alp*pi/180.0_wp
   bet  = bet*pi/180.0_wp
   gam  = gam*pi/180.0_wp

   end associate

end subroutine get_cell

!> read the coordinates from coord data group
subroutine get_coord(unit,line,err,mol)
   use mctc_econv
   use mctc_strings
   use tbdef_molecule
   implicit none
   integer,intent(in) :: unit            !< file handle
   character(len=:),allocatable :: line   !< string buffer
   integer,intent(out) :: err             !< status of unit
   type(tb_molecule),intent(inout) :: mol !< molecular structure information
   integer  :: narg
   integer  :: ncount,iat,icoord
   character(len=10) :: chdum
   real(wp) :: conv,ddum,xyz(3)
   logical  :: frac
   if (len(line).le.6) then
      frac = .false.
      conv = 1.0_wp ! default is bohr
   elseif (index(line(7:),'bohr').gt.0) then
      frac = .false.
      conv = 1.0_wp
   elseif (index(line(7:),'angs').gt.0) then
      frac = .false.
      conv = aatoau
   elseif (index(line(7:),'frac').gt.0) then
      frac = .true.
      conv = 1.0_wp
   else ! fall back to default
      iomsg = "Could not make sense out of unit '"//line(7:)//"'"
      frac = .false.
      conv = 1.0_wp
   endif

   ncount = 0
   do
      call getline(unit,line,err)
      if (err.eq.iostat_end) return
      if (index(line,flag).ne.0) return
      if (debug) print'(" ->",a)',line

      if (line.eq.'') cycle ! skip empty lines
      ncount = ncount + 1   ! but count non-empty lines first
      if (ncount.gt.mol%n) then
         iomsg = "Internal error while reading coord data group"
         return
      endif
      read(line,*,iostat=err) xyz(1), xyz(2), xyz(3), chdum
      if (err.ne.0) then
         iomsg = "not enough arguments in line in coord data group"
         return
      endif
      iat = chdum
      if (iat.eq.0) then
         iomsg = "invalid element input in line in coord data group"
         return
      endif
      mol%sym(ncount) = trim(chdum)
      mol%at(ncount) = iat
      ! in case of fractional coordinates we perform a gemv with the lattice
      if (frac) then
         mol%xyz(:,ncount) = matmul(mol%lattice,xyz)
      else
         mol%xyz(:,ncount) = conv*xyz
      endif
   enddo
end subroutine get_coord

end subroutine read_molecule_tmol


subroutine read_molecule_sdf(mol, unit, status, iomsg)
   use mctc_systools
   class(tb_molecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   integer :: error

   call read_molecule_molfile(mol, unit, status, iomsg)
   if (.not.status) return
   status = .false.

   error = 0
   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, '$$$$') == 1) exit
      call mol%info%push_back(line)
   enddo
   if (error /= 0) then
      iomsg = "failed while reading SDF key-value pairs"
      return
   endif

   status = .true.

end subroutine read_molecule_sdf


subroutine read_molecule_molfile(mol, unit, status, iomsg)
   use mctc_econv
   use mctc_systools
   class(tb_molecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   character(len=:), allocatable :: name
   integer :: i, iatom, jatom, ibond, btype, atomtype
   integer :: error, length, charge(2,15)
   integer :: number_of_atoms, number_of_bonds, number_of_atom_lists, &
      &       chiral_flag, number_of_stext_entries, i999
   integer :: list4(4), list12(12)
   real(wp) :: x, y, z
   character(len=2) :: sdf_dim
   character(len=3) :: symbol
   character(len=5) :: v2000
   integer, parameter :: ccc_to_charge(0:7) = [0, +3, +2, +1, 0, -1, -2, -3]

   status = .false.

   call getline(unit, name, error)
   call getline(unit, line, error)
   read(line, '(20x,a2)', iostat=error) sdf_dim
   if (error == 0 .and. (sdf_dim == '2D' .or. sdf_dim == '2d')) then
      iomsg = "two dimensional structures are not a valid input for this program"
      return
   endif
   call getline(unit, line, error)
   call getline(unit, line, error)
   read(line, '(3i3,3x,2i3,12x,i3,1x,a5)', iostat=error) &
      & number_of_atoms, number_of_bonds, number_of_atom_lists, &
      & chiral_flag, number_of_stext_entries, i999, v2000
   if (error /= 0) then
      iomsg = "could not read header of molfile"
      return
   endif

   call mol%allocate(number_of_atoms)
   allocate(mol%sdf(len(mol)), source=sdf_data())
   if (len(name) > 0) mol%name = name

   do iatom = 1, number_of_atoms
      call getline(unit, line, error)
      read(line, '(3f10.4,1x,a3,i2,11i3)', iostat=error) &
         & x, y, z, symbol, list12
      if (error /= 0) then
         iomsg = "could not coordinates from connection table"
         return
      endif
      atomtype = symbol
      mol%xyz(:, iatom) = [x, y, z] * aatoau
      mol%at(iatom) = atomtype
      mol%sym(iatom) = trim(symbol)
      mol%sdf(iatom)%isotope = list12(1)
      mol%sdf(iatom)%charge = ccc_to_charge(list12(2)) ! drop doublet radical
      mol%sdf(iatom)%hydrogens = list12(4)
      mol%sdf(iatom)%valence = list12(6)
   enddo

   call mol%bonds%allocate(size=number_of_bonds)
   do ibond = 1, number_of_bonds
      call getline(unit, line, error)
      read(line, '(7i3)', iostat=error) &
         & iatom, jatom, btype, list4
      if (error /= 0) then
         iomsg = "could not topology from connection table"
         return
      endif
      call mol%bonds%push_back([iatom, jatom])
   enddo

   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, 'M  END') == 1) exit
      if (index(line, 'M  CHG') == 1) then
         read(line(7:10), *) length
         read(line(11:), '(*(1x,i3,1x,i3))') (charge(:, i), i=1, length)
         do i = 1, length
            mol%sdf(charge(1, i))%charge = charge(2, i)
         enddo
      endif
   enddo
   if (error /= 0) then
      iomsg = "could not read connection table"
      return
   endif

   mol%chrg = sum(mol%sdf%charge)

   if (any(mol%sdf%hydrogens > 1)) then
      iomsg = "explicit hydrogen atoms are required"
      return
   endif

   status = .true.

end subroutine read_molecule_molfile


subroutine read_molecule_vasp(mol, unit, status, iomsg)
   use iso_fortran_env, wp => real64
   use mctc_econv
   use mctc_strings
   use mctc_systools
   use tbdef_molecule
   use pbc_tools
   logical, parameter :: debug = .false.
   type(tb_molecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   real(wp) :: lattice(3,3)
   integer, allocatable :: ncount(:)
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=80) :: args(90),args2(90)

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

   allocate( ncount(nn), source = 0 )
   ncheck=0
   do i=1,nn
      read(args2(i),*,iostat=err) ncount(i)
      iat = args(i)
      if (iat < 1 .or. ncount(i) < 1) then
         iomsg = 'unknown element.'
         return
      endif
   enddo

   call mol%allocate(sum(ncount))
   mol%pbc = .true.
   mol%npbc = 3
   k = 0
   do i = 1, nn
      iat = args(i)
      do j = 1, ncount(i)
         k = k+1
         mol%at(k) = iat
         mol%sym(k) = args(i)(1:2)
      enddo
   enddo

   call getline(unit,line,err)
   if (err.ne.0) then
      iomsg = "Could not read POSCAR"
      return
   endif
   if (debug) print'(">",a)',line
   line=adjustl(line)
   if (line(:1).eq.'s' .or. line(:1).eq.'S') then
      selective=.true.
      call getline(unit,line,err)
      if (debug) print'("->",a)',line
      if (err.ne.0) then
         iomsg = "Could not read POSCAR"
         return
      endif
      line=adjustl(line)
   endif

   cartesian=(line(:1).eq.'c' .or. line(:1).eq.'C' .or. &
      &       line(:1).eq.'k' .or. line(:1).eq.'K')
   do i = 1, len(mol)
      call getline(unit,line,err)
      if (err.ne.0) then
         iomsg = "Could not read geometry from POSCAR"
         return
      endif
      if (debug) print'("-->",a)',line
      read(line,*,iostat=err) coord
      if (err.ne.0) then
         iomsg = "Could not read geometry from POSCAR"
         return
      endif

      if (cartesian) then
         mol%xyz(:,i)=coord*scalar
      else
         mol%xyz(:,i)=matmul(lattice,coord)
      endif

   enddo

   mol%lattice = lattice
   ! save information about this POSCAR for later
   mol%vasp = vasp_info(scale=ddum, selective=selective, cartesian=cartesian)

   call dlat_to_cell(mol%lattice,mol%cellpar)
   call dlat_to_rlat(mol%lattice,mol%rec_lat)
   mol%volume = dlat_to_dvol(mol%lattice)

   call xyz_to_abc(mol%n,mol%lattice,mol%xyz,mol%abc,mol%pbc)

   status = .true.

end subroutine read_molecule_vasp


subroutine read_molecule_pdb(mol, unit, status, iomsg)
   use mctc_econv
   use mctc_systools
   use mctc_resize_arrays
   interface resize  ! add to overloaded resize interface
      module procedure :: resize_pdb_data
   end interface resize
   logical, parameter :: debug = .false.
   type(tb_molecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer, allocatable :: list(:)
   real(wp), allocatable :: xyz(:,:)
   type(pdb_data), allocatable :: pdb(:)
   character(len=2), allocatable :: sym(:)
   character(len=:), allocatable :: line
   character(len=2) :: a_charge
   integer :: iatom, jatom, iresidue, try, error, atom_type
   integer :: this_residue, last_residue
   real(wp) :: occ, temp, coords(3)
! ATOM   2461  HA3 GLY A 153     -10.977  -7.661   2.011  1.00  0.00           H
! TER    2462      GLY A 153
! a6----i5---xa4--aa3-xai4--axxxf8.3----f8.3----f8.3----f6.2--f6.2--xxxxxxa4--a2a2
! HETATM 2463  CHA HEM A 154       9.596 -13.100  10.368  1.00  0.00           C
   character(len=*), parameter :: pdb_format = &
      &  '(6x,i5,1x,a4,a1,a3,1x,a1,i4,a1,3x,3f8.3,2f6.2,6x,a4,a2,a2)'
   integer, parameter :: p_initial_size = 1000 ! this is going to be a proteine

   status = .false.

   allocate(list(p_initial_size), source=0)
   allocate(sym(p_initial_size), source='  ')
   allocate(xyz(3, p_initial_size), source=0.0_wp)
   allocate(pdb(p_initial_size), source=pdb_data())

   iatom = 0
   iresidue = 0

   error = 0
   do while(error == 0)
      call getline(unit, line, error)
      if (index(line, 'END') == 1) exit
      if (index(line, 'ATOM') == 1 .or. index(line, 'HETATM') == 1) then
         if (iatom >= size(xyz, 2)) call resize(xyz)
         if (iatom >= size(sym)) call resize(sym)
         if (iatom >= size(list)) call resize(list)
         if (iatom >= size(pdb)) call resize(pdb)
         iatom = iatom + 1
         pdb(iatom)%het = index(line, 'HETATM') == 1
         read(line, pdb_format) &
            & jatom, pdb(iatom)%name, pdb(iatom)%loc, pdb(iatom)%residue, &
            & pdb(iatom)%chains, this_residue, pdb(iatom)%code, &
            & coords, occ, temp, pdb(iatom)%segid, sym(iatom), a_charge
         xyz(:,iatom) = coords * aatoau
         atom_type = sym(iatom)
         if (atom_type == 0) then
            try = scan(pdb(iatom)%name, 'HCNOSPF')
            if (try > 0) sym(iatom) = pdb(iatom)%name(try:try)//' '
         endif
         if (this_residue /= last_residue) then
            iresidue = iresidue + 1
            last_residue = this_residue
         endif
         list(iatom) = iresidue
         read(a_charge, *, iostat=try) pdb(iatom)%charge
         if (try /= 0) pdb(iatom)%charge = 0
      endif
   enddo
   if (error /= 0) then
      iomsg = "could not read in coordinates, last line was: '"//line//"'"
      return
   endif

   call mol%allocate(iatom)
   mol%xyz = xyz(:,:iatom)
   mol%at = sym(:iatom)
   mol%sym = sym(:iatom)
   call mol%frag%allocate(list(:iatom))
   mol%pdb = pdb(:iatom)

   if (.not.all(mol%at > 0)) then
      iomsg = "invalid atom type found"
      return
   endif

   ! since PDB is used for biomolecules, this is a sensible check (prevents GIGO)
   if (.not.any(mol%at == 1)) then
      iomsg = "You get no calculation today, please add hydrogen atoms first"
      return
   endif

   status = .true.

end subroutine read_molecule_pdb

subroutine resize_pdb_data(var, n)
   type(pdb_data), allocatable, intent(inout) :: var(:)
   integer, intent(in), optional :: n
   type(pdb_data), allocatable :: tmp(:)
   integer :: length, current_length
   current_length = size(var)
   if (current_length > 0) then
      if (present(n)) then
         if (n <= current_length) return
         length = n
      else
         length = current_length + current_length/2 + 1
      endif
      allocate(tmp(length), source=pdb_data())
      tmp(:current_length) = var(:current_length)
      deallocate(var)
      call move_alloc(tmp, var)
   else
      if (present(n)) then
         length = n
      else
         length = 64
      endif
      allocate(var(length), source=pdb_data())
   endif
end subroutine resize_pdb_data


pure elemental subroutine symbol_to_number(number, symbol)
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

end subroutine symbol_to_number

end submodule molecule_reader
