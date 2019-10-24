submodule(tbdef_molecule) molecule_reader
   use tbdef_molecule
   implicit none

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
   case(p_ftype%sdf)
      call read_molecule_sdf(self, unit, status, iomsg=message)
   case(p_ftype%vasp)
      call read_molecule_vasp(self, unit, status, iomsg=message)
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

   rewind(unit) ! FIXME

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

      call elem(line,iat)
      if (debug) print'("->",g0)',iat
      if (iat > 0) then
         n = n+1
         mol%at(n) = iat
         mol%sym(n) = line(1:2)
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
      call elem(chdum,iat)
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
   use mctc_econv
   use mctc_systools
   class(tb_molecule), intent(out) :: mol
   integer, intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   character(len=:), allocatable :: line
   integer :: iatom, jatom, ibond, btype, atomtype
   integer :: error
   integer :: number_of_atoms, number_of_bonds, number_of_atom_lists, &
      &       chiral_flag, number_of_stext_entries, i999
   integer :: list4(4), list12(12)
   real(wp) :: x, y, z
   character(len=3) :: symbol
   character(len=5) :: v2000

   status = .false.

   call getline(unit, line, error)
   call getline(unit, line, error)
   call getline(unit, line, error)
   call getline(unit, line, error)
   read(line, '(3i3,3x,2i3,12x,i3,1x,a5)', iostat=error) &
      & number_of_atoms, number_of_bonds, number_of_atom_lists, &
      & chiral_flag, number_of_stext_entries, i999, v2000

   call mol%allocate(number_of_atoms)

   do iatom = 1, number_of_atoms
      call getline(unit, line, error)
      read(line, '(3f10.4,1x,a3,i2,11i3)', iostat=error) &
         & x, y, z, symbol, list12
      call elem(symbol, atomtype)
      mol%xyz(:, iatom) = [x, y, z] * aatoau
      mol%at(iatom) = atomtype
      mol%sym(iatom) = trim(symbol)
   enddo

   call mol%bonds%allocate(size=number_of_bonds)
   do ibond = 1, number_of_bonds
      call getline(unit, line, error)
      read(line, '(7i3)', iostat=error) &
         & iatom, jatom, btype, list4
      call mol%bonds%push_back([iatom, jatom])
   enddo

   do while(error /= 0)
      call getline(unit, line, error)
      if (index(line, '$$$$') == 1) exit
   enddo

   status = .true.

end subroutine read_molecule_sdf


subroutine read_molecule_vasp(mol, unit, status, iomsg)
   use iso_fortran_env, wp => real64
   use tbdef_molecule
   use pbc_tools
   logical, parameter :: debug = .false.
   type(tb_molecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer :: n

   status = .false.

   call get_atomnumber(unit,n)
   if (allocated(iomsg)) return

   if (n.lt.1) then
      iomsg = 'Found no atoms, cannot work without atoms!'
      return
   endif

   call mol%allocate(n)
   mol%npbc = 3 ! VASP is always 3D
   mol%pbc = .true.

   call get_coord(unit,mol%lattice,mol%n,mol%xyz,mol%at,mol%sym)
   if (allocated(iomsg)) return
   call dlat_to_cell(mol%lattice,mol%cellpar)
   call dlat_to_rlat(mol%lattice,mol%rec_lat)
   mol%volume = dlat_to_dvol(mol%lattice)

   call xyz_to_abc(mol%n,mol%lattice,mol%xyz,mol%abc,mol%pbc)

   status = .true.

contains

!> read the coordinates from POSCAR
subroutine get_coord(unit,lattice,n,xyz,at,sym)
   use mctc_econv
   use mctc_strings
   use mctc_systools

   implicit none

   integer, intent(in)  :: n
   real(wp),intent(out) :: xyz(3,n)
   real(wp),intent(out) :: lattice(3,3)
   integer, intent(out) :: at(n)
   character(len=2),intent(out) :: sym(n)
   integer, intent(in)  :: unit
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=80) :: args(90),args2(90)

   integer i,j,nn,ntype,ntype2,atnum,i_dummy1,i_dummy2,ncheck

   integer :: iat, inum, idum, err

   lattice=0

   rewind(unit) ! FIXME
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
   ncheck=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      call elem(args(i),iat)
      if (iat < 1 .or. inum < 1) then
         iomsg = 'unknown element.'
         return
      endif
      do j=1,inum
         ncheck=ncheck+1
         sym(ncheck) = args(i)(1:2)
         at(ncheck)=iat
      enddo
   enddo
   if (n.ne.ncheck) then
      iomsg = 'Error reading Number of Atoms'
      return
   endif

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
   do i=1,n
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
         xyz(:,i)=coord*scalar
      else
         xyz(:,i)=matmul(lattice,coord)
      endif

   enddo

end subroutine get_coord

!> read the coordinates from POSCAR
subroutine get_atomnumber(unit,n)
   use mctc_econv
   use mctc_strings
   use mctc_systools

   implicit none

   integer, intent(out) :: n
   integer, intent(in)  :: unit
   logical              :: selective=.false. ! Selective dynamics
   logical              :: cartesian=.true.  ! Cartesian or direct

   real(wp) :: ddum,latvec(3)
   real(wp) xx(10),scalar
   real(wp) :: coord(3)
   character(len=:),allocatable :: line
   character(len=80) :: args(90),args2(90)

   integer i,j,nn,ntype,ntype2,atnum

   integer :: iat, inum, idum, err

   rewind(unit) ! FIXME
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
   scalar = ddum*aatoau

   ! reading the lattice constants
   do i=1,3
      call getline(unit,line,err)
      if (err.ne.0) then
         iomsg = "Could not read lattice from POSCAR"
         return
      endif
      if (debug) print'("->",a)',line
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
   n=0
   do i=1,nn
      read(args2(i),*,iostat=err) inum
      call elem(args(i),iat)
      if (iat < 1 .or. inum < 1) then
         iomsg = 'unknown element'
         return
      endif
      do j=1,inum
         n=n+1
      enddo
   enddo

end subroutine get_atomnumber

end subroutine read_molecule_vasp


end submodule molecule_reader
