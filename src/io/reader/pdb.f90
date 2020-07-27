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

module xtb_io_reader_pdb
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_convert
   use xtb_mctc_resize
   use xtb_mctc_symbols, only : toNumber, toSymbol
   use xtb_mctc_systools
   use xtb_type_molecule
   use xtb_type_vendordata, only : pdb_data, resize
   implicit none
   private

   public :: readMoleculePDB


contains


subroutine readMoleculePDB(mol, unit, status, iomsg)
   use xtb_mctc_systools
   use xtb_mctc_resize
   logical, parameter :: debug = .false.
   type(TMolecule),intent(out) :: mol
   integer,intent(in) :: unit
   logical, intent(out) :: status
   character(len=:), allocatable, intent(out) :: iomsg
   integer, allocatable :: list(:)
   real(wp), allocatable :: xyz(:,:)
   type(pdb_data), allocatable :: pdb(:)
   character(len=4), allocatable :: sym(:)
   character(len=:), allocatable :: line
   character(len=2) :: a_charge
   integer :: iatom, jatom, iresidue, try, error, atom_type
   integer :: i, j
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
   allocate(sym(p_initial_size), source='    ')
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
            & pdb(iatom)%chains, pdb(iatom)%residue_number, pdb(iatom)%code, &
            & coords, occ, temp, pdb(iatom)%segid, sym(iatom), a_charge
         xyz(:,iatom) = coords * aatoau
         atom_type = toNumber(sym(iatom))
         if (atom_type == 0) then
            try = scan(pdb(iatom)%name, 'HCNOSPF')
            if (try > 0) sym(iatom) = pdb(iatom)%name(try:try)//' '
         endif         
         !read(a_charge, *, iostat=try) pdb(iatom)%charge
         call pdbchrg2int(a_charge,pdb(iatom)%charge)
         !if (try /= 0) pdb(iatom)%charge = 0
      endif
   enddo
   if (error /= 0) then
      iomsg = "could not read in coordinates, last line was: '"//line//"'"
      return
   endif
   j=0
   associate( rn => pdb(:iatom)%residue_number)
   do i = minval(rn),maxval(rn)
      if (any(i .eq. rn)) then
         j=j+1
         where (i .eq. rn) list = j
      endif
   enddo
   end associate

   call init(mol, sym(:iatom), xyz(:, :iatom))
   call mol%frag%allocate(list(:iatom))
   mol%pdb = pdb(:iatom)
   mol%chrg = sum(mol%pdb%charge)

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

end subroutine readMoleculePDB

subroutine pdbchrg2int(chrg_str,chrg_int)
   implicit none
   character(len=*) chrg_str
   integer n,k,ios, chrg_int
   n=len_trim(chrg_str)
   forall (k=1:n) chrg_str(k:k) = chrg_str(n-k+1:n-k+1)
   read(chrg_str,*, iostat=ios) chrg_int
   if (ios /= 0 ) chrg_int=0
end subroutine pdbchrg2int


end module xtb_io_reader_pdb
