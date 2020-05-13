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

!> TODO
module xtb_type_vendordata
   use xtb_mctc_accuracy, only : wp
   implicit none
   private

   public :: pdb_data, sdf_data, vasp_info, turbo_info, struc_info
   public :: resize


   !> atomic pdb data type
   !
   !  keeps information from PDB input that is currently not used by the
   !  caller program (like residues or chains) but is needed to write
   !  the PDB output eventually
   type :: pdb_data
! ATOM   2461  HA3 GLY A 153     -10.977  -7.661   2.011  1.00  0.00           H
! TER    2462      GLY A 153
! a6----i5---xa4--aa3-xai4--axxxf8.3----f8.3----f8.3----f6.2--f6.2--xxxxxxa4--a2a2
! HETATM 2463  CHA HEM A 154       9.596 -13.100  10.368  1.00  0.00           C
      logical :: het = .false.
      integer :: charge = 0
      integer :: residue_number = 0
      character(len=4) :: name = ' '
      character(len=1) :: loc = ' '
      character(len=3) :: residue = ' '
      character(len=1) :: chains = ' '
      character(len=1) :: code = ' '
      character(len=4) :: segid = ' '
   end type pdb_data

   !> sdf atomic data
   !
   !  we only support some entries, the rest is simply dropped.
   !  the format is: ddcccssshhhbbbvvvHHHrrriiimmmnnneee
   type :: sdf_data
      integer :: isotope = 0   !< d field
      integer :: charge = 0    !< c field
      integer :: hydrogens = 0 !< h field
      integer :: valence = 0   !< v field
   end type sdf_data

   !> vasp input data
   !
   !  contains specific vasp keywords that modify the appearance of the
   !  input file and can be used to reproduce it in the output
   type :: vasp_info
      real(wp) :: scale = 1.0_wp
      logical :: selective = .false.
      logical :: cartesian = .false.
   end type vasp_info

   !> Turbomole input data
   !
   !  Saves preferences for cartesian vs. direct coordinates and lattice vs. cell.
   !  Also saves units of input data groups.
   type :: turbo_info
      logical :: cartesian = .true.
      logical :: lattice = .true.
      logical :: angs_lattice = .false.
      logical :: angs_coord = .false.
   end type turbo_info


   !> structure input info
   !
   !  contains informations from different input file formats
   !  SDF: is the structure 2d or 3D?
   type :: struc_info
      logical :: two_dimensional = .false.
      logical :: add_hydrogen = .false.
   end type struc_info


   interface resize
      module procedure resize_pdb_data
   end interface


contains


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


end module xtb_type_vendordata
