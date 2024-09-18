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

!> molecular structure information
!
!  contains information about the molecular geometry, the atom types
!  the nuclear and total charge, atomic masses and all interatomic distances
!  In periodic calculations the lattice parameters, a Wigner--Seitz cell
!  as well as fractional coordinates are attached to the data type
!
!  Additionally the molecular structure can keep a topology information
!  usually containing the bonds of the system and a list of non-overlapping
!  fragments. Both data containers are optional and need not to be filled
!  by the provider of the structure data.
!
!  Vendor specific information can be stored along with the molecular structure
!  in auxilary data objects. Data objects like PDB or SDF atomic data should
!  be kept as light-weighted as possible and the user of the structure class
!  is not required to care about it existence
module xtb_type_molecule
   use mctc_io_structure, only : structure_type, new_structure
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_boundaryconditions, only : boundaryCondition
   use xtb_mctc_symbols, only : toNumber, toSymbol, symbolLength, getIdentity
   use xtb_type_wsc
   use xtb_type_topology
   use xtb_type_fragments
   use xtb_type_buffer
   use xtb_type_vendordata
   implicit none

   public :: TMolecule, new_molecule_api, init, assignment(=)
   public :: len, size

   private


   !> Molecular structure information
   type :: TMolecule

      !> Number of atoms
      integer  :: n = 0

      !> Number of unique species
      integer  :: nId = 0

      !> Total charge
      real(wp) :: chrg = 0.0_wp

      !> Number of unpaired electrons
      integer  :: uhf = 0

      !> Periodic dimensions
      logical  :: pbc(3) = .false.

      !> Periodicity of system
      integer  :: npbc = 0

      !> Boundary condition for this structure
      integer  :: boundaryCondition = boundaryCondition%cluster

      !> Element symbols
      character(len=symbolLength), allocatable :: sym(:)

      !> Ordinal numbers
      integer, allocatable :: at(:)

      !> Chemical identity
      integer, allocatable :: id(:)

      !> Cartesian coordinates in bohr
      real(wp),allocatable :: xyz(:,:)

      !> Fractional coordinates
      real(wp),allocatable :: abc(:,:)

      !> Interatomic distances
      real(wp),allocatable :: dist(:,:)

      !> Atomic masses in amu
      real(wp),allocatable :: atmass(:)

      !> Nuclear charges
      real(wp),allocatable :: z(:)

      !> Coordination number
      real(wp),allocatable :: cn(:)

      !> Cell parameters
      real(wp) :: cellpar(6) = 0.0_wp

      !> Direct lattice parameters
      real(wp) :: lattice(3,3) = 0.0_wp

      !> Reciprocal lattice parameters
      real(wp) :: rec_lat(3,3) = 0.0_wp

      !> Volume of unit cell
      real(wp) :: volume = 0.0_wp

      !> Wigner--Seitz cell
      type(tb_wsc) :: wsc

      !> File type of the input
      integer  :: ftype = 0

      character(len=:), allocatable :: name

      type(TTopology) :: bonds

      type(TFragments) :: frag

      !> PDB specific information about residues and chains
      type(pdb_data), allocatable :: pdb(:)

      !> SDF specific information about atom types
      type(sdf_data), allocatable :: sdf(:)

      !> Information on vendor specific data
      type(structure_info) :: info = structure_info()

   contains

      procedure :: allocate => allocate_molecule

      procedure :: deallocate => deallocate_molecule

      procedure :: calculate_distances => mol_calculate_distances

      procedure :: set_nuclear_charge => mol_set_nuclear_charge

      procedure :: set_atomic_masses => mol_set_atomic_masses

      procedure :: update

      procedure :: wrap_back

      procedure :: center_of_geometry

      procedure :: center_of_mass

      procedure :: shift_to_center_of_geometry

      procedure :: shift_to_center_of_mass

      procedure :: moments_of_inertia

      procedure :: align_to_principal_axes

      procedure ::  copy => copyMolecule

   end type TMolecule


   interface init
      module procedure :: initMolecule
      module procedure :: initMoleculeSymbols
      module procedure :: initMoleculeNumbers
   end interface


   interface TMolecule
      module procedure :: new_molecule_api
   end interface TMolecule


   interface len
      module procedure :: mol_length
   end interface


   interface assignment(=)
      module procedure :: structure_to_molecule
      module procedure :: molecule_to_structure
   end interface assignment(=)


contains


!> Copy molecule type
subroutine copyMolecule(self,mol0)

   class(TMolecule), intent(out) :: self
   type(TMolecule), intent(in) :: mol0

      Call init(self, mol0%at, mol0%sym, mol0%xyz, mol0%chrg, mol0%uhf, &
           & mol0%lattice, mol0%pbc)
end subroutine copyMolecule

!> Constructor for the molecular structure type
subroutine initMolecule &
     & (mol, at, sym, xyz, chrg, uhf, lattice, pbc)

   interface
      subroutine generate_wsc(mol,wsc)
         import :: TMolecule, tb_wsc
         type(TMolecule), intent(inout) :: mol
         type(tb_wsc),    intent(inout) :: wsc
      end subroutine generate_wsc
   end interface

   type(TMolecule), intent(out) :: mol
   integer, intent(in) :: at(:)
   character(len=*), intent(in) :: sym(:)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in), optional :: chrg
   integer, intent(in), optional :: uhf
   real(wp), intent(in), optional :: lattice(:, :)
   logical, intent(in), optional :: pbc(:)

   integer, allocatable :: id(:)
   character(len=symbolLength), allocatable :: sTmp(:)
   integer :: nAt, nId, iAt, iId

   nAt = min(size(at, dim=1), size(xyz, dim=2), size(sym, dim=1))

   call mol%allocate(nAt)
   
   mol%lattice = 0.0_wp
   if (present(lattice)) then
      if (size(lattice).ne.0) mol%lattice = lattice
   end if

   if (present(pbc)) then
      if (size(pbc) == 1) then
         mol%pbc = (/pbc, pbc, pbc/)
      else
         mol%pbc = pbc
      end if
      mol%npbc = count(pbc)
      if (any(pbc)) then
         mol%boundaryCondition = boundaryCondition%pbc3d
      else
         mol%boundaryCondition = boundaryCondition%cluster
      end if
   else
      if (present(lattice)) then
         mol%boundaryCondition = boundaryCondition%pbc3d
         mol%pbc = .true.
         mol%npbc = 3
      else
         mol%boundaryCondition = boundaryCondition%cluster
         mol%pbc = .false.
         mol%npbc = 0
      end if
   end if

   call getIdentity(mol%nId, mol%id, sym)
   mol%at(:) = at(:nAt)
   mol%sym(:) = sym(:nAt)

   mol%xyz(:, :) = xyz(:, :nAt)
   if (present(chrg)) then
      mol%chrg = chrg
   else
      mol%chrg = 0
   end if
   if (present(uhf)) then
      mol%uhf = uhf
   else
      mol%uhf = 0
   end if

   call mol%set_nuclear_charge
   call mol%set_atomic_masses

   call mol%update

   call generate_wsc(mol, mol%wsc)

end subroutine initMolecule


!> Constructor for the molecular structure type
subroutine initMoleculeNumbers &
      & (mol, at, xyz, chrg, uhf, lattice, pbc)
   type(TMolecule), intent(out) :: mol
   integer, intent(in) :: at(:)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in), optional :: chrg
   integer, intent(in), optional :: uhf
   real(wp), intent(in), optional :: lattice(3, 3)
   logical, intent(in), optional :: pbc(3)

   character(len=4), allocatable :: sym(:)
   integer :: nAt

   nAt = min(size(at, dim=1), size(xyz, dim=2))
   allocate(sym(nAt))
   sym(:) = toSymbol(at(:nAt))

   call init(mol, at, sym, xyz, chrg, uhf, lattice, pbc)

end subroutine initMoleculeNumbers


!> Constructor for the molecular structure type
subroutine initMoleculeSymbols &
      & (mol, sym, xyz, chrg, uhf, lattice, pbc)
   type(TMolecule), intent(out) :: mol
   character(len=*), intent(in) :: sym(:)
   real(wp), intent(in) :: xyz(:, :)
   real(wp), intent(in), optional :: chrg
   integer, intent(in), optional :: uhf
   real(wp), intent(in), optional :: lattice(3, 3)
   logical, intent(in), optional :: pbc(3)

   integer, allocatable :: at(:)
   integer :: nAt

   nAt = min(size(sym, dim=1), size(xyz, dim=2))
   allocate(at(nAt))
   at(:) = toNumber(sym(:nAt))

   call init(mol, at, sym, xyz, chrg, uhf, lattice, pbc)

end subroutine initMoleculeSymbols


!> Constructor for the molecular structure type, compatible with C input
type(TMolecule) function new_molecule_api &
      & (n, at, xyz, chrg, uhf, lattice, pbc) result(mol)
   use iso_c_binding
   use xtb_mctc_symbols, only : toSymbol
   integer(c_int), intent(in) :: n
   integer(c_int), intent(in) :: at(n)
   real(c_double), intent(in) :: xyz(3, n)
   real(c_double), intent(in), target, optional :: chrg
   integer(c_int), intent(in), target, optional :: uhf
   real(c_double), intent(in), target, optional :: lattice(3, 3)
   logical(c_bool), intent(in), target, optional :: pbc(3)

   integer :: i

   call mol%allocate(n)

   if (c_associated(c_loc(lattice))) then
      mol%lattice = lattice
   else
      mol%lattice = 0.0_wp
   end if

   if (c_associated(c_loc(pbc))) then
      mol%pbc = pbc
      mol%npbc = count(pbc)
      if (any(pbc)) then
         mol%boundaryCondition = boundaryCondition%pbc3d
      else
         mol%boundaryCondition = boundaryCondition%cluster
      end if
   else
      if (c_associated(c_loc(lattice))) then
         mol%boundaryCondition = boundaryCondition%pbc3d
         mol%pbc = .true.
         mol%npbc = 3
      else
         mol%boundaryCondition = boundaryCondition%cluster
         mol%pbc = .false.
         mol%npbc = 0
      end if
   end if

   mol%at = at
   do i = 1, n
      mol%sym(i) = toSymbol(at(i))
   end do

   call getIdentity(mol%nId, mol%id, mol%at)

   mol%xyz = xyz
   if (c_associated(c_loc(chrg))) then
      mol%chrg = chrg
   else
      mol%chrg = 0
   end if
   if (c_associated(c_loc(uhf))) then
      mol%uhf = uhf
   else
      mol%uhf = 0
   end if

   call mol%set_nuclear_charge
   call mol%set_atomic_masses

   call mol%update

end function new_molecule_api


subroutine structure_to_molecule(mol, struc)
   type(TMolecule), intent(inout) :: mol
   type(structure_type), intent(in) :: struc
   integer :: ibd

   call initMolecule(mol, struc%num(struc%id), struc%sym(struc%id), struc%xyz, &
      & chrg=struc%charge, uhf=struc%uhf, pbc=struc%periodic, lattice=struc%lattice)
   mol%info = struc%info
   if (allocated(struc%sdf)) then
      mol%sdf = struc%sdf
   end if
   if (allocated(struc%pdb)) then
      mol%pdb = struc%pdb
   end if
   if (struc%nbd > 0) then
      call mol%bonds%allocate(3, struc%nbd)
      do ibd = 1, struc%nbd
         call mol%bonds%push_back(struc%bond(:, ibd))
      end do
   end if
end subroutine structure_to_molecule

subroutine molecule_to_structure(struc, mol)
   type(structure_type), intent(inout) :: struc
   type(TMolecule), intent(in) :: mol
   integer :: ibd, idx(3)

   call new_structure(struc, mol%at, mol%sym, mol%xyz, charge=mol%chrg, uhf=mol%uhf, &
      & periodic=mol%pbc, lattice=mol%lattice, info=mol%info)
   if (allocated(mol%sdf)) then
      struc%sdf = mol%sdf
   end if
   if (allocated(mol%pdb)) then
      struc%pdb = mol%pdb
   end if
   if (len(mol%bonds) > 0) then
      allocate(struc%bond(3, len(mol%bonds)))
      struc%nbd = len(mol%bonds)
      do ibd = 1, len(mol%bonds)
         call mol%bonds%get_item(ibd, idx)
         struc%bond(:, ibd) = idx
      end do
   end if
end subroutine molecule_to_structure


!> obtain number of atoms for molecular structure
!
!  The molecular structure is assumed to be well-behaved in this respect
!  so there is no sanity check on the allocation status.
integer pure elemental function mol_length(self) result(length)
   class(TMolecule),intent(in) :: self !< molecular structure information
   length = self%n
end function mol_length


!> constructor for molecular structure
subroutine allocate_molecule(self,n)
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   integer,intent(in) :: n
   call self%deallocate
   self%n = n
   allocate( self%id(n),          source = 0 )
   allocate( self%at(n),          source = 0 )
   allocate( self%sym(n),         source = '    ' )
   allocate( self%xyz(3,n),       source = 0.0_wp )
   allocate( self%abc(3,n),       source = 0.0_wp )
   allocate( self%dist(n,n),      source = 0.0_wp )
   allocate( self%atmass(n),      source = 0.0_wp )
   allocate( self%z(n),           source = 0.0_wp )
   allocate( self%cn(n),          source = 0.0_wp )
end subroutine allocate_molecule

!> deconstructor for molecular structure
subroutine deallocate_molecule(self)
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   self%n = 0
   self%nId = 0
   self%pbc = .false.
   self%chrg = 0.0_wp
   self%uhf = 0
   self%lattice = 0.0_wp
   if (allocated(self%id))     deallocate(self%id)
   if (allocated(self%at))     deallocate(self%at)
   if (allocated(self%sym))    deallocate(self%sym)
   if (allocated(self%xyz))    deallocate(self%xyz)
   if (allocated(self%abc))    deallocate(self%abc)
   if (allocated(self%dist))   deallocate(self%dist)
   if (allocated(self%atmass)) deallocate(self%atmass)
   if (allocated(self%z))      deallocate(self%z)
   if (allocated(self%cn))     deallocate(self%cn)
   if (allocated(self%pdb))    deallocate(self%pdb)
   if (allocated(self%sdf))    deallocate(self%sdf)
   if (allocated(self%name))   deallocate(self%name)
   call self%bonds%deallocate
   call self%frag%deallocate
end subroutine deallocate_molecule

subroutine update(self)
   use xtb_mctc_accuracy, only : wp
   use xtb_pbc_tools

   implicit none
   class(TMolecule),intent(inout) :: self  
      !! molecular structure information
   
   !> For periodic calculations
   if (self%npbc > 0) then
      call dlat_to_cell(self%lattice,self%cellpar)
      call dlat_to_rlat(self%lattice,self%rec_lat)
      self%volume = dlat_to_dvol(self%lattice)

      call self%wrap_back
   endif
   
   call self%calculate_distances

end subroutine update

!> calculates all distances for molecular structures and minimum
!> image distances for periodic structures
subroutine mol_calculate_distances(self)
   use xtb_mctc_accuracy, only : wp
   use xtb_pbc_tools
   
   implicit none
   class(TMolecule),intent(inout) :: self 
      !! molecular structure information
   integer :: i,j
   
   if (self%npbc > 0) then
      do i = 1, self%n
         do j = 1, i-1
            self%dist(j,i) = minimum_image_distance(.false.,self%abc(:,i), &
               &              self%abc(:,j),self%lattice,self%pbc)
            self%dist(i,j) = self%dist(j,i)
         enddo
         self%dist(j,i) = minimum_image_distance(.true.,self%abc(:,i), &
            &              self%abc(:,i),self%lattice,self%pbc)
      enddo
   else
      do i = 1, self%n
         do j = 1, i-1
            self%dist(j,i) = sqrt(sum((self%xyz(:,j)-self%xyz(:,i))**2))
            self%dist(i,j) = self%dist(j,i)
         enddo
         self%dist(i,i) = 0.0_wp
      enddo
   endif

end subroutine mol_calculate_distances

!> get all nuclear charges
subroutine mol_set_nuclear_charge(self)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(inout) :: self  !< molecular structure information
   integer :: i
   do i = 1, self%n
      self%z(i) = real(self%at(i),wp) - real(ncore(self%at(i)))
      if (self%at(i) > 57 .and. self%at(i) < 72) self%z(i) = 3.0_wp
   enddo
contains

pure elemental integer function ncore(at)
  integer,intent(in) :: at
  if(at.le.2)then
     ncore=0
  elseif(at.le.10)then
     ncore=2
  elseif(at.le.18)then
     ncore=10
  elseif(at.le.29)then   !zn
     ncore=18
  elseif(at.le.36)then
     ncore=28
  elseif(at.le.47)then
     ncore=36
  elseif(at.le.54)then
     ncore=46
  elseif(at.le.71)then
     ncore=54
  elseif(at.le.79)then
     ncore=68
  elseif(at.le.86)then
     ncore=78
  elseif(at.le.103)then
     ncore=86
  endif
end function ncore
end subroutine mol_set_nuclear_charge

!> get all nuclear charges
subroutine mol_set_atomic_masses(self)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_param
   implicit none
   class(TMolecule),intent(inout) :: self  !< molecular structure information
   self%atmass = atomic_mass(self%at)
end subroutine mol_set_atomic_masses

!> wrap cartesian coordinates back into cell
!
!  This automatically done when calling @see xyz_to_abc, so we only have
!  to perform the transformation there and back again
subroutine wrap_back(self)
   use xtb_mctc_accuracy, only : wp
   use xtb_pbc_tools
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   call xyz_to_abc(self%n,self%lattice,self%xyz,self%abc,self%pbc)
   call abc_to_xyz(self%n,self%lattice,self%abc,self%xyz)
end subroutine wrap_back

pure function center_of_geometry(self) result(center)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(in) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: idir
   center = 0.0_wp

   do idir = 1, 3
      if (.not.self%pbc(idir)) &
         center(idir) = sum(self%xyz(idir,:))
   enddo
   center = center/real(self%n,wp)
end function center_of_geometry

pure subroutine shift_to_center_of_geometry(self)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_geometry()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_geometry

pure function molecular_mass(self) result(molmass)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(in) :: self !< molecular structure information
   real(wp) :: molmass
   molmass = sum(self%atmass)
end function molecular_mass

pure function center_of_mass(self) result(center)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(in) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: idir
   center = 0.0_wp
   do idir = 1, 3
      if (.not.self%pbc(idir)) &
         center(idir) = sum(self%atmass*self%xyz(idir,:))
   enddo
   center = center/sum(self%atmass)
end function center_of_mass

pure subroutine shift_to_center_of_mass(self)
   use xtb_mctc_accuracy, only : wp
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_mass()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_mass

pure function moments_of_inertia(self) result(moments)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : eigval3x3
   implicit none
   class(TMolecule),intent(in) :: self !< molecular structure information
   real(wp) :: moments(3)
   real(wp) :: center(3),atmass,inertia(3,3),work(9),tmp(3,3)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: iat,info
   ! currently not supported
   if (self%npbc > 0) then
      moments = -1.0_wp
      return
   endif
   center = self%center_of_mass()

   inertia(:, :) = 0.0_wp

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      inertia(1,1) = inertia(1,1) + atmass * (y2+z2)
      inertia(1,2) = inertia(1,2) - atmass * x*y
      inertia(2,2) = inertia(2,2) + atmass * (x2+z2)
      inertia(1,3) = inertia(1,3) - atmass * x*z
      inertia(2,3) = inertia(2,3) - atmass * y*z
      inertia(3,3) = inertia(3,3) + atmass * (x2+y2)
   enddo

   call eigval3x3(inertia, moments)

end function moments_of_inertia

pure function rotational_constants(self)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_la
   implicit none
   class(TMolecule),intent(in) :: self !< molecular structure information
   real(wp) :: moments(3)
   real(wp) :: rotational_constants(3)

   moments = self%moments_of_inertia()
   rotational_constants = 0.5_wp/[moments(3), moments(2), moments(1)]

end function rotational_constants

pure subroutine align_to_principal_axes(self,break_symmetry)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_math, only : eigvec3x3, matdet3x3
   use xtb_pbc_tools
   implicit none
   class(TMolecule),intent(inout) :: self !< molecular structure information
   logical, optional, intent(in)    :: break_symmetry
   real(wp) :: moments(3),det
   real(wp) :: center(3),atmass,inertia(3,3),axes(3,3)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: i,iat
   ! currently not supported
   if (self%npbc > 0) then
      return
   endif
   center = self%center_of_mass()
   inertia(:, :) = 0.0_wp
   if (present(break_symmetry)) then
      if (break_symmetry) then
         inertia(1,1) = 1.0e-10_wp
         inertia(1,2) = 2.0e-10_wp
         inertia(2,2) = 3.0e-10_wp
         inertia(1,3) = 4.0e-10_wp
         inertia(2,3) = 5.0e-10_wp
         inertia(3,3) = 6.0e-10_wp
      end if
   endif

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      inertia(1,1) = inertia(1,1) + atmass * (y2+z2)
      inertia(1,2) = inertia(1,2) - atmass * x*y
      inertia(2,2) = inertia(2,2) + atmass * (x2+z2)
      inertia(1,3) = inertia(1,3) - atmass * x*z
      inertia(2,3) = inertia(2,3) - atmass * y*z
      inertia(3,3) = inertia(3,3) + atmass * (x2+y2)
   enddo

   call eigvec3x3(inertia,moments,axes)

   det = matDet3x3(axes)
   if (det < 0) axes(:,1) = -axes(:,1)

   call coord_trafo(self%n,axes,self%xyz)

end subroutine align_to_principal_axes


end module xtb_type_molecule
