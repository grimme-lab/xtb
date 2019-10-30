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
module tbdef_molecule
   use iso_fortran_env, only : wp => real64
   use tbdef_wsc
   use tbdef_topology
   use tbdef_fragments
   use tbdef_buffer
   implicit none

   public :: tb_molecule
   public :: len, size

   private


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

   !> molecular structure information
   type :: tb_molecule
      integer  :: n = 0            !< number of atoms
      real(wp) :: chrg = 0.0_wp    !< total charge
      integer  :: uhf = 0          !< number of unpaired electrons
      logical  :: pbc(3) = .false. !< periodic dimensions
      integer  :: npbc = 0         !< periodicity of system
      character(len=2),allocatable :: sym(:) !< element symbols
      integer, allocatable :: at(:)          !< ordinal numbers
      real(wp),allocatable :: xyz(:,:)       !< cartesian coordinates in bohr
      real(wp),allocatable :: abc(:,:)       !< fractional coordinates
      real(wp),allocatable :: dist(:,:)      !< interatomic distances
      real(wp),allocatable :: atmass(:)      !< atomic masses in amu
      real(wp),allocatable :: z(:)           !< nuclear charges
      real(wp),allocatable :: cn(:)          !< coordination number
      real(wp) :: cellpar(6) = 0.0_wp        !< cell parameters
      real(wp) :: lattice(3,3) = 0.0_wp      !< direct lattice parameters
      real(wp) :: rec_lat(3,3) = 0.0_wp      !< reciprocal lattice parameters
      real(wp) :: volume = 0.0_wp            !< volume of unit cell
      type(tb_wsc) :: wsc                    !< Wigner--Seitz cell
      integer  :: ftype = 0                  !< file type of the input
      character(len=:), allocatable :: name
      type(tb_topology) :: bonds
      type(tb_fragments) :: frag
      !> PDB specific information about residues and chains
      type(pdb_data), allocatable :: pdb(:)
      !> SDF specific information about atom types
      type(sdf_data), allocatable :: sdf(:)
      !> VASP specific information about input type
      type(vasp_info) :: vasp = vasp_info()
      !> Turbomole specific information about input type
      type(turbo_info) :: turbo = turbo_info()
      !> raw input buffer
      type(tb_buffer) :: info
   contains
      procedure :: allocate => allocate_molecule
      procedure :: deallocate => deallocate_molecule
      procedure :: calculate_distances => mol_calculate_distances
      procedure :: set_nuclear_charge => mol_set_nuclear_charge
      procedure :: set_atomic_masses => mol_set_atomic_masses
      procedure :: write => write_molecule_generic
      procedure :: read => read_molecule_generic
      procedure :: update => mol_update
      procedure :: wrap_back => mol_wrap_back
      procedure :: center_of_geometry
      procedure :: center_of_mass
      procedure :: shift_to_center_of_geometry
      procedure :: shift_to_center_of_mass
      procedure :: moments_of_inertia
      procedure :: align_to_principal_axes
   end type tb_molecule


   interface
      module subroutine write_molecule_generic(self, unit, format, &
            &                                  energy, gnorm, number)
         class(tb_molecule), intent(in) :: self
         integer, intent(in) :: unit
         integer, intent(in), optional :: format
         real(wp), intent(in), optional :: energy
         real(wp), intent(in), optional :: gnorm
         integer, intent(in), optional :: number
      end subroutine write_molecule_generic

      module subroutine read_molecule_generic(self, unit, format)
         class(tb_molecule), intent(out) :: self
         integer, intent(in) :: unit
         integer, intent(in), optional :: format
      end subroutine read_molecule_generic
   end interface


   interface len
      module procedure :: mol_length
   end interface


contains


!> obtain number of atoms for molecular structure
!
!  The molecular structure is assumed to be well-behaved in this respect
!  so there is no sanity check on the allocation status.
integer pure elemental function mol_length(self) result(length)
   class(tb_molecule),intent(in) :: self !< molecular structure information
   length = self%n
end function mol_length


!> constructor for molecular structure
subroutine allocate_molecule(self,n)
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   integer,intent(in) :: n
   call self%deallocate
   self%n = n
   allocate( self%at(n),          source = 0 )
   allocate( self%sym(n),         source = '  ' )
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
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   self%n = 0
   self%pbc = .false.
   self%chrg = 0.0_wp
   self%uhf = 0
   self%lattice = 0.0_wp
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
   call self%info%deallocate
end subroutine deallocate_molecule

subroutine mol_update(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(tb_molecule),intent(inout) :: self  !< molecular structure information

   if (self%npbc > 0) then
      call dlat_to_cell(self%lattice,self%cellpar)
      call dlat_to_rlat(self%lattice,self%rec_lat)
      self%volume = dlat_to_dvol(self%lattice)

      call self%wrap_back
   endif

   call self%calculate_distances

end subroutine mol_update

!> calculates all distances for molecular structures and minimum
!  image distances for peridic structures
subroutine mol_calculate_distances(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
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
            self%dist(j,i) = norm2(self%xyz(:,j)-self%xyz(:,i))
            self%dist(i,j) = self%dist(j,i)
         enddo
         self%dist(i,i) = 0.0_wp
      enddo
   endif
end subroutine mol_calculate_distances

!> get all nuclear charges
subroutine mol_set_nuclear_charge(self)
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(inout) :: self  !< molecular structure information
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
  endif
end function ncore
end subroutine mol_set_nuclear_charge

!> get all nuclear charges
subroutine mol_set_atomic_masses(self)
   use iso_fortran_env, wp => real64
   use mctc_param
   implicit none
   class(tb_molecule),intent(inout) :: self  !< molecular structure information
   self%atmass = atomic_mass(self%at)
end subroutine mol_set_atomic_masses

!> wrap cartesian coordinates back into cell
! 
!  This automatically done when calling @see xyz_to_abc, so we only have
!  to perform the transformation there and back again
subroutine mol_wrap_back(self)
   use iso_fortran_env, wp => real64
   use pbc_tools
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   call xyz_to_abc(self%n,self%lattice,self%xyz,self%abc,self%pbc)
   call abc_to_xyz(self%n,self%lattice,self%abc,self%xyz)
end subroutine mol_wrap_back

pure function center_of_geometry(self) result(center)
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(in) :: self !< molecular structure information
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
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_geometry()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_geometry

pure function molecular_mass(self) result(molmass)
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(in) :: self !< molecular structure information
   real(wp) :: molmass
   molmass = sum(self%atmass)
end function molecular_mass

pure function center_of_mass(self) result(center)
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(in) :: self !< molecular structure information
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
   use iso_fortran_env, wp => real64
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   real(wp) :: center(3)
   integer  :: iat
   center = self%center_of_mass()
   do iat = 1, self%n
      self%xyz(:,iat) = self%xyz(:,iat) - center
   enddo
end subroutine shift_to_center_of_mass

pure function moments_of_inertia(self) result(moments)
   use iso_fortran_env, wp => real64
   use mctc_la
   implicit none
   class(tb_molecule),intent(in) :: self !< molecular structure information
   real(wp) :: moments(3)
   real(wp) :: center(3),atmass,t(6),work(9)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: iat,info
   ! currently not supported
   if (self%npbc > 0) then
      moments = -1.0_wp
      return
   endif
   center = self%center_of_mass()

   t = 0.0_wp

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      t(1) = t(1) + atmass * (y2+z2)
      t(2) = t(2) - atmass * x*y
      t(3) = t(3) + atmass * (x2+z2)
      t(4) = t(4) - atmass * x*z
      t(5) = t(5) - atmass * y*z
      t(6) = t(6) + atmass * (x2+y2)
   enddo

   call dspev('N','U',3,t,moments,work,3,work,info)
   if (info.ne.0) moments = -1.0_wp

end function moments_of_inertia

pure function rotational_constants(self)
   use iso_fortran_env, wp => real64
   use mctc_la
   implicit none
   class(tb_molecule),intent(in) :: self !< molecular structure information
   real(wp) :: moments(3)
   real(wp) :: rotational_constants(3)

   moments = self%moments_of_inertia()
   rotational_constants = 0.5_wp/[moments(3), moments(2), moments(1)]

end function rotational_constants

pure subroutine align_to_principal_axes(self,break_symmetry)
   use iso_fortran_env, wp => real64
   use mctc_la
   use pbc_tools
   implicit none
   class(tb_molecule),intent(inout) :: self !< molecular structure information
   logical, optional, intent(in)    :: break_symmetry
   real(wp) :: moments(3),det
   real(wp) :: center(3),atmass,t(6),work(9),axes(3,3)
   real(wp) :: x,x2,y,y2,z,z2
   integer  :: i,iat,info
   ! currently not supported
   if (self%npbc > 0) then
      return
   endif
   center = self%center_of_mass()
   t = 0.0_wp
   if (present(break_symmetry)) then
      if (break_symmetry) t = [(real(i,wp)*1.0e-10_wp,i=1,6)]
   endif

   do iat = 1, self%n
      atmass = self%atmass(iat)
      x = self%xyz(1,iat)-center(1); x2 = x**2
      y = self%xyz(2,iat)-center(2); y2 = y**2
      z = self%xyz(3,iat)-center(3); z2 = z**2
      t(1) = t(1) + atmass * (y2+z2)
      t(2) = t(2) - atmass * x*y
      t(3) = t(3) + atmass * (x2+z2)
      t(4) = t(4) - atmass * x*z
      t(5) = t(5) - atmass * y*z
      t(6) = t(6) + atmass * (x2+y2)
   enddo

   call dspev('V','U',3,t,moments,axes,3,work,info)
   if (info.ne.0) return

   det = mat_det_3x3(axes)
   if (det < 0) axes(:,1) = -axes(:,1)

   call coord_trafo(self%n,axes,self%xyz)

end subroutine align_to_principal_axes

end module tbdef_molecule
