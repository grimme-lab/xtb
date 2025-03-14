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

!> Topological data for force field type calculations
module xtb_gfnff_topology
   use xtb_mctc_accuracy, only : wp
   use xtb_type_dispersionmodel, only : TDispersionModel
   implicit none
   private

   public :: TGFFTopology, TPrintTopo


   !> Topology information for a given system
   type :: TGFFTopology

      !number of terms
      integer  :: nbond
      integer  :: nangl
      integer  :: ntors
      integer  :: nstors
      integer  :: nathbH
      integer  :: nathbAB
      integer  :: natxbAB
      integer  :: nbatm
      integer  :: nfrag
      integer  :: maxsystem   ! max. number of fragmentsfor hessian
      integer  :: bond_hb_nr  ! number of unique AH...B HB/bond terms
      integer  :: b_max      ! number of B atoms per unique AH bond

      !numbers that are rewritten, so must be stored for allocation
      integer  :: nbond_blist
      integer  :: nbond_vbond
      integer  :: nangl_alloc
      integer  :: ntors_alloc

      !file type read
      integer  :: read_file_type

      !lists
      integer,allocatable ::     hyb(:)   ! hybridization of every atom
      integer,allocatable ::  blist(:,:)   ! bonded atoms
      integer,allocatable ::  alist(:,:)   ! angles
      integer,allocatable ::  tlist(:,:)   ! torsions
      integer,allocatable :: b3list(:,:)   ! bond atm
      integer,allocatable :: sTorsl(:,:)
      !-----------------------------------------------
      integer,allocatable :: nr_hb(:)      ! Nr. of H bonds per O-H or N-H bond
      integer,allocatable :: bond_hb_AH(:,:) ! A, H atoms in bonds that are also part of HBs
      integer,allocatable :: bond_hb_B(:,:,:)  ! B atoms in bonds that are also part of HBs
      integer,allocatable :: bond_hb_Bn(:)   ! Nr. of B atoms for one AH bond pair
      !-----------------------------------------------
      integer,allocatable :: hbatABl(:,:)  ! AB atoms for HB
      integer,allocatable :: xbatABl(:,:)  ! AB atoms for XB
      integer,allocatable :: hbatHl (:,:)    ! H  atoms for HB
      integer,allocatable :: fraglist(:)   ! atoms in molecular fragments (for EEQ)
      integer,allocatable :: qpdb  (:)     ! atomic charge in residues from PDB file

      !potential parameters used in energy-gradient routine
      real(wp),allocatable:: vbond(:,:)    ! bonds
      real(wp),allocatable:: vangl(:,:)    ! angles
      real(wp),allocatable:: vtors(:,:)    ! torsions
      real(wp),allocatable:: chieeq(:)     ! atomic ENs for EEQ
      real(wp),allocatable:: gameeq(:)     ! atomic gamma for EEQ
      real(wp),allocatable:: alpeeq(:)     ! atomic alpha for EEQ, squared
      real(wp),allocatable:: alphanb(:,:,:)    ! non-bonded exponent for atom pairs
      real(wp),allocatable::    qa(:)      ! estimated atomic charges (fixed and obtained from topology EEQ)
      real(wp),allocatable::    xyze0(:,:) ! atom xyz, starting geom. (for Efield energy)
      real(wp),allocatable:: zetac6(:)     ! D4 scaling factor product
      real(wp),allocatable:: qfrag (:)     ! fragment charge (for EEQ)
      real(wp),allocatable:: hbbas (:)     ! HB donor atom basicity
      real(wp),allocatable:: hbaci (:)     ! HB acceptor atom acidity
      integer, allocatable:: hb_mapABH(:)  ! mapping of indices from all atoms to only AB and H separately
      logical, allocatable:: isABH(:)      ! logical set to true if the atom is part of a hydrogen bond
      integer :: hb_mapNAB                 ! number of AB atoms that are part of a hydrogen bond
      integer :: hb_mapNH                  ! number of H atoms that are part of a hydrogen bond

      integer, allocatable  :: ispinsyst(:,:)
      integer, allocatable  :: nspinsyst(:)
      integer               :: nsystem

      type(TDispersionModel) :: dispm



   contains

      procedure :: zero

   end type TGFFTopology

   ! logicals for GFN-FF topology list printout
   type :: TPrintTopo
     logical :: etot    = .false.
     logical :: gnorm   = .false.
     logical :: nb      = .false.
     logical :: bpair   = .false.
     logical :: alist   = .false.
     logical :: blist   = .false.
     logical :: tlist   = .false.
     logical :: vtors   = .false.
     logical :: vbond   = .false.
     logical :: vangl   = .false.
     logical :: hbbond  = .false.
     logical :: eeq     = .false.
     logical :: warning = .false.

   contains

   procedure :: any

   end type TPrintTopo

contains


subroutine zero(self)
   class(TGFFTopology), intent(out) :: self

   self%nbond = 0
   self%nangl = 0
   self%ntors = 0
   self%nathbH = 0
   self%nathbAB = 0
   self%natxbAB = 0
   self%nbatm = 0
   self%nfrag = 0
   self%maxsystem = 0
   self%bond_hb_nr = 0
   self%b_max = 0

   self%nbond_blist = 0
   self%nbond_vbond = 0
   self%nangl_alloc = 0
   self%ntors_alloc = 0

   self%read_file_type = 0

end subroutine zero

function any(self) result(tf)
  class(TPrintTopo), intent(in) :: self
  logical :: tf

  tf = self%nb.or.self%bpair.or.self%alist.or.self%blist.or. &
     & self%tlist.or.self%vtors.or.self%vbond.or.self%vangl.or. &
     & self%hbbond.or.self%eeq.or.self%etot.or.self%gnorm
end function any

end module xtb_gfnff_topology
