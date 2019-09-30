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

module tbdef_solvent
   use iso_fortran_env, wp => real64

   implicit none
   public :: tb_solvent
   public :: allocate_gbsa,deallocate_gbsa
   private

! ========================================================================
!  GBSA class: contains all molecule specific information
   type :: tb_solvent
!     number of atoms
      integer  :: nat
!     atom types
      integer, allocatable :: at(:)
!     number of pairs
      integer  :: ntpair
!     number of angular grid points
      integer :: nang 
!     angular grid
      real(wp), allocatable :: grida(:,:)
! ------------------------------------------------------------------------
!     van der Waals radii of the particles
      real(wp),allocatable :: vdwr(:)
!     greatest van der Waals radius
      real(wp) :: maxvdwr
! ------------------------------------------------------------------------
!     pair descreening approximation radii
      real(wp),allocatable :: rho(:)
!     offset van der Waals radii
      real(wp),allocatable :: svdw(:)
! ------------------------------------------------------------------------
!     Neighbor list:
!     cut-off radius for the Born radius NN list
      real(wp) :: lrcut
!     cut-off radius for the SASA NN list
      real(wp) :: srcut
!     number of neighbors for Born radii
      integer  :: nnrad
!     number of neighbors for SASA computation
      integer, allocatable :: nnsas(:)
!     neighbors of an atom for Born radii
      integer, allocatable :: nnlistr(:,:)
!     neighbors of an atom for SASA
      integer, allocatable :: nnlists(:,:)
!     all pairs indeces array
      integer, allocatable :: ppind(:,:)
!     all pairs vector differences and magnitudes array
      real(wp),allocatable :: ddpair(:,:)
! ------------------------------------------------------------------------
!     Atom specific surface data
      real(wp),allocatable :: vdwsa(:)
      real(wp),allocatable :: wrp(:)
      real(wp),allocatable :: trj2(:,:)
! ------------------------------------------------------------------------
!     Dielectric data
      real(wp) :: gborn
! ------------------------------------------------------------------------
!     Born radii
      real(wp),allocatable :: brad(:)
! ------------------------------------------------------------------------
!     Salt screening
      real(wp),allocatable :: ionscr(:)
      real(wp),allocatable :: discr(:)
! ------------------------------------------------------------------------
!     Atomic surfaces
      real(wp) :: gsasa
      real(wp) :: sasagam
      real(wp),allocatable :: gamsasa(:)
      real(wp),allocatable :: sasa(:)
! ------------------------------------------------------------------------
!     Hydrogen bond contribution
      real(wp) :: ghb
      real(wp),allocatable :: hbw(:)
! ------------------------------------------------------------------------
!  Gradient:
! ------------------------------------------------------------------------
!     Born radii gradient
      real(wp),allocatable :: brdr(:,:,:)
! ------------------------------------------------------------------------
!     Molecular Surface gradient
      real(wp),allocatable :: dsdr(:,:)
      real(wp),allocatable :: dsdrt(:,:,:)
! ------------------------------------------------------------------------
!     Hydrogen bond gradient
      real(wp),allocatable :: dhbdw(:)
! ------------------------------------------------------------------------
!     GB energy gradient
      real(wp),allocatable :: dbrdp(:)
   end type tb_solvent

contains

subroutine deallocate_gbsa(this)
   implicit none
   type(tb_solvent) :: this

   this%gborn=0.0_wp
   this%ghb=0.0_wp

   this%nat=0
   this%ntpair=0
   this%nnrad=0
   this%nang=0

   this%lrcut=0.0_wp
   this%srcut=0.0_wp
   this%maxvdwr=0.0_wp

   if (allocated(this%grida))   deallocate(this%grida)
   if (allocated(this%vdwr))    deallocate(this%vdwr)
   if (allocated(this%rho))     deallocate(this%rho)
   if (allocated(this%svdw))    deallocate(this%svdw)
   if (allocated(this%brad))    deallocate(this%brad)
   if (allocated(this%brdr))    deallocate(this%brdr)
   if (allocated(this%dbrdp))   deallocate(this%dbrdp)
   if (allocated(this%nnlistr)) deallocate(this%nnlistr)
   if (allocated(this%nnsas))   deallocate(this%nnsas)
   if (allocated(this%nnlists)) deallocate(this%nnlists)
   if (allocated(this%ppind))   deallocate(this%ppind)
   if (allocated(this%ddpair))  deallocate(this%ddpair)
   if (allocated(this%vdwsa))   deallocate(this%vdwsa)
   if (allocated(this%wrp))     deallocate(this%wrp)
   if (allocated(this%trj2))    deallocate(this%trj2)
   if (allocated(this%sasa))    deallocate(this%sasa)
   if (allocated(this%gamsasa)) deallocate(this%gamsasa)
   if (allocated(this%dsdr))    deallocate(this%dsdr)
   if (allocated(this%dsdrt))   deallocate(this%dsdrt)
   if (allocated(this%hbw))     deallocate(this%hbw)
   if (allocated(this%dhbdw))   deallocate(this%dhbdw)
   if (allocated(this%ionscr))  deallocate(this%ionscr)
   if (allocated(this%discr))   deallocate(this%discr)
end subroutine deallocate_gbsa

subroutine allocate_gbsa(this,n,nang)
   type(tb_solvent), intent(inout) :: this
   integer, intent(in) :: n
   integer, intent(in) :: nang

   call deallocate_gbsa(this)

   this%nat = n
   this%ntpair=n*(n-1)/2
   this%nang = nang

   ! initialize the vdw radii array
   allocate(this%vdwr(n), source = 0.0_wp)
   allocate(this%rho(n), source = 0.0_wp)
   allocate(this%svdw(n), source = 0.0_wp)
   allocate(this%at(n), source = 0)
   ! initialize Born radii
   allocate(this%brad(this%nat), source = 0.0_wp)
   allocate(this%brdr(3,this%nat,this%nat), source = 0.0_wp)
   allocate(this%dbrdp(this%nat), source = 0.0_wp)
   allocate(this%nnlistr(3,this%ntpair), source = 0)
   allocate(this%nnsas(this%nat), source = 0)
   allocate(this%nnlists(this%nat,this%nat), source = 0)
   allocate(this%ppind(2,this%ntpair), source = 0)
   allocate(this%ddpair(4,this%ntpair), source = 0.0_wp)
   ! initialize solvent-accessible atomic surface area computation (SASA)
   allocate(this%vdwsa(this%nat), source = 0.0_wp)
   allocate(this%wrp(this%nat), source = 0.0_wp)
   allocate(this%trj2(2,this%nat), source = 0.0_wp)
   allocate(this%sasa(this%nat), source = 0.0_wp)
   allocate(this%gamsasa(this%nat), source = 0.0_wp)
   allocate(this%dsdr(3,this%nat), source = 0.0_wp)
   allocate(this%dsdrt(3,this%nat,this%nat), source = 0.0_wp)
   allocate(this%grida(4,this%nang), source = 0.0_wp)
   ! initialize the hydrogen bonding contribution
   allocate(this%hbw(this%nat), source = 0.0_wp)
   allocate(this%dhbdw(this%nat), source = 0.0_wp)
   ! initialize the salt term
   allocate(this%ionscr(this%nat), source = 0.0_wp)
   allocate(this%discr(this%nat), source = 0.0_wp)

end subroutine allocate_gbsa

end module tbdef_solvent
