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

module xtb_type_solvent
   use xtb_mctc_accuracy, only : wp
   use xtb_solv_gbsa, only : TBorn

   implicit none
   public :: TBorn
   public :: allocate_gbsa,deallocate_gbsa
   private

contains

subroutine deallocate_gbsa(this)
   implicit none
   type(TBorn) :: this

   this%nat=0
   this%ntpair=0
   this%nnrad=0
   this%nang=0

   this%lrcut=0.0_wp
   this%srcut=0.0_wp

   if (allocated(this%angGrid)) deallocate(this%angGrid)
   if (allocated(this%angWeight)) deallocate(this%angWeight)
   if (allocated(this%vdwr))    deallocate(this%vdwr)
   if (allocated(this%rho))     deallocate(this%rho)
   if (allocated(this%svdw))    deallocate(this%svdw)
   if (allocated(this%brad))    deallocate(this%brad)
   if (allocated(this%brdr))    deallocate(this%brdr)
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
   if (allocated(this%at))      deallocate(this%at)
   if (allocated(this%hbmag))   deallocate(this%hbmag)
end subroutine deallocate_gbsa

subroutine allocate_gbsa(this,n,nang)
   type(TBorn), intent(inout) :: this
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
   allocate(this%angGrid(3, this%nAng), source = 0.0_wp)
   allocate(this%angWeight(this%nAng), source = 0.0_wp)
   ! initialize the hydrogen bonding contribution
   allocate(this%hbw(this%nat), source = 0.0_wp)
   allocate(this%hbmag(this%nat), source = 0.0_wp)
   allocate(this%dhbdw(this%nat), source = 0.0_wp)
   ! initialize the salt term
   allocate(this%ionscr(this%nat), source = 0.0_wp)
   allocate(this%discr(this%nat), source = 0.0_wp)

end subroutine allocate_gbsa

end module xtb_type_solvent
