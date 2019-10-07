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

!> defines the Wigner--Seitz-cell (wsc) data type
!
!  The Wigner--Seitz cell is used to define periodic boundary conditions
!  by the cyclic cluster model (CCM). This type is usually bound to
!  the molecule class but can in principle be used independently.
module tbdef_wsc
   use iso_fortran_env, only : wp => real64
   implicit none

   public :: tb_wsc

   private

!> definition of the the Wigner--Seitz-cell (wsc) data type
   type :: tb_wsc
      integer  :: n = 0       !< number of atoms in the WSC
      integer  :: cells = 0   !< number of cells used to generate WSC
      integer  :: rep(3) = 0  !< translations defining the number of cells
      real(wp) :: lattice(3,3) = 0.0_wp    !< lattice parameters
      integer, allocatable :: at(:,:)      !< define species
      integer, allocatable :: lattr(:,:,:,:) !< lattice translation
      real(wp),allocatable :: w(:,:)       !< define weights
      integer, allocatable :: itbl(:,:)    !< define index table
   contains
   procedure :: allocate => allocate_wsc
   procedure :: deallocate => deallocate_wsc
   procedure :: write => write_wsc
   end type tb_wsc

contains

!> constructor for Wigner--Seitz cell
subroutine allocate_wsc(self,n,rep)
   implicit none
   class(tb_wsc),intent(inout) :: self
   integer,intent(in) :: n      !< number of atoms
   integer,intent(in) :: rep(3) !< translations
   integer :: cells
   cells = product(2*rep+1)
   self%n       = n
   self%rep     = rep
   self%cells   = cells
   call self%deallocate
   allocate( self%at(n,n),            source = 0 )
   allocate( self%lattr(3,cells,n,n), source = 0 )
   allocate( self%w(n,n),             source = 0.0_wp )
   allocate( self%itbl(n,n),          source = 0 )
end subroutine allocate_wsc

!> @brief deconstructor for Wigner--Seitz cell
subroutine deallocate_wsc(self)
   implicit none
   class(tb_wsc),intent(inout) :: self
   if (allocated(self%at))   deallocate(self%at)
   if (allocated(self%lattr))deallocate(self%lattr)
   if (allocated(self%w))    deallocate(self%w)
   if (allocated(self%itbl)) deallocate(self%itbl)
end subroutine deallocate_wsc

subroutine write_wsc(self,iunit,comment)
   implicit none
   class(tb_wsc),   intent(in) :: self
   integer,         intent(in) :: iunit
   character(len=*),intent(in) :: comment
   character(len=*),parameter :: dfmt = '(1x,a,1x,"=",1x,g0)'

   write(iunit,'(72(">"))')
   write(iunit,'(1x,"*",1x,a)') "Writing 'tb_wsc' class"
   write(iunit,'(  "->",1x,a)') comment
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "status of the fields"
   write(iunit,dfmt) "integer :: n           ",self%n
   write(iunit,dfmt) "integer :: cells       ",self%cells
   write(iunit,dfmt) "integer :: rep(1)      ",self%rep(1)
   write(iunit,dfmt) "        &  rep(2)      ",self%rep(2)
   write(iunit,dfmt) "        &  rep(3)      ",self%rep(3)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "allocation status"
   write(iunit,dfmt) "allocated? at(:)       ",allocated(self%at)
   write(iunit,dfmt) "allocated? lattr(:,:,:,:)",allocated(self%lattr)
   write(iunit,dfmt) "allocated? w(:,:)      ",allocated(self%w)
   write(iunit,dfmt) "allocated? itbl(:,:)   ",allocated(self%itbl)
   write(iunit,'(72("-"))')
   write(iunit,'(1x,"*",1x,a)') "size of memory allocation"
   if (allocated(self%at)) then
   write(iunit,dfmt) "size(1) :: at(*,:)     ",size(self%at,1)
   write(iunit,dfmt) "size(2) :: at(:,*)     ",size(self%at,2)
   endif
   if (allocated(self%lattr)) then
   write(iunit,dfmt) "size(1) :: lattr(*,:,:,:)",size(self%lattr,1)
   write(iunit,dfmt) "size(2) :: lattr(:,*,:,:)",size(self%lattr,2)
   write(iunit,dfmt) "size(3) :: lattr(:,:,*,:)",size(self%lattr,3)
   write(iunit,dfmt) "size(4) :: lattr(:,:,:,*)",size(self%lattr,4)
   endif
   if (allocated(self%w)) then
   write(iunit,dfmt) "size(1) :: w(*,:)      ",size(self%w,1)
   write(iunit,dfmt) "size(2) :: w(:,*)      ",size(self%w,2)
   endif
   if (allocated(self%w)) then
   write(iunit,dfmt) "size(1) :: itbl(*,:)   ",size(self%itbl,1)
   write(iunit,dfmt) "size(2) :: itbl(:,*)   ",size(self%itbl,2)
   endif
   write(iunit,'(72("<"))')

end subroutine write_wsc

end module tbdef_wsc
