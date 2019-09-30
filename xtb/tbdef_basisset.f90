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

module tbdef_basisset
   use iso_fortran_env, only : wp => real64
   implicit none

   public :: tb_basisset

   private

   type :: tb_basisset
      integer  :: maxao = 0
      integer  :: n = 0
      integer  :: nbf = 0
      integer  :: nao = 0
      integer  :: nshell = 0
      integer, allocatable :: ash(:)
      integer, allocatable :: lsh(:)
      integer, allocatable :: ao2sh(:)
      integer, allocatable :: primcount(:)
      integer, allocatable :: caoshell(:,:)
      integer, allocatable :: saoshell(:,:)
      integer, allocatable :: lao(:)
      integer, allocatable :: nprim(:)
      integer, allocatable :: aoat(:)
      integer, allocatable :: fila(:,:)
      integer, allocatable :: lao2(:)
      integer, allocatable :: aoat2(:)
      integer, allocatable :: fila2(:,:)
      integer, allocatable :: valao(:)
      integer, allocatable :: valao2(:)
      real(wp),allocatable :: alp(:)
      real(wp),allocatable :: cont(:)
      real(wp),allocatable :: hdiag(:)
      real(wp),allocatable :: hdiag2(:)
      real(wp),allocatable :: aoexp(:)
   contains
   procedure :: allocate => allocate_basisset
   procedure :: deallocate => deallocate_basisset
   end type tb_basisset

contains

subroutine allocate_basisset(self,n,nbf,nao,nshell)
   implicit none
   class(tb_basisset),intent(inout) :: self
   integer,intent(in) :: n,nbf,nao,nshell
   self%n=n
   self%nbf=nbf
   self%nao=nao
   self%nshell=nshell
   call self%deallocate
   allocate( self%hdiag(nbf),     source = 0.0_wp )
   allocate( self%alp(6*nbf),     source = 0.0_wp )
   allocate( self%cont(6*nbf),    source = 0.0_wp )
   allocate( self%hdiag2(nao),    source = 0.0_wp )
   allocate( self%aoexp(nao),     source = 0.0_wp )
   allocate( self%ash(nao),       source = 0 )
   allocate( self%lsh(nao),       source = 0 )
   allocate( self%ao2sh(nao),     source = 0 )
   allocate( self%nprim(nbf),     source = 0 )
   allocate( self%primcount(nbf), source = 0 )
   allocate( self%caoshell(5,n),  source = 0 )
   allocate( self%saoshell(5,n),  source = 0 )
   allocate( self%fila(2,n),      source = 0 )
   allocate( self%fila2(2,n),     source = 0 )
   allocate( self%lao(nbf),       source = 0 )
   allocate( self%aoat(nbf),      source = 0 )
   allocate( self%valao(nbf),     source = 0 )
   allocate( self%lao2(nao),      source = 0 )
   allocate( self%aoat2(nao),     source = 0 )
   allocate( self%valao2(nbf),    source = 0 )
end subroutine allocate_basisset

subroutine deallocate_basisset(self)
   implicit none
   class(tb_basisset),intent(inout) :: self
   if(allocated(self%hdiag))     deallocate(self%hdiag)
   if(allocated(self%alp))       deallocate(self%alp)
   if(allocated(self%cont))      deallocate(self%cont)
   if(allocated(self%hdiag2))    deallocate(self%hdiag2)
   if(allocated(self%aoexp))     deallocate(self%aoexp)
   if(allocated(self%ash))       deallocate(self%ash)
   if(allocated(self%lsh))       deallocate(self%lsh)
   if(allocated(self%ao2sh))     deallocate(self%ao2sh)
   if(allocated(self%nprim))     deallocate(self%nprim)
   if(allocated(self%primcount)) deallocate(self%primcount)
   if(allocated(self%caoshell))  deallocate(self%caoshell)
   if(allocated(self%saoshell))  deallocate(self%saoshell)
   if(allocated(self%fila))      deallocate(self%fila)
   if(allocated(self%fila2))     deallocate(self%fila2)
   if(allocated(self%lao))       deallocate(self%lao)
   if(allocated(self%lao2))      deallocate(self%lao2)
   if(allocated(self%aoat))      deallocate(self%aoat2)
   if(allocated(self%valao))     deallocate(self%valao)
   if(allocated(self%valao2))    deallocate(self%valao2)
end subroutine deallocate_basisset

end module tbdef_basisset
