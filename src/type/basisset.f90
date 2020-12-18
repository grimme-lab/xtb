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

module xtb_type_basisset
   use xtb_mctc_accuracy, only : wp
   implicit none

   public :: TBasisset

   private

   type :: TBasisset
      integer  :: maxao = 0
      integer  :: n = 0
      integer  :: nbf = 0
      integer  :: nao = 0
      integer  :: nshell = 0
      ! ------------------------------------------------------------------------
      !  This variables describe the ATOMS
      integer, allocatable :: caoshell(:,:)  ! atom number -> basis function
      integer, allocatable :: saoshell(:,:)  ! atom number -> AO
      integer, allocatable :: fila(:,:)      ! atom number -> basis function
      integer, allocatable :: fila2(:,:)     ! atom number -> AO
      integer, allocatable :: shells(:,:)    ! atom number -> shell
      ! ------------------------------------------------------------------------
      !  This variables describe the indivdual SHELLS
      integer, allocatable :: lsh(:)   ! shell -> azimudal quantum number
      integer, allocatable :: ash(:)   ! shell -> atom number
      real(wp),allocatable :: zeta(:)  ! shell -> exponent
      real(wp),allocatable :: level(:) ! shell -> level energy
      real(wp),allocatable :: minalp(:)! shell -> most diffuse exponent
      integer, allocatable :: sh2bf(:,:) ! shell -> BF
      integer, allocatable :: sh2ao(:,:) ! shell -> AO
      integer, allocatable :: valsh(:) ! shell -> shell type
      ! ------------------------------------------------------------------------
      !  This variables describe the SPHERICAL atomic orbitals
      real(wp),allocatable :: aoexp(:)    ! AO -> exponent
      integer, allocatable :: ao2sh(:)    ! AO -> shell
      integer, allocatable :: lao2(:)     ! AO -> azimudal quantum number
      integer, allocatable :: aoat2(:)    ! AO -> atom number
      real(wp),allocatable :: hdiag2(:)   ! AO -> level energy
      integer, allocatable :: valao2(:)   ! AO -> AO type
      ! ------------------------------------------------------------------------
      !  This variables describe the CARTESIAN basis functions
      integer, allocatable :: primcount(:) ! BF -> primitive count
      real(wp),allocatable :: hdiag(:)     ! BF -> level energy
      integer, allocatable :: nprim(:)     ! BF -> primitive number
      integer, allocatable :: aoat(:)      ! BF -> atom
      integer, allocatable :: valao(:)     ! BF -> BF type
      integer, allocatable :: lao(:)       ! BF -> azimudal quantum number
      ! ------------------------------------------------------------------------
      !  This variables decribe the PRIMITIV basis functions
      real(wp),allocatable :: alp(:)  ! primitive -> primitive exponent
      real(wp),allocatable :: cont(:) ! primitive -> contraction coeffient
   contains
   procedure :: allocate => allocate_basisset
   procedure :: deallocate => deallocate_basisset
   end type TBasisset

contains

subroutine allocate_basisset(self,n,nbf,nao,nshell)
   implicit none
   class(TBasisset),intent(inout) :: self
   integer,intent(in) :: n,nbf,nao,nshell
   self%n=n
   self%nbf=nbf
   self%nao=nao
   self%nshell=nshell
   call self%deallocate
   allocate( self%shells(2,n),    source = 0 )
   allocate( self%sh2ao(2,nshell),source = 0 )
   allocate( self%sh2bf(2,nshell),source = 0 )
   allocate( self%minalp(nshell), source = 0.0_wp )
   allocate( self%level(nshell),  source = 0.0_wp )
   allocate( self%zeta(nshell),   source = 0.0_wp )
   allocate( self%valsh(nshell),  source = 0 )
   allocate( self%hdiag(nbf),     source = 0.0_wp )
   allocate( self%alp(9*nbf),     source = 0.0_wp )
   allocate( self%cont(9*nbf),    source = 0.0_wp )
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
   class(TBasisset),intent(inout) :: self
   if(allocated(self%shells))    deallocate(self%shells)
   if(allocated(self%sh2ao))     deallocate(self%sh2ao)
   if(allocated(self%sh2bf))     deallocate(self%sh2bf)
   if(allocated(self%minalp))    deallocate(self%minalp)
   if(allocated(self%valsh))     deallocate(self%valsh)
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
   if(allocated(self%zeta))      deallocate(self%zeta)
   if(allocated(self%level))     deallocate(self%level)
end subroutine deallocate_basisset

end module xtb_type_basisset
