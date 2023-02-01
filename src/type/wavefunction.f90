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

module xtb_type_wavefunction
   use xtb_mctc_accuracy, only : wp
   implicit none

   public :: TWavefunction

   private

   type :: TWavefunction
      integer :: n = 0
         !! Number of atoms 
      integer :: nel = 0
         !! Number of elctrons
      integer :: nopen = 0
         !! Number of unpaired electrons
      integer :: nao = 0
         !! Number of atomic orbitals
      integer :: nshell = 0
         !! Number of shells
      real(wp),allocatable :: P(:,:)    
         !! Density matrix
      real(wp),allocatable :: q(:)      
         !! Partial charges
      real(wp),allocatable :: qsh(:)    
         !! Shell charges
      real(wp),allocatable :: dipm(:,:) 
         !! Dipole moments
      real(wp),allocatable :: qp(:,:)   
         !! Quadrupole moments
      real(wp),allocatable :: wbo(:,:)  
         !! Wiberg bond orders
      integer :: ihomo = 0,ihomoa = 0,ihomob = 0 
         !! HOMO position
      real(wp) :: efa = 0.0_wp, efb = 0.0_wp
      real(wp),allocatable :: focc(:)   
         !! Fractional occupation
      real(wp),allocatable :: focca(:)  
         !! For alpha space
      real(wp),allocatable :: foccb(:)  
         !! For beta space
      real(wp),allocatable :: emo(:)    
         !! Orbital energies
      real(wp),allocatable :: C(:,:)    
         !! Molecular orbitals
   
   contains
   procedure :: allocate => allocate_wavefunction
   procedure :: deallocate => deallocate_wavefunction
   end type TWavefunction

contains

subroutine allocate_wavefunction(self,n,nshell,nao)
   class(TWavefunction),intent(inout) :: self
   integer,intent(in) :: n,nshell,nao
   self%n = n
   self%nshell = nshell
   self%nao = nao
   self%ihomo = 0
   self%ihomoa = 0
   self%ihomob = 0
   call self%deallocate
   allocate( self%P(nao,nao),  source = 0.0_wp )
   allocate( self%q(n),        source = 0.0_wp )
   allocate( self%qsh(nshell), source = 0.0_wp )
   allocate( self%dipm(3,n),   source = 0.0_wp )
   allocate( self%qp(6,n),     source = 0.0_wp )
   allocate( self%wbo(n,n),    source = 0.0_wp )
   allocate( self%focc(nao),   source = 0.0_wp )
   allocate( self%focca(nao),  source = 0.0_wp )
   allocate( self%foccb(nao),  source = 0.0_wp )
   allocate( self%emo(nao),    source = 0.0_wp )
   allocate( self%C(nao,nao),  source = 0.0_wp )
end subroutine allocate_wavefunction

subroutine deallocate_wavefunction(self)
   class(TWavefunction),intent(inout) :: self
   if(allocated(self%P))    deallocate(self%P)
   if(allocated(self%q))    deallocate(self%q)
   if(allocated(self%qsh))  deallocate(self%qsh)
   if(allocated(self%dipm)) deallocate(self%dipm)
   if(allocated(self%qp))   deallocate(self%qp)
   if(allocated(self%wbo))  deallocate(self%wbo)
   if(allocated(self%focc)) deallocate(self%focc)
   if(allocated(self%focca))deallocate(self%focca)
   if(allocated(self%foccb))deallocate(self%foccb)
   if(allocated(self%emo))  deallocate(self%emo)
   if(allocated(self%C))    deallocate(self%C)
end subroutine deallocate_wavefunction

end module xtb_type_wavefunction
