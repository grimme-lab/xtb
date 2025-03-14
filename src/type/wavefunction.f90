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
      
      !> number of atoms 
      integer :: n = 0
      
      !> number of electrons
      integer :: nel = 0
      
      !> number of unpaired electrons
      integer :: nopen = 0
      
      !> number of atomic orbitals
      integer :: nao = 0
         
      !> number of shells
      integer :: nshell = 0

      !> density matrix
      real(wp),allocatable :: P(:,:)    

      !> partial charges
      real(wp),allocatable :: q(:)      

      !> shell charges
      real(wp),allocatable :: qsh(:)    

      !> dipole moments
      real(wp),allocatable :: dipm(:,:) 
      
      !> quadrupole moments
      real(wp),allocatable :: qp(:,:)   
      
      !> wiberg bond orders
      real(wp),allocatable :: wbo(:,:)  
      
      !> HOMO position
      integer :: ihomo = 0,ihomoa = 0,ihomob = 0 

      !> fermi
      real(wp) :: efa = 0.0_wp, efb = 0.0_wp
      
      !> fractional occupation
      real(wp),allocatable :: focc(:)   
      
      !> alpha space
      real(wp),allocatable :: focca(:)  
      
      !> beta space
      real(wp),allocatable :: foccb(:)  
      
      !> orbital energies
      real(wp),allocatable :: emo(:)    
      
      !> molecular orbitals
      real(wp),allocatable :: C(:,:)    
   
   contains
   procedure :: allocate => allocate_wavefunction
   procedure :: deallocate => deallocate_wavefunction
   procedure :: print => print_wavefunction
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

!> print content of wavefunction (for debug)
subroutine print_wavefunction(self, level, unit)
   use iso_fortran_env, only : output_unit
   class(TWavefunction),intent(in) :: self
   integer,intent(in) :: level
   integer,intent(in),optional :: unit

   character(len=*), parameter :: fmt0 = '(3x,a,3x,i0)'
   character(len=30):: fmt1_n 
   character(len=30):: fmt1_nao 
   character(len=30):: fmt1_nshell
   integer :: out, i, ndim
   
   if (self%nao>20) then
      ndim=20
   else
      ndim=self%nao
   end if

   write(fmt1_n,'(a,i0,a)') '(3x,a,3x,',self%n,'(f6.2,1x))'
   write(fmt1_nao,'(a,i0,a)') '(3x,a,3x,',ndim,'(f6.2,1x))'
   write(fmt1_nshell,'(a,i0,a)') '(3x,a,3x,',self%nshell,'(f6.2,1x))'

   if (present(unit)) then
      out = unit
   else
      out = output_unit
   end if


   if (level >= 0 .and. level <= 2) then
      write(out, '(3x,a,/,2x,12("="))') 'Wavefunction'

      ! scalar values !
      write(out,fmt0) 'n      :', self%n
      write(out,fmt0) 'nel    :', self%nel
      write(out,fmt0) 'nopen  :', self%nopen
      write(out,fmt0) 'nao    :', self%nao
      write(out,fmt0) 'nshell :', self%nshell
      write(out,fmt0) 'ihomo  :', self%ihomo
      write(out,fmt0) 'ihomoa :', self%ihomoa
      write(out,fmt0) 'ihomob :', self%ihomob

      ! vector values !
      if (level>0) then
         write(out,fmt1_n) 'q      :', self%q
         write(out,fmt1_nshell) 'qsh    :', self%qsh
         write(out,fmt1_nao) 'focc   :', self%focc(:ndim)
         write(out,fmt1_nao) 'emo    :', self%emo(:ndim)

         ! matrix values !
         if (level>1) then
            
            write(out,'(/)')
            
            do i=1,3
               write(out,fmt1_n) 'dipm   :', self%dipm(i,:)
            end do
            
            write(out,'(/)')
            
            do i=1,6
               write(out,fmt1_n) 'qp     :', self%qp(i,:)
            end do
            
            write(out,'(/)')
            
            do i=1,ndim
               write(out,fmt1_nao) 'C      :', self%C(:ndim,i)
            end do

            write(out,'(/)')

            do i=1,ndim
               write(out,fmt1_nao) 'P      :', self%P(:ndim,i)
            end do
            
            write(out,'(/)')
            
         endif
      endif
   endif
end subroutine print_wavefunction

end module xtb_type_wavefunction
