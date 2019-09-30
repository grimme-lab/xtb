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

module ehtparam
   use iso_fortran_env, only : wp => real64
   implicit none
   private :: wp
   public

!   integer,parameter :: maxao = 50000
!   integer   :: lao,nprim,aoat,fila,lao2,aoat2,fila2,valao,valao2
!   real(wp)  :: alp,cont,hdiag,hdiag2
!   dimension :: lao(maxao),nprim(maxao),aoat(maxao),  &
!   &            fila(2,maxao/5),alp(maxao*5),cont(maxao*5),hdiag(maxao),  &
!   &            lao2(maxao),aoat2(maxao),fila2(2,maxao/5),valao(maxao),  &
!   &            valao2(maxao),hdiag2(maxao)

   integer, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)
!! ------------------------------------------------------------------------
!  This variables describe the ATOMS
   integer, allocatable :: caoshell(:,:)  ! atom number -> basis function
   integer, allocatable :: saoshell(:,:)  ! atom number -> AO
   integer, allocatable :: fila(:,:)      ! atom number -> basis function
   integer, allocatable :: fila2(:,:)     ! atom number -> AO
!! ------------------------------------------------------------------------
!  This variables describe the indivdual SHELLS
   integer, allocatable :: lsh(:) ! shell -> azimudal quantum number
   integer, allocatable :: ash(:) ! shell -> atom number
!! ------------------------------------------------------------------------
!  This variables describe the SPHERICAL atomic orbitals
   real(wp),allocatable :: aoexp(:)    ! AO -> exponent
   integer, allocatable :: ao2sh(:)    ! AO -> shell
   integer, allocatable :: lao2(:)     ! AO -> azimudal quantum number
   integer, allocatable :: aoat2(:)    ! AO -> atom number
   real(wp),allocatable :: hdiag2(:)   ! AO -> level energy
   integer, allocatable :: valao2(:)   ! AO -> AO type
!! ------------------------------------------------------------------------
!  This variables describe the CARTESIAN basis functions
   integer, allocatable :: primcount(:) ! BF -> primitive count
   real(wp),allocatable :: hdiag(:)     ! BF -> level energy
   integer, allocatable :: nprim(:)     ! BF -> primitive number
   integer, allocatable :: aoat(:)      ! BF -> atom
   integer, allocatable :: valao(:)     ! BF -> BF type
   integer, allocatable :: lao(:)       ! BF -> azimudal quantum number
!! ------------------------------------------------------------------------
!  This variables decribe the PRIMITIV basis functions
   real(wp),allocatable :: alp(:)  ! primitive -> primitive exponent
   real(wp),allocatable :: cont(:) ! primitive -> contraction coeffient

contains

subroutine init_ehtparam(n,nbf,nao,nshell)
   implicit none
   integer,intent(in) :: n,nbf,nao,nshell
   call clear_ehtparam
   allocate( hdiag(nbf),     source = 0.0_wp )
   allocate( alp(6*nbf),     source = 0.0_wp )
   allocate( cont(6*nbf),    source = 0.0_wp )
   allocate( hdiag2(nao),    source = 0.0_wp )
   allocate( aoexp(nao),     source = 0.0_wp )
   allocate( ash(nshell),    source = 0 )
   allocate( lsh(nshell),    source = 0 )
   allocate( ao2sh(nao),     source = 0 )
   allocate( nprim(nbf),     source = 0 )
   allocate( primcount(nbf), source = 0 )
   allocate( caoshell(5,n),  source = 0 )
   allocate( saoshell(5,n),  source = 0 )
   allocate( fila(2,n),      source = 0 )
   allocate( fila2(2,n),     source = 0 )
   allocate( lao(nbf),       source = 0 )
   allocate( aoat(nbf),      source = 0 )
   allocate( valao(nbf),     source = 0 )
   allocate( lao2(nao),      source = 0 )
   allocate( aoat2(nao),     source = 0 )
   allocate( valao2(nbf),    source = 0 )
end subroutine init_ehtparam

subroutine clear_ehtparam
   implicit none
   if(allocated(hdiag))     deallocate(hdiag)
   if(allocated(alp))       deallocate(alp)
   if(allocated(cont))      deallocate(cont)
   if(allocated(hdiag2))    deallocate(hdiag2)
   if(allocated(aoexp))     deallocate(aoexp)
   if(allocated(ash))       deallocate(ash)
   if(allocated(lsh))       deallocate(lsh)
   if(allocated(ao2sh))     deallocate(ao2sh)
   if(allocated(nprim))     deallocate(nprim)
   if(allocated(primcount)) deallocate(primcount)
   if(allocated(caoshell))  deallocate(caoshell)
   if(allocated(saoshell))  deallocate(saoshell)
   if(allocated(fila))      deallocate(fila)
   if(allocated(fila2))     deallocate(fila2)
   if(allocated(lao))       deallocate(lao)
   if(allocated(lao2))      deallocate(lao2)
   if(allocated(aoat))      deallocate(aoat)
   if(allocated(aoat2))      deallocate(aoat2)
   if(allocated(valao))     deallocate(valao)
   if(allocated(valao2))    deallocate(valao2)
end subroutine clear_ehtparam

subroutine import_basisset(basis)
   use tbdef_basisset
   implicit none
   type(tb_basisset),intent(in) :: basis
   integer :: n,nbf,nao,nshell

   n = basis%n
   nbf = basis%nbf
   nao = basis%nao
   nshell = basis%nshell

   call init_ehtparam(n,nbf,nao,nshell)

   hdiag = basis%hdiag
   alp = basis%alp
   cont = basis%cont
   hdiag2 = basis%hdiag2
   aoexp = basis%aoexp
   ash = basis%ash
   lsh = basis%lsh
   ao2sh = basis%ao2sh
   nprim = basis%nprim
   primcount = basis%primcount
   caoshell = basis%caoshell
   saoshell = basis%saoshell
   fila = basis%fila
   fila2 = basis%fila2
   lao = basis%lao
   aoat = basis%aoat
   valao = basis%valao
   lao2 = basis%lao2
   aoat2 = basis%aoat2
   valao2 = basis%valao2

end subroutine import_basisset

subroutine export_basisset(basis,n,nbf,nao,nshell)
   use tbdef_basisset
   implicit none
   type(tb_basisset),intent(inout) :: basis
   integer,intent(in) :: n,nbf,nao,nshell

   call basis%allocate(n,nbf,nao,nshell)

   basis%hdiag = hdiag
   basis%alp = alp
   basis%cont = cont
   basis%hdiag2 = hdiag2
   basis%aoexp = aoexp
   basis%ash = ash
   basis%lsh = lsh
   basis%ao2sh = ao2sh
   basis%nprim = nprim
   basis%primcount = primcount
   basis%caoshell = caoshell
   basis%saoshell = saoshell
   basis%fila = fila
   basis%fila2 = fila2
   basis%lao = lao
   basis%aoat = aoat
   basis%valao = valao
   basis%lao2 = lao2
   basis%aoat2 = aoat2
   basis%valao2 = valao2

end subroutine export_basisset



end module ehtparam
