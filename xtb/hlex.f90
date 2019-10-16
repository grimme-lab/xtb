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

subroutine hlex(nat,at,nbf,nao,ihomo,xyz,focc,s,cmo,eiga,basis)
   use iso_fortran_env, wp => real64
   use mctc_econv, only : autoev,evtoau
   use tbdef_basisset
   use dtrafo
   implicit none
   type(tb_basisset), intent(in) :: basis
   integer, intent(in)  :: nao,ihomo,nat,at(nat),nbf
   real(wp),intent(in)  :: cmo(nao,nao),eiga(nao),focc(nao),xyz(3,nat)
   real(wp),intent(in)  :: s(nao,nao)

   real(wp),allocatable :: cca(:)
   real(wp),allocatable :: dip2(:)
   real(wp),allocatable :: dip (:,:)

   integer  :: mo,n,i,j,k,idum,ii
   real(wp) :: dum,de,dtot(3),fl

   ! # of transformed MOs
   n=2

   allocate(cca(nao*nao))

   ! do only singly occ. ones
   cca=0
   ii =0
   do i=ihomo-1,ihomo
      ii=ii+1
      do j=1,nao
         cca(j+(ii-1)*nao)=cmo(j,i)
      enddo
   enddo

   ! boys
   !     dipole integrals
   allocate(dip2(nao*nao),dip(nbf*(nbf+1)/2,3))

   call Dints(nat,nbf,xyz,dip(1,1),dip(1,2),dip(1,3),basis)
   call cao2saop(nbf,nao,dip(1,1),basis)
   call cao2saop(nbf,nao,dip(1,2),basis)
   call cao2saop(nbf,nao,dip(1,3),basis)
   i=1
   j=2
   do mo=1,3
      call onetri(1,dip(1,mo),dip2,cca,nao,n)
      dtot(mo)=dip2(j+(i-1)*n)
   enddo
   dum=sqrt(dtot(1)**2+dtot(2)**2+dtot(3)**2)
   de=(eiga(ihomo)-eiga(ihomo-1))*evtoau
   write(*,*)
   write(*,*)'transition dipole moment (au) for excitation:',ihomo-1,ihomo
   write(*,*)'    X       Y       Z   '
   write(*,'(3f9.4,''  total (au/Debye): '',2f8.3)') &
   &      dtot(1),dtot(2),dtot(3),dum,dum*2.5418
   fl=sqrt(2.)*(2./3.)*dum*dum*de
   write(*,'('' dE (eV)             : '',f8.3)') de*autoev
   write(*,'('' oscillator strength : '',e12.5)') fl

end subroutine hlex
