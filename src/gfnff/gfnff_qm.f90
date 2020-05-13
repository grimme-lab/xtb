! This file is part of xtb.
!
! Copyright (C) 2019-2020 Stefan Grimme
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

! solve QM Hamiltonian A in overlap basis S (if ovlp=.true., ZDO otherwise)
! and return density matrix in A (and energy weighted density in S if ovlp=.true.)
! ndim is the dimension of the problem for nel electrons with nopen more alpha than beta

subroutine gfnffqmsolve(pr,A,S,ovlp,et,ndim,nopen,nel,eel,focc,e)
      use xtb_mctc_accuracy, only : wp
      use xtb_scc_core, only : dmat, fermismear, occu
      implicit none
      integer ndim      ! # basis
      integer nopen     ! # of open shells
      integer nel       ! # of electrons
      logical ovlp      ! in overlap basis?
      logical pr        !
      real(wp) et       ! electronic temp Fermi smear
      real(wp) eel      ! electronic energy = sum_occ nocc*eps
      real(wp) focc(ndim)! occupations
      real(wp) e(ndim)  ! eigenvalues
      real(wp) A(ndim,ndim)
      real(wp) S(ndim,ndim)

      integer ihomoa,ihomob,i,liwork,info,lwork
      real(wp) ga,gb,efa,efb,nfoda,nfodb
      real(wp),allocatable:: X(:,:)
      real(wp),allocatable:: focca(:),foccb(:),aux(:)
!     real*4  ,allocatable:: aux4(:),e4(:),A4(:,:)
      integer,allocatable :: iwork(:),ifail(:)


      allocate(focca(ndim),foccb(ndim))

! ZDO case
      if(.not.ovlp) then

      lwork  = 1 + 6*ndim + 2*ndim**2
!     allocate(aux4(lwork),e4(ndim),A4(ndim,ndim))
!     A4 = A
!     call ssyev ('V','U',ndim,A4,ndim,e4,aux4,lwork,info)
!     A = A4
!     e = e4
!     deallocate(aux4,A4,e4)
      allocate(aux(lwork))
      call dsyev ('V','U',ndim,A,ndim,e,aux,lwork,info)

      else
! with overlap case

      allocate (aux(1),iwork(1),ifail(ndim))
      call DSYGVD(1,'V','U',ndim,A,ndim,S,ndim,e,aux, &  !workspace query
     &               -1,IWORK,LIWORK,INFO)
      lwork=int(aux(1))
      liwork=iwork(1)
      deallocate(aux,iwork)

      allocate (aux(lwork),iwork(liwork))              !do it
      call DSYGVD(1,'V','U',ndim,A,ndim,S,ndim,e,aux, &
     &            LWORK,IWORK,LIWORK,INFO)

      endif
! DONE

! scale energy  so that gaps a roughly ok (in eV)

      e=e*0.1*27.2113957

      if(info.ne.0) then
         write(*,*) 'INFO',info
         stop 'diag error'
      endif

      if(et.gt.1.d-3) then
! Fermi smearing, convert restricted occ first to alpha/beta
        call occu(ndim,nel,nopen,ihomoa,ihomob,focca,foccb)
        if(ihomoa.le.ndim.and.ihomoa.gt.0) then
         call fermismear(.false.,ndim,ihomoa,et,e,focca,nfoda,efa,ga)
        else
          focca=0
        endif
        if(ihomob.le.ndim.and.ihomob.gt.0) then
          call fermismear(.false.,ndim,ihomob,et,e,foccb,nfodb,efb,gb)
        else
          foccb=0
        endif
        focc = focca + foccb
        if(ihomoa+1.le.ndim)then
        if(abs(focc(ihomoa)-focc(ihomoa+1)).lt.1.d-4) then ! a perfect birad is anit-aromatic
           focc = 0                                        ! and hence we break the sym
           do i=1,nel/2
              focc(i)=2.0d0
           enddo
           write(*,*) 'perfect biradical detected at FT-HMO level. Breaking the symmetry'
           write(*,*) 'because its assumed to be an anit-aromatic system like COT or CB.'
        endif
        endif
      else
        focc = 0
        do i=1,nel/2
           focc(i)=2.0d0
        enddo
        if(2*(nel/2).ne.nel) focc(nel/2+1)=1.0
      endif

      focca = focc * e
      eel = sum(focca)
! density matrix
      S = A
      call dmat(ndim,focc,S,A)

      if(ovlp) then
      allocate(X(ndim,ndim))
! energy weighted density on S
      call dmat(ndim,focca,S,X)
      S = X
      endif

      deallocate(focca,foccb)
      
end subroutine gfnffqmsolve

