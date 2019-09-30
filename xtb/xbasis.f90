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

module xbasis
   use iso_fortran_env, only : wp => real64
   implicit none

contains

subroutine xbasis_gfn1(n,at,basis,ok,diff)
   use tbdef_basisset
   use aoparam
   implicit none
   type(tb_basisset),intent(inout) :: basis
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   logical, intent(out) :: ok
   logical, intent(out) :: diff

   integer  :: elem
   integer  :: i,j,m,l,ibf,ipr,p,thisprim,thisprimR,idum,npq,npqR,pqn
   real(wp) :: a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi
   real(wp) :: aR(10),cR(10),ss
   real(wp) :: as(10),cs(10)
   real(wp) :: ap(10),cp(10)

!  if(n.gt.maxao/4) & ! error not necessary anymore
!  &   call raise('E','system too large. recompile code!',1)

   basis%hdiag(1:basis%nbf)=1.d+42
   ! note: Rydbergs are identified by valao(*)=0
   !       polarization by            valao(*)=-1
   basis%valao(1:basis%nbf)=1

   ibf=0
   ipr=0
   ok=.true.

   !  atom loop
   do i=1,n
      !     AO=shell loop
      basis%fila(1,i)=ibf+1
      do m=1,ao_n(at(i))
         !        principle QN
         npq=ao_pqn(m,at(i))
         l=ao_l(m,at(i))
! ========================================================================
!        H-He
! ========================================================================
         if(l.eq.0.and.at(i).le.2.and.npq.eq.1)then
            ! s
            ibf =ibf+1
            zeta=ao_exp(m,at(i))
            call setsto4(thisprim,npq,1,zeta,a,c)
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif

         if(l.eq.12.and.at(i).le.2)then
            ! diff sp
            ibf =ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum), &
            &           basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            split1=0
            if(at(i).eq.2)split1=2.25  ! set to SCS-CC2 value

            basis%hdiag(ibf)=ao_lev(m,at(i))-split1
            zeta=ao_exp(m,at(i))
            call setsto4(thisprim,2,2,zeta,a,c)
            do j=2,4
               ibf=ibf+1
               basis%valao(ibf)=0
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))+split1
            enddo
         endif

         if(l.eq.0.and.at(i).le.2.and.npq.eq.2)then
            ! diff s
            ibf =ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum), &
            &           basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif

         ! p polarization
         if(l.eq.1.and.at(i).le.2)then
            zeta=ao_exp(m,at(i))
            call setsto3(thisprim,npq,2,zeta,ap,cp)
            do j=2,4
               ibf=ibf+1
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=ap(p)
                  basis%cont(ipr)=cp(p)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%valao(ibf)=-1
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

! ========================================================================
!        general sp
! ========================================================================
         if(l.eq.0.and.at(i).gt.2)then
            !   s
            ibf=ibf+1
            zeta=ao_exp(m,at(i))
            call setsto6(thisprim,npq,1,zeta,as,cs)
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif
         !   p
         if(l.eq.1.and.at(i).gt.2)then
            zeta=ao_exp(m,at(i))
            call setsto6(thisprim,npq,2,zeta,ap,cp)
            do j=2,4
               ibf=ibf+1
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=ap(p)
                  basis%cont(ipr)=cp(p)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

! ========================================================================
! diffuse. important: the code assumes that the previous shells are
!          valence on which the functions are orthogonalized
!          i.e. prims of val. are in as,ap, cs,cp
! ========================================================================
         if(l.ge.12.and.at(i).gt.2)then
            !   s
            ibf=ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
               call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum), &
            &           basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
            !  p
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
               call setsto6(thisprimR,npq,2,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,2,zeta,aR,cR)
            endif
            call atovlp(1,thisprim,thisprimR,ap,aR,cp,cR,pp)
            do j=2,4
               ibf=ibf+1
               basis%valao(ibf)=0
               idum=ipr+1
               do p=1,thisprimR
                  ipr=ipr+1
                  basis%alp (ipr)=aR(p)
                  basis%cont(ipr)=cR(p)
               enddo
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=ap(p)
                  basis%cont(ipr)=-pp*cp(p)
               enddo
               basis%nprim(ibf)=thisprimR+thisprim
               call atovlp(1,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum), &
               &           basis%cont(idum),basis%cont(idum),ss)
               do p=1,basis%nprim(ibf)
                  basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
            !  d if its an spd Ryd shell
            if(l.eq.13)then
               zeta=ao_exp(m,at(i))
               if(npq.lt.5) then
                  call setsto3(thisprim,npq,3,zeta,a,c)
               else
                  call setsto4(thisprim,5  ,3,zeta,a,c)
               endif
               do j=5,10
                  ibf=ibf+1
                  basis%valao(ibf)=0
                  do p=1,thisprim
                     ipr=ipr+1
                     basis%alp (ipr)=a(p)
                     basis%cont(ipr)=c(p)
                     if(j.gt.7)basis%cont(ipr)=basis%cont(ipr)*sqrt(3.0d0)
                  enddo
                  basis%aoat (ibf)=i
                  basis%lao  (ibf)=j
                  basis%nprim(ibf)=thisprim
                  basis%hdiag(ibf)=ao_lev(m,at(i))
               enddo
            endif
         endif

         ! DZ s
         if(l.eq.11)then
            ibf=ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
               call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum), &
            &           basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
         endif

! ========================================================================
!  d
! ========================================================================
         if(l.eq.2)then
            zeta=ao_exp(m,at(i))
            !           valence
            if(ao_typ(m,at(i)).ne.-1)then
               call setsto4(thisprim,npq,3,zeta,a,c)
            else
               !           polarization
               call setsto4(thisprim,npq,3,zeta,a,c)
            endif
            do j=5,10
               ibf=ibf+1
               if(ao_typ(m,at(i)).eq.-1)basis%valao(ibf)=-1
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
                  if(j.gt.7)basis%cont(ipr)=basis%cont(ipr)*sqrt(3.0d0)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

! ========================================================================
!  f
! ========================================================================
         if(l.eq.3)then
            zeta=ao_exp(m,at(i))
            do j=11,20
               ibf=ibf+1
               !           valence
               call setsto4(thisprim,npq,4,zeta,a,c)
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
                  if(j.gt.13.and.j.lt.20)basis%cont(ipr)=basis%cont(ipr)*sqrt(5.0d0)
                  if(j.eq.20            )basis%cont(ipr)=basis%cont(ipr)*sqrt(15.0d0)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
         ! next shell
      enddo
      ! next atom
      basis%fila(2,i)=ibf
   enddo


   diff=.false.
   do i=1,ibf
      if(basis%valao(i).eq.0) diff=.true.
      if(basis%hdiag(i).gt.1.d+10)then
         write(*,*)'Hii not defined for',i,basis%aoat(i)
         ok=.false.
      endif
   enddo
   do i=1,ipr
      if(basis%alp(i).eq.0) then
         ok=.false.
         write(*,*)'alp=0 for',i
      endif
   enddo

   if(basis%nbf.ne.ibf) then
      write(*,*) ibf,basis%nbf
      error stop 'internal error'
   endif

   !     write(*,*) alp(1:ipr)
   !     write(*,*) cont(1:ipr)

   do i=1,n
      do j=1,ao_n(at(i))
         l = ao_l(j,at(i))
         if(l.eq.11)ao_l(j,at(i))=0
      enddo
   enddo

end subroutine xbasis_gfn1

subroutine xbasis_gfn2(n,at,basis,ok)
   use tbdef_basisset
   use aoparam
   implicit none
   type(tb_basisset),intent(inout) :: basis
   integer, intent(in)  :: n
   integer, intent(in)  :: at(n)
   logical, intent(out) :: ok

   integer  :: elem
   integer  :: i,j,m,l,ibf,ipr,p
   integer  :: thisprim,thisprimR,idum,npq,npqR,pqn
   integer  :: maingroup,info
   real(wp) :: a(10),c(10),zeta

!  if(n.gt.maxao/4) & ! error not necessary anymore
!  &   call raise('E','system too large. recompile code!',1)

   basis%hdiag(1:basis%nbf)=1.d+42
   basis%valao(1:basis%nbf)=1

   ibf=0
   ipr=0
   ok=.true.

! ========================================================================
! atom loop
! ========================================================================
   do i=1,n
      !     AO=shell loop
      basis%fila(1,i)=ibf+1
      do m=1,ao_n(at(i))
         !        principle QN
         npq=ao_pqn(m,at(i))
         l=ao_l(m,at(i))

! ========================================================================
!        H-He
! ========================================================================
         if(l.eq.0.and.at(i).le.2.and.npq.eq.1)then
            !  s
            ibf =ibf+1
            zeta=ao_exp(m,at(i))
            call setsto3(thisprim,npq,1,zeta,a,c)
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif

! ========================================================================
!        general spdf
! ========================================================================
         if(l.eq.0.and.at(i).gt.2)then
            !   s
            ibf=ibf+1
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
               call setsto6(thisprim,npq,1,zeta,a,c)
            else
               call setsto4(thisprim,npq,1,zeta,a,c)
            endif
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif
         !   p
         if(l.eq.1)then
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
               call setsto6(thisprim,npq,2,zeta,a,c)
            else
               call setsto4(thisprim,npq,2,zeta,a,c)
            endif
            do j=2,4
               ibf=ibf+1
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
         !   d
         if(l.eq.2)then
            zeta=ao_exp(m,at(i))
            call setsto3(thisprim,npq,3,zeta,a,c)
            do j=5,10
               ibf=ibf+1
               if(ao_typ(m,at(i)).eq.-1)basis%valao(ibf)=-1
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
                  if(j.gt.7)basis%cont(ipr)=basis%cont(ipr)*sqrt(3.0d0)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
         !    f
         if(l.eq.3)then
            zeta=ao_exp(m,at(i))
            do j=11,20
               ibf=ibf+1
               call setsto4(thisprim,npq,4,zeta,a,c)
               do p=1,thisprim
                  ipr=ipr+1
                  basis%alp (ipr)=a(p)
                  basis%cont(ipr)=c(p)
                  if(j.gt.13.and.j.lt.20)basis%cont(ipr)=basis%cont(ipr)*sqrt(5.0d0)
                  if(j.eq.20            )basis%cont(ipr)=basis%cont(ipr)*sqrt(15.0d0)
               enddo
               basis%aoat (ibf)=i
               basis%lao  (ibf)=j
               basis%nprim(ibf)=thisprim
               basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
         ! next shell
      enddo
      ! next atom
      basis%fila(2,i)=ibf
   enddo
! ========================================================================

   do i=1,ibf
      if(basis%hdiag(i).gt.1.d+10)then
         write(*,*)'Hii not defined for',i,basis%aoat(i)
         ok=.false.
      endif
   enddo
   do i=1,ipr
      if(basis%alp(i).eq.0) then
         ok=.false.
         write(*,*)'alp=0 for',i
      endif
   enddo

   if(basis%nbf.ne.ibf) then
      write(*,*) ibf,basis%nbf
      error stop 'internal error'
   endif

!  write(*,*) alp(1:ipr)
!  write(*,*) cont(1:ipr)

   do i=1,n
      do j=1,ao_n(at(i))
         l = ao_l(j,at(i))
         if(l.eq.11)ao_l(j,at(i))=0
      enddo
   enddo

end subroutine xbasis_gfn2

! ========================================================================
!  stuff
! ========================================================================
subroutine xbasis_cao2sao(n,at,basis)
   use tbdef_basisset
   use aoparam
   implicit none
   type(tb_basisset),intent(inout) :: basis
   integer,intent(in) :: n,at(n)
   integer :: i,ia,ishell,mi,j,k,kk,kkk,l,m
   integer, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

!  stuff for ovlp
!  PRIM INDEX
   k=0
   do i=1,basis%nbf
      basis%primcount(i)=k
      k=k+basis%nprim(i)
   enddo
!  AO/CAO SHELL INDEX
   kkk = 0
   k=0
   j=0
   kk=0
   do i=1,n
      ia=at(i)
!     aostart(i)=j+1
      do mi=1,ao_n(ia)
         kkk=kkk+1
         ishell=ao_l(mi,ia)
         basis%lsh(kkk)=ishell
         basis%ash(kkk)=i
         basis%caoshell(mi,i)=k
         basis%saoshell(mi,i)=j
         do m=1,llao2(ishell)
            kk=kk+1
            basis%ao2sh(kk)=kkk
            basis%aoexp(kk)=ao_exp(mi,ia)
         enddo
         if(ishell.eq.0) k=k+1
         if(ishell.eq.1) k=k+3
         if(ishell.eq.2) k=k+6
         if(ishell.eq.3) k=k+10
         if(ishell.eq.0) j=j+1
         if(ishell.eq.1) j=j+3
         if(ishell.eq.2) j=j+5
         if(ishell.eq.3) j=j+7
      enddo
!     aoend(i)=j
   enddo

   if(basis%nbf.ne.basis%nao)then
!     relocate fila,lao,hdiag
      call dtrafo2(basis)
   else
      basis%lao2  =basis%lao
      basis%aoat2 =basis%aoat
      basis%valao2=basis%valao
      basis%fila2 =basis%fila
      basis%hdiag2=basis%hdiag
   endif

end subroutine xbasis_cao2sao

subroutine dtrafo2(basis)
  use tbdef_basisset
  implicit none
  type(tb_basisset),intent(inout) :: basis
  integer i,j,be,en

  !     do i=1,n
  !        write(*,*) '6d ',fila(1:2,i)
  !     enddo

  j=0
  do i=1,basis%nbf
     if(basis%lao(i).ne.5)then
        j=j+1
        basis%lao2(j)=basis%lao(i)
        if(basis%lao(i).gt.4) basis%lao2(j)=basis%lao2(j)-1
        basis%aoat2 (j)=basis%aoat (i)
        basis%valao2(j)=basis%valao(i)
        basis%hdiag2(j)=basis%hdiag(i)
     endif
  enddo

  do j=1,basis%n
     be=100000
     en=-1
     do i=1,basis%nao
        if(basis%aoat2(i).eq.j.and.i.lt.be)be=i
        if(basis%aoat2(i).eq.j.and.i.gt.en)en=i
     enddo
     basis%fila2(1,j)=be
     basis%fila2(2,j)=en
  enddo

  !     do i=1,n
  !        write(*,*) '5d ',fila(1:2,i)
  !     enddo
  !     write(*,'(''6d '',40i3)') lao (1:nbf)
  !     write(*,'(''5d '',40i3)') lao2(1:nbf)

end subroutine dtrafo2


! ========================================================================
!  determine bf limits
! ========================================================================

subroutine xbasis0(n,at,basis)
   use tbdef_basisset
   use aoparam
   implicit none
   type(tb_basisset),intent(inout) :: basis
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer :: nbf
   integer :: nao
   integer :: nshell

   integer i,j,k,l

   basis%n = n

   call dim_basis(n,at,nshell,nao,nbf)

   call basis%allocate(n,nbf,nao,nshell)
   !call init_ehtparam(n,nbf,nao,nshell)
!  if (nbf.gt.maxao) &
!  &  call raise('E','TB basis too large for common storage',1)

end subroutine xbasis0

subroutine dim_basis(n,at,nshell,nao,nbf)
   use tbdef_basisset
   use aoparam
   implicit none
   integer,intent(in)  :: n
   integer,intent(in)  :: at(n)
   integer,intent(out) :: nshell
   integer,intent(out) :: nao
   integer,intent(out) :: nbf

   integer i,j,k,l

   nao=0
   nbf=0
   nshell=0

   do i=1,n
      k=0
      do j=1,ao_n(at(i))
         l = ao_l(j,at(i))
         k = k + 1
         nshell=nshell+1
!        if(l.eq.0) nbf=nbf+1
!        if(l.eq.1) nbf=nbf+3
!        if(l.eq.2) nbf=nbf+6
!        if(l.eq.3) nbf=nbf+10
!        if(l.eq.3) call raise('E','f-functions not implemented',1)
!        if(l.eq.4) call raise('E','g-functions not implemented',1)
!        if(l.eq.11)nbf=nbf+1
!        if(l.eq.12)nbf=nbf+4
!        if(l.eq.13)nbf=nbf+10
         select case(l)
         case(0) ! s
            nbf = nbf+1
            nao = nao+1
         case(1) ! p
            nbf = nbf+3
            nao = nao+3
         case(2) ! d
            nbf = nbf+6
            nao = nao+5
         case(3) ! f
            nbf = nbf+10
            nao = nao+7
            call raise('E','f-functions not implemented',1)
         case(4) ! g
            nbf = nbf+15
            nao = nao+9
            call raise('E','g-functions not implemented',1)
         case(11) ! diffuse s
            nbf = nbf+1
            nao = nao+1
         case(12) ! diffuse sp
            nbf = nbf+4
            nao = nao+4
         case(13) ! diffuse spd
            nbf = nbf+10
            nao = nao+ 9
         end select
      enddo
      if(k.eq.0) then
         write(*,*) 'no basis found for atom', i,' Z= ',at(i)
         call terminate(1)
      endif
   enddo

end subroutine dim_basis

subroutine xbasis_gfn0(n,at,basis,ok,diff)    !ppracht 10/2018
   use tbdef_basisset
      use aoparam
      implicit none
      type(tb_basisset),intent(inout) :: basis

      integer elem,n
      integer at(n)
      logical ok,diff
      !include 'ehtcommon.f'
      !include 'aoelementcommon.f'

      integer i,j,m,l,ibf,ipr,p,thisprim,thisprimR,idum,npq,npqR,pqn
      real(wp)  a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi
      real(wp)  aR(10),cR(10),ss
      real(wp)  as(10),cs(10)
      real(wp)  ap(10),cp(10)

      !if(n.gt.maxao/4) stop 'system too large. recompile code!'

      basis%hdiag(1:basis%nbf)=1.d+42
! note: Rydbergs are identified by valao(*)=0
!       polarization by            valao(*)=-1
      basis%valao(1:basis%nbf)=1

      ibf=0
      ipr=0
      ok=.true.

!     atom loop
      do i=1,n
!     AO=shell loop
      basis%fila(1,i)=ibf+1
      do m=1,ao_n(at(i))
!        principle QN
         npq=ao_pqn(m,at(i))
         l=ao_l(m,at(i))
!=================================
!        H-He
!=================================
         if(l.eq.0.and.at(i).le.2.and.npq.eq.1)then
!  s
            ibf =ibf+1
            zeta=ao_exp(m,at(i))
            call setsto3(thisprim,npq,1,zeta,a,c)  !GFN0
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif

         if(l.eq.12.and.at(i).le.2)then
! diff sp
            ibf =ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum),&
            &                                 basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            split1=0
            if(at(i).eq.2)split1=2.25  ! set to SCS-CC2 value

            basis%hdiag(ibf)=ao_lev(m,at(i))-split1
            zeta=ao_exp(m,at(i))
            call setsto4(thisprim,2,2,zeta,a,c)
            do j=2,4
            ibf=ibf+1
            basis%valao(ibf)=0
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))+split1
            enddo
         endif

         if(l.eq.0.and.at(i).le.2.and.npq.eq.2)then
! diff s
            ibf =ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            !call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call setsto2(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,thisprim,thisprimR,a,aR,c,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum),&
            &                                 basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif

! p polarization
         if(l.eq.1.and.at(i).le.2)then
            zeta=ao_exp(m,at(i))
            call setsto3(thisprim,npq,2,zeta,ap,cp)
            do j=2,4
            ibf=ibf+1
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=ap(p)
               basis%cont(ipr)=cp(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%valao(ibf)=-1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

!=================================
! general sp
!=================================
         if(l.eq.0.and.at(i).gt.2)then
!   s
            ibf=ibf+1
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprim,npq,1,zeta,as,cs)
            else
            call setsto4(thisprim,npq,1,zeta,as,cs) !GFN0
            endif
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
         endif
!   p
         if(l.eq.1.and.at(i).gt.2)then
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprim,npq,2,zeta,ap,cp)
            else
            call setsto3(thisprim,npq,2,zeta,ap,cp) !GFN0
            endif
            do j=2,4
            ibf=ibf+1
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=ap(p)
               basis%cont(ipr)=cp(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

!=================================
! diffuse. important: the code assumes that the previous shells are
!          valence on which the functions are orthogonalized
!          i.e. prims of val. are in as,ap, cs,cp
!=================================
         if(l.ge.12.and.at(i).gt.2)then
!   s
            ibf=ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum),&
            &                                 basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
!  p
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,2,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,2,zeta,aR,cR)
            endif
            call atovlp(1,thisprim,thisprimR,ap,aR,cp,cR,pp)
            do j=2,4
            ibf=ibf+1
            basis%valao(ibf)=0
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=ap(p)
               basis%cont(ipr)=-pp*cp(p)
            enddo
            basis%nprim(ibf)=thisprimR+thisprim
            call atovlp(1,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum),&
            &                                 basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
!  d if its an spd Ryd shell
            if(l.eq.13)then
            zeta=ao_exp(m,at(i))
            if(npq.lt.5) then
            call setsto3(thisprim,npq,3,zeta,a,c)
            else
            call setsto4(thisprim,5  ,3,zeta,a,c)
            endif
            do j=5,10
            ibf=ibf+1
            basis%valao(ibf)=0
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
               if(j.gt.7)basis%cont(ipr)=basis%cont(ipr)*sqrt(3.0d0)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
            endif
         endif

! DZ s
         if(l.eq.11)then
            ibf=ibf+1
            basis%valao(ibf)=0
            zeta=ao_exp(m,at(i))
            if(npq.gt.5) then
            call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,thisprim,thisprimR,as,aR,cs,cR,ss)
            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=1
            basis%nprim(ibf)=thisprimR+thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf),basis%alp(idum),basis%alp(idum),&
            &                                 basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo
         endif

!=================================
! d
!=================================
         if(l.eq.2)then
            zeta=ao_exp(m,at(i))
!           valence
            if(ao_typ(m,at(i)).ne.-1)then
            call setsto4(thisprim,npq,3,zeta,a,c)
            else
!           polarization
            call setsto4(thisprim,npq,3,zeta,a,c)
            endif
            do j=5,10
            ibf=ibf+1
            if(ao_typ(m,at(i)).eq.-1)basis%valao(ibf)=-1
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
               if(j.gt.7)basis%cont(ipr)=basis%cont(ipr)*sqrt(3.0d0)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif

!=================================
! f
!=================================
         if(l.eq.3)then
            zeta=ao_exp(m,at(i))
            do j=11,20
            ibf=ibf+1
!           valence
            call setsto4(thisprim,npq,4,zeta,a,c)
            do p=1,thisprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
               if(j.gt.13.and.j.lt.20)basis%cont(ipr)=basis%cont(ipr)*sqrt(5.0d0)
               if(j.eq.20            )basis%cont(ipr)=basis%cont(ipr)*sqrt(15.0d0)
            enddo
            basis%aoat (ibf)=i
            basis%lao  (ibf)=j
            basis%nprim(ibf)=thisprim
            basis%hdiag(ibf)=ao_lev(m,at(i))
            enddo
         endif
! next shell
         enddo
! next atom
         basis%fila(2,i)=ibf
      enddo


      diff=.false.
      do i=1,ibf
         if(basis%valao(i).eq.0) diff=.true.
         if(basis%hdiag(i).gt.1.d+10)then
           write(*,*)'Hii not defined for',i,basis%aoat(i)
           ok=.false.
         endif
      enddo
      do i=1,ipr
         if(basis%alp(i).eq.0) then
          ok=.false.
          write(*,*)'alp=0 for',i
         endif
      enddo

      if(basis%nbf.ne.ibf) then
      write(*,*) ibf,basis%nbf
      stop 'internal error'
      endif

!     write(*,*) alp(1:ipr)
!     write(*,*) cont(1:ipr)

      do i=1,n
         do j=1,ao_n(at(i))
            l = ao_l(j,at(i))
            if(l.eq.11)ao_l(j,at(i))=0
         enddo
      enddo

end subroutine xbasis_gfn0

! ------------------------------------------------------------------------
!  Helper functions

subroutine atovlp(l,npri,nprj,alpa,alpb,conta,contb,ss)
  implicit none
  integer l,npri,nprj
  real(wp) alpa(*),alpb(*)
  real(wp) conta(*),contb(*)
  real(wp) ss

  integer ii,jj
  real(wp) ab,s00,sss,pi,ab05
  data pi/3.1415926535897932384626433832795029_wp/

  SS=0.0_wp
  do ii=1,npri
     do jj=1,nprj
        ab =1./(alpa(ii)+alpb(jj))
        s00=(pi*ab)**1.50_wp
        if(l.eq.0)then
           sss=s00
        endif
        if(l.eq.1)then
           ab05=ab*0.5_wp
           sss=s00*ab05
        endif
        SS=SS+SSS*conta(ii)*contb(jj)
     enddo
  enddo

end subroutine atovlp


end module xbasis
