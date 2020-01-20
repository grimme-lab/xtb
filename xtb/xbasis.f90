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

   integer  :: elem,valao
   integer  :: i,j,m,l,iat,ati,ish,ibf,iao,ipr,p,nprim,thisprimR,idum,npq,npqR,pqn
   real(wp) :: a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi,level
   real(wp) :: aR(10),cR(10),ss
   real(wp) :: as(10),cs(10)
   real(wp) :: ap(10),cp(10)

   basis%hdiag(1:basis%nbf)=1.d+42

   ibf=0
   iao=0
   ipr=0
   ish=0
   ok=.true.

   atoms: do iat=1,n
      ati = at(iat)
      basis%shells(1,iat)=ish+1
      basis%fila  (1,iat)=ibf+1
      basis%fila2 (1,iat)=iao+1
      shells: do m=1,ao_n(ati)
         ish = ish+1
         ! principle QN
         npq=ao_pqn(m,ati)
         l=ao_l(m,ati)

         level = ao_lev(m,ati)
         zeta  = ao_exp(m,ati)

         basis%lsh(ish) = l
         if (l.eq.11) basis%lsh(ish) = 0
         basis%ash(ish) = iat
         basis%sh2bf(1,ish) = ibf
         basis%sh2ao(1,ish) = iao
         basis%caoshell(m,iat)=ibf
         basis%saoshell(m,iat)=iao

         ! add new shellwise information, for easier reference
         basis%level(ish) = level
         basis%zeta (ish) = zeta

         if (ao_typ(m,ati).eq.-1) then
            valao = -1
         else
            valao = 1
         endif

         ! H-He
         if(l.eq.0.and.ati.le.2.and.npq.eq.1)then
            ! s
            call setsto4(nprim,npq,1,zeta,a,c)
            basis%minalp(ish) = minval(a(:nprim))

            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 1
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = 1
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         if(l.eq.0.and.ati.le.2.and.npq.eq.2)then
            ! diff s
            call setsto3(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,nprim,thisprimR,a,aR,c,cR,ss)
            basis%minalp(ish) = min(minval(a(:nprim)),minval(aR(:thisprimR)))

            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 0
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum), &
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = 0
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! p polarization
         if(l.eq.1.and.ati.le.2)then
            call setsto3(nprim,npq,2,zeta,ap,cp)
            basis%minalp(ish) = minval(ap(:nprim))
            do j=2,4
               ibf=ibf+1
               basis%primcount(ibf) = ipr
               basis%aoat     (ibf) = iat
               basis%lao      (ibf) = j
               basis%valao    (ibf) = -1
               basis%nprim    (ibf) = nprim
               basis%hdiag    (ibf) = level

               do p=1,nprim
                  ipr=ipr+1
                  basis%alp (ipr)=ap(p)
                  basis%cont(ipr)=cp(p)
               enddo

               iao = iao+1
               basis%valao2(iao) = -1
               basis%aoat2 (iao) = iat
               basis%lao2  (iao) = j
               basis%hdiag2(iao) = level
               basis%aoexp (iao) = zeta
               basis%ao2sh (iao) = ish
            enddo
         endif

         ! general sp
         if(l.eq.0.and.ati.gt.2)then
            ! s
            call setsto6(nprim,npq,1,zeta,as,cs)
            basis%minalp(ish) = minval(as(:nprim))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 1
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=cs(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = 1
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif
         ! p
         if(l.eq.1.and.ati.gt.2)then
            call setsto6(nprim,npq,2,zeta,ap,cp)
            basis%minalp(ish) = minval(ap(:nprim))
            do j=2,4
               ibf=ibf+1
               basis%primcount(ibf) = ipr
               basis%valao    (ibf) = 1
               basis%aoat     (ibf) = iat
               basis%lao      (ibf) = j
               basis%nprim    (ibf) = nprim
               basis%hdiag    (ibf) = level

               do p=1,nprim
                  ipr=ipr+1
                  basis%alp (ipr)=ap(p)
                  basis%cont(ipr)=cp(p)
               enddo

               iao = iao+1
               basis%valao2(iao) = 1
               basis%aoat2 (iao) = iat
               basis%lao2  (iao) = j
               basis%hdiag2(iao) = level
               basis%aoexp (iao) = zeta
               basis%ao2sh (iao) = ish
            enddo
         endif

         ! DZ s
         if(l.eq.11)then
            if(npq.gt.5) then
               call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,nprim,thisprimR,as,aR,cs,cR,ss)
            basis%minalp(ish) = min(minval(as(:nprim)),minval(aR(:thisprimR)))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 0
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum), &
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = 0
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! d
         if(l.eq.2)then
            call set_d_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,4,zeta,level,valao)
         endif

         ! f
         if(l.eq.3)then
            call set_f_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,4,zeta,level,1)
         endif

         basis%valsh(ish) = basis%valao2(iao)

         basis%sh2bf(2,ish) = ibf-basis%sh2bf(1,ish)
         basis%sh2ao(2,ish) = iao-basis%sh2ao(1,ish)
      enddo shells
      basis%shells(2,iat)=ish
      basis%fila  (2,iat)=ibf
      basis%fila2 (2,iat)=iao
   enddo atoms

   ! note: Rydbergs are identified by valao(*)=0
   !       polarization by            valao(*)=-1
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

   if(basis%nao.ne.iao) then
      call raise('E','Internal error',1)
   endif

   do iat=1,n
      ati = at(iat)
      do j=1,ao_n(ati)
         l = ao_l(j,ati)
         if(l.eq.11)ao_l(j,ati)=0
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

   integer  :: elem,valao
   integer  :: i,iat,ati,j,m,l,ish,ibf,iao,ipr,p
   integer  :: nprim,thisprimR,idum,npq,npqR,pqn
   integer  :: info
   real(wp) :: a(10),c(10),zeta,level

   basis%hdiag(1:basis%nbf)=1.d+42

   iao=0
   ibf=0
   ipr=0
   ish=0
   ok=.true.

! ========================================================================
   atoms: do iat=1,n
      ati = at(iat)
! ========================================================================
      basis%shells(1,iat)=ish+1
      basis%fila  (1,iat)=ibf+1
      basis%fila2 (1,iat)=iao+1
      shells: do m=1,ao_n(ati)
         ish = ish+1
         ! principle QN
         npq=ao_pqn(m,ati)
         l=ao_l(m,ati)

         level = ao_lev(m,ati)
         zeta  = ao_exp(m,ati)

         basis%lsh(ish) = l
         if (l.eq.11) basis%lsh(ish) = 0
         basis%ash(ish) = iat
         basis%sh2bf(1,ish) = ibf
         basis%sh2ao(1,ish) = iao
         basis%caoshell(m,iat) = ibf
         basis%saoshell(m,iat) = iao

         ! add new shellwise information, for easier reference
         basis%level(ish) = level
         basis%zeta (ish) = zeta

         if (ao_typ(m,ati).eq.-1) then
            valao = -1
         else
            valao = 1
         endif

         ! H-He
         if(l.eq.0.and.ati.le.2.and.npq.eq.1)then
            call set_s_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,3,zeta,level,valao)
         endif
         ! general spdf
         if(l.eq.0.and.ati.gt.2)then
            if(npq.gt.5) then
               call set_s_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,6,zeta,level,valao)
            else
               call set_s_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,4,zeta,level,valao)
            endif
         endif
         if(l.eq.1)then
            if(npq.gt.5) then
               call set_p_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,6,zeta,level,valao)
            else
               call set_p_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,4,zeta,level,valao)
            endif
         endif
         if(l.eq.2)then
            call set_d_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,3,zeta,level,valao)
         endif
         if(l.eq.3)then
            call set_f_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,4,zeta,level,valao)
         endif

         basis%valsh(ish) = basis%valao2(iao)

         basis%sh2bf(2,ish) = ibf-basis%sh2bf(1,ish)
         basis%sh2ao(2,ish) = iao-basis%sh2ao(1,ish)
      enddo shells
      basis%shells(2,iat)=ish
      basis%fila  (2,iat)=ibf
      basis%fila2 (2,iat)=iao
   enddo atoms
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

   if(basis%nao.ne.iao) then
      call raise('E','Internal error',1)
   endif

end subroutine xbasis_gfn2

! ========================================================================
!> determine basisset limits
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

   call dim_basis(n,at,nshell,nao,nbf)

   call basis%allocate(n,nbf,nao,nshell)

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

   integer i,j,m,l,iat,ati,ish,ibf,iao,ipr,p,nprim,thisprimR,idum,npq,npqR,pqn,valao
   real(wp)  a(10),c(10),zeta,k1,k2,split1,pp,zqfR,zcnfR,qi,level
   real(wp)  aR(10),cR(10),ss
   real(wp)  as(10),cs(10)
   real(wp)  ap(10),cp(10)

   basis%hdiag(1:basis%nbf)=1.d+42
   ! note: Rydbergs are identified by valao(*)=0
   !       polarization by            valao(*)=-1

   ibf=0
   iao=0
   ipr=0
   ish=0
   ok=.true.

   atoms: do iat=1,n
      ati = at(iat)
      basis%shells(1,iat)=ish+1
      basis%fila  (1,iat)=ibf+1
      basis%fila2 (1,iat)=iao+1
      shells: do m=1,ao_n(ati)
         ish = ish+1
         ! principle QN
         npq=ao_pqn(m,ati)
         l=ao_l(m,ati)

         level = ao_lev(m,ati)
         zeta  = ao_exp(m,ati)

         basis%lsh(ish) = l
         if (l.eq.11) basis%lsh(ish) = 0
         basis%ash(ish) = iat
         basis%sh2bf(1,ish) = ibf
         basis%sh2ao(1,ish) = iao
         basis%caoshell(m,iat)=ibf
         basis%saoshell(m,iat)=iao

         ! add new shellwise information, for easier reference
         basis%level(ish) = level
         basis%zeta (ish) = zeta

         if (ao_typ(m,ati).eq.-1) then
            valao = -1
         else
            valao = 1
         endif

         ! H-He
         if(l.eq.0.and.ati.le.2.and.npq.eq.1)then
            !  s
            call setsto3(nprim,npq,1,zeta,a,c)  !GFN0
            basis%minalp(ish) = minval(a(:nprim))

            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 1
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=c(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = 1
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         if(l.eq.0.and.ati.le.2.and.npq.eq.2)then
            ! diff s
            call setsto2(thisprimR,npq,1,zeta,aR,cR)
            call atovlp(0,nprim,thisprimR,a,aR,c,cR,ss)
            basis%minalp(ish) = min(minval(a(:nprim)),minval(aR(:thisprimR)))
            ibf =ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 0
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=a(p)
               basis%cont(ipr)=-ss*c(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum),&
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = 0
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! p polarization
         if(l.eq.1.and.ati.le.2)then
            call set_p_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,3,zeta,level,-1)
         endif

         ! general sp
         if(l.eq.0.and.ati.gt.2)then
            ! s
            if(npq.gt.5) then
               call setsto6(nprim,npq,1,zeta,as,cs)
            else
               call setsto4(nprim,npq,1,zeta,as,cs) !GFN0
            endif
            basis%minalp(ish) = minval(as(:nprim))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 1
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = nprim
            basis%hdiag    (ibf) = level

            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=cs(p)
            enddo

            iao = iao+1
            basis%valao2(iao) = 1
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif
         ! p
         if(l.eq.1.and.ati.gt.2)then
            if(npq.gt.5) then
               call set_p_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,6,zeta,level,1)
            else
               call set_p_function(basis,iat,ish,iao,ibf,ipr, &
                  &                npq,l,3,zeta,level,1)
            endif
         endif

         ! DZ s
         if(l.eq.11)then
            if(npq.gt.5) then
               call setsto6(thisprimR,npq,1,zeta,aR,cR)
            else
               call setsto3(thisprimR,npq,1,zeta,aR,cR)
            endif
            call atovlp(0,nprim,thisprimR,as,aR,cs,cR,ss)
            basis%minalp(ish) = min(minval(as(:nprim)),minval(aR(:thisprimR)))

            ibf=ibf+1
            basis%primcount(ibf) = ipr
            basis%valao    (ibf) = 0
            basis%aoat     (ibf) = iat
            basis%lao      (ibf) = 1
            basis%nprim    (ibf) = thisprimR+nprim
            basis%hdiag    (ibf) = level

            idum=ipr+1
            do p=1,thisprimR
               ipr=ipr+1
               basis%alp (ipr)=aR(p)
               basis%cont(ipr)=cR(p)
            enddo
            do p=1,nprim
               ipr=ipr+1
               basis%alp (ipr)=as(p)
               basis%cont(ipr)=-ss*cs(p)
            enddo
            call atovlp(0,basis%nprim(ibf),basis%nprim(ibf), &
               &        basis%alp(idum),basis%alp(idum),&
               &        basis%cont(idum),basis%cont(idum),ss)
            do p=1,basis%nprim(ibf)
               basis%cont(idum-1+p)=basis%cont(idum-1+p)/sqrt(ss)
            enddo

            iao = iao+1
            basis%valao2(iao) = 0
            basis%aoat2 (iao) = iat
            basis%lao2  (iao) = 1
            basis%hdiag2(iao) = level
            basis%aoexp (iao) = zeta
            basis%ao2sh (iao) = ish
         endif

         ! d
         if(l.eq.2)then
            call set_d_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,4,zeta,level,valao)
         endif

         ! f
         if(l.eq.3)then
            call set_f_function(basis,iat,ish,iao,ibf,ipr, &
               &                npq,l,4,zeta,level,1)
         endif

         basis%valsh(ish) = basis%valao2(iao)

         basis%sh2bf(2,ish) = ibf-basis%sh2bf(1,ish)
         basis%sh2ao(2,ish) = iao-basis%sh2ao(1,ish)
      enddo shells
      basis%shells(2,iat)=ish
      basis%fila  (2,iat)=ibf
      basis%fila2 (2,iat)=iao
   enddo atoms

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

   if(basis%nao.ne.iao) then
      call raise('E','Internal error',1)
   endif

   do iat=1,n
      ati = at(iat)
      do j=1,ao_n(ati)
         l = ao_l(j,ati)
         if(l.eq.11)ao_l(j,ati)=0
      enddo
   enddo

end subroutine xbasis_gfn0

! ------------------------------------------------------------------------
!  Helper functions

subroutine atovlp(l,npri,nprj,alpa,alpb,conta,contb,ss)
  use mctc_constants, only: pi
  implicit none
  integer l,npri,nprj
  real(wp) alpa(*),alpb(*)
  real(wp) conta(*),contb(*)
  real(wp) ss

  integer ii,jj
  real(wp) ab,s00,sss,ab05

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

subroutine set_s_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   use tbdef_basisset
   implicit none
   type(tb_basisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: p
   real(wp) :: alp(10),cont(10)

   call expand_sto(nprim,npq,l,zeta,alp,cont)
   basis%minalp(ish) = minval(alp(:nprim))

   ibf = ibf+1
   basis%primcount(ibf) = ipr
   basis%valao    (ibf) = valao
   basis%aoat     (ibf) = iat
   basis%lao      (ibf) = 1
   basis%nprim    (ibf) = nprim
   basis%hdiag    (ibf) = level

   do p=1,nprim
      ipr = ipr+1
      basis%alp (ipr)=alp (p)
      basis%cont(ipr)=cont(p)
   enddo

   iao = iao+1
   basis%valao2(iao) = valao
   basis%aoat2 (iao) = iat
   basis%lao2  (iao) = 1
   basis%hdiag2(iao) = level
   basis%aoexp (iao) = zeta
   basis%ao2sh (iao) = ish
end subroutine set_s_function

subroutine set_p_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   use tbdef_basisset
   implicit none
   type(tb_basisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)

   call expand_sto(nprim,npq,l,zeta,alp,cont)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 2, 4

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)
      enddo

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_p_function

subroutine set_d_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   use tbdef_basisset
   implicit none
   type(tb_basisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)
   real(wp) :: trafo(5:10) = &
      & [1.0_wp, 1.0_wp, 1.0_wp, sqrt(3.0_wp), sqrt(3.0_wp), sqrt(3.0_wp)]

   call expand_sto(nprim,npq,l,zeta,alp,cont)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 5, 10

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)*trafo(j)
      enddo

      if (j .eq. 5) cycle

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j-1
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_d_function

subroutine set_f_function(basis,iat,ish,iao,ibf,ipr,npq,l,nprim,zeta,level,valao)
   use tbdef_basisset
   implicit none
   type(tb_basisset), intent(inout) :: basis
   integer, intent(in)    :: iat
   integer, intent(in)    :: ish
   integer, intent(inout) :: ibf
   integer, intent(inout) :: iao
   integer, intent(inout) :: ipr
   integer, intent(in)    :: npq
   integer, intent(in)    :: l
   integer, intent(in)    :: nprim
   integer, intent(in)    :: valao
   real(wp),intent(in)    :: zeta
   real(wp),intent(in)    :: level
   integer  :: j,p
   real(wp) :: alp(10),cont(10)
   real(wp) :: trafo(11:20) = &
      & [1.0_wp, 1.0_wp, 1.0_wp, sqrt(5.0_wp), sqrt(5.0_wp), &
      &  sqrt(5.0_wp), sqrt(5.0_wp), sqrt(5.0_wp), sqrt(5.0_wp), sqrt(15.0_wp)]

   call expand_sto(nprim,npq,l,zeta,alp,cont)
   basis%minalp(ish) = minval(alp(:nprim))

   do j = 11, 20

      ibf = ibf+1
      basis%primcount(ibf) = ipr
      basis%valao    (ibf) = valao
      basis%aoat     (ibf) = iat
      basis%lao      (ibf) = j
      basis%nprim    (ibf) = nprim
      basis%hdiag    (ibf) = level

      do p=1,nprim
         ipr = ipr+1
         basis%alp (ipr)=alp (p)
         basis%cont(ipr)=cont(p)*trafo(j)
      enddo

      if (j.ge.11 .and. j.le.13) cycle

      iao = iao+1
      basis%valao2(iao) = valao
      basis%aoat2 (iao) = iat
      basis%lao2  (iao) = j-3
      basis%hdiag2(iao) = level
      basis%aoexp (iao) = zeta
      basis%ao2sh (iao) = ish

   enddo
end subroutine set_f_function

subroutine expand_sto(nprim,npq,l,zeta,alp,cont)
   implicit none
   integer, intent(in)  :: nprim
   integer, intent(in)  :: npq
   integer, intent(in)  :: l
   real(wp),intent(in)  :: zeta
   real(wp),intent(out) :: alp(10)
   real(wp),intent(out) :: cont(10)
   integer  :: dum

   select case(nprim)
   case(2); call setsto2(dum,npq,l+1,zeta,alp,cont)
   case(3); call setsto3(dum,npq,l+1,zeta,alp,cont)
   case(4); call setsto4(dum,npq,l+1,zeta,alp,cont)
   case(6); call setsto6(dum,npq,l+1,zeta,alp,cont)
   end select

end subroutine expand_sto

end module xbasis
