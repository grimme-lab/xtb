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

module xtb_basic_geo
contains
!  .....................................................................

      Subroutine dthetadr(nat,xyz,j,k,i,theta,kijk,dei,dej,dek)
         use xtb_mctc_accuracy, only : wp

!  .....................................................................
!
!   Input
!   -----
!                              __              __
!    ra,rb          : distance ij and distnace ik
!
!    sinth          : sinus of angle j-i-k
!
!    costh          : cosinus of angle j-i-k
!
!
!   Output
!   ------
!
!    dei            : derivatives of the atom i
!
!    dej            : derivatives of the atom j
!
!    dek            : derivatives of the atom k
!
!  .....................................................................

      implicit none

      external vecnorm

      integer ic,nat,j,k,i

      real(wp) ran,rbn,theta,xyz(3,nat)
      real(wp) ra(3),rb(3),sinth,costh
      real(wp) dei(3),dej(3),dek(3),kijk
      real(wp) fac,vecnorm,eps

      parameter (eps=1.d-14)

! ... dtheta/dri, dtheta/drj and dtheta/drk

      sinth=sin(theta)
      costh=cos(theta)

      do ic=1,3
        ra(ic)=xyz(ic,j)-xyz(ic,i)
        rb(ic)=xyz(ic,k)-xyz(ic,i)
      end do

      ran=vecnorm(ra,3,0)
      rbn=vecnorm(rb,3,0)
      if (abs(sinth).le.eps) then
        fac=sign(kijk,costh)/(ran*rbn)
      else
        fac=1.d0/(sinth*ran*rbn)
      end if

        do ic=1,3
! ... dedri
          dei(ic)=fac*( ra(ic)+rb(ic)-(rbn*rbn*ra(ic)+ran*ran*rb(ic))* &
                        costh/(ran*rbn))
! ... dedrj
          dej(ic)=fac*(ra(ic)*costh*rbn/ran-rb(ic))
! ... dedrk
          dek(ic)=fac*(rb(ic)*costh*ran/rbn-ra(ic))
        end do

      End Subroutine dthetadr

      Subroutine dphidr(nat,xyz,i,j,k,l,phi, &
                        dphidri,dphidrj,dphidrk,dphidrl)
                     use xtb_mctc_accuracy, only : wp
!     the torsion derivatives

      implicit none

      external vecnorm

      integer  ic,i,j,k,l,nat

      real(wp) sinphi,cosphi,onenner,thab,thbc
      real(wp) ra(3),rb(3),rc(3),rab(3),rac(3),rbc(3),rbb(3)
      real(wp) raa(3),rba(3),rapba(3),rapbb(3),rbpca(3),rbpcb(3)
      real(wp) rapb(3),rbpc(3),na(3),nb(3),nan,nbn
      real(wp) dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3)
      real(wp) xyz(3,nat),phi,vecnorm,nenner,eps,vz

      parameter (eps=1.d-14)

      cosphi=cos(phi)
      sinphi=sin(phi)
      do ic=1,3
        ra(ic)=xyz(ic,j)-xyz(ic,i)
        rb(ic)=xyz(ic,k)-xyz(ic,j)
        rc(ic)=xyz(ic,l)-xyz(ic,k)

        rapb(ic)=ra(ic)+rb(ic)
        rbpc(ic)=rb(ic)+rc(ic)
      end do

      call crossprod(ra,rb,na)
      call crossprod(rb,rc,nb)
      nan=vecnorm(na,3,0)
      nbn=vecnorm(nb,3,0)

      nenner=nan*nbn*sinphi
      if (abs(nenner).lt.eps) then
        dphidri=0
        dphidrj=0
        dphidrk=0
        dphidrl=0
        onenner=1.0d0/(nan*nbn)
      else
        onenner=1.d0/nenner
      endif

      call crossprod(na,rb,rab)
      call crossprod(nb,ra,rba)
      call crossprod(na,rc,rac)
      call crossprod(nb,rb,rbb)
      call crossprod(nb,rc,rbc)
      call crossprod(na,ra,raa)

      call crossprod(rapb,na,rapba)
      call crossprod(rapb,nb,rapbb)
      call crossprod(rbpc,na,rbpca)
      call crossprod(rbpc,nb,rbpcb)

! ... dphidri
      do ic=1,3
        dphidri(ic)=onenner*(cosphi*nbn/nan*rab(ic)-rbb(ic))

! ... dphidrj
         dphidrj(ic)=onenner*(cosphi*(nbn/nan*rapba(ic) &
                                      +nan/nbn*rbc(ic)) &
                              -(rac(ic)+rapbb(ic)))
! ... dphidrk
        dphidrk(ic)=onenner*(cosphi*(nbn/nan*raa(ic) &
                                   +nan/nbn*rbpcb(ic)) &
                              -(rba(ic)+rbpca(ic)))
! ... dphidrl
         dphidrl(ic)=onenner*(cosphi*nan/nbn*rbb(ic)-rab(ic))
       end do

      End Subroutine dphidr

!  .....................................................................

      Subroutine domegadr &
                 (nat,xyz, &
                  i,j,k,l,omega, &
                  domegadri,domegadrj,domegadrk,domegadrl)
               use xtb_mctc_accuracy, only : wp

!     inversion derivatives
!  .....................................................................

      implicit none

      external vecnorm

      integer  ic,i,j,k,l,nat

      real(wp) omega,sinomega
      real(wp) xyz(3,nat),onenner,vecnorm,rnn,rvn
      real(wp) rn(3),rv(3),rd(3),re(3),rdme(3),rve(3)
      real(wp) rne(3),rdv(3),rdn(3)
      real(wp) rvdme(3),rndme(3),nenner
      real(wp) domegadri(3),domegadrj(3),domegadrk(3),domegadrl(3),eps

      parameter (eps=1.d-14)

      sinomega=sin(omega)

      do ic=1,3
        rv(ic)=xyz(ic,l)-xyz(ic,i)
        rd(ic)=xyz(ic,k)-xyz(ic,j)
        re(ic)=xyz(ic,i)-xyz(ic,j)

        rdme(ic)=rd(ic)-re(ic)
      end do

      call crossprod(re,rd,rn)
      rvn=vecnorm(rv,3,0)
      rnn=vecnorm(rn,3,0)

      call crossprod(rv,re,rve)
      call crossprod(rn,re,rne)
      call crossprod(rd,rv,rdv)
      call crossprod(rd,rn,rdn)
      call crossprod(rv,rdme,rvdme)
      call crossprod(rn,rdme,rndme)

      nenner=rnn*rvn*cos(omega)
      if (abs(nenner).gt.eps) then
        onenner=1.d0/nenner
        do ic=1,3
! ... domega/dri
          domegadri(ic)=onenner*(rdv(ic)-rn(ic)- &
                             sinomega*(rvn/rnn*rdn(ic)-rnn/rvn*rv(ic)))

! ... domega/drj
          domegadrj(ic)=onenner*(rvdme(ic)-sinomega*rvn/rnn*rndme(ic))

! ... domega/drk
          domegadrk(ic)=onenner*(rve(ic)-sinomega*rvn/rnn*rne(ic))

! ... domega/drl
          domegadrl(ic)=onenner*(rn(ic)-sinomega*rnn/rvn*rv(ic))
        end do
      else
        do ic=1,3
          domegadri(ic)=0.d0
          domegadrj(ic)=0.d0
          domegadrk(ic)=0.d0
          domegadrl(ic)=0.d0
        end do
      end if

      End Subroutine domegadr

!  .....................................................................

      real(wp) Function valijkl(nat,xyz,i,j,k,l)
         use xtb_mctc_accuracy, only : wp
      use xtb_mctc_constants, only: pi

!  .....................................................................

      implicit none

      external vecnorm,valijk

      integer ic,i,j,k,l,nat

      real(wp) xyz(3,nat)
      real(wp) eps,ra(3),rb(3),rc(3),na(3),nb(3)
      real(wp) rab,rbc,thab,thbc,valijk
      real(wp) vecnorm,nan,nbn,rcn,snanb,deter,test
      real(wp) raabs,rbabs,rcabs

      parameter (eps=1.0d-14)

      raabs=0.0d0
      rbabs=0.0d0
      rcabs=0.0d0
! ... get torsion coordinate
      do ic=1,3
        ra(ic)=xyz(ic,j)-xyz(ic,i)
        rb(ic)=xyz(ic,k)-xyz(ic,j)
        rc(ic)=xyz(ic,l)-xyz(ic,k)
        raabs=raabs+ra(ic)*ra(ic)
        rbabs=rbabs+rb(ic)*rb(ic)
        rcabs=rcabs+rc(ic)*rc(ic)
      end do
      raabs=sqrt(raabs)
      rbabs=sqrt(rbabs)
      rcabs=sqrt(rcabs)
! ... normalize vectors
      do ic=1,3
        ra(ic)=ra(ic)/raabs
        rb(ic)=rb(ic)/rbabs
        rc(ic)=rc(ic)/rcabs
      end do
! ... determinante of rb,ra,rc
      deter= ra(1)*(rb(2)*rc(3)-rb(3)*rc(2)) &
            -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1)) &
            +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

      call crossprod(-ra,rc,na)
      test=0.0d0
      do ic=1,3
        test=test+na(ic)*rb(ic)
      end do

      thab=valijk(nat,xyz,i,k,j)
      thbc=valijk(nat,xyz,j,l,k)
      call crossprod(ra,rb,na)
      call crossprod(rb,rc,nb)
      nan=vecnorm(na,3,1)
      nbn=vecnorm(nb,3,1)

      snanb=0.0d0
      do ic=1,3
        snanb=snanb+na(ic)*nb(ic)
      end do
      if (abs(abs(snanb)-1.d0).lt.eps) then
        snanb=sign(1.d0,snanb)
      end if

      valijkl=acos(snanb)
      if(test.lt.0.0d0) valijkl=-acos(snanb)

! the gradient dphir is only compatible with this subroutine
! if the statement below is commented out. If not, opt. and
! Hessian show large errors and imags. I don't understand
! this entirely but thats how it is.
! SG, Sat May 24 11:41:42 CEST 2014

!     if (deter.lt.0) then
!        valijkl=2.d0*pi-valijkl
!     end if

      End Function valijkl

!  .....................................................................

      real(wp) Function valijk(nat,xyz,j,k,i)
         use xtb_mctc_accuracy, only : wp

!  .....................................................................

      implicit none

      external vecnorm

      integer nat,j,k,i,ic

      real(wp) ra(3),rb(3),rab,eps
      real(wp) xyz(3,nat),vecnorm,ran,rbn

      parameter (eps=1.d-14)

      do ic=1,3
        ra(ic)=xyz(ic,j)-xyz(ic,i)
        rb(ic)=xyz(ic,k)-xyz(ic,i)
      end do

      ran=vecnorm(ra,3,1)
      rbn=vecnorm(rb,3,1)
      rab=0.d0
      do ic=1,3
        rab=rab+ra(ic)*rb(ic)
      end do

      if (abs(abs(rab)-1.d0).lt.eps) then
        rab=sign(1.d0,rab)
      end if
      valijk=acos(rab)

      End Function valijk

!  .....................................................................

      real(wp) Function omega (nat,xyz,i,j,k,l)
         use xtb_mctc_accuracy, only : wp

!   Calculates the inversion angle
!  .....................................................................

      implicit none

      external vecnorm

      integer ic,nat,i,j,k,l

      real(wp) xyz(3,nat)
      real(wp) rd(3),re(3),rn(3),rv(3),rnv
      real(wp) vecnorm,rkjn,rljn,rnn,rvn

      do ic=1,3
        re(ic)=xyz(ic,i)-xyz(ic,j)
        rd(ic)=xyz(ic,k)-xyz(ic,j)
        rv(ic)=xyz(ic,l)-xyz(ic,i)
      end do
      call crossprod(re,rd,rn)
      rnn=vecnorm(rn,3,1)
      rvn=vecnorm(rv,3,1)

      rnv=rn(1)*rv(1)+rn(2)*rv(2)+rn(3)*rv(3)
      omega=asin( rnv )

      End Function omega

!  .....................................................................

      Subroutine crossprod(ra,rb,rab)
         use xtb_mctc_accuracy, only : wp

      implicit none

      real(wp) ra(3),rb(3),rab(3)

      rab(1)=ra(2)*rb(3)-ra(3)*rb(2)
      rab(2)=ra(3)*rb(1)-ra(1)*rb(3)
      rab(3)=ra(1)*rb(2)-ra(2)*rb(1)

      End Subroutine crossprod

!  .....................................................................

      real(wp) Function  vecnorm (r,n,inorm)
         use xtb_mctc_accuracy, only : wp

      implicit none

      integer i,n,inorm

      real(wp) r(n),or,sp,rn

      sp=0.
      do i=1,n
        sp=sp+r(i)*r(i)
      end do
      rn=sqrt(sp)
      if (inorm.gt.0) then
        if (abs(rn).gt.1.d-14) then
          or=1.0d0/rn
          do i=1,n
            r(i)=or*r(i)
          end do
        end if
      end if
      vecnorm=rn

      End Function vecnorm

!  .....................................................................

      subroutine crprod(a,b,c)
         use xtb_mctc_accuracy, only : wp
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),c(3)

      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)

      return
      end subroutine crprod

!  .....................................................................

      subroutine vsc1(a,scale,tol)
         use xtb_mctc_accuracy, only : wp
      implicit integer (i-n)
      implicit double precision (a-h,o-z)
      dimension a(3)

      rlen = norm2(a)

      if (rlen.gt.tol) then
         do i=1,3
            a(i) = a(i)*scale/rlen
         end do
      endif

      return
      end subroutine vsc1

!  .....................................................................

      pure subroutine impsc(a,b,c)
         use xtb_mctc_accuracy, only : wp
      implicit none
      real(wp),intent(in)  :: a(3),b(3)
      real(wp),intent(out) :: c
      integer  :: i
      real(wp) :: rimp,al,bl

      rimp = 0.0d0

      do i=1,3
         rimp = rimp + a(i)*b(i)
      end do

      al = norm2(a)
      bl = norm2(b)

      if (al.gt.0.0d0.and.bl.gt.0.0d0) then
         c = rimp/(al*bl)
      else
         c = 0.0d0
      endif

      return
      end subroutine impsc

!  .....................................................................

      pure subroutine bangl(xyz,i,j,k,angle)
         use xtb_mctc_accuracy, only : wp
      implicit none
      real(wp),intent(in)  :: xyz(3,*)
      real(wp),intent(out) :: angle
      integer, intent(in)  :: i,j,k

      real(wp) d2ij,d2jk,d2ik,xy,temp

      d2ij = sum((xyz(:,i)-xyz(:,j))**2)
      d2jk = sum((xyz(:,j)-xyz(:,k))**2)
      d2ik = sum((xyz(:,i)-xyz(:,k))**2)
      xy = sqrt(d2ij*d2jk)
      temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
      if (temp .gt. 1.0d0)  temp= 1.0d0
      if (temp .lt. -1.0d0) temp=-1.0d0
      angle = acos( temp )

      end subroutine bangl
end module xtb_basic_geo
