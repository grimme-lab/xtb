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

module aespot
   use iso_fortran_env, only : wp => real64
   use intpack, only : olap,divpt,rhftce,prod,opab1,opab4,propa
   implicit none
   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

contains

! setdqlist:  precomputes the dipole and quadrupole potential terms for Fock matrix
! thr              : neglect matrix elements in potential
! nao              : # of spherical AOs (SAOs)
! dpint            : dipole integral matrix, dimension 3,nao*(nao+1)/2
! qpint            : quadrupole integral matrix, dimension 6,nao*(nao+1)/2
! ndp,nqp          : number of elements to be computed in Fock matrix with X-dip and X-qpole terms
! matdlst,matqlst  : index list, to which AO, the ndp/nqp potential terms refer to
subroutine setdqlist(nao,ndp,nqp,thr,dpint,qpint,matdlst,matqlst)
   implicit none
   integer, intent(in)    :: nao
   integer, intent(inout) :: ndp,nqp
   real(wp),intent(in)  ::  thr
   real(wp),intent(in)  :: dpint(3,nao,nao)
   real(wp),intent(in)  :: qpint(6,nao,nao)
   integer, intent(inout):: matqlst(2,nqp),matdlst(2,ndp)

   real(wp) tmp1,tmp2,tmp3,tmp4
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f
   ! stuff for potential

   integer i,j,k,l,m,ii,jj,ll,mm,kk,mq,md,ij

   ! INFO: this threshold must be slightly larger than max(0,thr2),
   !       where thr2 is the one used in screening in routine aesdqint
   thr2=thr*1.0d-2 ! we compare squared int-elements
   !      thr2=1.0d-20 ! conservative, keep all terms

   md=0
   mq=0
   ! set uo matrix lists
   ij=0
   do i=1,nao
      do j=1,i
         ij=ij+1
         tmp1=0.0_wp
         tmp2=0.0_wp
         kk=0
         do k=1,3
            tmp1=tmp1+dpint(k,j,i)*dpint(k,j,i)
            tmp2=tmp2-qpint(k,j,i)*qpint(k,j,i)
         enddo
         do k=1,6
            tmp2=tmp2+2.0_wp*qpint(k,j,i)*qpint(k,j,i)
         enddo
         if(tmp1.gt.thr2)then
            md=md+1
            matdlst(1,md)=int(i,2)
            matdlst(2,md)=int(j,2)
         endif
         if(tmp2.gt.thr2)then
            mq=mq+1
            matqlst(1,mq)=int(i,2)
            matqlst(2,mq)=int(j,2)
         endif
      enddo
   enddo
   ndp=md
   nqp=mq
end subroutine setdqlist

! scalecamm: scale all anisotropic CAMMs by element-specific parameters
! nat              : # of atoms
! at(nat)          : atom-to-element identifier
! dipm(3,nat)      : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)        : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
subroutine scalecamm(nat,at,dipm,qp)
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp), intent(inout):: dipm(3,nat),qp(6,nat)
   integer i,iat

   ! CAMM  scaling
   do i=1,nat
      qp(1:6,i)=qp(1:6,i)*(3./3.)
   enddo
end subroutine scalecamm

! unscalecamm: unscale all anisotropic CAMMs from element-specific parameters
! nat              : # of atoms
! at(nat)          : atom-to-element identifier
! dipm(3,nat)      : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)        : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
subroutine unscalecamm(nat,at,dipm,qp)
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp), intent(inout):: dipm(3,nat),qp(6,nat)
   real(wp) aesi
   integer i,iat

   ! CAMM  scaling
   do i=1,nat
      qp(1:6,i)=qp(1:6,i)*(3./3.)
   enddo
end subroutine unscalecamm




! mmpop:  compute the cumulative atomic dipole and quadrupole moments via Mulliken population analysis
! nat              : # of atoms
! nao              : # of spherical AOs (SAOs)
! aoat2(nao)       : SAO to atom intex
! s(nao,nao)       : overlap matrix
! xyz(3,nat)       : cartesian coordinates
! dpint            : dipole integral matrix, dimension 3,nao*(nao+1)/2
! qpint            : quadrupole integral matrix, dimension 6,nao*(nao+1)/2
! p(nao,nao)       : density matrix
! dipm(3,nat)      : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)        : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
subroutine mmompop(nat,nao,aoat2,xyz,p,s,dpint,qpint,dipm,qp)
   implicit none
   integer, intent(in) :: nao,nat,aoat2(nao)
   real(wp), intent(in) :: dpint(3,nao,nao),s(nao,nao)
   real(wp), intent(in) :: qpint(6,nao,nao),p(nao,nao)
   real(wp), intent(in) :: xyz(3,nat)
   real(wp), intent(out):: dipm(3,nat),qp(6,nat)

   real(wp) xk1,xl1,xk2,xl2,pij,tii,tjj
   real(wp) pqm,pdmk,pdml,ps,ra(3)

   integer i,j,k,l,ii,jj,ij,kl,kj,lin
   ! CAMM
   dipm=0.0_wp
   qp=0.0_wp
   ij=0
   do i=1,nao
      ii=aoat2(i)
      ra(1:3)=xyz(1:3,ii)
      do j=1,i-1
         ij=ij+1
         jj=aoat2(j)
         pij=p(j,i)
         ps=pij*s(j,i)
         kl=0
         !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
         !  when doing the Mulliken population, we switch to lin-compatible sorting
         !  i,e. xx,xy,yy,xz,yz,zz
         do k=1,3
            xk1=ra(k)
            xk2=xyz(k,jj)
            pdmk=pij*dpint(k,j,i)
            dipm(k,jj)=dipm(k,jj)+xk2*ps-pdmk
            dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
            ! off-diagonal
            do l=1,k-1
               kl=kl+1
               kj=k+l+1
               xl1=ra(l)
               xl2=xyz(l,jj)
               pdml=pij*dpint(l,j,i)
               pqm=pij*qpint(kj,j,i)
               tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
               tjj=pdmk*xl2+pdml*xk2-xl2*xk2*ps-pqm
               qp(kl,jj)=qp(kl,jj)+tjj
               qp(kl,ii)=qp(kl,ii)+tii
            enddo
            ! diagonal
            kl=kl+1
            pqm=pij*qpint(k,j,i)
            tii=2.0_wp*pdmk*xk1-xk1*xk1*ps-pqm
            tjj=2.0_wp*pdmk*xk2-xk2*xk2*ps-pqm
            qp(kl,jj)=qp(kl,jj)+tjj
            qp(kl,ii)=qp(kl,ii)+tii
         enddo
      enddo
      ij=ij+1
      pij=p(i,i)
      ps=pij*s(i,i)
      kl=0
      !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
      !  when doing the Mulliken population, we switch to lin-compatible sorting
      !  i,e. xx,xy,yy,xz,yz,zz
      do k=1,3
         xk1=ra(k)
         pdmk=pij*dpint(k,i,i)
         dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
         ! off-diagonal
         do l=1,k-1
            kl=kl+1
            kj=k+l+1 ! the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
            xl1=ra(l)
            pdml=pij*dpint(l,i,i)
            pqm=pij*qpint(kj,i,i)
            tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
            qp(kl,ii)=qp(kl,ii)+tii
         enddo
         !diagonal
         kl=kl+1
         pqm=pij*qpint(k,i,i)
         tii=2.0_wp*pdmk*xk1-xk1*xk1*ps-pqm
         qp(kl,ii)=qp(kl,ii)+tii
      enddo
   enddo
   ! remove trace
   do i=1,nat
      tii=qp(1,i)+qp(3,i)+qp(6,i)
      tii=0.50_wp*tii
      qp(1:6,i)=1.50_wp*qp(1:6,i)
      qp(1,i)=qp(1,i)-tii
      qp(3,i)=qp(3,i)-tii
      qp(6,i)=qp(6,i)-tii
   enddo
end subroutine mmompop


! distributed atomic multipole moment interactions: all interactions up to r**-3
! energy evaluation
! nat         : # of atoms
! xyz(3,nat)  : cartesian coordinates
! q(nat)      : atomic partial charges
! dipm(3,nat) : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! gab3,gab5   : damped R**-3 and R**-5 Coulomb laws, dimension: nat*(nat+1)/2
!               multiplication with numerator then leads to R**-2 and R**-3 decay, respectively
! e           : E_AES
subroutine aniso_electro(nat,at,xyz,q,dipm,qp,gab3,gab5,e,epol)
   use aoparam
   use lin_mod, only : lin
   implicit none
   integer nat,at(nat)
   real(wp) xyz(3,nat),q(nat)
   real(wp) e
   real(wp) qp1(6),rr(3),dp1(3),rij(3)
   real(wp) edd,e01,e02,e11,r2,tt,tt3,q1,qs2
   real(wp) ed,eq,epol
   ! stuff for potential
   real(wp) gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
   real(wp) dipm(3,nat),qp(6,nat)
   real(wp), parameter :: bohr=1.0_wp/0.52917726_wp

   integer i,j,k,l,m,ki,kj,kl
   e=0.0_wp
   epol=0.0_wp
   e01=0.0_wp
   e02=0.0_wp
   e11=0.0_wp
   do i=1,nat
      q1=q(i)
      rr(1:3)=xyz(1:3,i)
      dp1(1:3)=dipm(1:3,i)
      qp1(1:6)=qp(1:6,i)
      ! test: semilocal CT correction
      ! dipole
      tt=dp1(1)*dp1(1)+dp1(2)*dp1(2)+dp1(3)*dp1(3)
      ! qpole
      tt3=0.0_wp
      do k=1,3
         do l=1,3
            kl=lin(l,k)
            tt3=tt3+qp1(kl)*qp1(kl)
         enddo
      enddo
      epol=epol+dpolc(at(i))*tt+tt3*qpolc(at(i))
      ! ---
      do j=1,i-1             ! loop over all atoms
         kj=lin(j,i)
         rij(1:3)=xyz(1:3,j)-rr(1:3)
         r2=sum(rij*rij)
         ed=0.0_wp
         eq=0.0_wp
         edd=0.0_wp
         !           dipole - charge
         do k=1,3
            ed=ed+q(j)*dp1(k)*rij(k)
            ed=ed-dipm(k,j)*q1*rij(k)
            !              dip-dip & charge-qpole
            do l=1,3
               kl=lin(l,k)
               tt=rij(l)*rij(k)
               tt3=3.0_wp*tt
               eq=eq+q(j)*qp1(kl)*tt
               eq=eq+qp(kl,j)*q1*tt
               edd=edd-dipm(k,j)*dp1(l)*tt3
            enddo
            !              diagonal dip-dip term
            edd=edd+dipm(k,j)*dp1(k)*r2
         enddo
         e01=e01+ed*gab3(kj)
         e02=e02+eq*gab5(kj)
         e11=e11+edd*gab5(kj)
      enddo
   enddo
   e = e01 + e02 + e11
   !     write(*,'(''d,q,dd'',3f9.5)')  e01,e02,e11
   !      write(*,*) ' semilocal CT corr.: ',epol

end subroutine aniso_electro

! aniso-electro from Fock matrix elements
! nao              : # of spherical AOs (SAOs)
! s(nao,nao)       : overlap matrix
! aoat2(nao)       : SAO to atom intex
! dpint            : dipole integral matrix, dimension 3,nao*(nao+1)/2
! qpint            : quadrupole integral matrix, dimension 6,nao*(nao+1)/2
! p(nao,nao)       : density matrix
! vs(nat)          : overlap proportional potential
! vd(3,nat)        : dipint proportional potential
! vq(6,nat)        : quadrupole proportional potential
subroutine fockelectro(nat,nao,aoat2,p,s,dpint,qpint,vs,vd,vq,e)
   use lin_mod, only : lin
   implicit none
   integer, intent(in) :: nat,nao,aoat2(nao)
   real(wp), intent(in) :: dpint(3,nao,nao),s(nao,nao)
   real(wp), intent(in) :: qpint(6,nao,nao),p(nao,nao)
   real(wp), intent(in) :: vs(nat),vd(3,nat),vq(6,nat)
   real(wp), intent(out) :: e
   real(wp) eaes,pji,fji
   integer i,j,k,l,ii,jj,ij,kl,kj
   ! CAMM
   eaes=0.0_wp
   ij=0
   do i=1,nao
      ii=aoat2(i)
      do j=1,nao
         ij=lin(j,i)
         jj=aoat2(j)
         fji=0.0_wp
         pji=p(j,i)
         fji=fji+s(j,i)*(vs(ii)+vs(jj))
         do k=1,3
            fji=fji+dpint(k,j,i)*(vd(k,ii)+vd(k,jj))
         enddo
         do k=1,6
            fji=fji+qpint(k,j,i)*(vq(k,ii)+vq(k,jj))
         enddo
         eaes=eaes+pji*fji
      enddo
   enddo
   eaes=0.250_wp*eaes
   !      write(*,*) 'EAES',eaes
   e=eaes
end subroutine fockelectro


! set-up potential terms v, which are proportional to s, d, or q-integrals
! it is to be multiplied with Sji when stting up Fji (hence, termed vs)
! comes essentially at no cost, once cumulative atomic quadrupole moments are available.
! NEW: the CAMMs were already scaled by scalecamm, but the corresponding potential terms
!      including shift terms need to be scaled
! nat         : # of atoms
! at(nat)     : atom-to-element identifier
! xyz(3,nat)  : cartesian coordinates
! q(nat)      : atomic partial charges
! dipm(3,nat) : atomic dipole moments
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! gab3,gab5   : damped Coulomb laws, dimension: nat*(nat+1)/2
! vs(nat)     : s-proportional potential from all atoms acting on atom i
! vd(3,nat)   : dipint-proportional potential from all atoms acting on atom i
! vq(6,nat)   : qpole-int proportional potential from all atoms acting on atom i
subroutine setvsdq(nat,at,xyz,q,dipm,qp,gab3,gab5,vs,vd,vq)
   use aoparam
   use lin_mod, only : lin
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp), intent(in) ::  q(nat),dipm(3,nat)
   real(wp), intent(in) ::  xyz(3,nat),qp(6,nat)
   real(wp), intent(in) :: gab3(nat*(nat+1)/2)
   real(wp), intent(in) :: gab5(nat*(nat+1)/2)
   real(wp), intent(out) :: vs(nat),vd(3,nat),vq(6,nat)
   real(wp) ra(3),dra(3),rb(3),stmp,dum3a,dum5a,t1a,t2a,t3a,t4a,r2a
   real(wp) r2ab,t1b,t2b,t3b,t4b,dum3b,dum5b,dtmp(3),qtmp(6),g3,g5
   real(wp) qs1,qs2
   integer i,j,k,l1,l2,ll,m,mx,ki,kj,mm
   mm=0
   vs=0.0_wp
   vd=0.0_wp
   vq=0.0_wp
   ! set up overlap proportional potential
   do i=1,nat
      ra(1:3)=xyz(1:3,i)
      stmp=0.0_wp
      dtmp=0.0_wp
      qtmp=0.0_wp
      do j=1,nat
         mm=lin(j,i)
         g3=gab3(mm)
         g5=gab5(mm)
         rb(1:3)=xyz(1:3,j)
         dra(1:3)=ra(1:3)-rb(1:3)
         dum3a=0.0_wp ! collect gab3 dependent terms
         dum5a=0.0_wp ! collect gab5 dependent terms
         r2a=0.0_wp
         r2ab=0.0_wp
         t1a=0.0_wp
         t2a=0.0_wp
         t3a=0.0_wp
         t4a=0.0_wp
         ll=0
         do l1=1,3
            ! potential from dipoles
            r2a=r2a+ra(l1)*ra(l1)      ! R_C * R_C
            r2ab=r2ab+dra(l1)*dra(l1)  ! R_AC * R_AC
            t1a=t1a+ra(l1)*dra(l1)     ! R_C * R_AC  : for dip-q (q-shift) and dip-dip (q-shift)
            t2a=t2a+dipm(l1,j)*dra(l1) ! mu_A * R_AC : for q-dip and dip-dip (q-shift)
            t3a=t3a+ra(l1)*dipm(l1,j)  ! R_C * mu_A  : for diag. dip-dip (q-shift)
            t4a=t4a+dra(l1)*dra(l1)*ra(l1)*ra(l1) ! (R_C o R_AC)**"2(square of Hadamard product) :
            ! results from trace remove from q-pole (q-shift)
            do l2=1,3
               ll=lin(l1,l2)
               ! potential from quadrupoles
               dum5a=dum5a-qp(ll,j)*dra(l1)*dra(l2) &
                  & -1.50_wp*q(j)*dra(l1)*dra(l2)*ra(l1)*ra(l2)
               if(l2.ge.l1) cycle
               ki=l1+l2+1
               qtmp(ki)=qtmp(ki)-3.0_wp*q(j)*g5*dra(l2)*dra(l1)
            enddo
            qtmp(l1)=qtmp(l1)-1.50_wp*q(j)*g5*dra(l1)*dra(l1)
         enddo
         !
         ! set up S-dependent potential
         dum3a=-t1a*q(j)-t2a ! dip-q (q-shift) and q-dip
         dum5a=dum5a+t3a*r2ab-3.0_wp*t1a*t2a & !dip-dip (q-shift terms)
            & +0.50_wp*q(j)*r2a*r2ab !qpole-q (q-shift, trace removal)
         stmp=stmp+dum5a*g5+dum3a*g3
         do l1=1,3
            dum3a=dra(l1)*q(j)
            dum5a=3.0_wp*dra(l1)*t2a &           ! dipint-dip
               & -r2ab*dipm(l1,j) &            ! dipint-dip (diagonal)
               & -q(j)*r2ab*ra(l1) &           ! qpole-q (dipint-shift, trace removal)
               & +3.0_wp*q(j)*dra(l1)*t1a   ! qpole-q (dipint-shift)
            dtmp(l1)=dtmp(l1)+dum3a*g3+dum5a*g5
            qtmp(l1)=qtmp(l1)+0.50_wp*r2ab*q(j)*g5 !remove trace term
         enddo
      enddo
      vs(i)=stmp                       ! q terms
      vd(1:3,i)=dtmp(1:3)              ! dipints from atom i
      vq(1:6,i)=qtmp(1:6)              ! qpints from atom i
      ! --- CT correction terms
      qs1=dpolc(at(i))*2.0_wp
      qs2=qpolc(at(i))*6.0_wp ! qpole pot prefactors
      t3a=0.0_wp
      t2a=0.0_wp
      do l1=1,3
         ! potential from dipoles
         t3a=t3a+ra(l1)*dipm(l1,i)*qs1  ! R_C * mu_C  : for diag. dip-dip
         vd(l1,i)=vd(l1,i)-qs1*dipm(l1,i)
         do l2=1,l1-1
            ! potential from quadrupoles
            ll=lin(l1,l2)
            ki=l1+l2+1
            vq(ki,i)=vq(ki,i)-qp(ll,i)*qs2
            t3a=t3a-ra(l1)*ra(l2)*qp(ll,i)*qs2
            vd(l1,i)=vd(l1,i)+ra(l2)*qp(ll,i)*qs2
            vd(l2,i)=vd(l2,i)+ra(l1)*qp(ll,i)*qs2
         enddo
         ! diagonal
         ll=lin(l1,l1)
         vq(l1,i)=vq(l1,i)-qp(ll,i)*qs2*0.50_wp
         t3a=t3a-ra(l1)*ra(l1)*qp(ll,i)*qs2*0.50_wp
         vd(l1,i)=vd(l1,i)+ra(l1)*qp(ll,i)*qs2
         ! collect trace removal terms
         t2a=t2a+qp(ll,i)
      enddo
      vs(i)=vs(i)+t3a
      ! trace removal
      t2a=t2a*qpolc(at(i))
      do l1=1,3
         vq(l1,i)=vq(l1,i)+t2a
         vd(l1,i)=vd(l1,i)-2.0_wp*ra(l1)*t2a
         vs(i)=vs(i)+t2a*ra(l1)*ra(l1)
      enddo
      ! ---
   enddo

   !      call prmat(6,vs,nat,0,'vs')

end subroutine setvsdq

! set-up potential terms v used for nuclear gradients
! here, the shift terms are removed, since we use multipole derivatives w/ origin at the atoms
! nat         : # of atoms
! xyz(3,nat)  : cartesian coordinates
! q(nat)      : atomic partial charges
! dipm(3,nat) : atomic dipole moments
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! gab3,gab5   : damped Coulomb laws, dimension: nat*(nat+1)/2
! vs(nat)     : s-proportional potential from all atoms acting on atom i
! vd(3,nat)   : dipint-proportional potential from all atoms acting on atom i
! vq(6,nat)   : qpole-int proportional potential from all atoms acting on atom i
subroutine setdvsdq(nat,at,xyz,q,dipm,qp,gab3,gab5,vs,vd,vq)
   use aoparam
   use lin_mod, only : lin
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp), intent(in) ::  q(nat),dipm(3,nat)
   real(wp), intent(in) ::  xyz(3,nat),qp(6,nat)
   real(wp), intent(in) :: gab3(nat*(nat+1)/2)
   real(wp), intent(in) :: gab5(nat*(nat+1)/2)
   real(wp), intent(out) :: vs(nat),vd(3,nat),vq(6,nat)
   real(wp) ra(3),dra(3),rb(3),stmp,dum3a,dum5a,t1a,t2a,t3a,t4a,r2a
   real(wp) r2ab,t1b,t2b,t3b,t4b,dum3b,dum5b,dtmp(3),qtmp(6),g3,g5
   real(wp) qs1,qs2
   integer i,j,k,l1,l2,ll,m,mx,ki,kj,mm
   mm=0
   vs=0.0_wp
   vd=0.0_wp
   vq=0.0_wp
   ! set up overlap proportional potential
   do i=1,nat
      ra(1:3)=xyz(1:3,i)
      stmp=0.0_wp
      dtmp=0.0_wp
      qtmp=0.0_wp
      do j=1,nat
         mm=lin(j,i)
         g3=gab3(mm)
         g5=gab5(mm)
         rb(1:3)=xyz(1:3,j)
         dra(1:3)=ra(1:3)-rb(1:3)
         dum3a=0.0_wp ! collect gab3 dependent terms
         dum5a=0.0_wp ! collect gab5 dependent terms
         r2a=0.0_wp
         r2ab=0.0_wp
         t2a=0.0_wp
         ll=0
         do l1=1,3
            ! potential from dipoles
            r2ab=r2ab+dra(l1)*dra(l1)  ! R_AC * R_AC
            t2a=t2a+dipm(l1,j)*dra(l1) ! mu_A * R_AC : for q-dip and dip-dip (q-shift)
            do l2=1,3
               ll=lin(l1,l2)
               ! potential from quadrupoles
               dum5a=dum5a-qp(ll,j)*dra(l1)*dra(l2)
               if(l2.ge.l1) cycle
               ki=l1+l2+1
               qtmp(ki)=qtmp(ki)-3.0_wp*q(j)*g5*dra(l2)*dra(l1)
            enddo
            qtmp(l1)=qtmp(l1)-1.50_wp*q(j)*g5*dra(l1)*dra(l1)
         enddo
         dum3a=-t2a ! q-dip ! w/o shift terms
         stmp=stmp+dum3a*g3+dum5a*g5
         do l1=1,3
            dum3a=dra(l1)*q(j) ! w/o shift terms
            dum5a=3.0_wp*dra(l1)*t2a & ! dipint-dip
               & -r2ab*dipm(l1,j)  ! dipint-dip (diagonal)
            dtmp(l1)=dtmp(l1)+dum3a*g3+dum5a*g5
            qtmp(l1)=qtmp(l1)+0.50_wp*r2ab*q(j)*g5 !remove trace term
         enddo
      enddo
      vs(i)=stmp
      vd(1:3,i)=dtmp(1:3)
      vq(1:6,i)=qtmp(1:6)
      ! --- CT correction terms
      qs1=dpolc(at(i))*2.0_wp
      qs2=qpolc(at(i))*6.0_wp ! qpole pot prefactors
      t2a=0.0_wp
      do l1=1,3
         ! potential from dipoles
         vd(l1,i)=vd(l1,i)-qs1*dipm(l1,i)
         do l2=1,l1-1
            ! potential from quadrupoles
            ll=lin(l1,l2)
            ki=l1+l2+1
            vq(ki,i)=vq(ki,i)-qp(ll,i)*qs2
         enddo
         ll=lin(l1,l1)
         vq(l1,i)=vq(l1,i)-qp(ll,i)*qs2*0.50_wp
         ! collect trace removal terms
         t2a=t2a+qp(ll,i)
      enddo
      ! trace removal
      t2a=t2a*qpolc(at(i))
      do l1=1,3
         vq(l1,i)=vq(l1,i)+t2a
      enddo
      ! ---


   enddo

   !      call prmat(6,vs,nat,0,'vs')

end subroutine setdvsdq

! molmom: computes molecular multipole moments from CAMM
! n           : # of atoms
! xyz(3,n)    : cartesian coordinates
! q(n)        : atomic partial charges
! dipm(3,n)   : cumulative atomic dipole moments (x,y,z)
! qp(6,n)     : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
subroutine molmom(iunit,n,xyz,q,dipm,qp,dip,d3)
   use mctc_econv
   use lin_mod, only : lin
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   real(wp), intent(in)  :: xyz(3,n),q(n),dipm(3,n),qp(6,n)
   real(wp), intent(out) :: dip,d3(3)
   real(wp) rr1(3),rr2(3),tma(6),tmb(6),tmc(6),dum
   integer i,j,k,l
   rr1=0.0_wp
   rr2=0.0_wp
   write(iunit,'(a)')
   do i=1,n
      do j=1,3
         rr1(j)=rr1(j)+q(i)*xyz(j,i)
         rr2(j)=rr2(j)+dipm(j,i)
      enddo
   enddo
   d3(1:3)=rr1(1:3)+rr2(1:3)
   dip=sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
   write(iunit,'(a)',advance='yes')'molecular dipole:'
   write(iunit,'(a)',advance='no')'                 '
   write(iunit,'(a)',advance='yes') &
      & 'x           y           z       tot (Debye)'
   write(iunit,'(a,3f12.3)') ' q only: ',rr1(1:3)
   write(iunit,'(a,4f12.3)') '   full: ',d3(1:3),dip*autod

   tma=0.0_wp
   tmb=0.0_wp
   tmc=0.0_wp
   do i=1,n
      l=0
      do j=1,3
         do k=1,j
            l=lin(k,j)
            tma(l)=tma(l)+xyz(j,i)*xyz(k,i)*q(i)
            tmb(l)=tmb(l)+dipm(k,i)*xyz(j,i)+dipm(j,i)*xyz(k,i)
            tmc(l)=tmc(l)+qp(l,i)
         enddo
      enddo
   enddo
   ! remove traces and multiply with 3/2 in q and dip parts
   dum=tma(1)+tma(3)+tma(6)
   dum=0.50_wp*dum
   tma=1.50_wp*tma
   l=0
   do j=1,3
      l=l+j
      tma(l)=tma(l)-dum
   enddo
   dum=tmb(1)+tmb(3)+tmb(6)
   dum=0.50_wp*dum
   tmb=1.50_wp*tmb
   l=0
   do j=1,3
      l=l+j
      tmb(l)=tmb(l)-dum
   enddo
   write(iunit,'(a)',advance='yes')'molecular quadrupole (traceless):'
   write(iunit,'(a)',advance='no')'                '
   write(iunit,'(a)',advance='no')'xx          xy          yy          '
   write(iunit,'(a)',advance='yes')'xz          yz          zz'
   write(iunit,'(a,6f12.3)') ' q only: ',tma(1:6)
   write(iunit,'(a,6f12.3)') '  q+dip: ',tma(1:6)+tmb(1:6)
   write(iunit,'(a,6f12.3)') '   full: ',tma(1:6)+tmb(1:6)+tmc(1:6)

end subroutine molmom

! gradient evaluation from
! cumulative atomic multipole moment interactions: all interactions up to r**-3
! nat         : # of atoms
! xyz(3,nat)  : cartesian coordinates
! q(nat)      : atomic partial charges
! dipm(3,nat) : cumulative atomic dipole moments (x,y,z)
! qp(6,nat)   : traceless(!) cumulative atomic quadrupole moments (xx,xy,yy,xz,yz,zz)
! gab3,gab5   : damped R**-3 and R**-5 Coulomb laws, dimension: nat*(nat+1)/2
!               multiplication with numerator then leads to R**-2 and R**-3 decay, respectively
! radcn(nat)  : CN-depentent atomic radii
! dcn(3,i,j)  : derivative of radcn(j) w.r.t. cartesian directions of i
! exj         : exponent in gab, gab3, and gab5 - determines interpolation
! g           : nuclear gradient (3)
subroutine aniso_grad(nat,at,xyz,q,dipm,qp,kdmp3,kdmp5, &
      & radcn,draesdcn,dcn,gab3,gab5,g)
   use aoparam
   use lin_mod, only : lin
   !gab3 Hellmann-Feynman terms correct, shift terms to be tested yet
   implicit none
   integer, intent(in)   :: nat,at(nat)
   real(wp), intent(in)    :: xyz(3,nat),q(nat),dipm(3,nat),qp(6,nat)
   real(wp), intent(in)    :: gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
   real(wp), intent(in)    :: kdmp3,kdmp5,radcn(:),dcn(:,:,:)
   real(wp), intent(in)    :: draesdcn(:)
   real(wp), intent(inout) :: g(3,nat)
   real(wp) qp1(6),rr(3),dip(3),rij(3)
   real(wp) ed,eq,edd,e01,e02,e11,r2,tt,tt3,q1,dxi
   real(wp) tmp2,tmp3,rab,rabi,ddm2,ddm3a,ddm3b,qqa,qqb
   real(wp) dgab3,dgab5,damp1,damp2,ddamp,qs2

   integer iat,jat,k,l,m,ki,kj,kl,mm
   do iat=1,nat
      q1=q(iat)
      rr(1:3)=xyz(1:3,iat)
      dip(1:3)=dipm(1:3,iat)
      qp1(1:6)=qp(1:6,iat)
      tmp2=0.0_wp ! cumulate terms propto CN gradient -  to scale only quadratically
      do jat=1,nat  ! loop over other atoms
         if(iat.eq.jat) cycle
         kj=lin(jat,iat)
         rij(1:3)=xyz(1:3,jat)-rr(1:3)
         r2=sum(rij*rij)
         rabi=1.0_wp/sqrt(r2)
         call dzero(kdmp3,rabi,radcn(iat),radcn(jat),damp1,ddamp)
         dgab3=dgab(3.0_wp,rabi,damp1,ddamp)
         call dzero(kdmp5,rabi,radcn(iat),radcn(jat),damp2,ddamp)
         dgab5=dgab(5.0_wp,rabi,damp2,ddamp)
         ed=0.0_wp
         edd=0.0_wp
         eq=0.0_wp
         !           dipole - charge
         do k=1,3
            ed=ed+q(jat)*dip(k)*rij(k)
            ed=ed-dipm(k,jat)*q1*rij(k)
            tt=q1*dipm(k,jat)-q(jat)*dip(k)
            ! part of dip-q derivative
            g(k,iat)=g(k,iat)+gab3(kj)*tt
            !              dip-dip & charge-qpole
            ddm2=0.0_wp
            ddm3a=0.0_wp
            ddm3b=0.0_wp
            qqa=0.0_wp
            qqb=0.0_wp
            do l=1,3
               kl=lin(l,k)
               tt=rij(l)*rij(k)
               tt3=3.0_wp*tt
               eq=eq+q(jat)*qp1(kl)*tt
               eq=eq+qp(kl,jat)*q1*tt
               edd=edd-dipm(k,jat)*dip(l)*tt3
               ! extra d-d terms
               ddm2=ddm2+dipm(l,jat)*dip(l)
               ddm3a=ddm3a+dip(l)*rij(l)
               ddm3b=ddm3b+dipm(l,jat)*rij(l)
               ! extra q-qpole terms
               qqa=qqa+rij(l)*qp(kl,jat)
               qqb=qqb+rij(l)*qp1(kl)
            enddo
            edd=edd+dipm(k,jat)*dip(k)*r2
            g(k,iat)=g(k,iat)-2.0_wp*gab5(kj)*ddm2*rij(k)
            g(k,iat)=g(k,iat)+3.0_wp*gab5(kj)*ddm3a*dipm(k,jat)
            g(k,iat)=g(k,iat)+3.0_wp*gab5(kj)*ddm3b*dip(k)
            g(k,iat)=g(k,iat)-2.0_wp*gab5(kj)*qqa*q1
            g(k,iat)=g(k,iat)-2.0_wp*gab5(kj)*qqb*q(jat)
         enddo
         do k=1,3
            dxi=rij(k)*rabi
            g(k,iat)=g(k,iat)-ed*dgab3*dxi
            g(k,iat)=g(k,iat)-(eq+edd)*dgab5*dxi
         enddo
         ! collect terms for CN-dependent part
         rab=0.50_wp*(radcn(iat)+radcn(jat))
         tmp3=ed*kdmp3*gab3(kj)*(damp1/rab)*(rab*rabi)**kdmp3
         tmp2=tmp2+tmp3
         tmp3=(eq+edd)*kdmp5*gab5(kj)*(damp2/rab)*(rab*rabi)**kdmp5
         tmp2=tmp2+tmp3
      enddo
      ! CN-dependent part  - O(N^2)
      tmp2=3.0_wp*tmp2*draesdcn(iat)
      g = g + tmp2*dcn(:, :, iat)

   enddo
end subroutine aniso_grad

! check and print sparsity w.r.t. individual contribution
! to get an idea
subroutine checkspars(nao,ndp,nqp,nmat,matlist,mqlst,mdlst)
   use lin_mod, only : lin
   implicit none
   integer,intent(in) :: ndp,nqp,nmat,nao
   integer,intent(in) :: matlist(2,nao*(nao+1)/2)
   integer,intent(in) :: mqlst(2,nqp),mdlst(2,ndp)
   integer  :: i,j,m,k,ntot,mi,mj,ki,kj,mm,kk
   logical,allocatable ::  nzero(:)
   ! check overall sparsity
   allocate(nzero(nao*(nao+1)/2))
   nzero=.false.
   do k=1,ndp
      ki=mdlst(1,k)
      kj=mdlst(2,k)
      kk=lin(ki,kj)
      nzero(kk)=.true.
   enddo
   do k=1,nqp
      ki=mqlst(1,k)
      kj=mqlst(2,k)
      kk=lin(ki,kj)
      nzero(kk)=.true.
   enddo
   do k=1,nmat
      ki=matlist(1,k)
      kj=matlist(2,k)
      kk=lin(ki,kj)
      nzero(kk)=.true.
   enddo
   mm=nao*(nao+1)/2
   ntot=0
   do i=1,mm
      if(nzero(i)) ntot=ntot+1
   enddo
   write(*,'(a)',advance='yes') ' '
   write(*,'(a)')'% of non-zero elements in H:'
   write(*,'(''           by overlap:'',f6.2)') &
      & 100.*float(nmat)/float(mm)
   write(*,'(''      by dipole ints.:'',f6.2)') &
      & 100.*float(ndp)/float(mm)
   write(*,'(''  by quadrupole ints.:'',f6.2)') &
      & 100.*float(nqp)/float(mm)
   write(*,'(''                total:'',f6.2)') &
      & 100.*float(ntot)/float(mm)
   write(*,'(a)',advance='yes') ' '
   deallocate(nzero)
end subroutine checkspars

! zero-damped gab
subroutine mmomgabzero(nat,at,xyz,kdmp3,kdmp5,radcn,gab3,gab5)
   implicit none
   integer, intent(in) :: nat,at(nat)
   real(wp), intent(in)  ::  xyz(3,nat),radcn(nat)
   real(wp), intent(in)  ::  kdmp3,kdmp5
   real(wp), intent(out) :: gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
   real(wp) damp,ddamp

   real(wp) tmp1,tmp2,rr(3)
   integer i,j,k,l,lin

   !!!!!!! set up damped Coulomb operators for multipole interactions
   gab3=0.0_wp ! for r**-2 decaying q-dip term
   gab5=0.0_wp ! for r**-3 decaying terms (q-qpol,dip-dip)
   l=1
   gab3(1)=0.0_wp
   gab5(1)=0.0_wp
   do i=2,nat
      rr(1:3)=xyz(1:3,i)
      do j=1,i-1
         l=l+1
         tmp2=0.0_wp
         do k=1,3
            tmp1=xyz(k,j)-rr(k)
            tmp2=tmp2+tmp1**2
         enddo
         tmp1=1.0_wp/sqrt(tmp2)
         !           call dzero(2.0_wp,tmp1,at(i),at(j),damp,ddamp)
         call dzero(kdmp3,tmp1,radcn(i),radcn(j),damp,ddamp)
         gab3(l)=gab(3.0_wp,tmp1,damp)
         !           call dzero(3.0_wp,tmp1,at(i),at(j),damp,ddamp)
         call dzero(kdmp5,tmp1,radcn(i),radcn(j),damp,ddamp)
         gab5(l)=gab(5.0_wp,tmp1,damp)
      enddo
      l=l+1
      gab3(l)=0.0_wp
      gab5(l)=0.0_wp
   enddo

   !!!!!!! DEBUG
   !      gab3=0.0_wp
   !      gab5=0.0_wp
   !      gab3=100.0_wp*gab3
   !      gab5=100.0_wp*gab5
end subroutine mmomgabzero


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     setup CN-dependet atomic radii
!     n : # of atoms
! at(n) : atomic number array
! cn(n) : coordination number of atoms
! shift : global offset from cnval     (parameter)
! expo  : exponent scaling/ steepness  (parameter)
! rmax  : maximum radius               (parameter)
! radcn : CN-dependent radius
subroutine get_radcn(n,at,cn,shift,expo,rmax,radcn)
   use aoparam
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp), intent(in) :: cn(n),shift,expo,rmax
   real(wp), intent(out) :: radcn(n)
   real(wp) :: rco,t1,t2
   integer :: iat
   do iat=1,n
      rco=radaes(at(iat))             ! base radius of element
      t1 =cn(iat)-cnval(at(iat))-shift  ! CN - VALCN - SHIFT
      t2 =rco +(rmax-rco)/(1.0_wp+exp(-expo*t1))
      radcn(iat)=t2
   enddo
end subroutine get_radcn
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     derivative of CN-dependet atomic radii w.r.t. other atoms
!     n : # of atoms
! at(n) : atomic number array
! cn(n) : coordination number of atoms
! shift : global offset from cnval     (parameter)
! expo  : exponent scaling/ steepness  (parameter)
! rmax  : maximum radius               (parameter)
! dcn   : on input  : derivatives of CN(i) w.r.t. Cart. directions of j
!       : on output : derivatives of RADCN(j) w.r.t. Cart. directions of i, so we flip indices!!!
subroutine dradcn(n,at,cn,shift,expo,rmax,draesdcn)
   use aoparam
   implicit none
   integer, intent (in) :: n,at(:)
   real(wp), intent (in)  :: cn(:),shift,expo,rmax
   real(wp), intent (inout) :: draesdcn(:)
   real(wp) rco,t1,t2,t3,t4,tmp1(3),tmp2(3)
   integer iat
   do iat=1,n
      rco=radaes(at(iat))             ! base radius of element
      t1 =exp(-expo*(cn(iat)-cnval(at(iat))-shift))  ! CN - VALCN - SHIFT
      draesdcn(iat) = (rmax-rco)/(1.0_wp+2.0_wp*t1+t1*t1)*expo*t1
   enddo
end subroutine dradcn

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     zero-damping function and derivative
!     rscal         : scaling of radii
!     dex           : exponent in zero-damping function
!     rabinv        : inverse distance, i.e., 1/Rab
!     aradi,radj    : CN-dependent atomic number radii
!     damp          : zero damping function
!     ddamp         : derivative w.r.t. Rab
subroutine dzero(dex,rabinv,radi,radj,damp,ddamp)
   implicit none
   real(wp), intent (in) :: radi,radj
   real(wp), intent (in)  :: rabinv,dex
   real(wp), intent (out) :: damp,ddamp
   real(wp) rco,f1,f2
   !     rco=2.0/(1./radaes(ati)+1./radaes(atj))  ! unstable
   !     rco=sqrt(radaes(ati)*radaes(atj))        ! unstable
   rco=0.5*(radi+radj)
   ! zero-damping function and gradient w.r.t. Rab
   damp=1.0_wp/(1.0_wp+6.0_wp*(rco*rabinv)**dex)
   ddamp=-dex*rabinv*(damp*damp-damp)
end subroutine dzero

!!!   gab - computes the damped Coulomb type interaction with dex decay:
!     gab = damp * Rab**(-dex)
!     dex           : exponent defining the decay: Rab**(-dex)
!     rabinv        : inverse distance, i.e., 1/Rab
!     damp          : zero damping function
real(wp) function gab(dex,rabinv,damp)
   implicit none
   real(wp) dex,rabinv,damp
   ! compute r**dex decaying intermediate
   gab=damp*(rabinv**dex) ! LR-decay * damping
end function gab


!!!   dgab - computes the derivative w.r.t. Rab of a damped Coulomb type interaction with dex decay:
!     gab = damp * Rab**(-dex)
!     dex           : exponent defining the decay: Rab**(-dex)
!     rabinv        : inverse distance, i.e., 1/Rab
!     damp          : zero damping function
!     ddamp         : derivative of damping function
real(wp) function dgab(dex,rabinv,damp,ddamp)
   implicit none
   real(wp) dex,rabinv,damp,ddamp,tmp1
   ! compute r**dex decaying intermediate
   tmp1=-dex*rabinv*(rabinv**dex) ! LR-decay derivative
   dgab=tmp1*damp+ddamp*rabinv**dex
end function dgab

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

! this prints gnuplot scripts
subroutine printspars(nat,at,nao,nmat,ndp,nqp,xyz,s,dpint,qpint)
   use lin_mod, only : lin
   implicit none
   integer, intent(in)    :: nao,nat,at(nat)
   integer, intent(in) :: ndp,nqp,nmat
   real(wp), intent(in)  ::  s(nao,nao),xyz(3,nat)
   real(wp), intent(in)  :: dpint(3,nao,nao)
   real(wp), intent(in)  :: qpint(6,nao,nao)

   real(wp) tmp1,tmp2,tmp3,tmp4,imat
   real(wp) dd1(3),dd2(3),qq1(6),qq2(6)
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f
   ! stuff for potential
   real(wp)  rr(3),f1,f2,rco
   real(wp), allocatable :: x(:)

   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,mq,md,mi,mj,ij
   character(80) frmt

   allocate(x(nao))
   md=0
   mq=0
   x=0.0_wp
   kk=0
   imat=0

   !!! print gnuplot script for matrix sparsity
   frmt=''
   write(frmt,'(a,i0,a)') '(,',nao,'(f0.8,x))'
   mm=nao*(nao+1)/2
   ! dipole pot
   tmp4=100.0_wp*dble(ndp)/dble(mm)
   open(unit=46,file='dipplot',status='replace')
   write(46,'(a)',advance='no') 'set title "{/Symbol m}-int'
   write(46,'(a)',advance='yes') 'matrix potential" font ",18"'
   !      write(46,'(a,i0,a)') 'roportional potential (# of SAOs: ',nao,')"'
   write(46,'(a)') 'unset key'
   write(46,'(a)') 'set tic scale 0'
   write(46,'(a)') ' '
   write(46,'(a)') 'set size ratio 1'
   write(46,'(a)') 'set palette model RGB'
   write(46,'(a)') 'set palette defined (0 "white", 1 "blue")'
   write(46,'(a)') 'set format cb "10^{%L}"'
   write(46,'(a)') 'set cbrange [0.00000001:1]'
   write(46,'(a)') 'set logscale cb'
   write(46,'(a)') 'set cbtics offset 0.5'
   write(46,'(a)',advance='no') 'set cblabel "{/Symbol \263}" offset'
   write(46,'(a)',advance='yes') '-3.6,17.1 rotate by 0'
   write(46,'(a)') ' '
   write(46,'(a,i0,a)') 'set xrange [-0.5:',nao-1,'.5]'
   write(46,'(a,i0,a)') 'set yrange [-0.5:',nao-1,'.5]'
   write(46,'(a)') 'unset xtics'
   write(46,'(a)') 'unset ytics'
   write(46,'(a)',advance='no') 'set xlabel offset 0,2 "# of SAOs: '
   write(46,'(i0,a,f5.1,a)') nao,', non-zero: ',tmp4,'%" font ",16"'
   write(46,'(a)') 'unset ylabel'
   write(46,'(a)') 'set terminal postscript enhanced eps size 5,5'
   write(46,'(a)') 'set output "dipplot.eps"'
   !      write(46,'(a)') 'set terminal pdf size 5,5'
   !      write(46,'(a)') 'set output "dipplot.pdf"'
   write(46,'(a)') ' '
   write(46,'(a)') 'set view map'
   write(46,'(a)') "splot '-' matrix with image"
   do i=1,nao
      do j=1,nao
         k=lin(j,nao-i+1)
         x(j)=dpint(1,j,nao-i+1)*dpint(1,j,nao-i+1)+dpint(2,j,nao-i+1)*dpint(2,j,nao-i+1) &
            & +dpint(3,j,nao-i+1)*dpint(3,j,nao-i+1)
         x(j)=sqrt(x(j))
      enddo
      write(46,frmt) x(1:nao)
   enddo
   write(46,'(a)') 'e'
   write(46,'(a)') 'e'
   close(46)
   ! quadrupole pot
   tmp4=100.0_wp*dble(nqp)/dble(mm)
   open(unit=46,file='qplot',status='replace')
   write(46,'(a)',advance='no') 'set title "{/Symbol q}-int'
   write(46,'(a,i0,a)') ' matrix potential" font ",18"'
   write(46,'(a)') 'unset key'
   write(46,'(a)') 'set tic scale 0'
   write(46,'(a)') ' '
   write(46,'(a)') 'set size ratio 1'
   write(46,'(a)') 'set palette model RGB'
   write(46,'(a)') 'set palette defined (0 "white", 1 "blue")'
   write(46,'(a)') 'set format cb "10^{%L}"'
   write(46,'(a)') 'set cbrange [0.00000001:1]'
   write(46,'(a)') 'set logscale cb'
   write(46,'(a)') 'set cbtics offset 0.5'
   write(46,'(a)',advance='no') 'set cblabel "{/Symbol \263}" offset'
   write(46,'(a)',advance='yes') '-3.6,17.1 rotate by 0'
   write(46,'(a)') ' '
   write(46,'(a,i0,a)') 'set xrange [-0.5:',nao-1,'.5]'
   write(46,'(a,i0,a)') 'set yrange [-0.5:',nao-1,'.5]'
   write(46,'(a)') 'unset xtics'
   write(46,'(a)') 'unset ytics'
   write(46,'(a)',advance='no') 'set xlabel offset 0,2 "# of SAOs: '
   write(46,'(i0,a,f5.1,a)') nao,', non-zero: ',tmp4,'%" font ",16"'
   write(46,'(a)') 'unset ylabel'
   !      write(46,'(a)') 'set terminal pdf size 5,5'
   !      write(46,'(a)') 'set output "qplot.pdf"'
   write(46,'(a)') 'set terminal postscript enhanced eps size 5,5'
   write(46,'(a)') 'set output "qplot.eps"'
   write(46,'(a)') ' '
   write(46,'(a)') 'set view map'
   write(46,'(a)') "splot '-' matrix with image"
   do i=1,nao
      do j=1,nao
         k=lin(j,nao-i+1)
         x(j)=qpint(1,j,nao-i+1)*qpint(1,j,nao-i+1)+qpint(2,j,nao-i+1)*qpint(2,j,nao-i+1) &
            & +qpint(3,j,nao-i+1)*qpint(3,j,nao-i+1)+2.0_wp*qpint(4,j,nao-i+1)*qpint(4,j,nao-i+1) &
            & +2.0_wp*qpint(5,j,nao-i+1)*qpint(5,j,nao-i+1)+2.0_wp*qpint(6,j,nao-i+1)*qpint(6,j,nao-i+1)
         x(j)=sqrt(x(j))
      enddo
      write(46,frmt) x(1:nao)
   enddo
   write(46,'(a)') 'e'
   write(46,'(a)') 'e'
   close(46)
   ! total
   tmp4=100.0_wp*dble(nmat)/dble(mm)
   open(unit=46,file='splot',status='replace')
   write(46,'(a)',advance='no') 'set title "overlap matrix"'
   write(46,'(a,i0,a)') ' font ",18"'
   write(46,'(a)') 'unset key'
   write(46,'(a)') 'set tic scale 0'
   write(46,'(a)') ' '
   write(46,'(a)') 'set size ratio 1'
   write(46,'(a)') 'set palette model RGB'
   write(46,'(a)') 'set palette defined (0 "white", 1 "blue")'
   write(46,'(a)') 'set format cb "10^{%L}"'
   write(46,'(a)') 'set cbrange [0.00000001:1]'
   write(46,'(a)') 'set logscale cb'
   write(46,'(a)') 'set cbtics offset 0.5'
   write(46,'(a)',advance='no') 'set cblabel "{/Symbol \263}" offset'
   write(46,'(a)',advance='yes') '-3.6,17.1 rotate by 0'
   write(46,'(a)') ' '
   write(46,'(a,i0,a)') 'set xrange [-0.5:',nao-1,'.5]'
   write(46,'(a,i0,a)') 'set yrange [-0.5:',nao-1,'.5]'
   write(46,'(a)') 'unset xtics'
   write(46,'(a)') 'unset ytics'
   write(46,'(a)',advance='no') 'set xlabel offset 0,2 "# of SAOs: '
   write(46,'(i0,a,f5.1,a)') nao,', non-zero: ',tmp4,'%" font ",16"'
   write(46,'(a)') 'unset ylabel'
   !      write(46,'(a)') 'set terminal pdf size 5,5'
   !      write(46,'(a)') 'set output "splot.pdf"'
   write(46,'(a)') 'set terminal postscript enhanced eps size 5,5'
   write(46,'(a)') 'set output "splot.eps"'
   write(46,'(a)') ' '
   write(46,'(a)') 'set view map'
   write(46,'(a)') "splot '-' matrix with image"
   do i=1,nao
      write(46,frmt) abs(s(1:nao,nao-i+1))
   enddo
   write(46,'(a)') 'e'
   write(46,'(a)') 'e'
   close(46)


   deallocate(x)
end subroutine printspars

subroutine gfn2broyden_diff(n,istart,nbr,dipm,qp,q_in,dq)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: dipm(3,n),qp(6,n),q_in(nbr)
   real(wp), intent (inout) :: dq(nbr)
   integer i,j,k
   k=istart
   do i=1,n
      do j=1,3
         k=k+1
         dq(k)=dipm(j,i)-q_in(k)
      enddo
      do j=1,6
         k=k+1
         dq(k)=qp(j,i)-q_in(k)
      enddo
   enddo
   istart=k

end subroutine gfn2broyden_diff

subroutine gfn2broyden_save(n,istart,nbr,dipm,qp,q_in)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: dipm(3,n),qp(6,n)
   real(wp), intent (inout) :: q_in(nbr)
   integer i,j,k
   k=istart
   do i=1,n
      do j=1,3
         k=k+1
         q_in(k)=dipm(j,i)
      enddo
      do j=1,6
         k=k+1
         q_in(k)=qp(j,i)
      enddo
   enddo
   istart=k

end subroutine gfn2broyden_save



subroutine gfn2broyden_out(n,istart,nbr,q_in,dipm,qp)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: q_in(nbr)
   real(wp), intent (out) :: dipm(3,n),qp(6,n)
   integer i,j,k
   k=istart
   do i=1,n
      do j=1,3
         k=k+1
         dipm(j,i)=q_in(k)
      enddo
      do j=1,6
         k=k+1
         qp(j,i)=q_in(k)
      enddo
   enddo
   istart=k
end subroutine gfn2broyden_out

subroutine get_sdqint(mol, neighs, neighlist, nbf,nao,thr,intcut, &
      &               caoshell, saoshell, nprim, primcount, alp, cont, &
      &               sint, dpint, qpint, ndp, nqp)
   use tbdef_molecule
   use tbdef_neighbourlist
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   use lin_mod, only : lin
   type(tb_molecule), intent(in) :: mol
   type(tb_neighbourlist), intent(in) :: neighlist
   integer, intent(in) :: neighs(:)
   integer, intent(in) :: nao
   integer, intent(in) :: nbf
   real(wp), intent(in) :: thr
   real(wp), intent(in) :: intcut
   integer, intent(in) :: caoshell(:,:)
   integer, intent(in) :: saoshell(:,:)
   integer, intent(in) :: nprim(:)
   integer, intent(in) :: primcount(:)
   real(wp), intent(in) :: alp(:)
   real(wp), intent(in) :: cont(:)
   real(wp), intent(out) :: sint(nao,nao)
   real(wp), intent(out) :: dpint(3,nao,nao)
   real(wp), intent(out) :: qpint(6,nao,nao)
   integer, intent(out), optional :: ndp
   integer, intent(out), optional :: nqp

   integer :: i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij,ijao,img
   real(wp) :: tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) :: dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) :: skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,ab,est
   real(wp) :: ri(3),rj(3)
   real(wp) :: dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer :: ip,jp,iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax
   integer :: ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer, parameter :: itt(0:3) = (/0,1,4,10/)
   real(wp) :: saw(10)
   real(wp), parameter :: point(3) = 0.0_wp


   !call timing(t1,t3)
   thr2=(thr*1.0d-2)-thr*1.0d-12
   sint=0.0_wp
   dpint=0.0_wp
   qpint=0.0_wp

   kj=0
   !$omp parallel do default(none) &
   !$omp reduction(+:sint, dpint, qpint) shared(mol, neighs, neighlist, intcut) &
   !$omp shared(ao_n, ao_l, caoshell, saoshell, nprim, primcount, cont, alp) &
   !$omp private(ri, rj, r2, jat, ish, jsh, ati, atj, icao, jcao, ij, img, &
   !$omp&        naoi, naoj, jshmax, ishtyp, jshtyp, iptyp, jptyp, iao, jao, &
   !$omp&        ss, dd, qq, tmp)
   do iat=1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      do ij = 0, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat == jat) jshmax=ish
            do jsh=1,jshmax
               jshtyp=ao_l(jsh,atj)
               jcao=caoshell(jsh,jat)
               naoj=llao(jshtyp)
               jptyp=itt(jshtyp)
               ss=0.0_wp
               dd=0.0_wp
               qq=0.0_wp
               call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point, &
                  &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
               ! transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)
               do k=1,3
                  tmp(1:6,1:6)=dd(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  dd(k,1:6,1:6)=tmp(1:6,1:6)
               enddo
               do k=1,6
                  tmp(1:6,1:6)=qq(k,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  qq(k,1:6,1:6)=tmp(1:6,1:6)
               enddo
               do ii=1,llao2(ishtyp)
                  iao=ii+saoshell(ish,iat)
                  do jj=1,llao2(jshtyp)
                     jao=jj+saoshell(jsh,jat)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(:, jao, iao) = dpint(:, jao, iao) + dd(:, jj, ii)
                     qpint(:, jao, iao) = qpint(:, jao, iao) + qq(:, jj, ii)
                     if(jao >= iao .and. iat == jat) cycle
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     dpint(:, iao, jao) = dpint(:, iao, jao) + dd(:, jj, ii)
                     qpint(:, iao, jao) = qpint(:, iao, jao) + qq(:, jj, ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do

   if (present(ndp)) then
      ndp = 0
      do i=1,nao
         do j=1,i
            ij=ij+1
            tmp1=0.0_wp
            kk=0
            do k=1,3
               tmp1=tmp1+dpint(k,j,i)*dpint(k,j,i)
            enddo
            if(tmp1.gt.thr2) ndp=ndp+1
         enddo
      enddo
   endif
   if (present(nqp)) then
      nqp = 0
      do i=1,nao
         do j=1,i
            ij=ij+1
            tmp2=0.0_wp
            kk=0
            do k=1,3
               tmp2=tmp2-qpint(k,j,i)*qpint(k,j,i)
            enddo
            do k=1,6
               tmp2=tmp2+2.0_wp*qpint(k,j,i)*qpint(k,j,i)
            enddo
            if(tmp2.gt.thr2) nqp=nqp+1
         enddo
      enddo
   endif

end subroutine get_sdqint

!> computes the dipole and quadrupole integrals and performs screening to
!  determine, which contribute to potential
subroutine sdqint(nat,at,nbf,nao,xyz,thr,ndp,nqp,intcut,caoshell,saoshell, &
      &           nprim,primcount,alp,cont,sint,dpint,qpint)
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   use lin_mod, only : lin
   integer, intent(in)  :: nat
   integer, intent(in)  :: nao
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: xyz(3,nat)
   real(wp),intent(in)  :: thr
   real(wp),intent(in)  :: intcut
   integer, intent(in)  :: caoshell(:,:)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   real(wp),intent(out) :: sint(nao,nao)
   real(wp),intent(out) :: dpint(3,nao,nao)
   real(wp),intent(out) :: qpint(6,nao,nao)
   integer, intent(out) :: ndp,nqp


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,ab,est
   real(wp), parameter :: bohr=1.0_wp/0.52917726_wp

   real(wp) ri(3),rj(3),f1,f2
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer iat,jat,ati,atj,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp), parameter :: point(3) = 0.0_wp
   real(wp) :: saw(10)


   !call timing(t1,t3)
   thr2=(thr*1.0d-2)-thr*1.0d-12
   !     thr2=1.0d-20 ! conservative, keep all terms
   !     integrals
   sint=0.0_wp
   dpint=0.0_wp
   qpint=0.0_wp

   !$omp parallel do default(none) &
   !$omp reduction(+:sint, dpint, qpint) shared(nat, at, xyz, intcut) &
   !$omp shared(ao_n, ao_l, caoshell, saoshell, nprim, primcount, cont, alp) &
   !$omp private(ri, rj, r2, jat, ish, jsh, ati, atj, icao, jcao, &
   !$omp&        naoi, naoj, jshmax, ishtyp, jshtyp, iptyp, jptyp, iao, jao, &
   !$omp&        ss, dd, qq, tmp)
   do iat = 1, nat
      ri = xyz(:,iat)
      ati = at(iat)
      do jat = 1, iat
         rj = xyz(:,jat)
         atj = at(jat)
         r2 = sum((rj-ri)**2)
         ! ints < 1.d-9 for RAB > 40 Bohr
         if (r2 > 2000.0_wp) cycle
         do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat.eq.jat) jshmax = ish
            do jsh=1,jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               ss = 0.0_wp
               dd = 0.0_wp
               qq = 0.0_wp
               call get_multiints(icao, jcao, naoi, naoj, iptyp, jptyp, &
                  &               ri, rj, point, intcut, nprim, primcount, &
                  &               alp, cont, ss, dd, qq)
               !transform from CAO to SAO
               call dtrf2(ss, ishtyp, jshtyp)
               do k = 1, 3
                  tmp(1:6, 1:6) = dd(k, 1:6, 1:6)
                  call dtrf2(tmp, ishtyp, jshtyp)
                  dd(k, 1:6, 1:6)=tmp(1:6, 1:6)
               enddo
               do k = 1, 6
                  tmp(1:6, 1:6)=qq(k, 1:6, 1:6)
                  call dtrf2(tmp, ishtyp, jshtyp)
                  qq(k, 1:6, 1:6)=tmp(1:6, 1:6)
               enddo
               do ii = 1, llao2(ishtyp)
                  iao = ii + saoshell(ish, iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj + saoshell(jsh, jat)
                     sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                     dpint(:, jao, iao) = dpint(:, jao, iao) + dd(:, jj, ii)
                     qpint(:, jao, iao) = qpint(:, jao, iao) + qq(:, jj, ii)
                     if(jao >= iao .and. iat == jat) cycle
                     sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                     dpint(:, iao, jao) = dpint(:, iao, jao) + dd(:, jj, ii)
                     qpint(:, iao, jao) = qpint(:, iao, jao) + qq(:, jj, ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$omp end parallel do

   !if (present(ndp)) then
      ndp = 0
      do i=1,nao
         do j=1,i
            ij=ij+1
            tmp1=0.0_wp
            kk=0
            do k=1,3
               tmp1=tmp1+dpint(k,j,i)*dpint(k,j,i)
            enddo
            if(tmp1.gt.thr2) ndp=ndp+1
         enddo
      enddo
   !endif
   !if (present(nqp)) then
      nqp = 0
      do i=1,nao
         do j=1,i
            ij=ij+1
            tmp2=0.0_wp
            kk=0
            do k=1,3
               tmp2=tmp2-qpint(k,j,i)*qpint(k,j,i)
            enddo
            do k=1,6
               tmp2=tmp2+2.0_wp*qpint(k,j,i)*qpint(k,j,i)
            enddo
            if(tmp2.gt.thr2) nqp=nqp+1
         enddo
      enddo
   !endif

end subroutine sdqint

!> Computes the gradient of the overlap integral contributions to the Hamiltonian.
subroutine dsint(mol, neighs, neighlist, thr, caoshell, saoshell, &
      &          nprim, primcount, alp, cont, H, gradient, sigma)
   use tbdef_molecule
   use tbdef_neighbourlist
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   !> Integral cutoff according to prefactor from Gaussian product theorem.
   real(wp), intent(in) :: thr
   !> Map shell of atom to index in CAO space (lowest Cart. component is taken)
   integer, intent(in) :: caoshell(:,:)
   !> Map shell of atom to index in SAO space (lowest m_l component is taken)
   integer, intent(in) :: saoshell(:,:)
   !> Number of primitive basis functions for given CAO
   integer, intent(in) :: nprim(:)
   !> Index of first primitive (over entire system) of given CAO
   integer, intent(in) :: primcount(:)
   real(wp), intent(in) :: alp(:)
   real(wp), intent(in) :: cont(:)
   real(wp), intent(in) :: H(:, :)
   real(wp), intent(inout) :: gradient(:,:)
   real(wp), intent(inout) :: sigma(:,:)

   integer, parameter :: itt(0:3)=(/0,1,4,10/)
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, ixyz
   integer :: ii, jj, ij, jshmax, img, ip, jp, mli, mlj, iprim, jprim
   integer :: ishtyp, jshtyp, iptyp, jptyp, naoi, naoj
   real(wp) :: alpj, ci, cc, alpi, ab, est, orho, fact
   real(wp) :: r2, rij(3), ri(3), rj(3)
   real(wp) :: tmp(6,6), dumdum(3), dum, sdq(6,6), sdqg(3,6,6)
   real(wp) :: dG(3), dS(3,3)

   !$omp parallel do default(none) reduction(+:gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, H, caoshell, saoshell, nprim, primcount, &
   !$omp&       alp, cont, ao_n, ao_l, thr) &
   !$omp private(jat, ij, img, ati, atj, ishtyp, jshtyp, icao, jcao, naoi, naoj, &
   !$omp&        jshmax, iptyp, jptyp, iao, jao, ii, jj, ip, jp, iprim, jprim, &
   !$omp&        mli, mlj, ixyz, ri, rj, r2, rij, alpi, alpj, est, ci, cc, &
   !$omp&        dum, dumdum, sdq, sdqg, tmp, orho, fact, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         rij = ri - rj
         do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish,ati)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat.eq.jat) jshmax=ish
            do jsh = 1, jshmax
               jshtyp = ao_l(jsh,atj)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               ! we go through the primitives
               ! (because the screening is the same for all of them)
               sdqg = 0.0_wp
               sdq = 0.0_wp
               call get_grad_overlap(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                  &                  rj,thr,nprim,primcount, &
                  &                  alp,cont,sdq,sdqg)
               do ixyz = 1, 3
                  ! transform from CAO to SAO, only transform the gradient
                  tmp(1:6, 1:6)=sdqg(ixyz, 1:6, 1:6)
                  call dtrf2(tmp, ishtyp, jshtyp)
                  sdqg(ixyz, 1:6, 1:6)=tmp(1:6, 1:6)
               enddo
               do ii = 1, llao2(ishtyp)
                  iao = ii + saoshell(ish,iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj + saoshell(jsh,jat)
                     dG = sdqg(:, jj, ii) * 2*H(jao, iao)
                     dS = spread(dG, 1, 3) * spread(rij, 2, 3)
                     if (iat /= jat) then
                        sigma = sigma + dS
                        gradient(:, iat) = gradient(:, iat) + dG
                        gradient(:, jat) = gradient(:, jat) - dG
                     else
                        sigma = sigma + 0.5_wp * dS
                     endif
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo ! ish : loop over shells on iat
      enddo ! jat
   enddo ! iat
   !$omp end parallel do
end subroutine dsint

! ddqint: computes the gradient of the dipole/qpole integral contribution
! it starts with a numerical derivation of the multipole integrals w.r.t. nuclear displacements
! and collecting their contribution (w/ density matrix as input)
! thr         : integral cutoff according to prefactor from Gaussian product theorem
! nat         : # of atoms
! nao         : # of spherical AOs (SAOs)
! nbf         : # of Cartesian AOs (CAOs)
! at(nat)     : atomic numbers of atoms
! xyz(3,nat)  : cartesian coordinates
! caoshell    : map shell of atom to index in CAO space (lowest Cart. component is taken), dimension: (5,nat)
! saoshell     : map shell of atom to index in SAO space (lowest m_l component is taken), dimension: (5,nat)
! primcount   : index of first primitive (over entire system) of given CAO, dimension: nbf
! p(nao,nao)  : density matrix
subroutine ddqint(mol, neighs, neighlist, thr, caoshell, saoshell, &
      &           nprim, primcount, alp, cont, p, vs, vd, vq, H, &
      &           gradient, sigma)
   use tbdef_molecule
   use tbdef_neighbourlist
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   !> Molecular structure information.
   type(tb_molecule), intent(in) :: mol
   !> Static neighbourlist.
   type(tb_neighbourlist), intent(in) :: neighlist
   !> Number of neighbours for each atom.
   integer, intent(in) :: neighs(:)
   real(wp), intent(in) :: thr
   real(wp), intent(in) :: vd(:,:)
   real(wp), intent(in) :: vq(:,:)
   real(wp), intent(in) :: H(:,:)
   real(wp), intent(in) :: vs(:)
   integer, intent(in) :: caoshell(:,:)
   integer, intent(in) :: saoshell(:,:)
   integer, intent(in) :: nprim(:)
   integer, intent(in) :: primcount(:)
   real(wp), intent(in) :: alp(:)
   real(wp), intent(in) :: cont(:)
   real(wp), intent(in) :: p(:,:)
   real(wp), intent(inout) :: gradient(:,:)
   real(wp), intent(inout) :: sigma(:,:)

   integer, parameter :: itt(0:3)=(/0,1,4,10/)
   integer :: iat, jat, ati, atj, ish, jsh, icao, jcao, iao, jao, ixyz, k
   integer :: ii, jj, ij, jshmax, img, ip, jp, mli, mlj, iprim, jprim
   integer :: ishtyp, jshtyp, iptyp, jptyp, naoi, naoj
   real(wp) :: alpj, ci, cc, alpi, ab, est, orho, fact
   real(wp) :: r2, rij(3), ri(3), rj(3)
   real(wp) :: tmp(6,6), dumdum(3), dum, sdq(10,6,6), sdqg(3,19,6,6)
   real(wp) :: stmp(3), dtmp(3), qtmp(3)
   real(wp) :: dG(3), dS(3,3)

   !$omp parallel do default(none) reduction(+:gradient, sigma) &
   !$omp shared(mol, neighs, neighlist, H, P, vs, vd, vq, caoshell, saoshell, &
   !$omp&       nprim, primcount, alp, cont, ao_n, ao_l, thr) &
   !$omp private(jat, ij, img, ati, atj, ishtyp, jshtyp, icao, jcao, naoi, naoj, &
   !$omp&        jshmax, iptyp, jptyp, iao, jao, ii, jj, k, ixyz, &
   !$omp&        ri, rj, r2, rij, alpi, alpj, est, ci, cc, sdq, sdqg, tmp, &
   !$omp&        stmp, dtmp, qtmp, dG, dS)
   do iat = 1, len(mol)
      ati = mol%at(iat)
      ri = mol%xyz(:, iat)
      do ij = 1, neighs(iat)
         img = neighlist%ineigh(ij, iat)
         r2 = neighlist%dists2(ij, iat)
         rj = neighlist%coords(:, img)
         jat = neighlist%image(img)
         atj = mol%at(jat)
         rij = ri - rj
         do ish = 1, ao_n(ati)
            ishtyp = ao_l(ish, ati)
            icao = caoshell(ish, iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = ao_n(atj)
            if(iat.eq.jat) jshmax = ish
            do jsh = 1, jshmax ! jshells
               jshtyp = ao_l(jsh, atj)
               jcao = caoshell(jsh, jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               sdqg = 0.0_wp
               sdq = 0.0_wp
               call get_grad_multiint(icao, jcao, naoi, naoj, iptyp, jptyp, &
                  &                   ri, rj, thr, nprim, primcount, alp, cont, &
                  &                   sdq, sdqg)
               do k = 1, 19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                  do ixyz = 1, 3
                     tmp(1:6, 1:6) = sdqg(ixyz, k, 1:6, 1:6)
                     call dtrf2(tmp, ishtyp, jshtyp)
                     sdqg(ixyz, k, 1:6, 1:6) = tmp(1:6, 1:6)
                  enddo
               enddo
               do ii = 1, llao2(ishtyp)
                  iao = ii + saoshell(ish, iat)
                  do jj = 1, llao2(jshtyp)
                     jao = jj + saoshell(jsh, jat)
                     cc = p(jao, iao)
                     stmp = sdqg(:, 1, jj, ii)*(2*H(jao, iao) &
                        & +cc*(vs(iat)+vs(jat)))
                     dtmp=cc*(sdqg(:, 11, jj, ii)*vd(1, iat) &
                        &    +sdqg(:, 12, jj, ii)*vd(2, iat) &
                        &    +sdqg(:, 13, jj, ii)*vd(3, iat) &
                        &    +sdqg(:, 2, jj, ii)*vd(1, jat) &
                        &    +sdqg(:, 3, jj, ii)*vd(2, jat) &
                        &    +sdqg(:, 4, jj, ii)*vd(3, jat))
                     qtmp=cc*(sdqg(:, 14, jj, ii)*vq(1, iat) &
                        &    +sdqg(:, 15, jj, ii)*vq(2, iat) &
                        &    +sdqg(:, 16, jj, ii)*vq(3, iat) &
                        &    +sdqg(:, 17, jj, ii)*vq(4, iat) &
                        &    +sdqg(:, 18, jj, ii)*vq(5, iat) &
                        &    +sdqg(:, 19, jj, ii)*vq(6, iat) &
                        &    +sdqg(:, 5, jj, ii)*vq(1, jat) &
                        &    +sdqg(:, 6, jj, ii)*vq(2, jat) &
                        &    +sdqg(:, 7, jj, ii)*vq(3, jat) &
                        &    +sdqg(:, 8, jj, ii)*vq(4, jat) &
                        &    +sdqg(:, 9, jj, ii)*vq(5, jat) &
                        &    +sdqg(:, 10, jj, ii)*vq(6, jat))
                     dG = stmp + dtmp + qtmp
                     dS = spread(dG, 1, 3) * spread(rij, 2, 3)
                     if (iat /= jat) then
                        sigma = sigma + dS
                        gradient(:, iat) = gradient(:, iat) + dG
                        gradient(:, jat) = gradient(:, jat) - dG
                     else
                        sigma = sigma + 0.5_wp * dS
                     endif
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo ! ish : loop over shells on iat
      enddo ! jat
   enddo ! iat
   !$omp end parallel do
end subroutine ddqint

end module aespot
