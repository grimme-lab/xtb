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

module aespot
   use iso_fortran_env, only : wp => real64
   use intpack, only : olap,divpt,rhftce,prod,opab1,opab4,propa
   integer,private, parameter :: llao (0:3) = (/ 1, 3, 6,10/)
   integer,private, parameter :: llao2(0:3) = (/ 1, 3, 5, 7/)

contains

! aesdqint: computes the dipole and quadrupole integrals 
!           and performs screening to determine, which contribute to 
!           potential
! thr         : neglect matrix elements in potential
! intcut      : integral cutoff according to prefactor from Gaussian product theorem 
! nat         : # of atoms
! nao         : # of spherical AOs (SAOs)
! nbf         : # of Cartesian AOs (CAOs)
! xyz(3,nat)  : cartesian coordinates
! caoshell    : map shell of atom to index in CAO space (lowest Cart. component is taken), dimension: (5,nat)
! saoshell     : map shell of atom to index in SAO space (lowest m_l component is taken), dimension: (5,nat)
! primcount   : index of first primitive (over entire system) of given CAO, dimension: nbf
! sint        : overlap integral matrix, simension nao,nao
! dpint       : dipole integral matrix, dimension 3,nao*(nao+1)/2
! qpint       : quadrupole integral matrix, dimension 6,nao*(nao+1)/2
! ndp,nqp     : number of elements to be computed in Fock matrix, i.e., with X-dip and X-qpole potential terms

subroutine sdqint(nat,at,nbf,nao,xyz,thr,ndp,nqp,intcut,caoshell,saoshell, &
      &           nprim,primcount,alp,cont,sint,dpint,qpint)
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   use lin_mod, only : lin
   implicit none
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
   real(wp),intent(out) :: dpint(3,nao*(nao+1)/2)
   real(wp),intent(out) :: qpint(6,nao*(nao+1)/2)
   integer, intent(out) :: ndp,nqp


   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij
   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp), parameter :: bohr=1.0_wp/0.52917726_wp

   real(wp)  ra(3),rb(3),f1,f2,point(3)
   real(wp) dtmp(3),qtmp(6),ss(6,6),dd(3,6,6),qq(6,6,6),tmp(6,6)
   integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,jshmax
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)
   parameter(itt  =(/0,1,4,10/))
   real(wp) :: saw(10)


   !call timing(t1,t3)
   thr2=(thr*1.0d-2)-thr*1.0d-12
   !     thr2=1.0d-20 ! conservative, keep all terms
   !     integrals
   sint=0.0_wp
   dpint=0.0_wp
   qpint=0.0_wp
   ! --- Aufpunkt for moment operator
   point=0.0_wp

   ndp=0
   nqp=0
   kj=0
   !$OMP PARALLEL PRIVATE (iat,jat,iatyp,cc,ci,ra,rb,saw, &  
   !$omp&               rab2,jatyp,ish,ishtyp,icao,naoi,iptyp, &
   !$omp&               jsh,jshmax,jshtyp,jcao,naoj,jptyp,ss,dd,qq, &
   !$omp&               est,alpi,alpj,ab,iprim,jprim,ip,jp, &
   !$omp&               mli,mlj,tmp,tmp1,tmp2,iao,jao,ii,jj,k,ij) &
   !$omp&               reduction (+:ndp,nqp) shared (sint,dpint,qpint)
   !$OMP DO schedule(dynamic) 
   do iat=1,nat
      ra(1:3)=xyz(1:3,iat)
      iatyp=at(iat)
      do jat=1,iat
         rb(1:3)=xyz(1:3,jat)
         jatyp=at(jat)
         rab2=sum( (rb-ra)**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if(rab2.gt.2000) cycle
         do ish=1,ao_n(iatyp)
            ishtyp=ao_l(ish,iatyp)
            icao=caoshell(ish,iat)
            naoi=llao(ishtyp)
            iptyp=itt(ishtyp)
            jshmax=ao_n(jatyp)
            if(iat.eq.jat) jshmax=ish
            !              jshells
            do jsh=1,jshmax
               jshtyp=ao_l(jsh,jatyp)
               jcao=caoshell(jsh,jat)
               naoj=llao(jshtyp)
               jptyp=itt(jshtyp)
               ! we go through the primitives (because the screening is the same for all of them)
               ss=0.0_wp
               dd=0.0_wp
               qq=0.0_wp
               call get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ra,rb,point, &
                  &               intcut,nprim,primcount,alp,cont,ss,dd,qq)
               !transform from CAO to SAO
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
                     if(jao.gt.iao) cycle
                     ij=lin(iao,jao)
                     sint(iao,jao) = ss(jj,ii)
                     sint(jao,iao) = ss(jj,ii)
                     dpint(1:3,ij)=dd(1:3,jj,ii)
                     tmp1=sum(dd(1:3,jj,ii)**2)
                     if(tmp1.gt.thr2) ndp=ndp+1
                     qpint(1:6,ij)=qq(1:6,jj,ii)
                     tmp2=sum(qq(1:3,jj,ii)**2)+2*sum(qq(4:6,jj,ii)**2)
                     if(tmp2.gt.thr2) nqp=nqp+1
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   !$OMP END DO
   !$OMP END PARALLEL
   !! DEBUG
   !                                                      call timing(t2,t4)
   !                                     call prtime(6,t2-t1,t4-t3,'dqint2')
   !      open(unit=136,file='test2',status='replace')
   !      call prmat(136,dpint,3,mm,'dip2')
   !      call prmat(136,qpint,6,mm,'qp2')
   !      close(136)

end subroutine sdqint

pure subroutine get_multiints(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj,point,intcut, &
      &                       nprim,primcount,alp,cont,ss,dd,qq)
   use intgrad
   implicit none
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: iptyp
   integer, intent(in)  :: jptyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: point(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: ss(:,:)
   real(wp),intent(out) :: dd(:,:,:)
   real(wp),intent(out) :: qq(:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc
   real(wp) :: ab,est,saw(10)

   real(wp),parameter :: max_r2 = 2000.0_wp

   ss = 0.0_wp
   dd = 0.0_wp
   qq = 0.0_wp

   rij = rj - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   do ip=1,nprim(icao+1)
      iprim=ip+primcount(icao+1)
      alpi=alp(iprim) ! exponent the same for each l component
      do jp=1,nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         alpj=alp(jprim) ! exponent the same for each l component
         ab=1.0_wp/(alpi+alpj)
         est=rij2*alpi*alpj*ab
         if(est.gt.intcut) cycle
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli=1,naoi
            iprim=ip+primcount(icao+mli)
            ci=cont(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            do mlj=1,naoj
               jprim=jp+primcount(jcao+mlj)
               saw=0.0_wp
               ! prim-prim quadrupole and dipole integrals
               call build_sdq_ints(ri,rj,point,alpi,alpj, &
                  & iptyp+mli,jptyp+mlj,saw)
               cc=cont(jprim)*ci
               ! from primitive integrals fill CAO-CAO matrix for ish-jsh block
               !                             ! overlap
               ss(mlj,mli)=ss(mlj,mli)+saw(1)*cc
               ! dipole
               dd(:,mlj,mli)=dd(:,mlj,mli)+saw(2:4)*cc
               ! quadrupole
               qq(:,mlj,mli)=qq(:,mlj,mli)+saw(5:10)*cc
            enddo ! mlj
         enddo ! mli
      enddo ! jp
   enddo ! ip 

end subroutine get_multiints

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
   real(wp),intent(in)  :: dpint(3,nao*(nao+1)/2)
   real(wp),intent(in)  :: qpint(6,nao*(nao+1)/2)
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
            tmp1=tmp1+dpint(k,ij)*dpint(k,ij)
            tmp2=tmp2-qpint(k,ij)*qpint(k,ij)
         enddo
         do k=1,6
            tmp2=tmp2+2.0_wp*qpint(k,ij)*qpint(k,ij)
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
   real(wp), intent(in) :: dpint(3,nao*(nao+1)/2),s(nao,nao)
   real(wp), intent(in) :: qpint(6,nao*(nao+1)/2),p(nao,nao)
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
            pdmk=pij*dpint(k,ij)
            dipm(k,jj)=dipm(k,jj)+xk2*ps-pdmk
            dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
            ! off-diagonal 
            do l=1,k-1
               kl=kl+1
               kj=k+l+1 
               xl1=ra(l)
               xl2=xyz(l,jj)
               pdml=pij*dpint(l,ij) 
               pqm=pij*qpint(kj,ij)
               tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
               tjj=pdmk*xl2+pdml*xk2-xl2*xk2*ps-pqm
               qp(kl,jj)=qp(kl,jj)+tjj
               qp(kl,ii)=qp(kl,ii)+tii
            enddo
            ! diagonal 
            kl=kl+1
            pqm=pij*qpint(k,ij)
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
         pdmk=pij*dpint(k,ij)
         dipm(k,ii)=dipm(k,ii)+xk1*ps-pdmk
         ! off-diagonal 
         do l=1,k-1
            kl=kl+1
            kj=k+l+1 ! the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
            xl1=ra(l)
            pdml=pij*dpint(l,ij)
            pqm=pij*qpint(kj,ij)
            tii=pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
            qp(kl,ii)=qp(kl,ii)+tii
         enddo
         !diagonal 
         kl=kl+1
         pqm=pij*qpint(k,ij)
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
   real(wp), intent(in) :: dpint(3,nao*(nao+1)/2),s(nao,nao)
   real(wp), intent(in) :: qpint(6,nao*(nao+1)/2),p(nao,nao)
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
            fji=fji+dpint(k,ij)*(vd(k,ii)+vd(k,jj))
         enddo
         do k=1,6
            fji=fji+qpint(k,ij)*(vq(k,ii)+vq(k,jj))
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
      & radcn,dcn,gab3,gab5,g) 
   use aoparam
   use lin_mod, only : lin
   !gab3 Hellmann-Feynman terms correct, shift terms to be tested yet
   implicit none          
   integer, intent(in)   :: nat,at(nat)
   real(wp), intent(in)    :: xyz(3,nat),q(nat),dipm(3,nat),qp(6,nat)
   real(wp), intent(in)    :: gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
   real(wp), intent(in)    :: kdmp3,kdmp5,radcn(nat),dcn(3,nat,nat)    
   real(wp), intent(inout) :: g(3,nat)
   real(wp) qp1(6),rr(3),dip(3),rij(3)
   real(wp) ed,eq,edd,e01,e02,e11,r2,tt,tt3,q1,dxi
   real(wp) tmp2,tmp3,rab,rabi,ddm2,ddm3a,ddm3b,qqa,qqb
   real(wp) dgab3,dgab5,damp1,damp2,ddamp,qs2

   integer i,j,k,l,m,ki,kj,kl,mm
   do i=1,nat 
      q1=q(i)
      rr(1:3)=xyz(1:3,i)
      dip(1:3)=dipm(1:3,i)
      qp1(1:6)=qp(1:6,i)
      tmp2=0.0_wp ! cumulate terms propto CN gradient -  to scale only quadratically
      do j=1,nat            ! loop over other atoms            
         if(i.eq.j) cycle
         kj=lin(j,i)
         rij(1:3)=xyz(1:3,j)-rr(1:3)            
         r2=sum(rij*rij)
         rabi=1.0_wp/sqrt(r2)
         !           call dzero(2.0_wp,rabi,at(i),at(j),damp,ddamp)
         call dzero(kdmp3,rabi,radcn(i),radcn(j),damp1,ddamp)
         dgab3=dgab(3.0_wp,rabi,damp1,ddamp)
         !           call dzero(3.0_wp,rabi,at(i),at(j),damp,ddamp)
         call dzero(kdmp5,rabi,radcn(i),radcn(j),damp2,ddamp)
         dgab5=dgab(5.0_wp,rabi,damp2,ddamp)
         !!!         DEBUG
         !            dgab3=0.0_wp
         !            dgab5=0.0_wp 
         !            dgab3=dgab3*100.0_wp
         !            dgab5=dgab5*100.0_wp
         !!!
         ed=0.0_wp
         edd=0.0_wp
         eq=0.0_wp
         !           dipole - charge
         do k=1,3
            ed=ed+q(j)*dip(k)*rij(k)
            ed=ed-dipm(k,j)*q1*rij(k)
            tt=q1*dipm(k,j)-q(j)*dip(k)
            ! part of dip-q derivative
            g(k,i)=g(k,i)+gab3(kj)*tt 
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
               eq=eq+q(j)*qp1(kl)*tt
               eq=eq+qp(kl,j)*q1*tt
               edd=edd-dipm(k,j)*dip(l)*tt3
               ! extra d-d terms
               ddm2=ddm2+dipm(l,j)*dip(l)
               ddm3a=ddm3a+dip(l)*rij(l)
               ddm3b=ddm3b+dipm(l,j)*rij(l)
               ! extra q-qpole terms
               qqa=qqa+rij(l)*qp(kl,j)
               qqb=qqb+rij(l)*qp1(kl)
            enddo
            edd=edd+dipm(k,j)*dip(k)*r2
            g(k,i)=g(k,i)-2.0_wp*gab5(kj)*ddm2*rij(k)
            g(k,i)=g(k,i)+3.0_wp*gab5(kj)*ddm3a*dipm(k,j)
            g(k,i)=g(k,i)+3.0_wp*gab5(kj)*ddm3b*dip(k)
            g(k,i)=g(k,i)-2.0_wp*gab5(kj)*qqa*q1
            g(k,i)=g(k,i)-2.0_wp*gab5(kj)*qqb*q(j)
         enddo
         do k=1,3
            dxi=rij(k)*rabi 
            g(k,i)=g(k,i)-ed*dgab3*dxi
            g(k,i)=g(k,i)-(eq+edd)*dgab5*dxi
         enddo
         ! collect terms for CN-dependent part
         rab=0.50_wp*(radcn(i)+radcn(j)) 
         tmp3=ed*kdmp3*gab3(kj)*(damp1/rab)*(rab*rabi)**kdmp3
         tmp2=tmp2+tmp3
         tmp3=(eq+edd)*kdmp5*gab5(kj)*(damp2/rab)*(rab*rabi)**kdmp5
         tmp2=tmp2+tmp3
      enddo
      ! CN-dependent part  - O(N^2) 
      tmp2=3.0_wp*tmp2 
      do j=1,nat
         do k=1,3
            g(k,j)=g(k,j)-tmp2*dcn(k,j,i)
         enddo
      enddo

   enddo
end subroutine aniso_grad

subroutine ovlp2(thr,nat,nao,nbf,at,xyz,caoshell,saoshell,nprim,primcount,cont,alp,s)
   use mctc_constants, only : pi
   use aoparam
   use intgrad, only : dtrf2
   implicit none          
   integer, intent(in)  :: nat
   integer, intent(in)  :: nao
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: thr
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: caoshell(:,:)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   real(wp),intent(out) :: s(nao,nao)

   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2,dx,dy,dz,s00r,s00l,s00
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp)  ra(3),rb(3),f1,f2,dpl(3),dpr(3),qpl(6),qpr(6)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,li,lj,iprim,jprim
   integer itt(0:3)
   data itt  /0,1,4,10/
   real(wp) ss (6,6),stmp (6)

   thr2=thr

   s=0.0_wp
   kj=0
   mm=nbf*(nbf+1)/2
   step=1.0d-4
   step2=0.50_wp/step
!$OMP parallel default(none) &
!$omp&         private(iat,ra,iatyp,jat,rb,jatyp,dx,dy,dz,rab2,ish,ishtyp, &
!$omp&                 icao,naoi,iptyp,jshmax,jsh,jshtyp,jcao,naoj,ss,ip,iprim, &
!$omp&                 alpi,jp,jprim,ab,est,s00,li,ci,stmp,lj,cc,ii,iao,jj,jao) &
!$omp&         shared(s,cont,alp,saoshell,caoshell,thr2,primcount,nprim,ao_l, &
!$omp&                ao_n,itt,xyz,at,nat)
!$omp do schedule(dynamic)
   do iat=1,nat
      ra(1:3)=xyz(1:3,iat)
      iatyp=at(iat)
      !         write (*,*), 'i',iat, iatyp, ao_n(iatyp)
      do jat=1,iat
         rb(1:3)=xyz(1:3,jat)
         jatyp=at(jat)
         dx=ra(1)-rb(1)
         dy=ra(2)-rb(2)
         dz=ra(3)-rb(3)
         rab2=dx*dx+dy*dy+dz*dz           
         !c          ints < 1.d-9 for RAB > 40 Bohr            
         if(rab2.gt.2000) cycle
         !            write(*,*) 'j',jat,jatyp,ao_n(jatyp)
         do ish=1,ao_n(iatyp)
            ishtyp=ao_l(ish,iatyp)
            icao=caoshell(ish,iat)+1
            naoi=llao(ishtyp)
!           write(*,*) 'i',iat,iatyp,ish,ishtyp,icao,naoi
            iptyp=itt(ishtyp)
            jshmax=ao_n(jatyp)
            if(iat.eq.jat) jshmax=ish
            !              jshells
            do jsh=1,jshmax
               jshtyp=ao_l(jsh,jatyp)             
               jcao=caoshell(jsh,jat)+1
               naoj=llao(jshtyp)
!              write(*,*) 'j',jat,jatyp,jsh,jshtyp,jcao,naoj
               ! we go through the primitives (beacause the screening is the same for all of them) 
               ss=0.0_wp       
               do ip=1,nprim(icao)
                  iprim=ip+primcount(icao) 
                  alpi=alp(iprim)          ! exponent the same for each l component 
                  do jp=1,nprim(jcao)
                     jprim=jp+primcount(jcao)
                     ab=1.0_wp/(alpi+alp(jprim))
                     est=rab2*alpi*alp(jprim)*ab
                     if(est.gt.thr2) cycle  
                     s00=(pi*ab)**1.50_wp*exp(-est) 
                     ! now compute overlap for different components of i(e.g., px,py,pz)
                     do li=1,naoi
!                       iprim=ip+primcount(icao+li-1)
                        ci=cont(ip+primcount(icao+li-1))
!                       print*,allocated(alp),jprim,alp(jprim)
                        call sspd2(ra,rb,iptyp+li,jshtyp+1,alpi, &
                           & alp(jprim),ab,s00,stmp)         
                        do lj=1,naoj
!                          jprim=jp+primcount(jcao+lj-1)
                           cc=cont(jp+primcount(jcao+lj-1))*ci
                           ! fill CAO-CAO overlap matrix for ish-jsh pair
                           ss(lj,li)=ss(lj,li)+stmp(lj)*cc
                        enddo
                     enddo
                  enddo
               enddo
               !transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)
               do ii=1,llao2(ishtyp)
                  iao=ii+saoshell(ish,iat)
                  do jj=1,llao2(jshtyp)
                     jao=jj+saoshell(jsh,jat)
                     s(jao,iao)=ss(jj,ii)
                     s(iao,jao)=ss(jj,ii)
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
!$omp enddo
!$omp end parallel
end subroutine ovlp2




!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! calc. overlap integrals betweeen unormalized C
! cartesian gaussians centered at xyza and     C
! xyzb with alpha_a=alpi and beta_b=alpj C
! est is r_AB^2*alpha*beta/(alpha+beta)        C
! D is the result                              C
! S,P,D only !!!!!!!!!!!!!!!!!!!               C
! special task is to obtain all components of jC
! e.g. i=s,j=d gives 6 overlap integrals       C
! s.grimme, june 1998                          C
! recursion formulas from Obara and Saika,     C
! JCP 84 (1986) 3963                           C
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

pure subroutine sspd2(xyza,xyzb,iptyi,iptyj,alpi,alpj,ab,s00,d)
   implicit none
   real(wp),intent(in)  :: xyza(3)
   real(wp),intent(in)  :: xyzb(3)
   integer, intent(in)  :: iptyi
   integer, intent(in)  :: iptyj
   real(wp),intent(in)  :: alpi
   real(wp),intent(in)  :: alpj
   real(wp),intent(in)  :: ab
   real(wp),intent(in)  :: s00
   real(wp),intent(out) :: d(6)
   real(wp) :: p(3),ps(3),dp(3)
   real(wp) :: pip(3),pjp(3),deltapb(3)
   integer,parameter :: ii(5:10) = (/1,2,3,1,1,2/)
   integer,parameter :: jj(5:10) = (/1,2,3,2,3,3/)
   real(wp) :: dum
   real(wp) :: dp1
   real(wp) :: dp2
   real(wp) :: dp3
   real(wp) :: sp1 
   real(wp) :: sp2 
   real(wp) :: sp3 
   real(wp) :: ab05
   integer :: id
   real(wp) :: xps
   real(wp) :: sdxx
   real(wp) :: sdyy
   real(wp) :: sdzz
   real(wp) :: sdxy
   real(wp) :: sdxz
   real(wp) :: sdyz
   integer :: iii
   integer :: jjj
   integer :: kkk
   integer :: lll
   real(wp) :: fij
   real(wp) :: ds
   real(wp) :: fik
   real(wp) :: fjk
   integer :: m
   real(wp) :: ab05ds

   d=0

   if(iptyi.eq.1) then
      goto 100
   elseif(iptyi.ge.2.and.iptyi.le.4) then
      goto 200
   elseif(iptyi.ge.5.and.iptyi.le.10) then
      goto 200
   else
      return!stop 'programing error in sspd 1'
   endif

   !CCCCCCCCCCCCCCCCCCCCCCC
   !     S  --- S,P,D     C
   !CCCCCCCCCCCCCCCCCCCCCCC

   100 select case(iptyj)
   case(1)
      ! s-s
      d(1)=s00
   case(2)
      ! s-p
      dum=alpi*ab*s00
      d(1)=dum*(xyza(1)-xyzb(1))
      d(2)=dum*(xyza(2)-xyzb(2))
      d(3)=dum*(xyza(3)-xyzb(3))
   case(3)
      ! s-d
      p(1)=ab*(alpi*xyza(1)+alpj*xyzb(1))
      p(2)=ab*(alpi*xyza(2)+alpj*xyzb(2))
      p(3)=ab*(alpi*xyza(3)+alpj*xyzb(3))
      dum=0.5*ab*s00
      dp1=p(1)-xyzb(1)
      dp2=p(2)-xyzb(2)
      dp3=p(3)-xyzb(3)
      sp1 =s00*dp1            
      sp2 =s00*dp2           
      sp3 =s00*dp3           
      d(1)=    dp1       *sp1+dum       
      d(2)=    dp2       *sp2+dum        
      d(3)=    dp3       *sp3+dum       
      d(4)=    dp2       *sp1
      d(5)=    dp3       *sp1
      d(6)=    dp3       *sp2
   end select
   return

   !CCCCCCCCCCCCCCCCCCCCCCC
   !     P  --- S,P,D     C
   !CCCCCCCCCCCCCCCCCCCCCCC

   200 p(1)=ab*(alpi*xyza(1)+alpj*xyzb(1))
   p(2)=ab*(alpi*xyza(2)+alpj*xyzb(2))
   p(3)=ab*(alpi*xyza(3)+alpj*xyzb(3))
   ab05=ab*0.5

   id=iptyi-1
   select case(iptyj)
   case(1)
      ! p-s
      d(1)=(p(id)-xyza(id))*s00
   case(2)
      ! p-p
      xps=(p(id)-xyza(id))*s00         
      d(1)=(p(1)-xyzb(1))*xps
      d(2)=(p(2)-xyzb(2))*xps
      d(3)=(p(3)-xyzb(3))*xps
      d(id)=d(id)+ab05*s00
   case(3)
      ! p-d
      dum=ab05*s00
      dp1=p(1)-xyzb(1)
      dp2=p(2)-xyzb(2)
      dp3=p(3)-xyzb(3)
      sp1=s00*dp1
      sp2=s00*dp2
      sp3=s00*dp3
      sdxx=    dp1       *sp1+dum       
      sdyy=    dp2       *sp2+dum        
      sdzz=    dp3       *sp3+dum       
      sdxy=    dp2       *sp1
      sdxz=    dp3       *sp1
      sdyz=    dp3       *sp2
   end select

   select case(iptyi)
   case(1)
      return!stop 'programing error in sspd 2'
   case(2)
      dp1=p(1)-xyza(1)
      d(1)=dp1*sdxx+  ab*sp1 
      d(2)=dp1*sdyy
      d(3)=dp1*sdzz
      d(4)=dp1*sdxy+ab05*sp2
      d(5)=dp1*sdxz+ab05*sp3
      d(6)=dp1*sdyz
   case(3)
      dp1=p(2)-xyza(2)
      d(1)=dp1*sdxx
      d(2)=dp1*sdyy+  ab*sp2
      d(3)=dp1*sdzz
      d(4)=dp1*sdxy+ab05*sp1
      d(5)=dp1*sdxz
      d(6)=dp1*sdyz+ab05*sp3
   case(4)
      dp1=p(3)-xyza(3)
      d(1)=dp1*sdxx
      d(2)=dp1*sdyy
      d(3)=dp1*sdzz+  ab*sp3
      d(4)=dp1*sdxy
      d(5)=dp1*sdxz+ab05*sp1
      d(6)=dp1*sdyz+ab05*sp2
   end select
   return

   !CCCCCCCCCCCCCCCCCCCCCCC
   !     D  --- S,P,D     C
   !CCCCCCCCCCCCCCCCCCCCCCC

   300 p(1)=ab*(alpi*xyza(1)+alpj*xyzb(1))
   p(2)=ab*(alpi*xyza(2)+alpj*xyzb(2))
   p(3)=ab*(alpi*xyza(3)+alpj*xyzb(3))
   ab05=ab*0.5
   iii=ii(iptyi)
   jjj=jj(iptyi)   

   select case(iptyj)
   case(1)
      ! d-s
      fij=0.0_wp
      if(iii.eq.jjj) fij=ab05
      d(1)=(p(jjj)-xyza(jjj))*(p(iii)-xyza(iii))*s00+fij*s00
   case(2)
      ! d-p
      fij=0.0_wp
      if(iii.eq.jjj) fij=ab05
      ds=(p(jjj)-xyza(jjj))*(p(iii)-xyza(iii))*s00+fij*s00
      ps(1)=(p(1)-xyza(1))*s00
      ps(2)=(p(2)-xyza(2))*s00
      ps(3)=(p(3)-xyza(3))*s00
      do m=1,3
         fik=0.0_wp
         fjk=0.0_wp
         if(iii.eq.m) fik=ab05
         if(jjj.eq.m) fjk=ab05  
         d(m)=(p(m)-xyzb(m))*ds+fik*ps(jjj)+fjk*ps(iii)  
      enddo
   case(3)
      ! d-d
      fij=0.0_wp
      if(iii.eq.jjj) fij=ab05
      ds=(p(jjj)-xyza(jjj))*(p(iii)-xyza(iii))*s00+fij*s00
      ps(1)=(p(1)-xyza(1))*s00
      ps(2)=(p(2)-xyza(2))*s00
      ps(3)=(p(3)-xyza(3))*s00
      deltapb(1)=(p(1)-xyzb(1))
      deltapb(2)=(p(2)-xyzb(2))
      deltapb(3)=(p(3)-xyzb(3))
      do m=1,3
         fik=0.0_wp
         fjk=0.0_wp
         if(iii.eq.m) fik=ab05
         if(jjj.eq.m) fjk=ab05  
         dp(m)=deltapb(m)*ds+fik*ps(jjj)+fjk*ps(iii)  
      enddo
      do m=1,3
         pip(m)=deltapb(m)*ps(iii)
         if(iii.eq.m) pip(m)=pip(m)+ab05*s00
      enddo
      if(iii.ne.jjj) then
         do m=1,3
            pjp(m)=deltapb(m)*ps(jjj)
            if(jjj.eq.m) pjp(m)=pjp(m)+ab05*s00
         enddo
      else
         pjp(1)=pip(1)
         pjp(2)=pip(2)
         pjp(3)=pip(3)
      endif
      ab05ds=ab05*ds
      d(1)=deltapb(1)    *dp(1)+ab05ds 
      d(2)=deltapb(2)    *dp(2)+ab05ds 
      d(3)=deltapb(3)    *dp(3)+ab05ds
      d(4)=deltapb(2)    *dp(1) 
      d(5)=deltapb(3)    *dp(1) 
      d(6)=deltapb(3)    *dp(2) 
      do m=1,6
         kkk=ii(m+4)
         lll=jj(m+4)
         if(iii.eq.lll) d(m)=d(m)+ab05*pjp(kkk)
         if(jjj.eq.lll) d(m)=d(m)+ab05*pip(kkk)
      enddo
   end select

   return
end subroutine sspd2


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
   integer, intent (in) :: n,at(n)
   real(wp), intent (in)  :: cn(n),shift,expo,rmax
   real(wp), intent (out) :: radcn(n)
   real(wp) rco,t1,t2 
   integer i,j
   do i=1,n
      rco=radaes(at(i))             ! base radius of element
      t1 =cn(i)-cnval(at(i))-shift  ! CN - VALCN - SHIFT
      t2 =rco +(rmax-rco)/(1.0_wp+exp(-expo*t1))
      radcn(i)=t2
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
subroutine dradcn(n,at,cn,shift,expo,rmax,dcn)
   use aoparam
   implicit none
   integer, intent (in) :: n,at(n)
   real(wp), intent (in)  :: cn(n),shift,expo,rmax
   real(wp), intent (inout) :: dcn(3,n,n)
   real(wp) rco,t1,t2,t3,t4,tmp1,tmp2
   integer i,j,k
   do i=1,n
      rco=radaes(at(i))             ! base radius of element
      t1 =exp(-expo*(cn(i)-cnval(at(i))-shift))  ! CN - VALCN - SHIFT
      t2 =(rmax-rco)/(1.0_wp+2.0_wp*t1+t1*t1)
      t2 = t2*expo*t1
      do j=1,i
         rco=radaes(at(j))             ! base radius of element
         t3 =exp(-expo*(cn(j)-cnval(at(j))-shift))  ! CN - VALCN - SHIFT
         t4 =(rmax-rco)/(1.0_wp+2.0_wp*t3+t3*t3)
         t4 = t4*expo*t3
         do k=1,3
            tmp1=dcn(k,j,i)*t4
            tmp2=dcn(k,i,j)*t2
            dcn(k,i,j)=tmp1
            dcn(k,j,i)=tmp2
         enddo
      enddo
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
   real(wp), intent(in)  :: dpint(3,nao*(nao+1)/2)
   real(wp), intent(in)  :: qpint(6,nao*(nao+1)/2)

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
         x(j)=dpint(1,k)*dpint(1,k)+dpint(2,k)*dpint(2,k) &
            & +dpint(3,k)*dpint(3,k)
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
         x(j)=qpint(1,k)*qpint(1,k)+qpint(2,k)*qpint(2,k) &
            & +qpint(3,k)*qpint(3,k)+2.0_wp*qpint(4,k)*qpint(4,k) &
            & +2.0_wp*qpint(5,k)*qpint(5,k)+2.0_wp*qpint(6,k)*qpint(6,k)
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



! not used at the moment (maybe for GPUs?)  
! camm2: computes the CAMM directly from integrals (w/ density matrix as input)
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
subroutine camm2(thr,nat,nao,nbf,at,xyz,caoshell,saoshell,nprim,primcount,alp,cont,p,q,dipm,qp)
   use mctc_constants, only : pi
   use aoparam
   use intgrad, only : dtrf2
   implicit none          
   integer, intent(in)  :: nat
   integer, intent(in)  :: nao
   integer, intent(in)  :: nbf
   integer, intent(in)  :: at(nat)
   real(wp),intent(in)  :: thr
   real(wp),intent(in)  :: xyz(3,nat)
   integer, intent(in)  :: caoshell(:,:)
   integer, intent(in)  :: saoshell(:,:)
   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)
   real(wp),intent(in)  :: p(nao,nao)
   real(wp),intent(out) :: q(nat)
   real(wp),intent(out) :: dipm(3,nat)
   real(wp),intent(out) :: qp(6,nat)

   real(wp) tmp1,tmp2,tmp3,tmp4,step,step2
   real(wp) dx,dy,dz,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp), parameter :: bohr=1.0_wp/0.52917726_wp
   real(wp)  ra(3),rb(3),f1,f2,point(3) 
   real(wp) ss(6,6),stmp(6)
   real(wp) dtmp(3),qtmp(6),dd(3,6,6),qq(6,6,6)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   integer itt(0:3)
   data itt  /0,1,4,10/

   q=0.0_wp
   dipm=0.0_wp
   qp=0.0_wp
   thr2=thr!*1.0d4-thr*1.0d-6
   !      thr2=-1.0_wp ! conservative, keep all terms
   point=0.0_wp

   kj=0
   mm=nbf*(nbf+1)/2
   step=1.0d-4
   step2=0.50_wp/step
   do iat=1,nat
      ra(1:3)=xyz(1:3,iat)
      iatyp=at(iat)
      do jat=1,iat
         rb(1:3)=xyz(1:3,jat)
         jatyp=at(jat)
         dx=rb(1)-ra(1)
         dy=rb(2)-ra(2)
         dz=rb(3)-ra(3)
         rab2=dx*dx+dy*dy+dz*dz           
         !c          ints < 1.d-9 for RAB > 40 Bohr            
         if(rab2.gt.2000) cycle            
         do ish=1,ao_n(iatyp)
            ishtyp=ao_l(ish,iatyp)
            icao=caoshell(ish,iat)
            naoi=llao(ishtyp)
            iptyp=itt(ishtyp)
            jshmax=ao_n(jatyp)
            if(iat.eq.jat) jshmax=ish
            !              jshells
            do jsh=1,jshmax
               jshtyp=ao_l(jsh,jatyp)             
               jcao=caoshell(jsh,jat)
               naoj=llao(jshtyp)
               print*,naoj
               jptyp=itt(jshtyp)
               ! we go through the primitives (beacause the screening is the same for all of them) 
               ss=0.0_wp      
               dd=0.0_wp
               qq=0.0_wp 
               do ip=1,nprim(icao+1)
                  iprim=ip+primcount(icao+1) 
                  alpi=alp(iprim)          ! exponent the same for each l component 
                  do jp=1,nprim(jcao+1)
                     jprim=jp+primcount(jcao+1)
                     alpj=alp(jprim) ! exponent the same for each l component
                     ab=1.0_wp/(alpi+alpj)
                     est=rab2*alpi*alpj*ab
                     if(est.gt.thr2) cycle  
                     s00=(pi*ab)**1.50_wp*exp(-est) 
                     ! now compute integrals  for different components of i(e.g., px,py,pz)
                     do mli=1,naoi
                        iprim=ip+primcount(icao+mli)
                        ci=cont(iprim) ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
                        call sspd2(ra,rb,iptyp+mli,jshtyp+1,alpi, &
                           & alpj,ab,s00,stmp)        
                        do mlj=1,naoj
                           jprim=jp+primcount(jcao+mlj)
                           qtmp=0.0_wp
                           dtmp=0.0_wp
                           ! prim-prim qpole and dipole integrals
                           call propa(opab4,ra,rb,point,alpi, &
                              & alpj,iptyp+mli,jptyp+mlj,qtmp,6)
                           call propa(opab1,ra,rb,point,alpi, &
                              & alpj,iptyp+mli,jptyp+mlj,dtmp,3)

                           cc=cont(jprim)*ci
                           ! from primitive integrals fill CAO-CAO  matrix for ish-jsh block
                           ss(mlj,mli)=ss(mlj,mli)+stmp(mlj)*cc
                           ! dipole 
                           do k=1,3 
                              dd(k,mlj,mli)=dd(k,mlj,mli)-dtmp(k)*cc
                           enddo
                           ! quadrupole 
                           do k=1,6
                              qq(k,mlj,mli)=qq(k,mlj,mli)-qtmp(k)*cc
                           enddo
                        enddo ! mlj : Cartesian component of j prims
                     enddo    ! mli : Cartesian component of i prims
                  enddo ! jp : loop over j prims
               enddo    ! ip : loop over i prims

               !transform from CAO to SAO
               call dtrf2(ss,ishtyp,jshtyp)
               do k=1,3
                  call dtrf2(dd(k,:,:),ishtyp,jshtyp)
               enddo
               do k=1,6
                  call dtrf2(qq(k,:,:),ishtyp,jshtyp)
               enddo

               !!! directly contract w/ density matrix to compute atomic moments,
               !   w/o "shift terms" from lower moments (is done later atom-wise)

               ! special diagonal block case
               if(icao.eq.jcao) then 
                  do ii=1,llao2(ishtyp)
                     iao=ii+saoshell(ish,iat)
                     do jj=1,ii-1
                        jao=jj+saoshell(jsh,jat)
                        cc=2.0_wp*p(jao,iao)
                        q(jat)=q(jat)-cc*ss(jj,ii)
                        dipm(1:3,jat)=dipm(1:3,jat)-cc*dd(1:3,jj,ii)
                        qp(1:6,jat)=qp(1:6,jat)-cc*qq(1:6,jj,ii)
                     enddo
                     cc=p(iao,iao)
                     q(iat)=q(iat)-cc*ss(ii,ii)
                     dipm(1:3,iat)=dipm(1:3,iat)-cc*dd(1:3,ii,ii)
                     qp(1:6,iat)=qp(1:6,iat)-cc*qq(1:6,ii,ii)
                  enddo
               else
                  ! off-diagonal blocks
                  do ii=1,llao2(ishtyp)
                     iao=ii+saoshell(ish,iat)
                     do jj=1,llao2(jshtyp)
                        jao=jj+saoshell(jsh,jat)
                        cc=p(jao,iao)
                        q(jat)=q(jat)-cc*ss(jj,ii)
                        q(iat)=q(iat)-cc*ss(jj,ii)
                        dipm(1:3,jat)=dipm(1:3,jat)-cc*dd(1:3,jj,ii)
                        dipm(1:3,iat)=dipm(1:3,iat)-cc*dd(1:3,jj,ii)
                        qp(1:6,jat)=qp(1:6,jat)-cc*qq(1:6,jj,ii)
                        qp(1:6,iat)=qp(1:6,iat)-cc*qq(1:6,jj,ii)
                     enddo
                  enddo
               endif ! off-diag/dig?

            enddo ! jsh : loop over shells on jat
         enddo    ! ish : loop over shells on iat
      enddo ! jat
   enddo    ! iat

   ! now add "shift terms" from lower moments
   do i=1,nat  
      dipm(1:3,i)=dipm(1:3,i)-xyz(1:3,i)*q(i)
      do k=1,3
         do l=1,k-1
            kl=k+l+1 ! the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
            qp(kl,i)=qp(kl,i)-xyz(k,i)*dipm(l,i)
            qp(kl,i)=qp(kl,i)-xyz(l,i)*dipm(k,i)
            qp(kl,i)=qp(kl,i)-xyz(k,i)*xyz(l,i)*q(i)
         enddo
         qp(k,i)=qp(k,i)-2.0_wp*xyz(k,i)*dipm(k,i)
         qp(k,i)=qp(k,i)-xyz(k,i)*xyz(k,i)*q(i)
      enddo
      ! now resort quadrupole moments, scale, and remove trace
      qtmp(1)=qp(1,i)
      qtmp(2)=qp(4,i)
      qtmp(3)=qp(2,i)
      qtmp(4)=qp(5,i)
      qtmp(5)=qp(6,i)
      qtmp(6)=qp(3,i)
      tmp1=qtmp(1)+qtmp(3)+qtmp(6) 
      tmp1=tmp1*0.50_wp
      qtmp(1:6)=qtmp(1:6)*1.50_wp
      qp(1,i)=qtmp(1)-tmp1 
      qp(2,i)=qtmp(2)
      qp(3,i)=qtmp(3)-tmp1
      qp(4,i)=qtmp(4)
      qp(5,i)=qtmp(5)
      qp(6,i)=qtmp(6)-tmp1 
   enddo
end subroutine camm2


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
subroutine ddqint(intcut,nat,nao,nbf,at,xyz,caoshell,saoshell,nprim,primcount,alp,cont,p,vs,vd,vq,H,g)
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   implicit none
   integer, intent(in)    :: nat
   integer, intent(in)    :: nao
   integer, intent(in)    :: nbf
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: intcut
   real(wp),intent(in)    :: vd(3,nat)
   real(wp),intent(in)    :: vq(6,nat)
   real(wp),intent(in)    :: H(nao,nao)
   real(wp),intent(in)    :: vs(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   integer, intent(in)    :: caoshell(:,:)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   real(wp),intent(in)    :: p(nao,nao)
   real(wp),intent(inout) :: g(3,nat)

   ! OMP
   common /proc/ nproc
   integer nproc
   integer OMP_GET_THREAD_NUM

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4,thr2,f,ci,cc,cj,alpi,rij2,ab,est
   real(wp) f1,f2,point(3),tmp(6,6),rij(3),ri(3),rj(3)
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3,10),dum(10),sdq(10,6,6),sdqg(3,19,6,6)

   thr2=intcut
   point=0.0_wp
   ! call timing(t1,t3)
   !$OMP PARALLEL PRIVATE(iat,jat,ixyz,iatyp,cc,ci,rij2,jatyp,ish,ishtyp, &
   !$omp&               icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&               sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp,ri,rj,rij, &
   !$omp&               mli,mlj,dum,dumdum,tmp,stmp,dtmp,qtmp, &
   !$omp&               iao,jao,ii,jj,k ) reduction (+:g )
   !$OMP DO schedule(dynamic)
   do iat=1,nat
      ri = xyz(:,iat)
      iatyp=at(iat)
      do jat=1,iat-1
         rj = xyz(:,jat)
         !           if (jat.eq.iat) cycle
         jatyp=at(jat)
         rij = ri - rj
         rij2= sum( rij**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if (rij2.gt.2000) cycle
         do ish=1,ao_n(iatyp)
            ishtyp=ao_l(ish,iatyp)
            icao=caoshell(ish,iat)
            naoi=llao(ishtyp)
            iptyp=itt(ishtyp)
            jshmax=ao_n(jatyp)
            !              if(iat.eq.jat) jshmax=ish
            do jsh=1,jshmax ! jshells
               jshtyp=ao_l(jsh,jatyp)
               jcao=caoshell(jsh,jat)
               naoj=llao(jshtyp)
               jptyp=itt(jshtyp)
               sdqg=0;sdq=0
               call get_grad_multiint(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
                  &                   intcut,nprim,primcount,alp,cont,sdq,sdqg)
               do k=1,19 ! 1 S, 2-4 D, 5-10 Q, 11-13 D, 14-19 Q
                  do ixyz=1,3
                     ! transform from CAO to SAO
                     !                        call dtrf2(sdqg(ixyz,k,1:6,1:6),ishtyp,jshtyp)
                     tmp(1:6,1:6)=sdqg(ixyz,k,1:6,1:6)
                     call dtrf2(tmp,ishtyp,jshtyp)
                     sdqg(ixyz,k,1:6,1:6)=tmp(1:6,1:6)
                  enddo
               enddo
               do ii=1,llao2(ishtyp)
                  iao=ii+saoshell(ish,iat)
                  do jj=1,llao2(jshtyp)
                     jao=jj+saoshell(jsh,jat)
                     cc=p(jao,iao)
                     do ixyz=1,3
                        stmp=sdqg(ixyz,1,jj,ii)*(2*H(jao,iao) &
                           & +cc*(vs(iat)+vs(jat)))
                        dtmp=cc*sum(sdqg(ixyz,11:13,jj,ii)*vd(1:3,iat) &
                           & +sdqg(ixyz, 2:4, jj,ii)*vd(1:3,jat) )
                        qtmp=cc*sum( sdqg(ixyz,14:19,jj,ii)*vq(1:6,iat) &
                           & +sdqg(ixyz, 5:10,jj,ii)*vq(1:6,jat) )
                        g(ixyz,iat)=g(ixyz,iat)+stmp+dtmp+qtmp
                        g(ixyz,jat)=g(ixyz,jat)-stmp-dtmp-qtmp
                     enddo ! ixyz
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat
   !$OMP END DO
   !$OMP END PARALLEL
   !                                                      call timing(t2,t4)
   !                                     call prtime(6,t2-t1,t4-t3,'dqint5')
end subroutine ddqint

pure subroutine get_grad_multiint(icao,jcao,naoi,naoj,iptyp,jptyp,ri,rj, &
      &                           intcut,nprim,primcount,alp,cont,sdq,sdqg)
   use intgrad
   implicit none
   integer, intent(in)  :: icao
   integer, intent(in)  :: jcao
   integer, intent(in)  :: naoi
   integer, intent(in)  :: naoj
   integer, intent(in)  :: iptyp
   integer, intent(in)  :: jptyp
   real(wp),intent(in)  :: ri(3)
   real(wp),intent(in)  :: rj(3)
   real(wp),intent(in)  :: intcut
   real(wp),intent(out) :: sdq(:,:,:)
   real(wp),intent(out) :: sdqg(:,:,:,:)

   integer, intent(in)  :: nprim(:)
   integer, intent(in)  :: primcount(:)
   real(wp),intent(in)  :: alp(:)
   real(wp),intent(in)  :: cont(:)

   integer  :: ip,iprim,mli,jp,jprim,mlj
   real(wp) :: rij(3),rij2,alpi,alpj,ci,cj,cc
   real(wp) :: ab,est,saw(10),sawg(3,10)

   real(wp),parameter :: max_r2 = 2000.0_wp

   sdqg = 0.0_wp
   sdq  = 0.0_wp

   rij = ri - rj
   rij2 = rij(1)**2 + rij(2)**2 + rij(3)**2

   if(rij2.gt.max_r2) return

   ! we go through the primitives (because the screening is the same for all of them)
   do ip=1,nprim(icao+1)
      iprim=ip+primcount(icao+1)
      ! exponent the same for each l component
      alpi=alp(iprim)
      do jp=1,nprim(jcao+1)
         jprim=jp+primcount(jcao+1)
         ! exponent the same for each l component
         alpj=alp(jprim)
         est=alpi*alpj*rij2/(alpi+alpj)
         if(est.gt.intcut) cycle
         !--------------- compute gradient ----------
         ! now compute integrals  for different components of i(e.g., px,py,pz)
         do mli=1,naoi
            iprim=ip+primcount(icao+mli)
            ! coefficients NOT the same (contain CAO2SAO lin. comb. coefficients)
            ci=cont(iprim)
            do mlj=1,naoj
               jprim=jp+primcount(jcao+mlj)
               cc=cont(jprim)*ci
               saw=0;sawg=0
               call build_dsdq_ints(ri,rj,rj,alpi,alpj,iptyp+mli,jptyp+mlj,saw,sawg)
               sdq(:,mlj,mli) = sdq(:,mlj,mli) + saw*cc
               sdqg(:,:10,mlj,mli) = sdqg(:,:10,mlj,mli) + sawg*cc
            enddo ! mlj : Cartesian component of j prims
         enddo  ! mli : Cartesian component of i prims
      enddo ! jp : loop over j prims
   enddo  ! ip : loop over i prims
   do mli=1,naoi
      do mlj=1,naoj
         call shiftintg(sdqg(:,:,mlj,mli),sdq(:,mlj,mli),rij)
      enddo
   enddo
end subroutine get_grad_multiint

! dsint: computes the gradient of the dipole/qpole integral contribution 
! thr         : integral cutoff according to prefactor from Gaussian product theorem 
! nat         : # of atoms
! nao         : # of spherical AOs (SAOs)
! nbf         : # of Cartesian AOs (CAOs)
! at(nat)     : atomic numbers of atoms
! xyz(3,nat)  : cartesian coordinates
! caoshell    : map shell of atom to index in CAO space (lowest Cart. component is taken), dimension: (5,nat)
! saoshell     : map shell of atom to index in SAO space (lowest m_l component is taken), dimension: (5,nat)
! primcount   : index of first primitive (over entire system) of given CAO, dimension: nbf
subroutine dsint(thr,nat,nao,nbf,at,xyz,sqrab,caoshell,saoshell,nprim,primcount, &
      &          alp,cont,H,g)
   use mctc_constants, only : pi
   use aoparam
   use intgrad
   implicit none
   real(wp),intent(in)    :: thr
   integer, intent(in)    :: nat
   integer, intent(in)    :: nao
   integer, intent(in)    :: nbf
   integer, intent(in)    :: at(nat)
   real(wp),intent(in)    :: xyz(3,nat)
   real(wp),intent(in)    :: sqrab(nat*(nat+1)/2)
   integer, intent(in)    :: caoshell(:,:)
   integer, intent(in)    :: saoshell(:,:)
   integer, intent(in)    :: nprim(:)
   integer, intent(in)    :: primcount(:)
   real(wp),intent(in)    :: alp(:)
   real(wp),intent(in)    :: cont(:)
   real(wp),intent(in)    :: H(nao,nao)
   real(wp),intent(inout) :: g(3,nat)

   ! OMP
   common /proc/ nproc
   integer nproc
   integer OMP_GET_THREAD_NUM

   integer itt(0:3)
   parameter (itt=(/0,1,4,10/))
   real(wp) tmp1,tmp2,tmp3,step,step2,step3,s00r,s00l,s00,alpj
   real(wp) skj,r1,r2,tt,t1,t2,t3,t4
   real(wp) thr2,f,ci,cc,cj,alpi,rab2,ab,est
   real(wp) f1,f2,tmp(6,6),r(3)
   real(wp) stmp,ral(3,3),rar(3,3),rbl(3,3),pre
   real(wp) dtmp,qtmp,rbr(3,3),r2l(3),r2r(3),qqa(6,6,6,3)
   real(wp)  ss(6,6,3),ddc(3,6,6,3),qqc(6,6,6,3),dda(3,6,6,3)
   integer i,j,k,l,m,ii,jj,ll,mm,kk,ki,kj,kl,km,mi,mj,ij,lin,jshmax
   integer ip,jp,iat,jat,iatyp,jatyp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)

   !                                                      call timing(t1,t3)
   !$omp parallel private (iat,jat,ixyz,iatyp,cc,ci,rab2,jatyp,ish,ishtyp, &
   !$omp&               icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&               sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp, &
   !$omp&               mli,mlj,dum,dumdum,tmp,stmp, &
   !$omp&               iao,jao,ii,jj,k ) reduction (+:g )
   !$omp do schedule(dynamic)
   do iat=1,nat
      iatyp=at(iat)
      do jat=1,iat-1
         !           if (jat.eq.iat) cycle
         jatyp=at(jat)
         r=xyz(:,iat)-xyz(:,jat)
         rab2= sum( r**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if (rab2.gt.2000) cycle
         do ish=1,ao_n(iatyp)
            ishtyp=ao_l(ish,iatyp)
            icao=caoshell(ish,iat)
            naoi=llao(ishtyp)
            iptyp=itt(ishtyp)
            jshmax=ao_n(jatyp)
            !              if(iat.eq.jat) jshmax=ish
            do jsh=1,jshmax ! jshells
               jshtyp=ao_l(jsh,jatyp)
               jcao=caoshell(jsh,jat)
               naoj=llao(jshtyp)
               jptyp=itt(jshtyp)
               !                 we go through the primitives
               !                 (because the screening is the same for all of them)
               sdqg=0;sdq=0
               do ip=1,nprim(icao+1)
                  iprim=ip+primcount(icao+1)
                  !                    exponent the same for each l component
                  alpi=alp(iprim)
                  do jp=1,nprim(jcao+1)
                     jprim=jp+primcount(jcao+1)
                     !                       exponent the same for each l component
                     alpj=alp(jprim)
                     est=alpi*alpj*rab2/(alpi+alpj)
                     if(est.gt.thr) cycle
                     !--------------- compute gradient ----------
                     ! now compute integrals  for different components of i(e.g., px,py,pz)
                     do mli=1,naoi
                        iprim=ip+primcount(icao+mli)
                        !                          coefficients NOT the same
                        !                          (contain CAO2SAO lin. comb. coefficients)
                        ci=cont(iprim)
                        do mlj=1,naoj
                           jprim=jp+primcount(jcao+mlj)
                           cc=cont(jprim)*ci
                           dum=0;dumdum=0
                           call build_ds_ints(xyz(:,iat),xyz(:,jat), &
                              & alpi,alpj,iptyp+mli,jptyp+mlj, &
                              & dum,dumdum)
                           sdq(mlj,mli) = sdq(mlj,mli)+dum*cc
                           sdqg(:,mlj,mli) = sdqg(:,mlj,mli) &
                              & + dumdum(:)*cc
                        enddo ! mlj : Cartesian component of j prims
                     enddo  ! mli : Cartesian component of i prims
                  enddo ! jp : loop over j prims
               enddo  ! ip : loop over i prims
               do ixyz=1,3
                  !                    transform from CAO to SAO, only transform the gradient
                  tmp(1:6,1:6)=sdqg(ixyz,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  sdqg(ixyz,1:6,1:6)=tmp(1:6,1:6)
               enddo
               do ii=1,llao2(ishtyp)
                  iao=ii+saoshell(ish,iat)
                  do jj=1,llao2(jshtyp)
                     jao=jj+saoshell(jsh,jat)
                     do ixyz=1,3
                        stmp=sdqg(ixyz,jj,ii)*2*H(jao,iao)
                        g(ixyz,iat)=g(ixyz,iat)+stmp
                        g(ixyz,jat)=g(ixyz,jat)-stmp
                     enddo ! ixyz
                  enddo
               enddo
            enddo ! jsh : loop over shells on jat
         enddo  ! ish : loop over shells on iat
      enddo ! jat
   enddo  ! iat
   !$omp end do
   !$omp end parallel
   !                                                      call timing(t2,t4)
   !                                     call prtime(6,t2-t1,t4-t3,'dqint5')
end subroutine dsint

end module aespot
