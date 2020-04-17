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

module xtb_aespot
   use xtb_mctc_accuracy, only : wp
   use xtb_intpack, only : olap,divpt,rhftce,prod,opab1,opab4,propa
   use xtb_xtb_data
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
   thr2 = thr*1.0d-2 ! we compare squared int-elements
   !      thr2 = 1.0d-20 ! conservative, keep all terms

   md = 0
   mq = 0
   ! set uo matrix lists
   ij = 0
   do i = 1,nao
      do j = 1,i
         ij = ij+1
         tmp1 = 0.0_wp
         tmp2 = 0.0_wp
         kk = 0
         do k = 1,3
            tmp1 = tmp1+dpint(k,i,j)*dpint(k,i,j)
            tmp2 = tmp2-qpint(k,i,j)*qpint(k,i,j)
         enddo
         do k = 1,6
            tmp2 = tmp2+2.0_wp*qpint(k,i,j)*qpint(k,i,j)
         enddo
         if(tmp1.gt.thr2)then
            md = md+1
            matdlst(1,md) = int(i,2)
            matdlst(2,md) = int(j,2)
         endif
         if(tmp2.gt.thr2)then
            mq = mq+1
            matqlst(1,mq) = int(i,2)
            matqlst(2,mq) = int(j,2)
         endif
      enddo
   enddo
   ndp = md
   nqp = mq
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
   do i = 1,nat
      qp(1:6,i) = qp(1:6,i)*(3./3.)
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
   do i = 1,nat
      qp(1:6,i) = qp(1:6,i)*(3./3.)
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
   dipm = 0.0_wp
   qp = 0.0_wp
   ij = 0
   do i = 1,nao
      ii = aoat2(i)
      ra(1:3) = xyz(1:3,ii)
      do j = 1,i-1
         ij = ij+1
         jj = aoat2(j)
         pij = p(j,i)
         ps = pij*s(j,i)
         kl = 0
         !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
         !  when doing the Mulliken population, we switch to lin-compatible sorting
         !  i,e. xx,xy,yy,xz,yz,zz
         do k = 1,3
            xk1 = ra(k)
            xk2 = xyz(k,jj)
            pdmk = pij*dpint(k,i,j)
            dipm(k,jj) = dipm(k,jj)+xk2*ps-pdmk
            dipm(k,ii) = dipm(k,ii)+xk1*ps-pdmk
            ! off-diagonal
            do l = 1,k-1
               kl = kl+1
               kj = k+l+1
               xl1 = ra(l)
               xl2 = xyz(l,jj)
               pdml = pij*dpint(l,i,j)
               pqm = pij*qpint(kj,i,j)
               tii = pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
               tjj = pdmk*xl2+pdml*xk2-xl2*xk2*ps-pqm
               qp(kl,jj) = qp(kl,jj)+tjj
               qp(kl,ii) = qp(kl,ii)+tii
            enddo
            ! diagonal
            kl = kl+1
            pqm = pij*qpint(k,i,j)
            tii = 2.0_wp*pdmk*xk1-xk1*xk1*ps-pqm
            tjj = 2.0_wp*pdmk*xk2-xk2*xk2*ps-pqm
            qp(kl,jj) = qp(kl,jj)+tjj
            qp(kl,ii) = qp(kl,ii)+tii
         enddo
      enddo
      ij = ij+1
      pij = p(i,i)
      ps = pij*s(i,i)
      kl = 0
      !  the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
      !  when doing the Mulliken population, we switch to lin-compatible sorting
      !  i,e. xx,xy,yy,xz,yz,zz
      do k = 1,3
         xk1 = ra(k)
         pdmk = pij*dpint(k,i,i)
         dipm(k,ii) = dipm(k,ii)+xk1*ps-pdmk
         ! off-diagonal
         do l = 1,k-1
            kl = kl+1
            kj = k+l+1 ! the qpint is stored as xx,yy,zz,xy,xz,yz (from integral routine)
            xl1 = ra(l)
            pdml = pij*dpint(l,i,i)
            pqm = pij*qpint(kj,i,i)
            tii = pdmk*xl1+pdml*xk1-xl1*xk1*ps-pqm
            qp(kl,ii) = qp(kl,ii)+tii
         enddo
         !diagonal
         kl = kl+1
         pqm = pij*qpint(k,i,i)
         tii = 2.0_wp*pdmk*xk1-xk1*xk1*ps-pqm
         qp(kl,ii) = qp(kl,ii)+tii
      enddo
   enddo
   ! remove trace
   do i = 1,nat
      tii = qp(1,i)+qp(3,i)+qp(6,i)
      tii = 0.50_wp*tii
      qp(1:6,i) = 1.50_wp*qp(1:6,i)
      qp(1,i) = qp(1,i)-tii
      qp(3,i) = qp(3,i)-tii
      qp(6,i) = qp(6,i)-tii
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
subroutine aniso_electro(aesData,nat,at,xyz,q,dipm,qp,gab3,gab5,e,epol)
   use xtb_lin, only : lin
   implicit none
   class(TMultipoleData), intent(in) :: aesData
   integer nat,at(nat)
   real(wp) xyz(3,nat),q(nat)
   real(wp) e
   real(wp) qp1(6),rr(3),dp1(3),rij(3)
   real(wp) edd,e01,e02,e11,r2,tt,tt3,q1,qs2
   real(wp) ed,eq,epol
   ! stuff for potential
   real(wp) gab3(nat*(nat+1)/2),gab5(nat*(nat+1)/2)
   real(wp) dipm(3,nat),qp(6,nat)
   real(wp), parameter :: bohr = 1.0_wp/0.52917726_wp

   integer i,j,k,l,m,ki,kj,kl
   e = 0.0_wp
   epol = 0.0_wp
   e01 = 0.0_wp
   e02 = 0.0_wp
   e11 = 0.0_wp
   do i = 1,nat
      q1 = q(i)
      rr(1:3) = xyz(1:3,i)
      dp1(1:3) = dipm(1:3,i)
      qp1(1:6) = qp(1:6,i)
      ! test: semilocal CT correction
      ! dipole
      tt = dp1(1)*dp1(1)+dp1(2)*dp1(2)+dp1(3)*dp1(3)
      ! qpole
      tt3 = 0.0_wp
      do k = 1,3
         do l = 1,3
            kl = lin(l,k)
            tt3 = tt3+qp1(kl)*qp1(kl)
         enddo
      enddo
      epol = epol+aesData%dipKernel(at(i))*tt+tt3*aesData%quadKernel(at(i))
      ! ---
      do j = 1,i-1             ! loop over all atoms
         kj = lin(j,i)
         rij(1:3) = xyz(1:3,j)-rr(1:3)
         r2 = sum(rij*rij)
         ed = 0.0_wp
         eq = 0.0_wp
         edd = 0.0_wp
         !           dipole - charge
         do k = 1,3
            ed = ed+q(j)*dp1(k)*rij(k)
            ed = ed-dipm(k,j)*q1*rij(k)
            !              dip-dip & charge-qpole
            do l = 1,3
               kl = lin(l,k)
               tt = rij(l)*rij(k)
               tt3 = 3.0_wp*tt
               eq = eq+q(j)*qp1(kl)*tt
               eq = eq+qp(kl,j)*q1*tt
               edd = edd-dipm(k,j)*dp1(l)*tt3
            enddo
            !              diagonal dip-dip term
            edd = edd+dipm(k,j)*dp1(k)*r2
         enddo
         e01 = e01+ed*gab3(kj)
         e02 = e02+eq*gab5(kj)
         e11 = e11+edd*gab5(kj)
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
   use xtb_lin, only : lin
   implicit none
   integer, intent(in) :: nat,nao,aoat2(nao)
   real(wp), intent(in) :: dpint(3,nao,nao),s(nao,nao)
   real(wp), intent(in) :: qpint(6,nao,nao),p(nao,nao)
   real(wp), intent(in) :: vs(nat),vd(3,nat),vq(6,nat)
   real(wp), intent(out) :: e
   real(wp) eaes,pji,fji
   integer i,j,k,l,ii,jj,ij,kl,kj
   ! CAMM
   eaes = 0.0_wp
   ij = 0
   do i = 1,nao
      ii = aoat2(i)
      do j = 1,nao
         ij = lin(j,i)
         jj = aoat2(j)
         fji = 0.0_wp
         pji = p(j,i)
         fji = fji+s(j,i)*(vs(ii)+vs(jj))
         do k = 1,3
            fji = fji+dpint(k,i,j)*(vd(k,ii)+vd(k,jj))
         enddo
         do k = 1,6
            fji = fji+qpint(k,i,j)*(vq(k,ii)+vq(k,jj))
         enddo
         eaes = eaes+pji*fji
      enddo
   enddo
   eaes = 0.250_wp*eaes
   !      write(*,*) 'EAES',eaes
   e = eaes
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
subroutine setvsdq(aesData,nat,at,xyz,q,dipm,qp,gab3,gab5,vs,vd,vq)
   use xtb_lin, only : lin
   implicit none
   class(TMultipoleData), intent(in) :: aesData
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
   mm = 0
   vs = 0.0_wp
   vd = 0.0_wp
   vq = 0.0_wp
   ! set up overlap proportional potential
   do i = 1,nat
      ra(1:3) = xyz(1:3,i)
      stmp = 0.0_wp
      dtmp = 0.0_wp
      qtmp = 0.0_wp
      do j = 1,nat
         mm = lin(j,i)
         g3 = gab3(mm)
         g5 = gab5(mm)
         rb(1:3) = xyz(1:3,j)
         dra(1:3) = ra(1:3)-rb(1:3)
         dum3a = 0.0_wp ! collect gab3 dependent terms
         dum5a = 0.0_wp ! collect gab5 dependent terms
         r2a = 0.0_wp
         r2ab = 0.0_wp
         t1a = 0.0_wp
         t2a = 0.0_wp
         t3a = 0.0_wp
         t4a = 0.0_wp
         ll = 0
         do l1 = 1,3
            ! potential from dipoles
            r2a = r2a+ra(l1)*ra(l1)      ! R_C * R_C
            r2ab = r2ab+dra(l1)*dra(l1)  ! R_AC * R_AC
            t1a = t1a+ra(l1)*dra(l1)     ! R_C * R_AC  : for dip-q (q-shift) and dip-dip (q-shift)
            t2a = t2a+dipm(l1,j)*dra(l1) ! mu_A * R_AC : for q-dip and dip-dip (q-shift)
            t3a = t3a+ra(l1)*dipm(l1,j)  ! R_C * mu_A  : for diag. dip-dip (q-shift)
            t4a = t4a+dra(l1)*dra(l1)*ra(l1)*ra(l1) ! (R_C o R_AC)**"2(square of Hadamard product) :
            ! results from trace remove from q-pole (q-shift)
            do l2 = 1,3
               ll = lin(l1,l2)
               ! potential from quadrupoles
               dum5a = dum5a-qp(ll,j)*dra(l1)*dra(l2) &
                  & -1.50_wp*q(j)*dra(l1)*dra(l2)*ra(l1)*ra(l2)
               if(l2.ge.l1) cycle
               ki = l1+l2+1
               qtmp(ki) = qtmp(ki)-3.0_wp*q(j)*g5*dra(l2)*dra(l1)
            enddo
            qtmp(l1) = qtmp(l1)-1.50_wp*q(j)*g5*dra(l1)*dra(l1)
         enddo
         !
         ! set up S-dependent potential
         dum3a = -t1a*q(j)-t2a ! dip-q (q-shift) and q-dip
         dum5a = dum5a+t3a*r2ab-3.0_wp*t1a*t2a & !dip-dip (q-shift terms)
            & +0.50_wp*q(j)*r2a*r2ab !qpole-q (q-shift, trace removal)
         stmp = stmp+dum5a*g5+dum3a*g3
         do l1 = 1,3
            dum3a = dra(l1)*q(j)
            dum5a = 3.0_wp*dra(l1)*t2a &           ! dipint-dip
               & -r2ab*dipm(l1,j) &            ! dipint-dip (diagonal)
               & -q(j)*r2ab*ra(l1) &           ! qpole-q (dipint-shift, trace removal)
               & +3.0_wp*q(j)*dra(l1)*t1a   ! qpole-q (dipint-shift)
            dtmp(l1) = dtmp(l1)+dum3a*g3+dum5a*g5
            qtmp(l1) = qtmp(l1)+0.50_wp*r2ab*q(j)*g5 !remove trace term
         enddo
      enddo
      vs(i) = stmp                       ! q terms
      vd(1:3,i) = dtmp(1:3)              ! dipints from atom i
      vq(1:6,i) = qtmp(1:6)              ! qpints from atom i
      ! --- CT correction terms
      qs1 = aesData%dipKernel(at(i))*2.0_wp
      qs2 = aesData%quadKernel(at(i))*6.0_wp ! qpole pot prefactors
      t3a = 0.0_wp
      t2a = 0.0_wp
      do l1 = 1,3
         ! potential from dipoles
         t3a = t3a+ra(l1)*dipm(l1,i)*qs1  ! R_C * mu_C  : for diag. dip-dip
         vd(l1,i) = vd(l1,i)-qs1*dipm(l1,i)
         do l2 = 1,l1-1
            ! potential from quadrupoles
            ll = lin(l1,l2)
            ki = l1+l2+1
            vq(ki,i) = vq(ki,i)-qp(ll,i)*qs2
            t3a = t3a-ra(l1)*ra(l2)*qp(ll,i)*qs2
            vd(l1,i) = vd(l1,i)+ra(l2)*qp(ll,i)*qs2
            vd(l2,i) = vd(l2,i)+ra(l1)*qp(ll,i)*qs2
         enddo
         ! diagonal
         ll = lin(l1,l1)
         vq(l1,i) = vq(l1,i)-qp(ll,i)*qs2*0.50_wp
         t3a = t3a-ra(l1)*ra(l1)*qp(ll,i)*qs2*0.50_wp
         vd(l1,i) = vd(l1,i)+ra(l1)*qp(ll,i)*qs2
         ! collect trace removal terms
         t2a = t2a+qp(ll,i)
      enddo
      vs(i) = vs(i)+t3a
      ! trace removal
      t2a = t2a*aesData%quadKernel(at(i))
      do l1 = 1,3
         vq(l1,i) = vq(l1,i)+t2a
         vd(l1,i) = vd(l1,i)-2.0_wp*ra(l1)*t2a
         vs(i) = vs(i)+t2a*ra(l1)*ra(l1)
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
subroutine setdvsdq(aesData,nat,at,xyz,q,dipm,qp,gab3,gab5,vs,vd,vq)
   use xtb_lin, only : lin
   implicit none
   class(TMultipoleData), intent(in) :: aesData
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
   mm = 0
   vs = 0.0_wp
   vd = 0.0_wp
   vq = 0.0_wp
   ! set up overlap proportional potential
   do i = 1,nat
      ra(1:3) = xyz(1:3,i)
      stmp = 0.0_wp
      dtmp = 0.0_wp
      qtmp = 0.0_wp
      do j = 1,nat
         mm = lin(j,i)
         g3 = gab3(mm)
         g5 = gab5(mm)
         rb(1:3) = xyz(1:3,j)
         dra(1:3) = ra(1:3)-rb(1:3)
         dum3a = 0.0_wp ! collect gab3 dependent terms
         dum5a = 0.0_wp ! collect gab5 dependent terms
         r2a = 0.0_wp
         r2ab = 0.0_wp
         t2a = 0.0_wp
         ll = 0
         do l1 = 1,3
            ! potential from dipoles
            r2ab = r2ab+dra(l1)*dra(l1)  ! R_AC * R_AC
            t2a = t2a+dipm(l1,j)*dra(l1) ! mu_A * R_AC : for q-dip and dip-dip (q-shift)
            do l2 = 1,3
               ll = lin(l1,l2)
               ! potential from quadrupoles
               dum5a = dum5a-qp(ll,j)*dra(l1)*dra(l2)
               if(l2.ge.l1) cycle
               ki = l1+l2+1
               qtmp(ki) = qtmp(ki)-3.0_wp*q(j)*g5*dra(l2)*dra(l1)
            enddo
            qtmp(l1) = qtmp(l1)-1.50_wp*q(j)*g5*dra(l1)*dra(l1)
         enddo
         dum3a = -t2a ! q-dip ! w/o shift terms
         stmp = stmp+dum3a*g3+dum5a*g5
         do l1 = 1,3
            dum3a = dra(l1)*q(j) ! w/o shift terms
            dum5a = 3.0_wp*dra(l1)*t2a & ! dipint-dip
               & -r2ab*dipm(l1,j)  ! dipint-dip (diagonal)
            dtmp(l1) = dtmp(l1)+dum3a*g3+dum5a*g5
            qtmp(l1) = qtmp(l1)+0.50_wp*r2ab*q(j)*g5 !remove trace term
         enddo
      enddo
      vs(i) = stmp
      vd(1:3,i) = dtmp(1:3)
      vq(1:6,i) = qtmp(1:6)
      ! --- CT correction terms
      qs1 = aesData%dipKernel(at(i))*2.0_wp
      qs2 = aesData%quadKernel(at(i))*6.0_wp ! qpole pot prefactors
      t2a = 0.0_wp
      do l1 = 1,3
         ! potential from dipoles
         vd(l1,i) = vd(l1,i)-qs1*dipm(l1,i)
         do l2 = 1,l1-1
            ! potential from quadrupoles
            ll = lin(l1,l2)
            ki = l1+l2+1
            vq(ki,i) = vq(ki,i)-qp(ll,i)*qs2
         enddo
         ll = lin(l1,l1)
         vq(l1,i) = vq(l1,i)-qp(ll,i)*qs2*0.50_wp
         ! collect trace removal terms
         t2a = t2a+qp(ll,i)
      enddo
      ! trace removal
      t2a = t2a*aesData%quadKernel(at(i))
      do l1 = 1,3
         vq(l1,i) = vq(l1,i)+t2a
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
   use xtb_mctc_convert
   use xtb_lin, only : lin
   implicit none
   integer, intent(in) :: iunit
   integer, intent(in) :: n
   real(wp), intent(in)  :: xyz(3,n),q(n),dipm(3,n),qp(6,n)
   real(wp), intent(out) :: dip,d3(3)
   real(wp) rr1(3),rr2(3),tma(6),tmb(6),tmc(6),dum
   integer i,j,k,l
   rr1 = 0.0_wp
   rr2 = 0.0_wp
   write(iunit,'(a)')
   do i = 1,n
      do j = 1,3
         rr1(j) = rr1(j)+q(i)*xyz(j,i)
         rr2(j) = rr2(j)+dipm(j,i)
      enddo
   enddo
   d3(1:3) = rr1(1:3)+rr2(1:3)
   dip = sqrt(d3(1)**2+d3(2)**2+d3(3)**2)
   write(iunit,'(a)',advance='yes')'molecular dipole:'
   write(iunit,'(a)',advance='no')'                 '
   write(iunit,'(a)',advance='yes') &
      & 'x           y           z       tot (Debye)'
   write(iunit,'(a,3f12.3)') ' q only: ',rr1(1:3)
   write(iunit,'(a,4f12.3)') '   full: ',d3(1:3),dip*autod

   tma = 0.0_wp
   tmb = 0.0_wp
   tmc = 0.0_wp
   do i = 1,n
      l = 0
      do j = 1,3
         do k = 1,j
            l = lin(k,j)
            tma(l) = tma(l)+xyz(j,i)*xyz(k,i)*q(i)
            tmb(l) = tmb(l)+dipm(k,i)*xyz(j,i)+dipm(j,i)*xyz(k,i)
            tmc(l) = tmc(l)+qp(l,i)
         enddo
      enddo
   enddo
   ! remove traces and multiply with 3/2 in q and dip parts
   dum = tma(1)+tma(3)+tma(6)
   dum = 0.50_wp*dum
   tma = 1.50_wp*tma
   l = 0
   do j = 1,3
      l = l+j
      tma(l) = tma(l)-dum
   enddo
   dum = tmb(1)+tmb(3)+tmb(6)
   dum = 0.50_wp*dum
   tmb = 1.50_wp*tmb
   l = 0
   do j = 1,3
      l = l+j
      tmb(l) = tmb(l)-dum
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
   use xtb_lin, only : lin
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
   do i = 1,nat
      q1 = q(i)
      rr(1:3) = xyz(1:3,i)
      dip(1:3) = dipm(1:3,i)
      qp1(1:6) = qp(1:6,i)
      tmp2 = 0.0_wp ! cumulate terms propto CN gradient -  to scale only quadratically
      do j = 1,nat            ! loop over other atoms
         if(i.eq.j) cycle
         kj = lin(j,i)
         rij(1:3) = xyz(1:3,j)-rr(1:3)
         r2 = sum(rij*rij)
         rabi = 1.0_wp/sqrt(r2)
         !           call dzero(2.0_wp,rabi,at(i),at(j),damp,ddamp)
         call dzero(kdmp3,rabi,radcn(i),radcn(j),damp1,ddamp)
         dgab3 = dgab(3.0_wp,rabi,damp1,ddamp)
         !           call dzero(3.0_wp,rabi,at(i),at(j),damp,ddamp)
         call dzero(kdmp5,rabi,radcn(i),radcn(j),damp2,ddamp)
         dgab5 = dgab(5.0_wp,rabi,damp2,ddamp)
         !!!         DEBUG
         !            dgab3 = 0.0_wp
         !            dgab5 = 0.0_wp
         !            dgab3 = dgab3*100.0_wp
         !            dgab5 = dgab5*100.0_wp
         !!!
         ed = 0.0_wp
         edd = 0.0_wp
         eq = 0.0_wp
         !           dipole - charge
         do k = 1,3
            ed = ed+q(j)*dip(k)*rij(k)
            ed = ed-dipm(k,j)*q1*rij(k)
            tt = q1*dipm(k,j)-q(j)*dip(k)
            ! part of dip-q derivative
            g(k,i) = g(k,i)+gab3(kj)*tt
            !              dip-dip & charge-qpole
            ddm2 = 0.0_wp
            ddm3a = 0.0_wp
            ddm3b = 0.0_wp
            qqa = 0.0_wp
            qqb = 0.0_wp
            do l = 1,3
               kl = lin(l,k)
               tt = rij(l)*rij(k)
               tt3 = 3.0_wp*tt
               eq = eq+q(j)*qp1(kl)*tt
               eq = eq+qp(kl,j)*q1*tt
               edd = edd-dipm(k,j)*dip(l)*tt3
               ! extra d-d terms
               ddm2 = ddm2+dipm(l,j)*dip(l)
               ddm3a = ddm3a+dip(l)*rij(l)
               ddm3b = ddm3b+dipm(l,j)*rij(l)
               ! extra q-qpole terms
               qqa = qqa+rij(l)*qp(kl,j)
               qqb = qqb+rij(l)*qp1(kl)
            enddo
            edd = edd+dipm(k,j)*dip(k)*r2
            g(k,i) = g(k,i)-2.0_wp*gab5(kj)*ddm2*rij(k)
            g(k,i) = g(k,i)+3.0_wp*gab5(kj)*ddm3a*dipm(k,j)
            g(k,i) = g(k,i)+3.0_wp*gab5(kj)*ddm3b*dip(k)
            g(k,i) = g(k,i)-2.0_wp*gab5(kj)*qqa*q1
            g(k,i) = g(k,i)-2.0_wp*gab5(kj)*qqb*q(j)
         enddo
         do k = 1,3
            dxi = rij(k)*rabi
            g(k,i) = g(k,i)-ed*dgab3*dxi
            g(k,i) = g(k,i)-(eq+edd)*dgab5*dxi
         enddo
         ! collect terms for CN-dependent part
         rab = 0.50_wp*(radcn(i)+radcn(j))
         tmp3 = ed*kdmp3*gab3(kj)*(damp1/rab)*(rab*rabi)**kdmp3
         tmp2 = tmp2+tmp3
         tmp3 = (eq+edd)*kdmp5*gab5(kj)*(damp2/rab)*(rab*rabi)**kdmp5
         tmp2 = tmp2+tmp3
      enddo
      ! CN-dependent part  - O(N^2)
      tmp2 = 3.0_wp*tmp2
      g(:,:) = g-tmp2*dcn(:,:,i)

   enddo
end subroutine aniso_grad


! check and print sparsity w.r.t. individual contribution
! to get an idea
subroutine checkspars(nao,ndp,nqp,nmat,matlist,mqlst,mdlst)
   use xtb_lin, only : lin
   implicit none
   integer,intent(in) :: ndp,nqp,nmat,nao
   integer,intent(in) :: matlist(2,nao*(nao+1)/2)
   integer,intent(in) :: mqlst(2,nqp),mdlst(2,ndp)
   integer  :: i,j,m,k,ntot,mi,mj,ki,kj,mm,kk
   logical,allocatable ::  nzero(:)
   ! check overall sparsity
   allocate(nzero(nao*(nao+1)/2))
   nzero = .false.
   do k = 1,ndp
      ki = mdlst(1,k)
      kj = mdlst(2,k)
      kk = lin(ki,kj)
      nzero(kk) = .true.
   enddo
   do k = 1,nqp
      ki = mqlst(1,k)
      kj = mqlst(2,k)
      kk = lin(ki,kj)
      nzero(kk) = .true.
   enddo
   do k = 1,nmat
      ki = matlist(1,k)
      kj = matlist(2,k)
      kk = lin(ki,kj)
      nzero(kk) = .true.
   enddo
   mm = nao*(nao+1)/2
   ntot = 0
   do i = 1,mm
      if(nzero(i)) ntot = ntot+1
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
   gab3 = 0.0_wp ! for r**-2 decaying q-dip term
   gab5 = 0.0_wp ! for r**-3 decaying terms (q-qpol,dip-dip)
   l = 1
   gab3(1) = 0.0_wp
   gab5(1) = 0.0_wp
   do i = 2,nat
      rr(1:3) = xyz(1:3,i)
      do j = 1,i-1
         l = l+1
         tmp2 = 0.0_wp
         do k = 1,3
            tmp1 = xyz(k,j)-rr(k)
            tmp2 = tmp2+tmp1**2
         enddo
         tmp1 = 1.0_wp/sqrt(tmp2)
         !           call dzero(2.0_wp,tmp1,at(i),at(j),damp,ddamp)
         call dzero(kdmp3,tmp1,radcn(i),radcn(j),damp,ddamp)
         gab3(l) = gab(3.0_wp,tmp1,damp)
         !           call dzero(3.0_wp,tmp1,at(i),at(j),damp,ddamp)
         call dzero(kdmp5,tmp1,radcn(i),radcn(j),damp,ddamp)
         gab5(l) = gab(5.0_wp,tmp1,damp)
      enddo
      l = l+1
      gab3(l) = 0.0_wp
      gab5(l) = 0.0_wp
   enddo

   !!!!!!! DEBUG
   !      gab3 = 0.0_wp
   !      gab5 = 0.0_wp
   !      gab3 = 100.0_wp*gab3
   !      gab5 = 100.0_wp*gab5
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
subroutine get_radcn(aesData,n,at,cn,shift,expo,rmax,radcn)
   implicit none
   class(TMultipoleData), intent(in) :: aesData
   integer, intent (in) :: n,at(n)
   real(wp), intent (in)  :: cn(n),shift,expo,rmax
   real(wp), intent (out) :: radcn(n)
   real(wp) rco,t1,t2
   integer i,j
   do i = 1,n
      rco = aesData%multiRad(at(i))             ! base radius of element
      t1 =cn(i)-aesData%valenceCN(at(i))-shift  ! CN - VALCN - SHIFT
      t2 =rco +(rmax-rco)/(1.0_wp+exp(-expo*t1))
      radcn(i) = t2
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
! dcn   : on input  : derivatives of CN(j) w.r.t. Cart. directions of i
!       : on output : derivatives of RADCN(j) w.r.t. Cart. directions of i
subroutine dradcn(aesData,n,at,cn,shift,expo,rmax,dcn)
   implicit none
   class(TMultipoleData), intent(in) :: aesData
   integer, intent (in) :: n,at(n)
   real(wp), intent (in)  :: cn(n),shift,expo,rmax
   real(wp), intent (inout) :: dcn(3,n,n)
   real(wp) rco,t1,t2,t3,t4,tmp1,tmp2
   integer i,j,k
   do i = 1,n
      rco = aesData%multiRad(at(i))             ! base radius of element
      t1 =exp(-expo*(cn(i)-aesData%valenceCN(at(i))-shift))  ! CN - VALCN - SHIFT
      t2 =(rmax-rco)/(1.0_wp+2.0_wp*t1+t1*t1)
      t2 = t2*expo*t1
      dcn(:,:,i) = dcn(:,:,i)*t2
   end do
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
   !     rco = 2.0/(1./aesData%multiRad(ati)+1./aesData%multiRad(atj))  ! unstable
   !     rco = sqrt(aesData%multiRad(ati)*aesData%multiRad(atj))        ! unstable
   rco = 0.5*(radi+radj)
   ! zero-damping function and gradient w.r.t. Rab
   damp = 1.0_wp/(1.0_wp+6.0_wp*(rco*rabinv)**dex)
   ddamp = -dex*rabinv*(damp*damp-damp)
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
   gab = damp*(rabinv**dex) ! LR-decay * damping
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
   tmp1 = -dex*rabinv*(rabinv**dex) ! LR-decay derivative
   dgab = tmp1*damp+ddamp*rabinv**dex
end function dgab

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine gfn2broyden_diff(n,istart,nbr,dipm,qp,q_in,dq)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: dipm(3,n),qp(6,n),q_in(nbr)
   real(wp), intent (inout) :: dq(nbr)
   integer i,j,k
   k = istart
   do i = 1,n
      do j = 1,3
         k = k+1
         dq(k) = dipm(j,i)-q_in(k)
      enddo
      do j = 1,6
         k = k+1
         dq(k) = qp(j,i)-q_in(k)
      enddo
   enddo
   istart = k

end subroutine gfn2broyden_diff

subroutine gfn2broyden_save(n,istart,nbr,dipm,qp,q_in)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: dipm(3,n),qp(6,n)
   real(wp), intent (inout) :: q_in(nbr)
   integer i,j,k
   k = istart
   do i = 1,n
      do j = 1,3
         k = k+1
         q_in(k) = dipm(j,i)
      enddo
      do j = 1,6
         k = k+1
         q_in(k) = qp(j,i)
      enddo
   enddo
   istart = k

end subroutine gfn2broyden_save



subroutine gfn2broyden_out(n,istart,nbr,q_in,dipm,qp)
   implicit none
   integer, intent (in) :: n,nbr
   integer, intent (inout) ::  istart
   real(wp), intent (in)  :: q_in(nbr)
   real(wp), intent (out) :: dipm(3,n),qp(6,n)
   integer i,j,k
   k = istart
   do i = 1,n
      do j = 1,3
         k = k+1
         dipm(j,i) = q_in(k)
      enddo
      do j = 1,6
         k = k+1
         qp(j,i) = q_in(k)
      enddo
   enddo
   istart = k
end subroutine gfn2broyden_out


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
subroutine dsint(nShell,hData,thr,nat,nao,nbf,at,xyz,sqrab,caoshell,saoshell,nprim,primcount, &
      &          alp,cont,H,g)
   use xtb_mctc_constants, only : pi
   use xtb_intgrad
   implicit none
   integer, intent(in) :: nShell(:)
   type(THamiltonianData), intent(in) :: hData
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
   integer ip,jp,iat,jat,izp,jzp,ish,jsh,icao,jcao,iao,jao,ixyz
   integer ishtyp,jshtyp,iptyp,jptyp,naoi,naoj,mli,mlj,iprim,jprim
   real(wp) :: dumdum(3),dum,sdq(6,6),sdqg(3,6,6)

   !                                                      call timing(t1,t3)
   !$omp parallel private (iat,jat,ixyz,izp,cc,ci,rab2,jzp,ish,ishtyp, &
   !$omp&               icao,naoi,iptyp,jsh,jshmax,jshtyp,jcao,naoj,jptyp, &
   !$omp&               sdq,sdqg,est,alpi,alpj,ab,iprim,jprim,ip,jp, &
   !$omp&               mli,mlj,dum,dumdum,tmp,stmp, &
   !$omp&               iao,jao,ii,jj,k ) reduction (+:g )
   !$omp do schedule(dynamic)
   do iat = 1,nat
      izp = at(iat)
      do jat = 1,iat-1
         !           if (jat.eq.iat) cycle
         jzp = at(jat)
         r = xyz(:,iat)-xyz(:,jat)
         rab2 =  sum( r**2 )
         !           ints < 1.d-9 for RAB > 40 Bohr
         if (rab2.gt.2000) cycle
         do ish = 1,nShell(izp)
            ishtyp = hData%angShell(ish,izp)
            icao = caoshell(ish,iat)
            naoi = llao(ishtyp)
            iptyp = itt(ishtyp)
            jshmax = nShell(jzp)
            !              if(iat.eq.jat) jshmax = ish
            do jsh = 1,jshmax ! jshells
               jshtyp = hData%angShell(jsh,jzp)
               jcao = caoshell(jsh,jat)
               naoj = llao(jshtyp)
               jptyp = itt(jshtyp)
               !                 we go through the primitives
               !                 (because the screening is the same for all of them)
               sdqg = 0;sdq = 0
               do ip = 1,nprim(icao+1)
                  iprim = ip+primcount(icao+1)
                  !                    exponent the same for each l component
                  alpi = alp(iprim)
                  do jp = 1,nprim(jcao+1)
                     jprim = jp+primcount(jcao+1)
                     !                       exponent the same for each l component
                     alpj = alp(jprim)
                     est = alpi*alpj*rab2/(alpi+alpj)
                     if(est.gt.thr) cycle
                     !--------------- compute gradient ----------
                     ! now compute integrals  for different components of i(e.g., px,py,pz)
                     do mli = 1,naoi
                        iprim = ip+primcount(icao+mli)
                        !                          coefficients NOT the same
                        !                          (contain CAO2SAO lin. comb. coefficients)
                        ci = cont(iprim)
                        do mlj = 1,naoj
                           jprim = jp+primcount(jcao+mlj)
                           cc = cont(jprim)*ci
                           dum = 0;dumdum = 0
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
               do ixyz = 1,3
                  !                    transform from CAO to SAO, only transform the gradient
                  tmp(1:6,1:6) = sdqg(ixyz,1:6,1:6)
                  call dtrf2(tmp,ishtyp,jshtyp)
                  sdqg(ixyz,1:6,1:6) = tmp(1:6,1:6)
               enddo
               do ii = 1,llao2(ishtyp)
                  iao = ii+saoshell(ish,iat)
                  do jj = 1,llao2(jshtyp)
                     jao = jj+saoshell(jsh,jat)
                     do ixyz = 1,3
                        stmp = sdqg(ixyz,jj,ii)*2*H(jao,iao)
                        g(ixyz,iat) = g(ixyz,iat)+stmp
                        g(ixyz,jat) = g(ixyz,jat)-stmp
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

end module xtb_aespot
