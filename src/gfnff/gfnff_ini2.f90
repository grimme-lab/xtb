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
module xtb_gfnff_ini2
   use xtb_gfnff_data, only : TGFFData
   use xtb_gfnff_neighbourlist, only : TGFFNeighbourList
   use xtb_gfnff_topology, only : TGFFTopology
   use xtb_type_environment, only : TEnvironment
   implicit none
   private
   public :: gfnff_neigh, getnb, nbondmat
   public :: pairsbond, pilist, nofs, xatom, ctype, amide, amideH, alphaCO
   public :: ringsatom, ringsbond, ringsbend, ringstors, ringstorl
   public :: chktors, chkrng, hbonds, getring36, ssort, goedeckera, qheavy
   public :: gfnff_hbset, gfnff_hbset0, bond_hbset, bond_hbset0
   public :: bond_hb_AHB_set, bond_hb_AHB_set1, bond_hb_AHB_set0

contains

subroutine gfnff_neigh(env,makeneighbor,natoms,at,xyz,rab,fq,f_in,f2_in,lintr,mchar,hyb,itag,nbm,nbf,param,topo)
      use xtb_gfnff_param
      implicit none
      character(len=*), parameter :: source = 'gfnff_ini2_neigh'
      type(TEnvironment), intent(inout) :: env
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(inout) :: topo
      logical makeneighbor
      integer at(natoms),natoms
      integer hyb (natoms)
      integer itag(natoms)
      integer nbm(20,natoms)                 ! needed for ring assignment (done without metals)
      integer nbf(20,natoms)                 ! full needed for fragment assignment
      real*8  rab   (natoms*(natoms+1)/2)
      real*8  xyz(3,natoms)
      real*8  mchar(natoms)
      real*8  fq
      real*8  f_in,f2_in               ! radius scaling for atoms/metal atoms recpectively
      real*8  lintr                    ! threshold for linearity

      logical etacoord,da,strange_iat,metal_iat
      integer,allocatable :: nbdum(:,:)
      real*8 ,allocatable :: cn(:),rtmp(:)
      integer iat,i,j,k,ni,ii,jj,kk,ll,lin,ati,nb20i,nbdiff,hc_crit,nbmdiff,nnf,nni,nh,nm
      integer ai,aj,nn,im,ncm,l,no
      real*8 r,pi,a1,f,f1,phi,f2,rco,fat(86)
      data pi/3.1415926535897932384626433832795029d0/
      data fat   / 86 * 1.0d0 /

!     special hacks
      fat( 1)=1.02
      fat( 4)=1.03
      fat( 5)=1.02
      fat( 8)=1.02
      fat( 9)=1.05
      fat(10)=1.10
      fat(11)=1.01
      fat(12)=1.02
      fat(15)=0.97
      fat(18)=1.10
      fat(19)=1.02
      fat(20)=1.02
      fat(38)=1.02
      fat(34)=0.99
      fat(50)=1.01
      fat(51)=0.99
      fat(52)=0.95
      fat(53)=0.98
      fat(56)=1.02
      fat(76)=1.02
      fat(82)=1.06
      fat(83)=0.95

      allocate(cn(natoms),rtmp(natoms*(natoms+1)/2),nbdum(20,natoms))

! determine the neighbor list
      if(makeneighbor) then

        topo%nb =0  ! without highly coordinates atoms
        nbm=0  ! without any metal
        nbf=0  ! full

        do i=1,natoms
           cn(i)=dble(param%normcn(at(i)))
        enddo
        call gfnffrab(natoms,at,cn,rtmp) ! guess RAB based on "normal" CN
        do i=1,natoms
           ai=at(i)
           f1=fq
           if(param%metal(ai) > 0) f1 = f1 * 2.0d0
           do j=1,i-1
              f2=fq
              aj=at(j)
              if(param%metal(aj) > 0) f2 = f2 * 2.0d0
              k=lin(j,i)
              rco=rtmp(k)
              rtmp(k)=rtmp(k)-topo%qa(i)*f1-topo%qa(j)*f2 ! change radius of atom i and j with charge
!             element specials
              rtmp(k)=rtmp(k)*fat(ai)*fat(aj)
           enddo
        enddo

        call getnb(natoms,at,rtmp,rab,mchar,1,f_in,f2_in,nbdum,nbf,param) ! full
        call getnb(natoms,at,rtmp,rab,mchar,2,f_in,f2_in,nbf  ,topo%nb,param) ! no highly coordinates atoms
        call getnb(natoms,at,rtmp,rab,mchar,3,f_in,f2_in,nbf  ,nbm,param) ! no metals and unusually coordinated stuff

! take the input
      else

        nbf = topo%nb
        nbm = topo%nb

      endif
! done

      itag = 0 ! save special hyb info

! tag atoms in nb(19,i) if they belong to a cluster (which avoids the ring search)
      do i=1,natoms
         if(nbf(20,i).eq.0.and.param%group(at(i)).ne.8)then
            write(env%unit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
            write(env%unit,'(''  warning: no bond partners for atom'',i4)')i
            write(env%unit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
         endif
         if(at(i).lt.11.and.nbf(20,i).gt.2)then
            do k=1,nbf(20,i)
               kk=nbf(k,i)
               if(param%metal(at(kk)).ne.0.or.topo%nb(20,kk).gt.4) then
                  topo%nb (19,i)=1
                  nbf(19,i)=1
                  nbm(19,i)=1
               endif
            enddo
         endif
!        write(env%unit,*) i,(topo%nb(j,i),j=1,topo%nb(20,i))
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hybridization states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do i=1,natoms
         ati  = at(i)
!        important: determine cases where atom is pi bonded to a metal and thus
!        the hyb must be obtained from the reduced (wo metals) neighbor list
         etacoord=.false.
         if(ati.le.10)then
            if(ati.eq.6.and.nbf(20,i).ge.4.and.nbm(20,i).eq.3) etacoord=.true.  ! CP case
            if(ati.eq.6.and.nbf(20,i).eq.3.and.nbm(20,i).eq.2) etacoord=.true.  ! alkyne case
            nm=0
            do k=1,nbf(20,i)  ! how many metals ? and which
               kk=nbf(k,i)
               if(param%metal(at(kk)).ne.0) then
                  nm=nm+1
                  im=kk
               endif
            enddo
            if(nm.eq.0) then
               etacoord=.false.  ! etacoord makes no sense without metals!
            elseif(nm.eq.1)then  ! distinguish M-CR2-R i.e. not an eta coord.
               ncm=0
               do k=1,nbf(20,i)  !
                  if(nbf(k,i).ne.im)then ! all neighbors that are not the metal im
                     kk=nbf(k,i)
                     do l=1,nbf(20,kk)
                        if(nbf(l,kk).eq.im) ncm=ncm+1 ! ncm=1 is alkyne, =2 is cp
                     enddo
                  endif
               enddo
               if(ncm.eq.0) etacoord=.false.
            endif
         endif
         if(etacoord)then
          itag(i)=-1
          nbdum(1:20,i)=nbm(1:20,i)
         else
          nbdum(1:20,i)=nbf(1:20,i) ! take full set of neighbors by default
         endif
      enddo

      do i=1,natoms
         ati  = at(i)
         hyb(i)=0    ! don't know it
         nbdiff =nbf(20,i)-topo%nb (20,i)
         nbmdiff=nbf(20,i)-nbm(20,i)
         nb20i=nbdum(20,i)
         nh=0
         no=0
         do j=1,nb20i
            if(at(nbdum(j,i)).eq.1) nh=nh+1
            if(at(nbdum(j,i)).eq.8) no=no+1
         enddo
! H
         if(param%group(ati).eq.1) then
            if(nb20i.eq.2)               hyb(i)=1 ! bridging H
            if(nb20i.gt.2)               hyb(i)=3 ! M+ tetra coord
            if(nb20i.gt.4)               hyb(i)=0 ! M+ HC
         endif
! Be
         if(param%group(ati).eq.2) then
            if(nb20i.eq.2)               hyb(i)=1 ! bridging M
            if(nb20i.gt.2)               hyb(i)=3 ! M+ tetra coord
            if(nb20i.gt.4)               hyb(i)=0 !
         endif
! B
         if(param%group(ati).eq.3) then
            if(nb20i.gt.4)                               hyb(i)=3
            if(nb20i.gt.4.and.ati.gt.10.and.nbdiff.eq.0) hyb(i)=5
            if(nb20i.eq.4)                               hyb(i)=3
            if(nb20i.eq.3)                               hyb(i)=2
            if(nb20i.eq.2)                               hyb(i)=1
         endif
! C
         if(param%group(ati).eq.4) then
            if(nb20i.ge.4)                               hyb(i)=3
            if(nb20i.gt.4.and.ati.gt.10.and.nbdiff.eq.0) hyb(i)=5
            if(nb20i.eq.3)                               hyb(i)=2
            if(nb20i.eq.2) then
              call bangl(xyz,nbdum(1,i),i,nbdum(2,i),phi)
              if(phi*180./pi.lt.150.0)then                         ! geometry dep. setup! GEODEP
                                                         hyb(i)=2  ! otherwise, carbenes will not be recognized
                                                        itag(i)=1  ! tag for Hueckel and HB routines
              else
                                                         hyb(i)=1  ! linear triple bond etc
              endif
              if(topo%qa(i).lt.-0.4)                          then
                                                         hyb(i)=2
                                                        itag(i)=0  ! tag for Hueckel and HB routines
              endif
            endif
            if(nb20i.eq.1)                               hyb(i)=1  ! CO
         endif
! N
         if(param%group(ati).eq.5) then
            if(nb20i.ge.4)                               hyb(i)=3
            if(nb20i.gt.4.and.ati.gt.10.and.nbdiff.eq.0) hyb(i)=5
            if(nb20i.eq.3)                               hyb(i)=3
            if(nb20i.eq.3.and.ati.eq.7)  then
               kk=0
               ll=0
               nn=0
               do j=1,3
                  jj=nbdum(j,i)
                  if(at(jj).eq. 8.and.topo%nb(20,jj).eq.1) kk=kk+1 ! check for NO2 or R2-N=O
                  if(at(jj).eq. 5.and.topo%nb(20,jj).eq.4) ll=ll+1 ! check for B-N, if the CN(B)=4 the N is loosely bound and sp2
                  if(at(jj).eq.16.and.topo%nb(20,jj).eq.4) nn=nn+1 ! check for N-SO2-
               enddo
               if(nn.eq.1.and.ll.eq.0.and.kk.eq.0)       hyb(i)=3
               if(ll.eq.1.and.nn.eq.0)                   hyb(i)=2
               if(kk.ge.1) then
                                                         hyb(i)=2
                                                        itag(i)=1  ! tag for Hueckel with no el. for the N in NO2
               endif
               if(nbmdiff.gt.0.and.nn.eq.0)              hyb(i)=2  ! pyridin N coord. to heavy atom
            endif
            if(nb20i.eq.2) then
                                                         hyb(i)=2
               call bangl(xyz,nbdum(1,i),i,nbdum(2,i),phi)
               jj=nbdum(1,i)
               kk=nbdum(2,i)
               if(nbdum(20,jj).eq.1.and.at(jj).eq.6)     hyb(i)=1  ! R-N=C
               if(nbdum(20,kk).eq.1.and.at(kk).eq.6)     hyb(i)=1  ! R-N=C
               if(nbdum(20,jj).eq.1.and.at(jj).eq.7)     hyb(i)=1  ! R-N=N in e.g. diazomethane
               if(nbdum(20,kk).eq.1.and.at(kk).eq.7)     hyb(i)=1  ! R-N=N in e.g. diazomethane
               if(nbdum(1,i).gt.0.and.param%metal(at(nbdum(1,i))).gt.0) hyb(i)=1 ! M-NC-R in e.g. nitriles
               if(nbdum(2,i).gt.0.and.param%metal(at(nbdum(2,i))).gt.0) hyb(i)=1 ! M-NC-R in e.g. nitriles
               if(at(jj).eq.7.and.at(kk).eq.7.and. &
     &          nbdum(20,jj).le.2.and.nbdum(20,kk).le.2) hyb(i)=1  ! N=N=N
               if(phi*180./pi.gt.lintr)                 hyb(i)=1  ! geometry dep. setup! GEODEP
            endif
            if(nb20i.eq.1)                               hyb(i)=1
         endif
! O
         if(param%group(ati).eq.6) then
            if(nb20i.ge.3)                               hyb(i)=3
            if(nb20i.gt.3.and.ati.gt.10.and.nbdiff.eq.0) hyb(i)=5
            if(nb20i.eq.2)                               hyb(i)=3
            if(nb20i.eq.2.and.nbmdiff.gt.0) then
               call nn_nearest_noM(i,natoms,at,topo%nb,rab,j,param) ! CN of closest non-M atom
                                        if(j.eq.3)       hyb(i)=2 ! M-O-X konj
                                        if(j.eq.4)       hyb(i)=3 ! M-O-X non
            endif
            if(nb20i.eq.1)                               hyb(i)=2
            if(nb20i.eq.1.and.nbdiff.eq.0) then
            if(topo%nb(20,topo%nb(1,i)).eq.1)                      hyb(i)=1 ! CO
            endif
         endif
! F
         if(param%group(ati).eq.7) then
            if(nb20i.eq.2)                               hyb(i)=1
            if(nb20i.gt.2.and.ati.gt.10)                 hyb(i)=5
         endif
! Ne
         if(param%group(ati).eq.8) then
                                                         hyb(i)=0
            if(nb20i.gt.0.and.ati.gt.2)                  hyb(i)=5
         endif
! done with main groups
         if(param%group(ati).le.0) then ! TMs
            nni=nb20i
            if(nh.ne.0.and.nh.ne.nni) nni=nni-nh ! don't count Hs
            if(nni.le.2)                               hyb(i)=1
            if(nni.le.2.and.param%group(ati).le.-6)          hyb(i)=2
            if(nni.eq.3)                               hyb(i)=2
            if(nni.eq.4.and.param%group(ati).gt.-7)          hyb(i)=3  ! early TM, tetrahedral
            if(nni.eq.4.and.param%group(ati).le.-7)          hyb(i)=3  ! late TM, square planar
            if(nni.eq.5.and.param%group(ati).eq.-3)          hyb(i)=3  ! Sc-La are tetrahedral CN=5
         endif
      enddo

      topo%nb = nbdum ! list is complete but hyb determination is based only on reduced (without metals) list

      deallocate(nbdum)

      j = 0
      do i=1,natoms
         if(topo%nb(20,i).gt.12) j = j +1
         do k=1,topo%nb(20,i)
            kk=topo%nb(k,i)
            if(at(kk).eq.6.and.at(i).eq.6.and.itag(i).eq.1.and.itag(kk).eq.1) then ! check the very special situation of
               itag(i) =0                                                          ! two carbene C bonded which is an arine
               itag(kk)=0
            endif
         enddo
      enddo
      if(dble(j)/dble(natoms).gt.0.3) then
         call env%error(' too many atoms with extreme high CN', source)
      end if

      end subroutine gfnff_neigh

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! fill neighbor list
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine getnb(n,at,rad,r,mchar,icase,f,f2,nbf,nb,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer n,at(n),nbf(20,n),nb(20,n)
      real*8 rad(n*(n+1)/2),r(n*(n+1)/2),mchar(n),f,f2

      integer i,j,k,nn,icase,hc_crit,nnfi,nnfj,lin
      integer tag(n*(n+1)/2)
      real*8 rco,fm

      nb  = 0 ! resulting array (nbf is full from first call)
      tag = 0
      do i=1,n
         nnfi=nbf(20,i)                  ! full CN of i, only valid for icase > 1
         do j=1,i-1
            nnfj=nbf(20,j)               ! full CN of i
            fm=1.0d0
!           full case
            if(icase.eq.1)then
               if(param%metal(at(i)).eq.2) fm=fm*f2 !change radius for metal atoms
               if(param%metal(at(j)).eq.2) fm=fm*f2
               if(param%metal(at(i)).eq.1) fm=fm*(f2+0.025)
               if(param%metal(at(j)).eq.1) fm=fm*(f2+0.025)
            endif
!           no HC atoms
            if(icase.eq.2)then
               hc_crit = 6
               if(param%group(at(i)).le.2) hc_crit = 4
               if(nnfi.gt.hc_crit) cycle
               hc_crit = 6
               if(param%group(at(j)).le.2) hc_crit = 4
               if(nnfj.gt.hc_crit) cycle
            endif
!           no metals and unusually coordinated stuff
            if(icase.eq.3)then
               if(mchar(i).gt.0.25 .or. param%metal(at(i)).gt.0) cycle   ! metal case TMonly ?? TODO
               if(mchar(j).gt.0.25 .or. param%metal(at(j)).gt.0) cycle   ! metal case
               if(nnfi.gt.param%normcn(at(i)).and.at(i).gt.10)   cycle   ! HC case
               if(nnfj.gt.param%normcn(at(j)).and.at(j).gt.10)   cycle   ! HC case
            endif
            k=lin(j,i)
            rco=rad(k) !(rad(i)+rad(j))/0.5291670d0
!               R         est. R0
            if(r(k).lt. fm * f * rco) tag(k)=1
         enddo
      enddo

      do i=1,n
         nn = 0
         do j=1,n
            if(tag(lin(j,i)).eq.1.and.i.ne.j) then
               nn=nn+1
               nb(nn,i)=j
            endif
         enddo
         nb(20,i) = nn
      enddo

      end subroutine getnb 

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! find the CN of nearest non metal of atom i
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nn_nearest_noM(ii,n,at,nb,r,nn,param)
      implicit none
      type(TGFFData), intent(in) :: param
      integer ii,n,at(n),nn,nb(20,n)
      real*8 r(n*(n+1)/2)

      integer jmin,j,jj,lin
      real*8 rmin

      nn=0
      rmin=1.d+42
      jmin=0
      do j=1,nb(20,ii)
         jj=nb(j,ii)
         if(param%metal(at(jj)).ne.0) cycle
         if(r(lin(jj,ii)).lt.rmin)then
            rmin=r(lin(jj,ii))
            jmin=jj
         endif
      enddo

      if(jmin.gt.0) nn=nb(20,jmin)

      end subroutine nn_nearest_noM

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! find smallest ring in which atom i is located
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ringsatom(n,i,c,s,rings)
      implicit none
      integer n,i,k,l,c(10,20,n),s(20,n),rings,rings1

      rings=99
      do k=1,s(20,i)    ! all rings of atom i
         if(s(k,i).lt.rings) rings=s(k,i)
      enddo

      end subroutine ringsatom

!ccccccccccccccccccccccccccccccccccccccccccccccccc
! find smallest ring in which bond i-j is located
!ccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ringsbond(n,i,j,c,s,rings)
      implicit none
      integer n,i,j,k,l,c(10,20,n),s(20,n),rings,rings1,rings2

      rings1=99
      rings2=99
      do k=1,s(20,i)    ! all rings of atom i
         do l=1,s(k,i)  ! all atoms of ring k
            if(c(l,k,i).eq.j.and.s(k,i).lt.rings1)then
               rings1=s(k,i)
            endif
         enddo
      enddo
      do k=1,s(20,j)    ! all rings of atom i
         do l=1,s(k,j)  ! all atoms of ring k
            if(c(l,k,j).eq.i.and.s(k,j).lt.rings2)then
               rings2=s(k,j)
            endif
         enddo
      enddo
      continue
      rings=min(rings1,rings2)
      if(rings.eq.99) rings=0

      end subroutine ringsbond

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
! find smallest ring in which angle i-j-k is located
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ringsbend(n,i,j,k,c,s,rings)
      implicit none
      integer n,i,j,k,rings
      integer c(10,20,n),s(20,n)
      integer itest,rings1,rings2,rings3,m,l

      if(s(20,i).eq.0.or.s(20,j).eq.0.or.s(20,k).eq.0)then
         rings=0
         return
      endif

      rings1=99
      rings2=99
      rings3=99

      do m=1,s(20,i)    ! all rings of atom i
         itest=0
         do l=1,s(m,i)  ! all atoms of ring m
            if(c(l,m,i).eq.j.or.c(l,m,i).eq.k) itest=itest+1
         enddo
         if(itest.eq.2.and.s(m,i).lt.rings1) rings1=s(m,i)
      enddo
      do m=1,s(20,j)    ! all rings of atom j
         itest=0
         do l=1,s(m,j)  ! all atoms of ring m
            if(c(l,m,j).eq.i.or.c(l,m,j).eq.k) itest=itest+1
         enddo
         if(itest.eq.2.and.s(m,j).lt.rings2) rings2=s(m,j)
      enddo
      do m=1,s(20,k)    ! all rings of atom k
         itest=0
         do l=1,s(m,k)  ! all atoms of ring m
            if(c(l,m,k).eq.i.or.c(l,m,k).eq.j) itest=itest+1
         enddo
         if(itest.eq.2.and.s(m,j).lt.rings3) rings3=s(m,k)
      enddo

      rings=min(rings1,rings2,rings3)
      if(rings.eq.99) rings=0

      end subroutine ringsbend

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! find smallest torsion in which angle i-j-k-l is located
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ringstors(n,i,j,k,l,c,s,rings)
      implicit none
      integer n,i,j,k,l,rings
      integer c(10,20,n),s(20,n)
      integer itest,rings1,rings2,rings3,rings4,m,a
!     integer ht1,ht2,ht3,ht4

      if(s(20,i).eq.0.or.s(20,j).eq.0.or.s(20,k).eq.0.or.s(20,l).eq.0)then
         rings=0
         return
      endif

      rings1=99
      rings2=99
      rings3=99
      rings4=99

      do m=1,s(20,i)    ! all rings of atom i
         itest=0
         do a=1,s(m,i)  ! all atoms of ring m
            if(c(a,m,i).eq.j.or.c(a,m,i).eq.k.or.c(a,m,i).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,i).lt.rings1) then
            rings1=s(m,i)
!           ht1=c(m,19,i)
         endif
      enddo
      do m=1,s(20,j)    ! all rings of atom j
         itest=0
         do a=1,s(m,j)  ! all atoms of ring m
            if(c(a,m,j).eq.i.or.c(a,m,j).eq.k.or.c(a,m,j).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,j).lt.rings2) then
            rings2=s(m,j)
!           ht2=c(m,19,j)
         endif
      enddo
      do m=1,s(20,k)    ! all rings of atom k
         itest=0
         do a=1,s(m,k)  ! all atoms of ring m
            if(c(a,m,k).eq.i.or.c(a,m,k).eq.j.or.c(a,m,k).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,k).lt.rings3) then
            rings3=s(m,k)
!           ht3=c(m,19,k)
         endif
      enddo
      do m=1,s(20,l)    ! all rings of atom k
         itest=0
         do a=1,s(m,l)  ! all atoms of ring m
            if(c(a,m,l).eq.i.or.c(a,m,l).eq.k.or.c(a,m,l).eq.j) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,l).lt.rings4) then
            rings4=s(m,l)
!           ht4=c(m,19,l)
         endif
      enddo

      rings=min(rings1,rings2,rings3,rings4)
      if(rings.eq.99) then
         rings=0
      endif

      end subroutine ringstors

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! find smallest torsion in which angle i-j-k-l is located
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine ringstorl(n,i,j,k,l,c,s,ringl)
      implicit none
      integer n,i,j,k,l,ringl
      integer c(10,20,n),s(20,n)
      integer itest,rings1,rings2,rings3,rings4,m,a

      if(s(20,i).eq.0.or.s(20,j).eq.0.or.s(20,k).eq.0.or.s(20,l).eq.0)then
         ringl=0
         return
      endif

      rings1=-99
      rings2=-99
      rings3=-99
      rings4=-99

      do m=1,s(20,i)    ! all rings of atom i
         itest=0
         do a=1,s(m,i)  ! all atoms of ring m
            if(c(a,m,i).eq.j.or.c(a,m,i).eq.k.or.c(a,m,i).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,i).gt.rings1) then
            rings1=s(m,i)
         endif
      enddo
      do m=1,s(20,j)    ! all rings of atom j
         itest=0
         do a=1,s(m,j)  ! all atoms of ring m
            if(c(a,m,j).eq.i.or.c(a,m,j).eq.k.or.c(a,m,j).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,j).gt.rings2) then
            rings2=s(m,j)
         endif
      enddo
      do m=1,s(20,k)    ! all rings of atom k
         itest=0
         do a=1,s(m,k)  ! all atoms of ring m
            if(c(a,m,k).eq.i.or.c(a,m,k).eq.j.or.c(a,m,k).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,k).gt.rings3) then
            rings3=s(m,k)
         endif
      enddo
      do m=1,s(20,l)    ! all rings of atom k
         itest=0
         do a=1,s(m,l)  ! all atoms of ring m
            if(c(a,m,l).eq.i.or.c(a,m,l).eq.k.or.c(a,m,l).eq.j) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,l).gt.rings4) then
            rings4=s(m,l)
         endif
      enddo

      ringl=max(rings1,rings2,rings3,rings4)
      if(ringl.eq.-99) then
         ringl=0
      endif

      end subroutine ringstorl

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      logical function chktors(n,xyz,i,j,k,l)  ! true if dihedral angle is bad i.e. near 180
      implicit none
      integer n,i,j,k,l
      real*8 xyz(3,n),phi

      chktors=.true.

      call bangl(xyz,j,i,k,phi)
!     write(env%unit,*) phi*180./3.1415926d0
      if(phi*180./3.1415926d0.gt.170.0d0) return
      call bangl(xyz,i,j,l,phi)
!     write(env%unit,*) phi*180./3.1415926d0
      if(phi*180./3.1415926d0.gt.170.0d0) return

      chktors=.false.

      end function chktors

      logical function chkrng(nn,n,c)
      implicit none
      integer n,idum(nn),nn,c(10),i,j
      chkrng=.true.
      idum=0
      do i=1,n
         idum(c(i))=idum(c(i))+1
      enddo
      j=0
      do i=1,nn
         if(idum(i).eq.1) j=j+1
      enddo
      if(j.ne.n) chkrng=.false.
      end function chkrng

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gfnff_hbset(n,at,xyz,sqrab,topo,nlist,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      type(TGFFNeighbourList), intent(inout) :: nlist
      integer n
      integer at(n)
      real(wp) sqrab(n*(n+1)/2)
      real(wp) xyz(3,n)
      real(wp), intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh
      real(wp) rab,rmsd
      logical ijnonbond

      rmsd = sqrt(sum((xyz-nlist%hbrefgeo)**2))/dble(n)

      if(rmsd.lt.1.d-6 .or. rmsd.gt. 0.3d0) then ! update list if first call or substantial move occured

      nlist%nhb1=0
      nlist%nhb2=0
      do ix=1,topo%nathbAB
         i=topo%hbatABl(1,ix)
         j=topo%hbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr1)cycle
         ijnonbond=topo%bpair(ij).ne.1
         do k=1,topo%nathbH
            nh=topo%hbatHl(k)
            inh=lin(i,nh)
            jnh=lin(j,nh)
            if(topo%bpair(inh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               nlist%nhb2=nlist%nhb2+1
               nlist%hblist2(1,nlist%nhb2)=i
               nlist%hblist2(2,nlist%nhb2)=j
               nlist%hblist2(3,nlist%nhb2)=nh
            elseif(topo%bpair(jnh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               nlist%nhb2=nlist%nhb2+1
               nlist%hblist2(1,nlist%nhb2)=j
               nlist%hblist2(2,nlist%nhb2)=i
               nlist%hblist2(3,nlist%nhb2)=nh
            elseif(rab+sqrab(inh)+sqrab(jnh).lt.hbthr2) then
               nlist%nhb1=nlist%nhb1+1
               nlist%hblist1(1,nlist%nhb1)=i
               nlist%hblist1(2,nlist%nhb1)=j
               nlist%hblist1(3,nlist%nhb1)=nh
            endif
         enddo
      enddo

      nlist%nxb =0
      do ix=1,topo%natxbAB
         i =topo%xbatABl(1,ix)
         j =topo%xbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr2)cycle
         nlist%nxb=nlist%nxb+1
         nlist%hblist3(1,nlist%nxb)=i
         nlist%hblist3(2,nlist%nxb)=j
         nlist%hblist3(3,nlist%nxb)=topo%xbatABl(3,ix)
      enddo

      nlist%hbrefgeo = xyz

      endif  ! else do nothing

      end subroutine gfnff_hbset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bond_hbset(n,at,xyz,sqrab,bond_hbn,bond_hbl,topo,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in) :: n
      integer,intent(in) :: at(n)
      integer,intent(in) :: bond_hbn
      integer,intent(out) :: bond_hbl(3,bond_hbn)
      real(wp),intent(in) :: sqrab(n*(n+1)/2)
      real(wp),intent(in) :: xyz(3,n)
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh
      integer bond_nr
      real(wp) rab
      logical ijnonbond


      bond_nr=0
      bond_hbl=0
      do ix=1,topo%nathbAB
         i=topo%hbatABl(1,ix)
         j=topo%hbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr1)cycle
         ijnonbond=topo%bpair(ij).ne.1
         do k=1,topo%nathbH
            nh=topo%hbatHl(k)
            inh=lin(i,nh)
            jnh=lin(j,nh)
            if(topo%bpair(inh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               bond_nr=bond_nr+1
               bond_hbl(1,bond_nr)=i
               bond_hbl(2,bond_nr)=j
               bond_hbl(3,bond_nr)=nh
            elseif(topo%bpair(jnh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               bond_nr=bond_nr+1
               bond_hbl(1,bond_nr)=j
               bond_hbl(2,bond_nr)=i
               bond_hbl(3,bond_nr)=nh
            endif
         enddo
      enddo

end subroutine bond_hbset

subroutine bond_hbset0(n,at,xyz,sqrab,bond_hbn,topo,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in) :: n
      integer,intent(in) :: at(n)
      integer,intent(out) :: bond_hbn
      real(wp),intent(in) :: sqrab(n*(n+1)/2)
      real(wp),intent(in) :: xyz(3,n)
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh
      real(wp) rab
      logical ijnonbond


      bond_hbn=0
      do ix=1,topo%nathbAB
         i=topo%hbatABl(1,ix)
         j=topo%hbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr1)cycle
         ijnonbond=topo%bpair(ij).ne.1
         do k=1,topo%nathbH
            nh=topo%hbatHl(k)
            inh=lin(i,nh)
            jnh=lin(j,nh)
            if(topo%bpair(inh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               bond_hbn=bond_hbn+1
            elseif(topo%bpair(jnh).eq.1.and.ijnonbond)then ! exclude cases where A and B are bonded
               bond_hbn=bond_hbn+1
            endif
         enddo
      enddo

end subroutine bond_hbset0

subroutine bond_hb_AHB_set(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,topo)
      use xtb_gfnff_param
      implicit none
      !Dummy
      type(TGFFTopology), intent(inout) :: topo
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(3,bond_hbn)
      integer,intent(in)  :: tot_AHB_nr
      integer,intent(inout) :: lin_AHB(0:tot_AHB_nr)
      !Stack
      integer :: i,j
      integer :: ii,jj
      integer :: ia,ja
      integer :: lin
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB
      integer :: tot_count
      integer :: AH_count
      integer :: B_count
      integer :: lin_diff

      tot_count=0
      AH_count=0
      B_count=0
      lin_diff=0

      do i=1,numbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ia=at(ii)
         ja=at(jj)
         if (ia.eq.1) then
           hbH=ii
           hbA=jj
         else if (ja.eq.1) then
           hbH=jj
           hbA=ii
         else
           cycle
         end if
         if (at(hbA).eq.7.or.at(hbA).eq.8) then
            do j = 1, bond_hbn
               Bat = bond_hbl(2,j)
               atB = at(Bat)
               Aat  = bond_hbl(1,j)
               Hat  = bond_hbl(3,j)
               if (hbA.eq.Aat.and.hbH.eq.Hat) then
                  if (atB.eq.7.or.atB.eq.8) then
                     tot_count = tot_count + 1
                     lin_AHB(tot_count) = lin(hbA,hbH)
                     lin_diff=lin_AHB(tot_count)-lin_AHB(tot_count-1)
                     if (lin_diff.eq.0) then
                        B_count = B_count + 1
                     end if
                     !Next AH pair
                     if (lin_diff.ne.0) then
                        AH_count = AH_count + 1
                        topo%bond_hb_AH(1,AH_count) = hbA
                        topo%bond_hb_AH(2,AH_count) = hbH
                        !Reset B count
                        B_count = 1
                     end if
                     topo%bond_hb_Bn(AH_count) = B_count
                     topo%bond_hb_B(B_count,AH_count) = Bat
                  end if
               else
                 cycle
               end if
            end do
         end if
      end do

end subroutine bond_hb_AHB_set

subroutine bond_hb_AHB_set1(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,AH_count,bmax,topo)
      use xtb_gfnff_param
      implicit none
      !Dummy
      type(TGFFTopology), intent(inout) :: topo
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(3,bond_hbn)
      integer,intent(in)  :: tot_AHB_nr
      integer,intent(inout) :: lin_AHB(0:tot_AHB_nr)
      integer,intent(out) :: AH_count
      integer,intent(out) :: bmax
      !Stack
      integer :: i,j
      integer :: ii,jj
      integer :: ia,ja
      integer :: lin
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB
      integer :: tot_count
      integer :: B_count
      integer :: lin_diff

      tot_count=0
      AH_count=0
      B_count=1
      bmax=1
      lin_diff=0

      do i=1,numbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ia=at(ii)
         ja=at(jj)
         if (ia.eq.1) then
           hbH=ii
           hbA=jj
         else if (ja.eq.1) then
           hbH=jj
           hbA=ii
         else
           cycle
         end if
         if (at(hbA).eq.7.or.at(hbA).eq.8) then
            do j = 1, bond_hbn
               Bat = bond_hbl(2,j)
               atB = at(Bat)
               Aat  = bond_hbl(1,j)
               Hat  = bond_hbl(3,j)
               if (hbA.eq.Aat.and.hbH.eq.Hat) then
                 if (atB.eq.7.or.atB.eq.8) then
                     tot_count = tot_count + 1
                     lin_AHB(tot_count) = lin(hbA,hbH)
                     lin_diff=lin_AHB(tot_count)-lin_AHB(tot_count-1)
                     if (lin_diff.eq.0) B_count = B_count + 1
                     if (lin_diff.ne.0) then
                        AH_count = AH_count + 1
                        B_count = 1
                     end if
                     if (B_count.gt.bmax) bmax = B_count
                  end if
               else
                 cycle
               end if
            end do
            topo%nr_hb(i) = B_count
         end if
      end do

end subroutine bond_hb_AHB_set1

subroutine bond_hb_AHB_set0(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,topo)
      use xtb_gfnff_param
      implicit none
      !Dummy
      type(TGFFTopology), intent(in) :: topo
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(3,bond_hbn)
      integer,intent(out) :: tot_AHB_nr
      !Stack
      integer :: i,j
      integer :: ii,jj
      integer :: ia,ja
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB

      tot_AHB_nr=0

      do i=1,numbond
         ii=topo%blist(1,i)
         jj=topo%blist(2,i)
         ia=at(ii)
         ja=at(jj)
         if (ia.eq.1) then
           hbH=ii
           hbA=jj
         else if (ja.eq.1) then
           hbH=jj
           hbA=ii
         else
           cycle
         end if
         if (at(hbA).eq.7.or.at(hbA).eq.8) then
            do j = 1, bond_hbn
               Bat = bond_hbl(2,j)
               atB = at(Bat)
               Aat  = bond_hbl(1,j)
               Hat  = bond_hbl(3,j)
               if (hbA.eq.Aat.and.hbH.eq.Hat) then
                 if (atB.eq.7.or.atB.eq.8) then
                     tot_AHB_nr = tot_AHB_nr + 1
                  end if
               else
                 cycle
               end if
            end do
         end if
      end do

end subroutine bond_hb_AHB_set0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gfnff_hbset0(n,at,xyz,sqrab,topo,nhb1,nhb2,nxb,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      integer, intent(out) :: nhb1
      integer, intent(out) :: nhb2
      integer, intent(out) :: nxb
      integer n
      integer at(n)
      real(wp) sqrab(n*(n+1)/2)
      real(wp) xyz(3,n)
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh
      logical ijnonbond
      real(wp) rab

      nhb1=0
      nhb2=0
      do ix=1,topo%nathbAB
         i=topo%hbatABl(1,ix)
         j=topo%hbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr1)cycle
         ijnonbond=topo%bpair(ij).ne.1
         do k=1,topo%nathbH
            nh=topo%hbatHl(k)
            inh=lin(i,nh)
            jnh=lin(j,nh)
            if(topo%bpair(inh).eq.1.and.ijnonbond)then
               nhb2=nhb2+1
            elseif(topo%bpair(jnh).eq.1.and.ijnonbond)then
               nhb2=nhb2+1
            elseif(rab+sqrab(inh)+sqrab(jnh).lt.hbthr2) then
               nhb1=nhb1+1
            endif
         enddo
      enddo

      nxb =0
      do ix=1,topo%natxbAB
         i =topo%xbatABl(1,ix)
         j =topo%xbatABl(2,ix)
         ij=j+i*(i-1)/2
         rab=sqrab(ij)
         if(rab.gt.hbthr2)cycle
         nxb=nxb+1
      enddo

      end subroutine gfnff_hbset0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! HB strength
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hbonds(i,j,ci,cj,param,topo)
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      type(TGFFData), intent(in) :: param
      integer i,j
      integer ati,atj
      real*8 ci(2),cj(2)
      ci(1)=topo%hbbas(i)
      cj(1)=topo%hbbas(j)
      ci(2)=topo%hbaci(i)
      cj(2)=topo%hbaci(j)
      end subroutine hbonds

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ring analysis routine, don't touch ;-)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine getring36(n,at,nbin,a0_in,cout,irout)
      implicit none
      integer cout(10,20),irout(20)  ! output: atomlist, ringsize, # of rings in irout(20)
      integer at(n)
      integer n,nbin(20,n),a0,i,nb(20,n),a0_in
      integer    i1,i2,i3,i4,i5,i6,i7,i8,i9,i10
      integer n0,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10
      integer    a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
      integer maxr
      parameter (maxr=500)
      integer list(n),m,mm,nn,c(10),cdum(10,maxr),iring
      integer adum1(0:n),adum2(0:n),kk,j,idum(maxr),same(maxr),k
      real*8  w(n),av,sd

      cout = 0
      irout= 0
!     if(n.le.2.or.at(a0_in).gt.10.or.nbin(19,a0_in).eq.1) return
      if(n.le.2.                  .or.nbin(19,a0_in).eq.1) return

      nn=nbin(20,a0_in)

      cdum=0
      kk=0
      do m=1,nn
!     if(nb(m,a0_in).eq.1) cycle ! check (comment out)
      nb=nbin
      if(nb(m,a0_in).eq.1) cycle ! check (comment out)
      do i=1,n
         if(nb(20,i).eq.1)nb(20,i)=0
      enddo

      do mm=1,nn
         w(mm)=dble(mm)
         list(mm)=mm
      enddo
      w(m)   =0.0d0
      call ssort(nn,w,list)
      do mm=1,nn
         nb(mm,a0_in)=nbin(list(mm),a0_in)
      enddo

      iring =0
      c     =0

      a0=a0_in
      n0=nb(20,a0)

      do i1=1,n0
         a1=nb(i1,a0)
         if(a1.eq.a0) cycle
         n1=nb(20,a1)
         do i2=1,n1
            a2=nb(i2,a1)
            if(a2.eq.a1) cycle
            n2=nb(20,a2)
            do i3=1,n2
               a3=nb(i3,a2)
               n3=nb(20,a3)
               if(a3.eq.a2) cycle
               c(1)=a1
               c(2)=a2
               c(3)=a3
               if(a3.eq.a0.and.chkrng(n,3,c))then
                iring=3
                if(kk.eq.maxr) goto 99
                kk=kk+1
                cdum(1:iring,kk)=c(1:iring)
                idum(kk)=iring
               endif
               do i4=1,n3
                  a4=nb(i4,a3)
                  n4=nb(20,a4)
                  if(a4.eq.a3) cycle
                  c(4)=a4
                  if(a4.eq.a0.and.chkrng(n,4,c))then
                   iring=4
                   if(kk.eq.maxr) goto 99
                   kk=kk+1
                   cdum(1:iring,kk)=c(1:iring)
                   idum(kk)=iring
                  endif
                  do i5=1,n4
                     a5=nb(i5,a4)
                     n5=nb(20,a5)
                     if(a5.eq.a4) cycle
                     c(5)=a5
                     if(a5.eq.a0.and.chkrng(n,5,c))then
                      iring=5
                      if(kk.eq.maxr) goto 99
                      kk=kk+1
                      cdum(1:iring,kk)=c(1:iring)
                      idum(kk)=iring
                     endif
                     do i6=1,n5
                        a6=nb(i6,a5)
                        n6=nb(20,a6)
                        if(a6.eq.a5) cycle
                        c(6)=a6
                        if(a6.eq.a0.and.chkrng(n,6,c))then
                         iring=6
                         if(kk.eq.maxr) goto 99
                         kk=kk+1
                         cdum(1:iring,kk)=c(1:iring)
                         idum(kk)=iring
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

 99   continue

      enddo

! compare
      same=0
      do i=1,kk
         do j=i+1,kk
            if(idum(i).ne.idum(j)) cycle ! different ring size
            if(same(j).eq.1      ) cycle ! already double
            adum1=0
            adum2=0
            do m=1,10
               i1=cdum(m,i)
               i2=cdum(m,j)
               adum1(i1)=1
               adum2(i2)=1
            enddo
            if(sum(abs(adum1-adum2)).ne.0) then
                  same(j)=0
            else
                  same(j)=1
            endif
         enddo
      enddo

      m=0
      do i=1,kk
         if(same(i).eq.0) then
            m=m+1
            irout(m)=idum(i)     ! number of atoms in ring m
            nn=idum(i)
            cout(1:nn,m)=cdum(1:nn,i)
            if(m.gt.19) then
               m=19
               goto 999
            endif
         endif
      enddo
999   irout(20)=m  ! number of rings for this atom

      return
      end subroutine getring36

      subroutine ssort(n,edum,ind)
      implicit none
      integer n,ii,k,j,m,i,sc1
      real*8 edum(n),pp
      integer ind(n)

      do 140   ii = 2, n
         i = ii - 1
         k = i
         pp= edum(i)
         do 120   j = ii, n
            if (edum(j) .gt. pp) go to 120
            k = j
            pp= edum(j)
  120    continue
         if (k .eq. i) go to 140
         edum(k) = edum(i)
         edum(i) = pp
         sc1=ind(i)
         ind(i)=ind(k)
         ind(k)=sc1
  140 continue

      end subroutine ssort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! neighbor only version of EEQ model
! included up to 1,4 interactions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine goedeckera(env,n,at,nb,pair,q,es,topo)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : mctc_sytrf, mctc_sytrs
   implicit none
   character(len=*), parameter :: source = 'gfnff_ini2_goedeckera'
   type(TEnvironment), intent(inout) :: env
   type(TGFFTopology), intent(in) :: topo
   integer, intent(in)  :: n          ! number of atoms
   integer, intent(in)  :: at(n)      ! ordinal numbers
   integer, intent(in)  :: nb(20,n)   ! neighbors
   real(wp),intent(in)  :: pair(n*(n+1)/2)
   real(wp),intent(out) :: q(n)       ! output charges
   real(wp),intent(out) :: es         ! ES energy

!  local variables
   logical :: exitRun
   integer  :: m,i,j,k,l,ii,jj,kk
   integer  :: ij,lj
   integer,allocatable :: ipiv(:)

   real(wp) :: gammij,sief1,sief2
   real(wp) :: r2,r0
   real(wp) :: rij
   real(wp) :: tsqrt2pi,bohr
   real(wp) :: tmp
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: x(:)

!  parameter
   parameter (tsqrt2pi = 0.797884560802866_wp)

   m=n+topo%nfrag ! # atoms frag constrain
   allocate(A(m,m),x(m),ipiv(m))

!  call prmati(6,pair,n,0,'pair')

   A = 0

!  setup RHS
   do i=1,n
      x(i)    =topo%chieeq(i) ! EN of atom
      A(i,i)  =topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i))
   enddo

!  setup A matrix
   do i=1,n
      do j=1,i-1
         ij = i*(i-1)/2+j
         rij=pair(ij)
         r2 = rij*rij
         gammij=1.d0/sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
         tmp = erf(gammij*rij)/rij
         A(j,i) = tmp
         A(i,j) = tmp
      enddo
   enddo

!  fragment charge constrain
   do i=1,topo%nfrag
      x(n+i)=topo%qfrag(i)
      do j=1,n
         if(topo%fraglist(j).eq.i) then
            A(n+i,j)=1
            A(j,n+i)=1
         endif
      enddo
   enddo
!  call prmat(6,A,m,m,'A ini')

   call mctc_sytrf(env, a, ipiv)
   call mctc_sytrs(env, a, x, ipiv)

   call env%check(exitRun)
   if (exitRun) then
      call env%error('Solving linear equations failed', source)
      return
   end if

   q(1:n) = x(1:n)

   if(n.eq.1) q(1)=topo%qfrag(1)

!  energy
      es = 0.0_wp
      do i=1,n
      ii = i*(i-1)/2
      do j=1,i-1
         ij = ii+j
         rij=pair(ij)
         gammij=1.d0/sqrt(topo%alpeeq(i)+topo%alpeeq(j)) ! squared above
         tmp = erf(gammij*rij)/rij
         es = es + q(i)*q(j)*tmp/rij
      enddo
      es = es - q(i)* topo%chieeq(i) &
     &        + q(i)*q(i)*0.5d0*(topo%gameeq(i)+tsqrt2pi/sqrt(topo%alpeeq(i)))
      enddo

end subroutine goedeckera

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! condense charges to heavy atoms based on topology
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine qheavy(n,at,nb,q)
      implicit none
      integer,intent(in)   ::  n,nb(20,n),at(n)
      real*8,intent(inout) ::  q(n)

      integer i,j,k
      real*8 qtmp(n)

      qtmp = q
      do i=1,n
         if(at(i).ne.1) cycle
         qtmp(i)=0
         do j=1,nb(20,i)
            k=nb(j,i)
            qtmp(k)=qtmp(k)+q(i)/dble(nb(20,i))  ! could be a bridging H
         enddo
      enddo

      q = qtmp

      end subroutine qheavy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine number of cov. bonds between atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine nbondmat(n,nb,pair)
      implicit none
!     Dummy
      integer,intent(in)  :: n
      integer,intent(in)  :: nb(20,n)
      integer,intent(out) ::  pair(n*(n+1)/2)
!     Stack
      integer i,ni,newi,j,newatom,tag,d,i1,ni1,iii,ii,jj,k,lin
      integer,allocatable :: list(:,:),nlist(:,:),nnn(:),nn(:)
      logical da

      allocate(nnn(n),nn(n),list(5*n,n),nlist(5*n,n))

      nn(1:n)=nb(20,1:n)

      pair=0
      list=0
      do i=1,n
         ni=nn(i)
         list(1:ni,i)=nb(1:ni,i)
      enddo

      nlist=list

      pair=0
      do i=1,n
         do j=1,nb(20,i)
            k=nb(j,i)
            pair(lin(k,i))=1
         enddo
      enddo

!     one bond, tag=1
      tag=1
!     call pairsbond(n,nn,list,pair,tag)

!     determine up to 3 bonds in between
      do d=1,2

!     loop over atoms
      do i=1,n
         ni=nn(i)
         newi=ni
!        all neighbors of i
         do ii=1,ni
            i1=list(ii,i)
            ni1=nb(20,i1)
!           all neighbors of neighbors of i
            do iii=1,ni1
               newatom=nb(iii,i1)
               da=.false.
               do j=1,newi
                  if(newatom.eq.list(j,i))da=.true.
               enddo
               if(.not.da)then
                 newi=newi+1
                 nlist(newi,i)=newatom
               endif
            enddo
         enddo
         nnn(i)=newi
      enddo

      list=nlist
      nn  =nnn

!     one bond more
      tag=tag+1
      call pairsbond(n,nn,list,pair,tag)

      enddo
      do i=1,n
         do j=1,i
            if(i.ne.j.and.pair(lin(j,i)).eq.0) pair(lin(j,i))=5
         enddo
      enddo

      end subroutine nbondmat

      subroutine pairsbond(n,nn,list,pair,tag)
      implicit none
      integer n,nn(n),list(5*n,n),tag
      integer i,j,k,ni,nj,ii,jj,ij
      integer pair(n*(n+1)/2)
      logical dai,daj

      do i=1,n
         ni=nn(i)
         ij= i*(i-1)/2
         do j=1,i-1
            k = ij+j
            nj=nn(j)
            dai=.false.
            daj=.false.
            do ii=1,ni
               if(list(ii,i).eq.j)daj=.true.
            enddo
            do jj=1,nj
               if(list(jj,j).eq.i)dai=.true.
            enddo
            if(dai.and.daj.and.pair(k).eq.0) then
               pair(k)=tag
            endif
         enddo
      enddo

      end subroutine pairsbond

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function pilist(ati)
      integer ati
      pilist=.false.
!     if(ati.eq.5.or.ati.eq.6.or.ati.eq.7.or.ati.eq.8.or.ati.eq.9.or.ati.eq.16) pilist=.true.
      if(ati.eq.5.or.ati.eq.6.or.ati.eq.7.or.ati.eq.8.or.ati.eq.9.or.ati.eq.16.or.ati.eq.17) pilist=.true.
      end function pilist

      logical function nofs(ati)
      integer ati
      nofs=.false.
!     if(ati.eq.7.or.ati.eq.8.or.ati.eq.9.or.ati.eq.16) nofs=.true.
      if(ati.eq.7.or.ati.eq.8.or.ati.eq.9.or.ati.eq.16.or.ati.eq.17) nofs=.true.
      end function nofs

      logical function xatom(ati)
      integer ati
      xatom=.false.
      if(ati.eq.17.or.ati.eq.35.or.ati.eq.53.or.&  ! X in A-X...B
     &   ati.eq.16.or.ati.eq.34.or.ati.eq.52.or.&
     &   ati.eq.15.or.ati.eq.33.or.ati.eq.51) xatom=.true.
      end function xatom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer function ctype(n,at,nb,pi,a)
      integer n,a,at(n),nb(20,n),pi(n)
      integer i,no,j

      ctype = 0 ! don't know

      no=0
      do i=1,nb(20,a)
         j=nb(i,a)
         if(at(j).eq.8.and.pi(j).ne.0) no = no + 1
      enddo

      if(no.eq.1.and.pi(a).ne.0) ctype = 1 ! a C=O carbon

      end function ctype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function alphaCO(n,at,hyb,nb,pi,a,b)
      integer n,a,b,at(n),hyb(n),nb(20,n),pi(n)
      integer i,j,no,nc

      alphaCO = .false. 
      if(pi(a).ne.0.and.hyb(b).eq.3.and.at(a).eq.6.and.at(b).eq.6) then  
         no = 0
         do i=1,nb(20,a)
            j=nb(i,a)
            if(at(j).eq. 8.and.pi(j).ne.0.and.nb(20,j).eq.1) no = no +1 ! a pi =O on the C?
         enddo
         if(no.eq.1) then
            alphaCO = .true.
            return
         endif
      endif
      if(pi(b).ne.0.and.hyb(a).eq.3.and.at(b).eq.6.and.at(a).eq.6) then  
         no = 0
         do i=1,nb(20,b)
            j=nb(i,b)
            if(at(j).eq. 8.and.pi(j).ne.0.and.nb(20,j).eq.1) no = no +1 ! a pi =O on the C?
         enddo
         if(no.eq.1) then
            alphaCO = .true.
            return
         endif
      endif
 
      end function alphaCO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function amide(n,at,hyb,nb,pi,a)
      integer n,a,at(n),hyb(n),nb(20,n),pi(n)
      integer i,j,no,nc,ic

      amide = .false. ! don't know
      if(pi(a).eq.0.or.hyb(a).ne.3.or.at(a).ne.7) return

      nc=0
      no=0
      do i=1,nb(20,a)
         j=nb(i,a)
         if(at(j).eq.6.and.pi(j).ne.0) then  ! a pi C on N?
            nc = nc + 1
            ic = j
         endif
      enddo

      if(nc.eq.1)then
         do i=1,nb(20,ic)
            j=nb(i,ic)
            if(at(j).eq. 8.and.pi(j).ne.0.and.nb(20,j).eq.1) no = no +1 ! a pi =O on the C?
!           if(at(j).eq.16.and.pi(j).ne.0.and.nb(20,j).eq.1) no = no +1 ! a pi =S on the C?
         enddo
      endif

      if( no .eq. 1 ) amide = .true.

      end function amide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function amideH(n,at,hyb,nb,pi,a)
      integer n,a,at(n),hyb(n),nb(20,n),pi(n)
      integer i,j,nc,nn
      !logical amide

      amideH = .false. ! don't know
      if(nb(20,a).ne.1               )  return
      nn=nb(1,a)       ! the N
      if(.not.amide(n,at,hyb,nb,pi,nn)) return

      nc=0
      do i=1,nb(20,nn)
         j=nb(i,nn)
         if(at(j).eq.6.and.hyb(j).eq.3) then  ! a sp3 C on N?
            nc = nc + 1
         endif
      enddo

      if( nc .eq. 1 ) amideH = .true.

      end function amideH

end module xtb_gfnff_ini2
