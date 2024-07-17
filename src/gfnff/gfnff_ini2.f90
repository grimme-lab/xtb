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
   use xtb_gfnff_neighbor
   use xtb_type_molecule, only : TMolecule
   implicit none
   private
   public :: gfnff_neigh, getnb, nbondmat, nbondmat_pbc
   public :: pairsbond, pilist, nofs, xatom, ctype, amide, amideH, alphaCO
   public :: ringsatom, ringsbond, ringsbend, ringstors, ringstorl
   public :: chktors, chkrng, hbonds, getring36, ssort, goedeckera, goedeckera_PBC, qheavy
   public :: gfnff_hbset, gfnff_hbset0, bond_hbset, bond_hbset0
   public :: bond_hb_AHB_set, bond_hb_AHB_set1, bond_hb_AHB_set0

contains

subroutine gfnff_neigh(env,makeneighbor,natoms,at,xyz,rab,fq,f_in,f2_in,lintr, &
                      & mchar,hyb,itag,param,topo,mol,neigh,nb_call)
      use xtb_gfnff_param
      implicit none
      character(len=*), parameter :: source = 'gfnff_ini2_neigh'
      type(TEnvironment), intent(inout) :: env
      type(TGFFData), intent(in) :: param
      type(TGFFTopology), intent(inout) :: topo
      type(TMolecule), intent(in) :: mol
      type(TNeigh), intent(inout) :: neigh ! contains nb, nbf and nbm
      logical, intent(in) :: makeneighbor, nb_call
      integer at(natoms),natoms
      integer hyb (natoms)
      integer itag(natoms)
      real*8  rab   (natoms*(natoms+1)/2)
      real*8  xyz(3,natoms)
      real*8  mchar(natoms)
      real*8  fq
      real*8  f_in,f2_in               ! radius scaling for atoms/metal atoms recpectively
      real*8  lintr                    ! threshold for linearity

      logical etacoord,da,strange_iat,metal_iat
      integer,allocatable :: nbdum(:,:,:), nbdum2(:,:), locarr(:,:)
      real*8 ,allocatable :: cn(:),rtmp(:)
      integer iat,i,j,k,ni,ii,jj,kk,ll,lin,ati,nb20i,nbdiff,hc_crit,nbmdiff,nnf,nni,nh,nm
      integer ai,aj,nn,im,ncm,l,no, iTr, iTr2, numnbf, numnbm, numnb, idx, idxdum, idxdum2, numctr
      real*8 r,pi,a1,f,f1,phi,f2,rco,fat(103)
      data pi/3.1415926535897932384626433832795029d0/
      data fat   / 103 * 1.0d0 /

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

      allocate(cn(natoms),rtmp(natoms*(natoms+1)/2),nbdum2(20,natoms))
      rtmp = 0.0

! determine the neighbor list
      if(makeneighbor) then

        do i=1,natoms
           cn(i)=dble(param%normcn(at(i)))
        enddo
        call gfnffrab(natoms,at,cn,rtmp) ! guess RAB based on "normal" CN
        do i=1,natoms
           ai=at(i)
           f1=fq
           if(param%metal(ai) > 0) f1 = f1 * 2.0d0
           do j=1,i
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

        call neigh%get_nb(mol, rab, rtmp, mchar, 1, f_in, f2_in, param) ! nbf
        ! neigh%nb only used for hyb states, then overwritten with nbf
        call neigh%get_nb(mol, rab, rtmp, mchar, 2, f_in, f2_in, param) ! nb 
        call neigh%get_nb(mol, rab, rtmp, mchar, 3, f_in, f2_in, param) ! nbm

      ! take the input
      else

         neigh%nbf = neigh%nb
         neigh%nbm = neigh%nb 

      endif
! done

      itag = 0 ! save special hyb info
      numctr = neigh%numctr ! number of central cells considered (e.g. 1 for molec case)

! tag atoms in nb(19,i) if they belong to a cluster (which avoids the ring search)
      do i=1,natoms
         if(sum(neigh%nbf(neigh%numnb,i,:)).eq.0.and.param%group(at(i)).ne.8)then
            write(env%unit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
            write(env%unit,'(''  warning: no bond partners for atom'',i4)')i
            write(env%unit,'(''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'')')
         endif
         if(at(i).lt.11.and.sum(neigh%nbf(neigh%numnb,i,:)).gt.2)then
           do iTr=1, numctr 
             do k=1,neigh%nbf(neigh%numnb,i,iTr)
                 kk=neigh%nbf(k,i,iTr)
                 if(param%metal(at(kk)).ne.0.or.sum(neigh%nb(neigh%numnb,kk,:)).gt.4) then
                    neigh%nb(neigh%numnb-1,i,1) =1  ! ring search is limited to unit cell.
                    neigh%nbf(neigh%numnb-1,i,1)=1  ! Assumption: If the conditions are true
                    neigh%nbm(neigh%numnb-1,i,1)=1  ! in one cell they are true in all cells
                 endif
             enddo
           enddo
         endif
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! hybridization states
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(.not. allocated(nbdum)) &
        & allocate(nbdum(neigh%numnb, mol%n, numctr), source=0)
      do i=1,natoms
         ati  = at(i)
         numnbf=sum(neigh%nbf(neigh%numnb,i,:))
         numnbm=sum(neigh%nbm(neigh%numnb,i,:))
!        important: determine cases where atom is pi bonded to a metal and thus
!        the hyb must be obtained from the reduced (wo metals) neighbor list
         etacoord=.false.
         if(ati.le.10)then
            if(ati.eq.6.and.numnbf.ge.4.and.numnbm.eq.3) etacoord=.true.  ! CP case
            if(ati.eq.6.and.numnbf.eq.3.and.numnbm.eq.2) etacoord=.true.  ! alkyne case
            nm=0
           do iTr=1, numctr 
            do k=1, neigh%nbf(neigh%numnb,i,iTr)  ! how many metals ? and which
               kk=neigh%nbf(k,i,iTr)
               if(param%metal(at(kk)).ne.0) then
                  nm=nm+1
                  im=kk
               endif
            enddo
           enddo
           if(nm.eq.0) then
               etacoord=.false.  ! etacoord makes no sense without metals!
            elseif(nm.eq.1)then  ! distinguish M-CR2-R i.e. not an eta coord.
              ncm=0
              do iTr=1, numctr
                do k=1,numnbf  !
                  if(neigh%nbf(k,i,iTr).ne.im)then ! all neighbors that are not the metal im
                    kk=neigh%nbf(k,i,iTr)
                    do l=1,sum(neigh%nbf(neigh%numnb,kk,:))
                      if(neigh%nbf(l,kk,iTr).eq.im) ncm=ncm+1 ! ncm=1 is alkyne, =2 is cp
                    enddo
                  endif
                enddo
              enddo
              if(ncm.eq.0) etacoord=.false.
            endif
         endif
         if(etacoord)then
          itag(i)=-1
          nbdum(:,i,:)=neigh%nbm(:,i,:)
         else
          nbdum(:,i,:)=neigh%nbf(:,i,:) ! take full set of neighbors by default
        endif
      enddo

      do i=1,natoms
         ati  = at(i)
         hyb(i)=0    ! don't know it 
         numnbm=sum(neigh%nbm(neigh%numnb,i,:))
         numnbf=sum(neigh%nbf(neigh%numnb,i,:))
         numnb =sum(neigh%nb (neigh%numnb,i,:))
         nbdiff =numnbf-numnb
         nbmdiff=numnbf-numnbm
         nb20i =sum(nbdum(neigh%numnb,i,:))
         ! go over all nb's of i, count number nb's that are H or O 
         nh=0
         no=0
         do iTr=1, numctr
           do j=1,nb20i
             if(nbdum(j,i,iTr).eq.0) cycle
             if(at(nbdum(j,i,iTr)).eq.1) nh=nh+1
             if(at(nbdum(j,i,iTr)).eq.8) no=no+1
           enddo
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
              ! get the right neighbors for bangl call
              call neigh%nbLoc(natoms, nbdum, i, locarr)
              if (size(locarr,dim=2).eq.1) then
                idxdum  = locarr(1,1)
                idxdum2 = locarr(2,1)
                iTr = locarr(neigh%numnb,1)
                iTr2 = locarr(neigh%numnb,1)
                deallocate(locarr)
              elseif(size(locarr,dim=2).eq.2) then
                idxdum  = locarr(1,1)
                idxdum2 = locarr(1,2)
                iTr = locarr(neigh%numnb,1)
                iTr2 = locarr(neigh%numnb,2)
                deallocate(locarr)
              else
                call env%error(' Hybridization failed. Neighbors could not be located.', source)
              endif
              call banglPBC(1,xyz,idxdum,i,idxdum2,iTr,iTr2,neigh,phi)
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
               do iTr=1, numctr
                 do j=1,3
                   jj= nbdum(j,i,iTr)
                   if(jj.eq.0) exit ! if there is no 1st nb there is no nb at all
                   if(at(jj).eq. 8.and.sum(neigh%nb(neigh%numnb,jj,:)).eq.1) kk=kk+1 ! check for NO2 or R2-N=O
                   if(at(jj).eq. 5.and.sum(neigh%nb(neigh%numnb,jj,:)).eq.4) ll=ll+1 ! check for B-N, if the CN(B)=4 the N is loosely bound and sp2
                   if(at(jj).eq.16.and.sum(neigh%nb(neigh%numnb,jj,:)).eq.4) nn=nn+1 ! check for N-SO2-
                 enddo
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
               ! get the two neighbors for bangl call 
               !>> idxdum and idxdum2 are the indices of the two neighbors
               call neigh%nbLoc(natoms, nbdum, i, locarr)
               if (size(locarr,dim=2).eq.1) then
                 idxdum  = locarr(1,1)
                 idxdum2 = locarr(2,1)
                 iTr = locarr(neigh%numnb,1)
                 iTr2 = locarr(neigh%numnb,1)
                 deallocate(locarr)
               elseif(size(locarr,dim=2).eq.2) then
                 idxdum  = locarr(1,1)
                 idxdum2 = locarr(1,2)
                 iTr = locarr(neigh%numnb,1)
                 iTr2 = locarr(neigh%numnb,2)
                 deallocate(locarr)
               else
                 call env%error(' Hybridization failed. Neighbors could not be located.', source)
               endif
               call banglPBC(1,xyz,idxdum,i,idxdum2,iTr,iTr2,neigh,phi)
               jj=idxdum
               kk=idxdum2
               if(sum(nbdum(neigh%numnb,jj,:)).eq.1.and.at(jj).eq.6)     hyb(i)=1  ! R-N=C
               if(sum(nbdum(neigh%numnb,kk,:)).eq.1.and.at(kk).eq.6)     hyb(i)=1  ! R-N=C
               if(sum(nbdum(neigh%numnb,jj,:)).eq.1.and.at(jj).eq.7)     hyb(i)=1  ! R-N=N in e.g. diazomethane
               if(sum(nbdum(neigh%numnb,kk,:)).eq.1.and.at(kk).eq.7)     hyb(i)=1  ! R-N=N in e.g. diazomethane
               if(idxdum.gt.0.and.param%metal(at(idxdum)).gt.0) hyb(i)=1 ! M-NC-R in e.g. nitriles
               if(idxdum2.gt.0.and.param%metal(at(idxdum2)).gt.0) hyb(i)=1 ! M-NC-R in e.g. nitriles
               if(at(jj).eq.7.and.at(kk).eq.7.and. &
     &         sum(nbdum(neigh%numnb,jj,:)).le.2.and.sum(nbdum(neigh%numnb,kk,:)).le.2) hyb(i)=1  ! N=N=N
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
              call nn_nearest_noM(i,natoms,at,mol%xyz,neigh,rab,j,param) ! CN of closest non-M atom
                                        if(j.eq.3)       hyb(i)=2 ! M-O-X konj
                                        if(j.eq.4)       hyb(i)=3 ! M-O-X non
            endif
            if(nb20i.eq.1)                               hyb(i)=2
            if(nb20i.eq.1.and.nbdiff.eq.0) then
            call neigh%nbLoc(natoms,neigh%nb,i,locarr)     
            iTr = locarr(neigh%numnb,1)
            deallocate(locarr)
            if(sum(neigh%nb(neigh%numnb,neigh%nb(1,i,iTr),:)).eq.1)      hyb(i)=1 ! CO
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

      neigh%nb = nbdum ! list is complete but hyb determination is based only on reduced (without metals) list
      deallocate(nbdum)

      j = 0
      do i=1,natoms
        numnb = sum(neigh%nb(neigh%numnb,i,:))
        if(numnb.gt.12) j = j +1
        do iTr=1, neigh%numctr
          do k=1, neigh%nb(neigh%numnb,i,iTr)
            kk=neigh%nb(k,i,iTr)
            if(at(kk).eq.6.and.at(i).eq.6.and.itag(i).eq.1.and.itag(kk).eq.1) then ! check the very special situation of
              itag(i) =0                                                           ! two carbene C bonded which is an arine
              itag(kk)=0
            endif
          enddo
        enddo
      enddo
      if(dble(j)/dble(natoms).gt.0.3.and.nb_call) then
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

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! find the CN of the non metal atom that is closest to atom i
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine nn_nearest_noM(ii,n,at,xyz,neigh,r,nn,param)
      use xtb_mctc_accuracy, only : wp
      implicit none
      type(TGFFData), intent(in) :: param
      type(TNeigh), intent(in) :: neigh
      integer, intent(in) :: ii,n,at(n)
      real(wp), intent(in) :: xyz(3,n)
      integer, intent(inout) :: nn
      real*8, intent(in) :: r(n*(n+1)/2)

      integer jmin,j,jj,lin, numnb, iTr
      real(wp) :: dist
      real*8 rmin

      numnb=neigh%numnb
      nn=0
      rmin=1.d+42
      jmin=0
      dist=0.0
      do iTr=1, neigh%numctr
        do j=1,neigh%nb(numnb,ii,iTr)
          jj=neigh%nb(j,ii,iTr)
          if(param%metal(at(jj)).ne.0) cycle
          dist = NORM2(xyz(:,ii)-(xyz(:,jj)+neigh%transVec(:,iTr)))
          if(dist.lt.rmin)then
            rmin=dist
            jmin=jj
          endif
        enddo
      enddo

      if(jmin.gt.0) nn=sum(neigh%nb(numnb,jmin,:))

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
      do k=1,s(20,j)    ! all rings of atom j
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
         endif
      enddo
      do m=1,s(20,j)    ! all rings of atom j
         itest=0
         do a=1,s(m,j)  ! all atoms of ring m
            if(c(a,m,j).eq.i.or.c(a,m,j).eq.k.or.c(a,m,j).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,j).lt.rings2) then
            rings2=s(m,j)
         endif
      enddo
      do m=1,s(20,k)    ! all rings of atom k
         itest=0
         do a=1,s(m,k)  ! all atoms of ring m
            if(c(a,m,k).eq.i.or.c(a,m,k).eq.j.or.c(a,m,k).eq.l) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,k).lt.rings3) then
            rings3=s(m,k)
         endif
      enddo
      do m=1,s(20,l)    ! all rings of atom k
         itest=0
         do a=1,s(m,l)  ! all atoms of ring m
            if(c(a,m,l).eq.i.or.c(a,m,l).eq.k.or.c(a,m,l).eq.j) itest=itest+1
         enddo
         if(itest.eq.3.and.s(m,l).lt.rings4) then
            rings4=s(m,l)
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

      logical function chktors(n,xyz,i,j,k,l,iTrj,iTrk,iTrl,neigh)  ! true if dihedral angle is bad i.e. near 180
      use xtb_gfnff_neighbor
      implicit none
      type(TNeigh), intent(in) :: neigh
      integer n,i,j,k,l,iTrj,iTrk,iTrl
      real*8 xyz(3,n),phi

      chktors=.true.

      call banglPBC(1,xyz,j,i,k,iTrj,iTrk,neigh,phi)
      if(phi*180./3.1415926d0.gt.170.0d0) return
      call banglPBC(2,xyz,i,j,l,iTrj,iTrl,neigh,phi)
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

subroutine gfnff_hbset(n,at,xyz,topo,neigh,nlist,hbthr1,hbthr2)
      use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(inout) :: neigh 
      type(TGFFNeighbourList), intent(inout) :: nlist
      integer n
      integer at(n)
      real(wp) xyz(3,n)
      real(wp), intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh, iTri,iTrj,iTrDum
      real(wp) rab,rmsd, rih,rjh
      logical ijnonbond,free


      rmsd = sqrt(sum((xyz-nlist%hbrefgeo)**2))/dble(n)

      if(rmsd.lt.1.d-6 .or. rmsd.gt. 0.3d0) then ! update list if first call or substantial move occured
      nlist%nhb1=0
      nlist%nhb2=0
      ! loop over hb-relevant AB atoms
      do ix=1,topo%nathbAB  
        i=topo%hbatABl(1,ix) 
        j=topo%hbatABl(2,ix)
        do iTri=1, neigh%nTrans ! go through i shifts
          do iTrj=1, neigh%nTrans ! go through j shifts
            ! get adjustet iTr -> for use of neigh% distances and bpair with two shifted atoms
            iTrDum=neigh%fTrSum(neigh%iTrNeg(iTri),iTrj) 
            if(iTrDum.gt.neigh%nTrans.or.iTrDum.lt.-1.or.iTrDum.eq.0) cycle ! cycle nonsense 
            rab=NORM2((xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
            if(rab.gt.hbthr1) cycle
            ! check if ij bonded
            if(iTrDum.le.neigh%numctr.and.iTrDum.gt.0) then
              ijnonbond=neigh%bpair(j,i,iTrDum).ne.1
            else
              ! i and j are not in neighboring cells
              ijnonbond=.true. 
            endif
            ! loop over relevant H atoms
            do k=1,topo%nathbH
              free=.true. ! tripplet not assigned yet
              nh  =topo%hbatHl(1,k) ! nh always in central cell
              ! distances for non-cov bonded case
              rih=NORM2(xyz(:,nh)-(xyz(:,i)+neigh%transVec(:,iTri)))**2
              rjh=NORM2(xyz(:,nh)-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
              ! check if i is the bonded A    
              if(iTri.le.neigh%numctr) then ! nh is not shifted so bpair works without adjustment 
                if(neigh%bpair(i,nh,iTri).eq.1.and.ijnonbond) then
                  nlist%nhb2=nlist%nhb2+1
                  nlist%hblist2(1,nlist%nhb2)=i
                  nlist%hblist2(2,nlist%nhb2)=j
                  nlist%hblist2(3,nlist%nhb2)=nh
                  nlist%hblist2(4,nlist%nhb2)=iTri
                  nlist%hblist2(5,nlist%nhb2)=iTrj
                  free=.false. ! not available for nhb1 !!!
                endif 
              endif  
              ! check if j is the bonded A
              if(iTrj.le.neigh%numctr.and.free) then
                if(neigh%bpair(j,nh,iTrj).eq.1.and.ijnonbond) then
                  nlist%nhb2=nlist%nhb2+1
                  nlist%hblist2(1,nlist%nhb2)=j
                  nlist%hblist2(2,nlist%nhb2)=i
                  nlist%hblist2(3,nlist%nhb2)=nh
                  nlist%hblist2(4,nlist%nhb2)=iTrj
                  nlist%hblist2(5,nlist%nhb2)=iTri
                  free=.false. ! not available for nhb1 !!!
                endif
              endif  
              ! check for non-cov bonded A  
              if(rab+rih+rjh.lt.hbthr2.and.free) then ! sum of rAB,rAH,rBH is below threshold
                nlist%nhb1=nlist%nhb1+1
                nlist%hblist1(1,nlist%nhb1)=i 
                nlist%hblist1(2,nlist%nhb1)=j 
                nlist%hblist1(3,nlist%nhb1)=nh
                nlist%hblist1(4,nlist%nhb1)=iTri
                nlist%hblist1(5,nlist%nhb1)=iTrj
              endif
            enddo ! k: relevant H atoms
          enddo ! iTrj
        enddo ! iTri
      enddo ! ix: relevant AB atoms 

      ! for nxb list only i is not shifted
      nlist%nxb =0
      do ix=1,topo%natxbAB
          i =topo%xbatABl(1,ix)   ! A
          j =topo%xbatABl(2,ix)   ! B
          iTrj=topo%xbatABl(4,ix) ! iTrB
          if(iTrj.gt.neigh%nTrans.or.iTrj.lt.-1.or.iTrj.eq.0) cycle ! cycle nonsense 
          rab=NORM2(xyz(:,j)-xyz(:,i)+neigh%transVec(:,iTrj))**2
          if(rab.gt.hbthr2)cycle
          nlist%nxb=nlist%nxb+1
          nlist%hblist3(1,nlist%nxb)=i                  ! A
          nlist%hblist3(2,nlist%nxb)=j                  ! B
          nlist%hblist3(3,nlist%nxb)=topo%xbatABl(3,ix) ! X
          nlist%hblist3(4,nlist%nxb)=iTrj
          nlist%hblist3(5,nlist%nxb)=topo%xbatABl(5,ix) ! iTrX
        !enddo
      enddo

      nlist%hbrefgeo = xyz

      endif  ! else do nothing


      end subroutine gfnff_hbset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bond_hbset(n,at,xyz,npbc,bond_hbn,bond_hbl,topo,neigh,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(inout) :: neigh
      integer,intent(in) :: n
      integer,intent(in) :: at(n)
      integer,intent(in) :: bond_hbn
      integer,intent(out) :: bond_hbl(6,bond_hbn)  ! A-H...B Tripletts
      real(wp),intent(in) :: xyz(3,n)
      integer, intent(in) :: npbc
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,iTri,iTrj,iTrH,iTrDum
      integer bond_nr
      real(wp) dum, rab
      logical ijnonbond

      bond_nr=0
      bond_hbl=0
      do ix=1,topo%nathbAB  ! loop over hb-relevant AB atoms
        !get i and j from AB list, it is not decided who is A and who is B yet
        i=topo%hbatABl(1,ix)
        j=topo%hbatABl(2,ix)
        do iTri=1,neigh%nTrans
          do iTrj=1,neigh%nTrans
            iTrDum=neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
            rab=NORM2((xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
            if(rab.gt.hbthr1) cycle
            ! check if ij bonded
            if(iTrDum.le.neigh%numctr.and.iTrDum.gt.0) then
                ijnonbond=neigh%bpair(j,i,iTrDum).ne.1
            else
                ! i and j are not in neighboring cells
                ijnonbond=.true.
            endif
            ! loop over relevant H atoms
            do k=1,topo%nathbH  
              nh  =topo%hbatHl(1,k)  ! ALWAYS in central cell
              !i is the A (only possible if A is adjacent to central cell, since H always there)
              if(iTri.le.neigh%numctr) then ! check if a in or adjacent to central cell
                if(neigh%bpair(i,nh,iTri).eq.1.and.ijnonbond) then
                  bond_nr=bond_nr+1
                  bond_hbl(1,bond_nr)=i     ! 1=A
                  bond_hbl(2,bond_nr)=j     ! 2=B
                  bond_hbl(3,bond_nr)=nh    ! 3=H
                  bond_hbl(4,bond_nr)=iTri  ! 4=iTrA
                  bond_hbl(5,bond_nr)=iTrj  ! 5=iTrB
                  bond_hbl(6,bond_nr)=1     ! 6=iTrH
                endif
              endif  
              ! j is the A
              if(iTrj.le.neigh%numctr) then
                if(neigh%bpair(j,nh,iTrj).eq.1.and.ijnonbond)then
                  bond_nr=bond_nr+1
                  bond_hbl(1,bond_nr)=j     ! 1=A
                  bond_hbl(2,bond_nr)=i     ! 2=B
                  bond_hbl(3,bond_nr)=nh    ! 3=H
                  bond_hbl(4,bond_nr)=iTrj  ! 4=iTrA
                  bond_hbl(5,bond_nr)=iTri  ! 5=iTrB
                  bond_hbl(6,bond_nr)=1     ! 6=iTrH
                endif
              endif  
            enddo
          enddo
        enddo
      enddo        

end subroutine bond_hbset


subroutine bond_hbset0(n,at,xyz,npbc,bond_hbn,topo,neigh,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(inout) :: neigh
      integer,intent(in) :: n
      integer,intent(in) :: at(n)
      integer,intent(out) :: bond_hbn
      real(wp),intent(in) :: xyz(3,n)
      integer, intent(in) :: npbc
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ix,iTri,iTrj,iTrH,iTrDum
      !integer bond_nr
      real(wp) rab
      logical ijnonbond
      !
      real(wp) :: vec1, vec2

      ! get size of bond_hbl, for more comments see bond_hbset
      bond_hbn=0 
      do ix=1,topo%nathbAB  ! loop over hb-relevant AB atoms
        !get i and j from AB list, it is not decided who is A and who is B yet
        i=topo%hbatABl(1,ix)
        j=topo%hbatABl(2,ix)
        do iTri=1,neigh%nTrans
          do iTrj=1,neigh%nTrans
            iTrDum=neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
            if(iTrDum.eq.-1.or.iTrDum.gt.neigh%numctr) then
              rab=NORM2((xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
              if(rab.gt.hbthr1) cycle
              ijnonbond=.true. ! i and j are not in neighboring or same cell for iTrDum=-1
            else
              ! check ij distance
              rab=NORM2(xyz(:,i)-(xyz(:,j)+neigh%transVec(:,iTrDum)))**2
              if(rab.gt.hbthr1) cycle
              ! check if ij bonded
              ijnonbond=neigh%bpair(j,i,iTrDum).ne.1
            endif
            ! loop over relevant H atoms
            do k=1,topo%nathbH  
              nh  =topo%hbatHl(1,k)  ! ALWAYS in central cell
              ! i is the A
              if(iTri.le.neigh%numctr) then
                if(neigh%bpair(i,nh,iTri).eq.1.and.ijnonbond) then
                  bond_hbn = bond_hbn + 1
                endif
              endif  
              ! j is the A
              if(iTrj.le.neigh%numctr) then
                if(neigh%bpair(j,nh,iTrj).eq.1.and.ijnonbond)then!
                  bond_hbn = bond_hbn + 1
                endif
              endif  
            enddo
          enddo
        enddo
      enddo       

end subroutine bond_hbset0


subroutine bond_hb_AHB_set(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,topo,neigh)
      use xtb_gfnff_param
      implicit none
      !Dummy
      type(TNeigh), intent(inout) :: neigh
      type(TGFFTopology), intent(inout) :: topo
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(6,bond_hbn)
      integer,intent(in)  :: tot_AHB_nr
      integer,intent(inout) :: lin_AHB(4,0:tot_AHB_nr)
      !Stack
      integer :: i,j
      integer :: ii,jj,iTr,iTrA,iTrH,iTrB
      integer :: ia,ja
      integer :: lin
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB
      integer :: tot_count,t
      integer :: AH_count
      integer :: B_count
      integer :: lin_diff
      logical :: iTr_same, AH_diff, ldum(4)

      tot_count=0
      AH_count=0
      B_count=0
      lin_diff=0

      do i=1,numbond
         jj=neigh%blist(1,i)
         ii=neigh%blist(2,i)
         iTr=neigh%blist(3,i)
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
               iTrB= bond_hbl(5,j)
               atB = at(Bat)
               Aat  = bond_hbl(1,j)
               iTrA = bond_hbl(4,j)
               Hat  = bond_hbl(3,j)
               iTrH = bond_hbl(6,j) ! always 1
               if ((hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrA.and.iTrH.eq.1).or.&
                   &hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrH.and.iTrA.eq.1) then
                  if (atB.eq.7.or.atB.eq.8) then
                     tot_count = tot_count + 1
                     t=tot_count
                     !lin_AHB(1,tot_count) = lin(hbA,hbH)
                     lin_AHB(1,tot_count) = hbA
                     lin_AHB(2,tot_count) = hbH
                     lin_AHB(3,tot_count) = iTrA
                     lin_AHB(4,tot_count) = iTrH
                     if (lin_AHB(1,tot_count)-lin_AHB(1,tot_count-1).eq.0.and.&
                        &lin_AHB(2,tot_count)-lin_AHB(2,tot_count-1).eq.0) then
                       lin_diff=0
                     else
                       lin_diff=1
                     endif
                     iTr_same=iTrA.eq.lin_AHB(3,tot_count-1).and.iTrH.eq.lin_AHB(4,tot_count-1)
                     if (lin_diff.eq.0.and.iTr_same) B_count = B_count + 1
                     if (lin_diff.ne.0.or..not.iTr_same) then
                        AH_count = AH_count + 1
                        topo%bond_hb_AH(1,AH_count) = hbA
                        topo%bond_hb_AH(2,AH_count) = hbH
                        topo%bond_hb_AH(3,AH_count) = iTrA
                        topo%bond_hb_AH(4,AH_count) = iTrH
                        !Reset B count
                        B_count = 1
                     end if
                     topo%bond_hb_Bn(AH_count) = B_count
                     topo%bond_hb_B(1,B_count,AH_count) = Bat
                     topo%bond_hb_B(2,B_count,AH_count) = iTrB
                  end if
               else
                 cycle
               end if
            end do
         end if
      end do

end subroutine bond_hb_AHB_set

subroutine bond_hb_AHB_set1(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,lin_AHB,AH_count,bmax,topo,neigh)
      use xtb_gfnff_param
      use xtb_gfnff_neighbor
      implicit none
      !Dummy
      type(TNeigh), intent(inout) :: neigh
      type(TGFFTopology), intent(inout) :: topo
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(6,bond_hbn)
      integer,intent(in)  :: tot_AHB_nr
      integer,intent(inout) :: lin_AHB(4,0:tot_AHB_nr)
      integer,intent(out) :: AH_count
      integer,intent(out) :: bmax
      !Stack
      integer :: i,j
      integer :: ii,jj,iTr,iTrA,iTrH,iTrB
      integer :: ia,ja
      integer :: lin
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB
      integer :: tot_count,t
      integer :: B_count
      integer :: lin_diff
      logical :: iTr_same, AH_diff, ldum(4)

      tot_count=0
      AH_count=0
      B_count=1
      bmax=1
      lin_diff=0

      do i=1,numbond ! loop over blist to get same order as blist
         jj=neigh%blist(1,i)
         ii=neigh%blist(2,i)
         iTr=neigh%blist(3,i)
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
               iTrA = bond_hbl(4,j)
               Hat  = bond_hbl(3,j)
               iTrH = bond_hbl(6,j) ! always 1
               if ((hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrA.and.iTrH.eq.1).or.&
                   &(hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrH.and.iTrA.eq.1)) then
                 if (atB.eq.7.or.atB.eq.8) then
                     tot_count = tot_count + 1
                     lin_AHB(1,tot_count) = hbA
                     lin_AHB(2,tot_count) = hbH
                     lin_AHB(3,tot_count) = iTrA
                     lin_AHB(4,tot_count) = iTrH
                     topo%isABH(Bat)=.true.
                     if (lin_AHB(1,tot_count)-lin_AHB(1,tot_count-1).eq.0.and.&
                        &lin_AHB(2,tot_count)-lin_AHB(2,tot_count-1).eq.0) then
                       lin_diff=0
                     else
                       lin_diff=1
                     endif
                     iTr_same=iTrA.eq.lin_AHB(3,tot_count-1).and.iTrH.eq.lin_AHB(4,tot_count-1)
                     if (lin_diff.eq.0.and.iTr_same) B_count = B_count + 1
                     if (lin_diff.ne.0.or..not.iTr_same) then
                        AH_count = AH_count + 1
                        B_count = 1
                     end if
                     if (B_count.gt.bmax) bmax = B_count
                  end if
               else
                 cycle
               end if
            end do
            topo%isABH(hbA)=.true.
            topo%isABH(hbH)=.true.
            neigh%nr_hb(i) = B_count
         end if
      end do

end subroutine bond_hb_AHB_set1

subroutine bond_hb_AHB_set0(n,at,numbond,bond_hbn,bond_hbl,tot_AHB_nr,neigh)
      use xtb_gfnff_param
      implicit none
      !Dummy
      !type(TGFFTopology), intent(in) :: topo
      type(TNeigh), intent(inout) :: neigh
      integer,intent(in)  :: n
      integer,intent(in)  :: numbond
      integer,intent(in)  :: at(n)
      integer,intent(in)  :: bond_hbn
      integer,intent(in)  :: bond_hbl(6,bond_hbn)
      integer,intent(out) :: tot_AHB_nr
      !Stack
      integer :: i,j,k
      integer :: ii,jj,iTr,iTrA,iTrH,iTrB
      integer :: ia,ja
      integer :: hbH,hbA
      integer :: Hat,Aat
      integer :: Bat,atB

      tot_AHB_nr=0
      do i=1,numbond ! loop over blist to get same order as blist, more important in set1
         jj= neigh%blist(1,i)
         ii= neigh%blist(2,i)
         iTr=neigh%blist(3,i)
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
               iTrA = bond_hbl(4,j)
               Hat  = bond_hbl(3,j)
               iTrH = bond_hbl(6,j) ! always 1
               if ((hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrA.and.iTrH.eq.1).or.&
                   &hbA.eq.Aat.and.hbH.eq.Hat.and.iTr.eq.iTrH.and.iTrA.eq.1) then
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

subroutine gfnff_hbset0(n,at,xyz,topo,nhb1,nhb2,nxb,neigh,nlist,hbthr1,hbthr2)
use xtb_mctc_accuracy, only : wp
      use xtb_gfnff_param
      implicit none
      type(TGFFTopology), intent(in) :: topo
      integer, intent(out) :: nhb1
      integer, intent(out) :: nhb2
      integer, intent(out) :: nxb
      type(TNeigh), intent(inout) :: neigh 
      type(TGFFNeighbourList), intent(inout) :: nlist
      integer n
      integer at(n)
      real(wp) xyz(3,n)
      real(wp),intent(in) :: hbthr1, hbthr2

      integer i,j,k,nh,ia,ix,lin,ij,inh,jnh 
      integer :: iTri,iTrj,iTrDum !iTrA,iTrH,iTrDum2,iTrAH,sw,nsw,shift
      logical ijnonbond, free
      real(wp) rab,rih,rjh
      nhb1=0
      nhb2=0
      ! loop over hb-relevant AB atoms
      do ix=1,topo%nathbAB  
        i=topo%hbatABl(1,ix) 
        j=topo%hbatABl(2,ix)
        do iTri=1, neigh%nTrans ! go through i shifts
          do iTrj=1, neigh%nTrans ! go through j shifts
            ! get adjustet iTr -> for use of neigh% distances and bpair with two shifted atoms
            iTrDum=neigh%fTrSum(neigh%iTrNeg(iTri),iTrj)
            if(iTrDum.gt.neigh%nTrans.or.iTrDum.lt.-1.or.iTrDum.eq.0) cycle 
            if(iTrDum.eq.-1) then
              rab=NORM2((xyz(:,i)+neigh%transVec(:,iTri))-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
              if(rab.gt.hbthr1) cycle
              ijnonbond=.true. ! i and j are not in neighboring or same cell for iTrDum=-1
            else
              ! check
              rab=NORM2(xyz(:,i)-(xyz(:,j)+neigh%transVec(:,iTrDum)))**2
              if(rab.gt.hbthr1) cycle
              ! check if ij bonded
              if(iTrDum.le.neigh%numctr) then
                ijnonbond=neigh%bpair(j,i,iTrDum).ne.1
              else
                ! i and j are not in neighboring cells
                ijnonbond=.true.
              endif
            endif
            ! loop over relevant H atoms
            do k=1,topo%nathbH  
              free=.true.
              nh  =topo%hbatHl(1,k) ! nh always in central cell
              ! distances for non-cov bonded case
              rih=NORM2(xyz(:,nh)-(xyz(:,i)+neigh%transVec(:,iTri)))**2
              rjh=NORM2(xyz(:,nh)-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
              ! check if i is the bonded A    
              if(iTri.le.neigh%numctr) then ! nh is not shifted so bpair works without adjustment 
                if(neigh%bpair(i,nh,iTri).eq.1.and.ijnonbond) then
                  nhb2=nhb2+1
                  free=.false.
                endif  
              endif  
              ! check if j is the bonded A
              if(iTrj.le.neigh%numctr.and.free) then
                if(neigh%bpair(j,nh,iTrj).eq.1.and.ijnonbond) then
                  nhb2=nhb2+1
                  free=.false.
                endif  
              endif  
              ! check for non-cov bonded A  
              if(rab+rih+rjh.lt.hbthr2.and.free) then ! sum of rAB,rAH,rBH is below threshold
                nhb1=nhb1+1
              endif
            enddo ! k: relevant H atoms
          enddo ! iTrj
        enddo ! iTri
      enddo ! ix: relevant AB atoms 
      nxb =0
      do ix=1,topo%natxbAB
        !do iTrj=1, neigh%numctr
          i =topo%xbatABl(1,ix)
          j =topo%xbatABl(2,ix)
          iTrj=topo%xbatABl(4,ix)
          rab=NORM2(xyz(:,i)-(xyz(:,j)+neigh%transVec(:,iTrj)))**2
          if(rab.gt.hbthr2)cycle
          nxb=nxb+1
        !enddo
      enddo

! the actual size can be larger, so make it save
      if(neigh%numctr.gt.1)then
        nhb1=(nhb1*27)
        nhb2=(nhb2*27)
      else
        nhb1=(nhb1*6)
        nhb2=(nhb2*6)
      endif
      if (nxb.gt.1000) then
        nxb =(nxb *25)
      else
        nxb =(nxb *10)
      endif        

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

subroutine getring36(n,at,numnb,numctr,nbin,a0_in,cout,irout)
      implicit none
      integer, intent(out) :: cout(10,20),irout(20)  ! output: atomlist, ringsize, # of rings in irout(20)
      integer, intent(in) :: n,at(n),numnb,numctr,nbin(numnb,n),a0_in
      integer :: i,nb(numnb,n), a0
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
      if(n.le.2.                  .or.nbin(numnb-1,a0_in).eq.1) return

      nn=nbin(numnb,a0_in)

      cdum=0
      kk=0
      do m=1,nn
      nb=nbin
      ! ring search only in unit cell -> adjusted nbin to include neighbors from other cells
      if(nb(m,a0_in).eq.1) cycle
      do i=1,n
         if(nb(numnb,i).eq.1)nb(numnb,i)=0
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
      n0=nb(numnb,a0)

      do i1=1,n0
         a1=nb(i1,a0)
         if(a1.eq.a0) cycle
         n1=nb(numnb,a1)
         do i2=1,n1
            a2=nb(i2,a1)
            if(a2.eq.a1) cycle
            n2=nb(numnb,a2)
            do i3=1,n2
               a3=nb(i3,a2)
               n3=nb(numnb,a3)
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
                  n4=nb(numnb,a4)
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
                     n5=nb(numnb,a5)
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
                        n6=nb(numnb,a6)
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

subroutine goedeckera(env,n,at,pair,q,es,topo)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : mctc_sytrf, mctc_sytrs
   implicit none
   character(len=*), parameter :: source = 'gfnff_ini2_goedeckera'
   type(TEnvironment), intent(inout) :: env
   type(TGFFTopology), intent(in) :: topo
   integer, intent(in)  :: n          ! number of atoms
   integer, intent(in)  :: at(n)      ! ordinal numbers
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
         tmp = erf(gammij*rij)/rij  ! apart from diagonal(=0), if ij non-bonded rij=1.0d+12 
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


!> version of EEQ
subroutine goedeckera_PBC(env,mol,pair,topo,q,es)
   use xtb_mctc_accuracy, only : wp
   use xtb_mctc_lapack, only : mctc_sytrf, mctc_sytrs
   implicit none
   character(len=*), parameter :: source = 'gfnff_ini2_goedeckera'
   type(TEnvironment), intent(inout) :: env
   type(TMolecule), intent(in) :: mol 
   real(wp),intent(in)  :: pair(mol%n*(mol%n+1)/2) !estimate bond distances =1.0d+12 if no bond
   type(TGFFTopology), intent(in) :: topo
   real(wp),intent(out) :: q(mol%n)       ! output charges
   real(wp),intent(out) :: es         ! ES energy

!  local variables
   logical :: exitRun
   integer  :: m,i,j,k,l,ii,jj,kk,n
   integer  :: ij,lj
   integer,allocatable :: ipiv(:)

   real(wp) :: gammij,sief1,sief2
   real(wp) :: r2,r0
   real(wp) :: rij
   real(wp) :: tsqrt2pi,bohr
   real(wp) :: tmp
   real(wp),allocatable :: A (:,:)
   real(wp),allocatable :: x(:)
   real(wp), allocatable :: rTrans(:,:)
   real(wp), allocatable :: gTrans(:,:)
   real(wp) :: vec(3)
   integer :: iRp,iT1,iT2,iT3,iG1,iG2,iG3
   integer, parameter :: ewaldCutD(3) = 2
   integer, parameter :: ewaldCutR(3) = 2
   real(wp), parameter :: sqrtpi = 1.772453850905516_wp 
   real(wp) :: cf !convergence factor                   
!  parameter
   parameter (tsqrt2pi = 0.797884560802866_wp)
   
   n=mol%n

   m=n+topo%nfrag ! # atoms frag constrain
   allocate(A(m,m),x(m),ipiv(m))

   iRp = 0                                                                  
   allocate(gTrans(3, product(2*ewaldCutR+1)-1))                            
   do iG1 = -ewaldCutR(1), ewaldCutR(1)                                     
     do iG2 = -ewaldCutR(2), ewaldCutR(2)                                  
       do iG3 = -ewaldCutR(3), ewaldCutR(3)                               
         if (iG1 == 0 .and. iG2 == 0 .and. iG3 == 0) cycle               
         iRp = iRp + 1                                                   
         vec(:) = [iG1, iG2, iG3]                                        
         gTrans(:, iRp) = matmul(mol%rec_lat, vec)                       
       end do                                                             
     end do                                                                
   end do                                                                   
                                                                                
   iRp = 0                                                                  
   allocate(rTrans(3, product(2*ewaldCutD+1)))                              
   do iT1 = -ewaldCutD(1), ewaldCutD(1)                                     
     do iT2 = -ewaldCutD(2), ewaldCutD(2)                                  
       do iT3 = -ewaldCutD(3), ewaldCutD(3)                               
         iRp = iRp + 1                                                   
         vec(:) = [iT1, iT2, iT3]                                        
         rTrans(:, iRp) = matmul(mol%lattice, vec)                       
       end do                                                             
     end do                                                                
   end do                                                                   
   
   ! cf, aka ewald parameter              
   cf = sqrtpi/mol%volume**(1.0_wp/3.0_wp)

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
         tmp = erf(gammij*rij)/rij! apart from diagonal(=0), if ij non-bonded rij=1.0d+12 
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

end subroutine goedeckera_PBC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! condense charges to heavy atoms based on topology
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine qheavy(n,at,numnb,numctr,nb,q)
      implicit none
      integer, intent(in)  :: numnb, numctr
      integer,intent(in)   ::  n,nb(numnb,n,numctr),at(n)
      real*8,intent(inout) ::  q(n)

      integer i,j,k,iTr
      real*8 qtmp(n)
      qtmp = q
      do i=1,n
         if(at(i).ne.1) cycle
         qtmp(i)=0
         do iTr=1, numctr
           do j=1,nb(numnb,i,iTr)
             k=nb(j,i,iTr)
             qtmp(k)=qtmp(k)+q(i)/dble(sum(nb(numnb,i,:)))  ! could be a bridging H
           enddo
         enddo
      enddo

      q = qtmp

      end subroutine qheavy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! determine number of cov. bonds between atoms
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine nbondmat(n,numnb,numctr,nb,pair)
      implicit none
!     Dummy
      integer,intent(in)  :: n, numnb, numctr
      integer,intent(in)  :: nb(numnb,n,numctr)
      integer,intent(out) ::  pair(n*(n+1)/2)
!     Stack
      integer i,ni,newi,j,newatom,tag,d,i1,ni1,iii,ii,jj,k,lin
      integer,allocatable :: list(:,:),nlist(:,:),nnn(:),nn(:)
      logical da

      allocate(nnn(n),nn(n),list(5*n,n),nlist(5*n,n))

      nn(1:n)=nb(numnb,1:n,1)

      pair=0
      list=0
      do i=1,n
         ni=nn(i)
         list(1:ni,i)=nb(1:ni,i,1)
      enddo

      nlist=list

      pair=0
      do i=1,n
         do j=1,nb(numnb,i,1)
            k=nb(j,i,1)
            pair(lin(k,i))=1
         enddo
      enddo

!     one bond, tag=1
      tag=1

!     determine up to 3 bonds in between
      do d=1,2

!     loop over atoms
      do i=1,n
         ni=nn(i)
         newi=ni
!        all neighbors of i
         do ii=1,ni
            i1=list(ii,i)
            ni1=nb(numnb,i1,1)
!           all neighbors of neighbors of i
            do iii=1,ni1
               newatom=nb(iii,i1,1)
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


      subroutine nbondmat_pbc(n,numnb,numctr,nb,iTrNeg,neigh,pair)
      implicit none
!     
      type(TNeigh), intent(in) :: neigh ! for locating neighbor
      integer,intent(in)  :: n, numnb, numctr
      integer,intent(in)  :: nb(numnb,n,numctr)
      integer,intent(in)  :: iTrNeg(numctr)
      integer,allocatable,intent(out) ::  pair(:,:,:)
      integer :: nnbi, nbi(2,numnb), cval
      integer :: i,inew,j,inb,ixnb,xnb,iTr,iTrnew,sumiTr,k,iTr2,l
      integer  :: nbr(numnb,n,numctr) ! reduced neighbor list no unpaired bonds
      logical :: hasnb ! true if atom i and k have a paired bond (both have each other as nb)
      integer :: tmpp(3,10*n),nt
      ! using temporary list tmpp to match setup of nbondmat for mindless03
      !  in that case the sum of single bonds in bpair connecting atoms A and B
      !  does not necessarily add up to bpair(A,B) 
      !  if setup should be consistent (excluding half bonds) delete tmpp_usage

      ! only paired bonds should be considered for setting up bpair
      ! e.g. metals have one sided bonds where the bond partner does not
      !  have the metal as a neighbor 
      ! Therefore a list set up with only paired bonds
      tmpp = 0 ! tmpp_usage
      nt=0     ! tmpp_usage
      nbr = nb
      do i=1, n
        do iTr=1,neigh%numctr
          do j=1, nb(neigh%numnb,i,iTr)
            ! k are all the neighbors of i in cell iTr
            k=nb(j,i,iTr)
            hasnb=.false.
            do iTr2=1,neigh%numctr
              do l=1, nb(neigh%numnb,k,iTr2)
                ! check if neighbor k of atom i has i as a neighbor aswell
                if(nb(l,k,iTr2).eq.i) then
                  hasnb=.true.
                endif
              enddo
            enddo 
            ! if hasnb is still false then the bond is not paired and should be deleted
            if (.not.hasnb) then
              ! now delete neighbor and move other neighbors
              do l=j, neigh%numnb-2 ! fails if an atom has numnb-1 neighbors !
                nbr(l,i,iTr)=nb(l+1,i,iTr)
              enddo
              ! reduce neighbors by one
              nbr(numnb,i,iTr) = nbr(numnb,i,iTr) - 1
              ! save the neighbor  
              nt = nt + 1       ! tmpp_usage
              tmpp(1,nt) = k    ! tmpp_usage
              tmpp(2,nt) = i    ! tmpp_usage
              tmpp(3,nt) = iTr  ! tmpp_usage
            endif
          enddo
        enddo
      enddo

      allocate(pair(n,n,numctr), source=0)
      do i=1, n
        ! all first neighbors of i
        do iTr=1, numctr
          do inb=1, nbr(numnb, i, iTr)
            j = nbr(inb,i,iTr)
            pair(j,i,iTr) = 1
          enddo
        enddo
! Can only detect bond paths that lie in at most two cells
! To get num bonds between atoms i&j, which are connected over atoms that lie in
! different cells than i AND j, then first nb() would have to be expanded
! to track cells that i lies in.
        !2nd-3th neighbors
        do cval=1,2
          call countf(n,numctr,numnb,pair,i,cval,nnbi,nbi)
          ! go over i's xth neighbors (x=cval)
          do ixnb=1, nnbi
            ! go over the xth neighbors neighbors
            inew  = nbi(1,ixnb)
            iTrnew= nbi(2,ixnb)
            if (iTrnew.eq.1) then
              !
              do iTr=1, numctr
                do inb=1, nbr(numnb, inew, iTr)
                  j = nbr(inb,inew,iTr)
                  if (pair(j,i,iTr).ne.0) cycle ! take shortest path only
                  if (j.eq.i.AND.iTr.eq.1) cycle ! dont count to-self-bonds
                  pair(j,i,iTr) = cval+1
                enddo
              enddo
            else
                ! idx of neighbors in same cell are the same in every cell
                ! therefore first only consider nb(numnb,inew,iTr=1) 
                do inb=1, nbr(numnb, inew, 1) 
                j = nbr(inb,inew,1)
                  if (pair(j,i,iTrnew).ne.0) cycle ! take shortest path only
                  ! to-self-bonds not possible since iTr.ne.1
                  pair(j,i,iTrnew) = cval+1
                enddo
                ! now consider neighbors in other cells (iTr=2,...)
                do iTr=2, numctr
                  do inb=1, nbr(numnb, inew, iTr) 
                    j = nbr(inb,inew,iTr)       ! obtain nb index j with iTr
                    sumiTr=neigh%fTrSum(iTr,iTrnew)  ! but for pair get correct iTr=sumiTr here
                    if (sumiTr.eq.-1.or.(sumiTr.eq.1.and.j.eq.i)) cycle
                                            ! sumiTr=-1 is outside of 27 centr cells
                                            ! sumiTr= 1 would be to-self-bond
                    if (sumiTr.gt.27) cycle ! 
                    if (pair(j,i,sumiTr).ne.0) cycle ! take shortest path only
                    
                    pair(j,i,sumiTr) = cval+1
                  enddo
                enddo
            endif
          enddo ! ixnb 
        enddo ! cval
      enddo ! i
      ! set pair to 5 for atoms that are even further apart (unless atoms i and j are the same)
      do i=1, n
        do j=1, n
          do iTr=1, numctr
            if(pair(j,i,iTr).eq.0.and.j.ne.i) pair(j,i,iTr)=5
          enddo
        enddo
      enddo
      
      ! tmpp_usage: finaly overwrite the deleted bonds to be bonds again using the saved indices
      do l=1, 10*n
        if (tmpp(1,l).eq.0) exit
        k=tmpp(1,l)
        i=tmpp(2,l)
        iTr=tmpp(3,l)
        pair(k,i,iTr)=1
        pair(i,k,neigh%iTrNeg(iTr))=1
      enddo
      

      end subroutine nbondmat_pbc

      subroutine pairsbond_pbc(n,nn,list,pair,tag)
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

      end subroutine pairsbond_pbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine countf(n,numctr,numnb,pair,i,cval,nnbi,nbi)
        ! n atoms, numctr central cells, search for cval e.g. 1st, 2nd beighbors
        integer, intent(in) :: n, numctr, numnb, i, cval
        ! number of cov bonds between atoms
        integer, intent(in) :: pair(n,n,numctr)
        ! nnbi number of neighbors, nbi idx and cell of nb
        integer, intent(inout) :: nnbi, nbi(2,numnb) 
        integer :: k, l, m
        !
        nnbi = 0
        nbi = 0
          do l=1, n
            do m=1, numctr
              if(pair(l,i,m).eq.cval) then
                ! number of 1st, 2nd, ... neighbors depending on cval
                nnbi = nnbi + 1  
                nbi(1, nnbi) = l  ! idx of nb,  j
                nbi(2, nnbi) = m  ! which cell, iTr
              endif
            enddo
        enddo
      end subroutine countf  

      logical function pilist(ati)
      integer ati
      pilist=.false.
      if(ati.eq.5.or.ati.eq.6.or.ati.eq.7.or.ati.eq.8.or.ati.eq.9.or.ati.eq.16.or.ati.eq.17) pilist=.true.
      end function pilist

      logical function nofs(ati)
      integer ati
      nofs=.false.
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

      integer function ctype(n,at,numnb,numctr,nb,pi,a)
      integer n,a,at(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
      integer i,no,j, iTr

      ctype = 0 ! don't know

      no=0
      do iTr=1, numctr
        do i=1,nb(numnb,a,iTr)
          j=nb(i,a,iTr)
          if(at(j).eq.8.and.pi(j).ne.0) no = no + 1
        enddo
      enddo

      if(no.eq.1.and.pi(a).ne.0) ctype = 1 ! a C=O carbon

      end function ctype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function alphaCO(n,at,hyb,numnb,numctr,nb,pi,a,b)
      integer n,a,b,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
      integer i,j,no,nc,iTr

      alphaCO = .false.
      if(pi(a).ne.0.and.hyb(b).eq.3.and.at(a).eq.6.and.at(b).eq.6) then
         no = 0
         do iTr=1, numctr
           do i=1,nb(numnb,a,iTr)
             j=nb(i,a,iTr)
             if(at(j).eq. 8.and.pi(j).ne.0.and.sum(nb(numnb,j,:)).eq.1) no = no +1 ! a pi =O on the C?
           enddo
         enddo
         if(no.eq.1) then
            alphaCO = .true.
            return
         endif
      endif
      if(pi(b).ne.0.and.hyb(a).eq.3.and.at(b).eq.6.and.at(a).eq.6) then
         no = 0
         do iTr=1, numctr
           do i=1,nb(numnb,b,iTr)
             j=nb(i,b,iTr)
             if(at(j).eq. 8.and.pi(j).ne.0.and.sum(nb(numnb,j,:)).eq.1) no = no +1 ! a pi =O on the C?
           enddo
         enddo
         if(no.eq.1) then
            alphaCO = .true.
            return
         endif
      endif

      end function alphaCO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function amide(n,at,hyb,numnb,numctr,nb,pi,a)
      integer n,a,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
      integer i,j,no,nc,ic,iTr

      amide = .false. ! don't know
      if(pi(a).eq.0.or.hyb(a).ne.3.or.at(a).ne.7) return

      nc=0
      no=0
      do iTr=1, numctr
        do i=1,nb(numnb,a,iTr)
          j=nb(i,a,iTr)
          if(at(j).eq.6.and.pi(j).ne.0) then  ! a pi C on N?
            nc = nc + 1
            ic = j
          endif
        enddo
      enddo

      if(nc.eq.1)then
         do iTr=1, numctr      
           do i=1,nb(numnb,ic,iTr)
             j=nb(i,ic,iTr)
             if(at(j).eq. 8.and.pi(j).ne.0.and.nb(numnb,j,iTr).eq.1) no = no +1 ! a pi =O on the C?
           enddo
         enddo
      endif

      if( no .eq. 1 ) amide = .true.

      end function amide

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      logical function amideH(n,at,hyb,numnb,numctr,nb,pi,a,neigh)
      type(TNeigh), intent(in) :: neigh ! for locating neighbor
      integer n,a,at(n),hyb(n),numnb,numctr,nb(numnb,n,numctr),pi(n)
      integer, allocatable :: locarr(:,:)
      integer i,j,nc,nn,iTr
      !logical amide

      amideH = .false. ! don't know
      if(sum(nb(numnb,a,:)).ne.1               )  return
      call neigh%nbLoc(n,nb,a,locarr) ! locarr gives iTr of cell with the neighbor
      if(size(locarr, dim=2).gt.1) write(*,*) 'WARNING: Neighbors in more cells than expected! source: ini2, amideH'
      nn=nb(1,a,locarr(numnb,1))       ! the N
      deallocate(locarr)
      if(.not.amide(n,at,hyb,numnb,numctr,nb,pi,nn)) return

      nc=0
      do iTr=1, numctr
      do i=1,nb(numnb,nn,iTr)
         j=nb(i,nn,iTr)
         if(at(j).eq.6.and.hyb(j).eq.3) then  ! a sp3 C on N?
            nc = nc + 1
         endif
      enddo
      enddo

      if( nc .eq. 1 ) amideH = .true.

      end function amideH

end module xtb_gfnff_ini2
