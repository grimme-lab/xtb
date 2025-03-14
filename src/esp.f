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

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module esp
      contains
      subroutine espplot(n,nmo,nbf,at,xyz,z,occ,C,basis)
      use xtb_mctc_systools
      use xtb_type_basisset
      use xtb_intpack
      use xtb_printout, only : writecosmofile
      use xtb_setparam, only : set, get_namespace
      implicit none
      type(TBasisset), intent(in) :: basis
      integer nproc,n,nbf,nmo,at(n)
      real*8 xyz(3,n),z(n),occ(nmo),C(nbf,nmo)
      integer,external :: omp_get_num_threads

      real*8,allocatable :: pa(:,:),espe(:), atom_weight(:,:)
      real*8,allocatable :: ee(:,:),efact(:),gama(:),dd(:,:)
      real*8,allocatable :: P(:)
      integer,allocatable:: indp(:,:),nnn(:)
      real*8 point(3),rx,ry,rz,r,intcut,t0,t1,w0,w1,pthr,av
      character(len=:),allocatable :: grid_file
      character(len=:),allocatable :: line
      integer :: err

      integer i,j,k,l,m,done,mprim,npp,ind1,ind2,np
      logical ex

      nproc = 1
      !$ nproc = omp_get_num_threads()
      intcut=10. ! loose primitive cut-off, =20 makes a factor of 2 compared to =10
                 ! errors are around 1e-4 Eh for Etoposid
      pthr=1.d-4 ! -4 def

      if (allocated(set%esp_gridfile)) then
          grid_file = set%esp_gridfile
      else
          grid_file = get_namespace('esp_coord')
      endif

      write(*,*)
      write(*,*)'ESP plotting routine'
      inquire(file=grid_file,exist=ex)
      if(ex) write(*,*)
     &'reading <esp_coord> for xyz triples (in Bohr) to be used'
      write(*,*)'writing <xtb_esp.dat> and <xtb_esp_profile.dat>'

!     call setfmc

      if(.not.ex) call surfac(grid_file,n,xyz,at) ! generate surface coords

      open(unit=83,file=grid_file)
      open(unit=84,file='xtb_esp.dat')

      write(*,*)'computing CAO density maxtrix ...'
      allocate(P(nbf*(nbf+1)/2))
      P=0
      do i=1,nmo
         if(occ(i).gt.1.d-6) then
            m=0
            do k=1,nbf
               do l=1,k-1
                  m=m+1
                  P(m)=P(m)+occ(i)*C(k,i)*C(l,i)*2.0 ! include factor of 2 from contraction loop
               enddo
               m=m+1
               P(m)=P(m)+occ(i)*C(k,i)*C(k,i)
            enddo
         endif
      enddo

!     compute primitive pair data
      write(*,*)'computing primitive pair data ...'
      call preints(n,nbf,xyz,intcut,mprim,npp,pthr,P,basis)
      allocate(efact(npp),gama(npp),ee(3,npp),dd(35,npp),
     &         nnn(npp),indp(mprim,mprim))
      rewind(103)
      indp = 0
      do i=1,npp
         read(103)ind1,ind2,nnn(i),efact(i),gama(i),ee(1:3,i),dd(1:35,i)
         indp(ind1,ind2)=i
      enddo
      close(103,status='delete')

!     write(*,*)
!    &' Cartesian coordinates                     ESP (Hartree)'
      done=0
 10   continue
      read(83,*,end=100) point
      done=done+1
      goto 10
100   continue

      np=done
      allocate(pa(3,np),espe(np))
      allocate(atom_weight(2,np))
      rewind(83)
      do i=1,np
         call getline(83,line,err)
         if(err.ne.0)call raise('E',"Could not read points for ESP!")
         read(line,*,iostat=err) pa(1:3,i), atom_weight(1:2,i)
         if (err.ne.0) then
            read(line,*,iostat=err) pa(1:3,i)
            if (err.ne.0) call raise('E',"Could not read coordinates!")
         endif
      enddo
      write(*,*) np,' points read.'
      call timing(t0,w0)
      point=0
      call espints(n,nbf,xyz,intcut,point,pthr,P,
     &             mprim,npp,nnn,indp,efact,gama,ee,dd,espe(1),basis)
      call timing(t1,w1)
      write(*,'('' estimated wall time (s) '',f6.1)') np*(w1-w0)/nproc
      write(*,*)'computing ESP ...'

!$OMP  PARALLEL PRIVATE(i,j,rx,ry,rz,r,point)
!$OMP& SHARED(espe,pa,basis)
!$OMP& DEFAULT(SHARED)
!$OMP DO
      do i=1,np
!     ESP at point(3)
      espe(i)=0
      point(1:3)=pa(1:3,i)
!     core
      do j=1,n
         rx=xyz(1,j)-point(1)
         ry=xyz(2,j)-point(2)
         rz=xyz(3,j)-point(3)
         r=dsqrt(rx**2+ry**2+rz**2)
         espe(i)=espe(i)+z(j)/(1.0d-14+r)
      enddo
!     ints and contraction with P
      call espints(n,nbf,xyz,intcut,point,pthr,P,
     &             mprim,npp,nnn,indp,efact,gama,ee,dd,espe(i),basis)
      enddo
!$OMP END DO
!$OMP END PARALLEL

      av=0
      do i=1,np
         write(84,'(3E18.10,F14.8)') pa(1:3,i),espe(i)
!        write(*, '(3E14.6,F14.8 )') pa(1:3,i),espe(i)
         av=av+espe(i)
      enddo
      call writecosmofile(np,pa,espe,get_namespace('xtb_esp.cosmo')
     & ,n,at,xyz,atom_weight)
      write(*,'(''maximum/minimum/av ESP value :'',3F12.6)')
     &          maxval(espe),minval(espe),av/np

      call histo(1,np,espe,'xtb_esp_profile.dat')

      if(.not.ex) then
          close(83,status='delete')
      else
          close(83)
      endif
      close(84)

      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine espints(n,nbf,xyz,intcut,point,pthr,P,
     &                   mprim,npp,nnn,indp,efact,gama,ee,dd,espe,basis)
      use xtb_type_basisset
      use xtb_intpack
      implicit none
      type(TBasisset), intent(in) :: basis
      integer n,nbf
      integer mprim,npp,indp(mprim,mprim),nnn(npp)
      real*8  P (nbf*(nbf+1)/2)
      real*8  xyz(3,n)
      real*8  point(3)
      real*8  efact(npp),gama(npp),ee(3,npp),dd(35,npp)
      real*8  espe,intcut,pthr

      integer i,j,k,l,iprimcount,jprimcount,ppp
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real*8  xyza(3),xyzb(3),rab,est,gam,arg
      real*8  tt1,ttt,thr,cce,pthr2

      pthr2=pthr*5.d-3

      k=0
      iprimcount=0
      do i=1,nbf
! aufpunkt i
         xyza(1:3)=xyz(1:3,basis%aoat(i))
! #prims
         npri=basis%nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=basis%nprim(j)
            if(abs(P(k)).lt.pthr) goto 42   ! loose P mat neglect threshold
! aufpunkt j
              xyzb(1:3)=xyz(1:3,basis%aoat(j))
! prim loop
              tt1=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    ppp=indp(iii,jjj)
! cutoff
                    if(ppp.gt.0)then
!                   call propa(opaa0,xyza,xyzb,point,
!    &                         alp(iii),alp(jjj),
!    &                         lao(i),lao(j),ttt,1)
                    cce=basis%cont(iii)*basis%cont(jjj)*efact(ppp)*P(K)
                    if(abs(cce).gt.pthr2)then
                    call propa1(opaa0,point,nnn(ppp),
     &                         gama(ppp),ee(1,ppp),dd(1,ppp),
     &                         ttt)
                    tt1=tt1+ttt*cce
                    endif
                    endif
                 enddo
              enddo
            espe=espe-tt1
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine preints(n,nbf,xyz,intcut,nprimp,ppair,pthr,P,basis)
      use xtb_type_basisset
      use xtb_intpack
      implicit none
      type(TBasisset), intent(in) :: basis
      integer n,nbf,nprimp,ppair
      real*8  P(nbf*(nbf+1)/2)
      real*8  xyz(3,n),intcut,pthr

      integer i,j,k,l,iprimcount,jprimcount
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real*8  xyza(3),xyzb(3),rab,est,gama,arg,cc

      open(unit=103,file='xtb_TMP',form='unformatted')

      k=0
      ppair=0
      iprimcount=0
      do i=1,nbf
! aufpunkt i
         xyza(1:3)=xyz(1:3,basis%aoat(i))
! #prims
         npri=basis%nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=basis%nprim(j)
            if(abs(P(k)).lt.pthr) goto 42   ! loose P mat neglect threshold
! aufpunkt j
              xyzb(1:3)=xyz(1:3,basis%aoat(j))
              rab=(xyza(1)-xyzb(1))**2
     &           +(xyza(2)-xyzb(2))**2
     &           +(xyza(3)-xyzb(3))**2
! prim loop
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    gama=1.0d0/(basis%alp(iii)+basis%alp(jjj))
                    est=rab*basis%alp(iii)*basis%alp(jjj)*gama
                    cc=basis%cont(iii)*basis%cont(jjj)
! cutoff
                    if(est.lt.intcut.and.abs(cc*P(k)).gt.1.d-5)then
                    ppair=ppair+1
                    call propa0(xyza,xyzb,
     &                         basis%alp(iii),basis%alp(jjj),
     &                         basis%lao(i),basis%lao(j),iii,jjj)
                    endif
                 enddo
              enddo
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo
      nprimp=iprimcount

      write(*,*) 'number of prim pairs',ppair

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine propa0(a,b,etaij4,etakl4,iff1,iff2,ind1,ind2)
      use xtb_intpack
      implicit integer(i-n)
      implicit real*8(a-h,o-z)
! aufpunkte
      real*8 a(3),b(3)
! local
      common /abfunc/ ra(3),rb(3),ga,gb,ia,ib
      dimension d(3),dd(84),ddd(35)
      dimension e(3),aa(20),bb(20)
      dimension lin(84),min(84),nin(84)
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,
     8 1,1,4,2/

! --- a,b are centres of gaussians
      ra = a
      rb = b
! --- ga,gb are their exponents
      ga=etaij4
      gb=etakl4
      aa = 0
      bb = 0
      aa(iff1)=1.0d0
      bb(iff2)=1.0d0
      ia=iff1
      ib=iff2
! --- ia,ib are the canonical indices for monom
! --- apply product theorem
      cij=0.0d0
      ckl=0.0d0
      call divpt(a,etaij4,b,etakl4,cij,ckl,e,gama,efact)

! --- calculate cartesian prefactor for first gaussian
      call rhftce(aa,a,e,iff1)
! --- calculate cartesian prefactor for second gaussian
      call rhftce(bb,b,e,iff2)
! --- form their product
      call prod(aa,bb,dd,iff1,iff2)

      ddd(1:35)=dd(1:35)
      call lprimprod(iff1,iff2,nnn)  ! l-type of pair
      write(103) ind1,ind2,nnn,efact,gama,e,ddd

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! no f-functions !

      subroutine propa1(aname,c,nnn,gama,e,dd,va)
      use xtb_intpack
      implicit integer(i-n)
      implicit real*8(a-h,o-z)
      external aname
! aufpunkte,ref point,intarray
      real*8 c(3),va,dd(35),e(3)  ! dd(84) for f-fct
! local
      dimension lin(84),min(84),nin(84),d(3)
      data lin/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,4,0,0,3,3,1,0,1,0
     1 ,2,2,0,2,1,1,5,0,0,3,3,2,2,0,0,4,4,1,0,0,1,1,3,1,2,2,1,6,0,0,3,
     2 3,0,5,5,1,0,0,1,4,4,2,0,2,0,3,3,1,2,2,1,4,1,1,2/,min/0,0,1,0,0,
     3 2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,0,4,0,1,0,3,3,0,1,2,0,2,1,2,1,0,5,
     4 0,2,0,3,0,3,2,1,0,4,4,1,0,1,1,3,2,1,2,0,6,0,3,0,3,1,0,0,1,5,5,2,0
     5,0,2,4,4,2,1,3,1,3,2,1,4,1,2/,nin/0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,
     6 1,2,2,1,0,0,4,0,1,0,1,3,3,0,2,2,1,1,2,0,0,5,0,2,0,3,2,3,0,1,0,1,
     7 4,4,3,1,1,1,2,2,0,0,6,0,3,3,0,1,5,5,1,0,0,2,4,4,0,2,1,2,2,3,1,3,
     8 1,1,4,2/

      va =0
!     e is center of product gaussian with exponent gama, c is reference point
      d = e - c
!     aname represents an external function
      do i=1,nnn
         if (dabs(dd(i))-1.d-8 .gt. 0.d0) then
            call aname(lin(i),min(i),nin(i),gama,v,d)
            va=dd(i)*v+va
         end if
      end do
      return

      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine histo(n1,n2,aval,ftmp)
      implicit none
      integer n1,n2
      real*8 aval(*)
      character*(*) ftmp

      real*8 mina,maxa,dx,val,bin1,bin2,dum,du
      integer i,j,k,nbin,bin(10)

      nbin=10
      mina= 1.d+42
      maxa=-1.d+42
      do i=n1,n2
         if(aval(i).lt.mina)mina=aval(i)
         if(aval(i).gt.maxa)maxa=aval(i)
      enddo
      mina=mina-1.d-2
      maxa=maxa+1.d-2
      dx=(maxa-mina)/(nbin)

      bin=0
      do i=n1,n2
         val=aval(i)
         bin1=mina
         do k=1,nbin
            bin2=bin1+dx
            if(val.ge.bin1.and.val.lt.bin2)then
            bin(k)=bin(k)+1
            endif
            bin1=bin1+dx
         enddo
      enddo

      dum=0
      do i=1,nbin
         dum=dum+dx*bin(i)
      enddo

      open(unit=11,file=ftmp)
      bin1=mina
      do i=1,nbin
         du=bin1+0.5*dx
         write(11,*) du,bin(i)/dum
         bin1=bin1+dx
      enddo
!     write(11,*)n2-n1+1
!     do i=n1,n2
!        write(11,*) aval(i)
!     enddo
      close(11)

      end subroutine histo

      subroutine lprimprod(iff1,iff2,nnn)
      implicit none
      integer iff1,iff2,nnn

      if(iff1.gt.4.or.iff2.gt.4) goto 120
      if(iff1.gt.1.or.iff2.gt.1) goto 130
!           s - s
      nnn=1
      return
  130 continue
      if(iff1.ne.1.and.iff2.ne.1) goto 131
!           s - p
      nnn=4
      return
  131 continue
!           p - p
      nnn=10
      return
  120 continue
      if(iff1.ne.1.and.iff2.ne.1) goto 121
!           s - d
      nnn=10
      return
  121 continue
      if(iff1.gt.4.and.iff2.gt.4) goto 122
!           p - d
      nnn=20
      return
  122 continue
!           d - d
      nnn=35

      end subroutine lprimprod
      end module esp
