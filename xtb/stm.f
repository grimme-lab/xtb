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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c simple STM picture simulation routine using PT according to
c Tersoff & Hamann, PRL 50 (1983), 1998.
c assumes that the molecule is lying mainly in xy plane
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module stm
      contains
      subroutine stmpic(n,nmo,nbf,at,xyz,C,efermi,emo,basis)
      use iso_fortran_env, wp => real64
      use tbdef_basisset
      use setparam
      use scc_core, only : dmat
      use esp
      implicit none
      type(tb_basisset), intent(in) :: basis
      integer nproc,n,nbf,nmo,at(n)
      real(wp) xyz(3,n),emo(nmo),C(nbf,nmo),efermi

      real(wp),allocatable :: pa(:,:),espe(:)
      real(wp),allocatable :: ee(:,:),efact(:),gama(:),dd(:,:)
      real(wp),allocatable :: P(:),eocc(:),ptmp(:,:)
      integer,allocatable:: indp(:,:),nnn(:)
      real(wp) point(3),rx,ry,rz,r,intcut,t0,t1,w0,w1,pthr,av,scal,alp
      real(wp) stepr,dx,dy,dz,xx1,yy1,zz1,bord,dum(3),targ,zmin,tmin
      real(wp) sqrt2pi,cut,bohr,bordxy,pot,de,stepz
      parameter (sqrt2pi = sqrt(2*3.141592653589790d0))
      parameter (bohr=0.52917726d0)

      integer i,j,k,l,m,done,mprim,npp,ind1,ind2,np,m1,m2,m3
      !$ integer :: omp_get_num_threads
      nproc = 1
      !$ nproc = omp_get_num_threads()

      write(*,*)
      write(*,*)'STM picture simulation routine'
      write(*,*)'in constant current mode'
      write(*,'('' tip DOS energy broadening (eV):'',f9.1)')stm_alp
      write(*,'('' tip potential (V)             :'',f9.6)')stm_pot
      write(*,'('' constant current value (a.u.) :'',f9.6)')stm_targ
      write(*,'('' grid width (Bohr)             :'',f9.6)')stm_grid
      write(*,'('' int and Pmat neglect stm_thr  :'',f9.6)')stm_thr

c prog parameters (should be made a user input)
      alp =stm_alp  !1.500 ! width of tip DOS energy broadening (eV)
      targ=stm_targ !1.d-4 ! constant current value (a.u.)
      stepr=stm_grid!0.5d0! grid width (Bohr), half that value along Z
      pot  =stm_pot ! potential of tip vs. molecule, negative values let e flow from mol to tip
                    ! ie occ space of mol is probed
c fix
      bord=4.5d0                    ! Z-limit (Bohr), max Z of molecule added
      intcut=20.-3.0*log10(stm_thr) ! primitive cut-off, =20 makes a factor of 2 compared to =10
      pthr=1.d-5*stm_thr            ! dmat neglect threshold
      write(*,'('' incut                         :'',f9.6)')intcut
      write(*,'('' Pthr                          :'',e9.3)')pthr

      allocate(P(nbf*(nbf+1)/2),ptmp(nbf,nbf),eocc(nmo))

c start calc
      scal=1.0d0/(alp*sqrt2pi)
      do i=1,nmo
         de=efermi-emo(i)-pot
         eocc(i)=scal*exp(-0.5*(de/alp)**2)  ! assume Gaussian DOS for tip instead of delta function
                                             ! as in the original work (eq 11)
      enddo

      write(*,*)'computing MO energy weighted density maxtrix ...'
      call dmat(nbf,eocc,C,ptmp)
      m=0
      do k=1,nbf
         do l=1,k
            m=m+1
            P(m)=ptmp(l,k)
         enddo
      enddo
      deallocate(ptmp,eocc)

c     compute primitive pair data
      write(*,*)'computing primitive pair data ...'
      call preints(n,nbf,xyz,intcut,mprim,npp,pthr,P,basis)
      allocate(efact(npp),gama(npp),ee(3,npp),dd(35,npp),
     .         nnn(npp),indp(mprim,mprim))
      rewind(103)
      indp = 0
      do i=1,npp
         read(103)ind1,ind2,nnn(i),efact(i),gama(i),ee(1:3,i),dd(1:35,i)
         indp(ind1,ind2)=i
      enddo
      close(103,status='delete')

c     R space limits
      bordxy=6.0
      stepz=stepr*0.5 ! better resol in Z
      dx=maxval(xyz(1,1:n))-minval(xyz(1,1:n))+bordxy
      dy=maxval(xyz(2,1:n))-minval(xyz(2,1:n))+bordxy
      dz=maxval(xyz(3,1:n))+bord
      xx1=minval(xyz(1,1:n))-bordxy*0.5
      yy1=minval(xyz(2,1:n))-bordxy*0.5
      zz1=0
      m1=1+int(dx/stepr)
      m2=1+int(dy/stepr)
      m3=1+int(dz/stepz)
      np=m1*m2*m3
      write(*,'('' # xyz range:'',3f8.2)') dx,dy,dz
      write(*,'('' # points   :'',i8  )') np
      allocate(pa(3,np),espe(np))

cccccccccccccc
c R grid
cccccccccccccc
      l=0
      dx=xx1
      do i=1,m1
         dum(1)=dx
         dy=yy1
         do j=1,m2
            dum(2)=dy
            dz=zz1
            do k=1,m3
               dum(3)=dz
               l=l+1
               pa(1:3,l)=dum(1:3)
               dz=dz+stepz
            enddo
            dy=dy+stepr
         enddo
         dx=dx+stepr
      enddo

      write(*,*)'computing density ...'

!$OMP  PARALLEL PRIVATE(r,i,point)
!$OMP& SHARED(espe,pa)
!$OMP& DEFAULT(SHARED)
!$OMP DO
      do i=1,np
c     at point(3)
         point(1:3)=pa(1:3,i)
c     if(mod(i,np/10).eq.0.and.nproc.eq.1) write(*,*) 100*i/np,' % done'
c     ints and contraction with P
         call densints(n,nbf,xyz,intcut,point,pthr,P,
     .                 mprim,npp,nnn,indp,efact,gama,ee,dd,cut,r,basis)
         espe(i)=r ! prefacto
      enddo
!$OMP END DO
!$OMP END PARALLEL

c makes no sense because values not treated (flaged 99) ar in
c     av=0
c     do i=1,np
c        av=av+espe(i)
c     enddo
c     write(*,'('' maximum/minimum/av current value :'',2F12.6,F9.4)')
c    .          maxval(espe),minval(espe),av/np

      write(*,*)
     .' writing STM picture in xyz (Angstroem!) to file <stm.dat>'
c output in gnuplot format, constant current mode
      open(unit=42,file='stm.dat')
      l=0
      dx=xx1
      do i=1,m1
         dy=yy1
         do j=1,m2
            dz=zz1
            tmin=1.d+42
            zmin=0
            do k=1,m3
               l=l+1
               if(abs(espe(l)-targ).lt.tmin)then ! Z point closest to target current
                  tmin=abs(espe(l)-targ)
                  zmin=dz
               endif
               dz=dz+stepz
            enddo
            write(42,'(3F14.8)') dx*bohr,dy*bohr,zmin*bohr
            dy=dy+stepr
         enddo
         write(42,*)
         dx=dx+stepr
      enddo
      close(42)

      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c density routine in analogy to esp plot
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine densints(n,nbf,xyz,intcut,point,pthr,P,
     &                    mprim,npp,nnn,indp,efact,gama,ee,dd,cut,espe,
     &                    basis)
      use iso_fortran_env, wp => real64
      use tbdef_basisset
      use intpack
      use esp
      implicit none
      type(tb_basisset), intent(in) :: basis
      integer n,nbf
      integer mprim,npp,indp(mprim,mprim),nnn(npp)
      real(wp) P (nbf*(nbf+1)/2)
      real(wp) xyz(3,n)
      real(wp) point(3)
      real(wp) efact(npp),gama(npp),ee(3,npp),dd(35,npp)
      real(wp) espe,intcut,pthr,cut

      integer i,j,k,l,iprimcount,jprimcount,ppp
      integer npri,nprj,ii,iii,jj,jjj,ll,m,li,lj,mm,nn
      real(wp) xyza(3),xyzb(3),rab,est,gam,arg
      real(wp) tt1,ttt,thr,cce,pthr2

c check if point is close to any atom and if yes dens is too high and exit
      do j=1,n
         rab=(xyz(1,j)-point(1))**2+(xyz(2,j)-point(2))**2
     .      +(xyz(3,j)-point(3))**2
         if(rab.lt.9.)then
            espe=99.
            return
         endif
      enddo

      pthr2=pthr*1.d-1
      espe=0

      k=0
      iprimcount=0
      do i=1,nbf
c aufpunkt i
         xyza(1:3)=xyz(1:3,basis%aoat(i))
c #prims
         npri=basis%nprim(i)
         jprimcount=0
         do j=1,i
            k=k+1
            nprj=basis%nprim(j)
            if(abs(P(k)).lt.pthr) goto 42   ! loose P mat neglect threshold
c aufpunkt j
              xyzb(1:3)=xyz(1:3,basis%aoat(j))
c prim loop
              tt1=0.0d0
              do ii=1,npri
                 iii=iprimcount+ii
                 do jj=1,nprj
                    jjj=jprimcount+jj
                    ppp=indp(iii,jjj)
c cutoff
                    if(ppp.gt.0)then
                    cce=basis%cont(iii)*basis%cont(jjj)*efact(ppp)*P(K)
                    if(abs(cce).gt.pthr2)then
                    call propa1(opac3,point,nnn(ppp),
     .                         gama(ppp),ee(1,ppp),dd(1,ppp),
     .                         ttt)
                    tt1=tt1+ttt*cce
                    endif
                    endif
                 enddo
              enddo
            espe=espe+tt1
            if(espe.gt.cut) return     ! value too large i.e. we're looking for constant I in space
                                       ! and dens is positive
 42         jprimcount=jprimcount+nprj
         enddo
         iprimcount=iprimcount+npri
      enddo

      end


      end module stm
