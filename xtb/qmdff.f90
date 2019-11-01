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

module qmdff
   use iso_fortran_env, wp => real64
   implicit none
!  integer, private :: ndim
!  parameter(ndim=50000)
   integer, private,parameter   :: max_elem = 94
   integer, private,parameter   :: ntterm = 4
   integer, private :: nbond   = 0
   integer, private :: nbond12 = 0
   integer, private :: nangl   = 0
   integer, private :: ntors   = 0
   integer, private :: nhb     = 0
   integer, private :: nnci    = 0
   integer, private,allocatable :: bond (:,:)
   integer, private,allocatable :: angl (:,:)
   integer, private,allocatable :: tors (:,:)
   integer, private,allocatable :: hb   (:,:)
   integer, private,allocatable :: nci  (:,:)
   real(wp),private,allocatable :: vbond(:,:)
   real(wp),private,allocatable :: vangl(:,:)
   real(wp),private,allocatable :: vtors(:,:)
   real(wp),private,allocatable :: vhb  (:,:)
   real(wp),private :: scalehb(max_elem) = 0.0_wp
   real(wp),private :: scalexb(max_elem) = 0.0_wp
   real(wp),private :: morsethr = 99.0_wp
   real(wp),private,allocatable :: qff(:)
   real(wp),private,allocatable :: c6ff(:,:)                    
   real(wp),private :: eps1(6)=(/0.00_wp,0.00_wp,0.85_wp,1.00_wp,1.00_wp,0.00_wp/)
   real(wp),private :: eps2(6)=(/0.00_wp,0.00_wp,0.50_wp,0.50_wp,1.00_wp,1.00_wp/)
   real(wp),private :: zabff (max_elem,max_elem) = 0.0_wp
   real(wp),private :: r094ff(max_elem,max_elem) = 0.0_wp
   real(wp),private :: sr42ff(max_elem,max_elem) = 0.0_wp
   real(wp),private :: r0abff(max_elem,max_elem) = 0.0_wp

   real(wp),private,parameter   :: pi  = 3.14159265358979323846264338327950_wp
   real(wp),private,parameter   :: pi2 = 6.28318530717958623199592693708837_wp
   real(wp),private,parameter   :: pi6 = 961.389193575304212170602974626298_wp
   real(wp),private,parameter   :: spi = 1.77245385090551599275151910313925_wp
contains

subroutine ff_ini(n,at,xyz,cn,s6)
   use aoparam
   use dftd4, only : r2r4 => r4r2, rcov
   implicit none
   integer, intent(in) :: n
   integer, intent(in) :: at(n)
   real(wp),intent(in) :: xyz(3,n)
   real(wp),intent(in) :: cn(n)
   real(wp),intent(in) :: s6

   real(wp) :: e
   real(wp) :: g(3,n)
   real(wp) :: vz(max_elem)
   real(wp) :: s8 = 2.70_wp
   real(wp) :: a1 = 0.45_wp
   real(wp) :: a2 = 4.00_wp
   real(wp) :: c6
   integer  :: i,j,i1,i2,iz1,iz2
   logical  :: exist

!  include 'ffparam.f90'

   ! clear heap from previous runs
   if (allocated(bond))  deallocate(bond)
   if (allocated(angl))  deallocate(angl)
   if (allocated(tors))  deallocate(tors)
   if (allocated(hb))    deallocate(hb)
   if (allocated(nci))   deallocate(nci)
   if (allocated(vbond)) deallocate(vbond)
   if (allocated(vangl)) deallocate(vangl)
   if (allocated(vtors)) deallocate(vtors)
   if (allocated(vhb))   deallocate(vhb)
   if (allocated(qff))   deallocate(qff)
   if (allocated(c6ff))  deallocate(c6ff)

   inquire(file='solvent',exist=exist)
   if(.not.exist) then
      call raise('E','FF run requested but solvent file does not exist!',1)
   endif

   write(output_unit,'("initializing FF ..")')

!  if(n.gt.5000) call raise('E','too many atoms for FF',1) ! fixed this for you

   scalehb=0
   scalehb(7  )=0.8
   scalehb(8  )=0.3
   scalehb(9  )=0.1
   scalehb(15 )=2.0
   scalehb(16 )=2.0
   scalehb(17 )=2.0
   scalehb(34 )=2.0
   scalehb(35 )=2.0
   scalexb=0
   scalexb(17 )=0.30
   scalexb(35 )=0.60
   scalexb(53 )=0.80
   scalexb(85 )=1.00
   eps1(1)=0
   eps1(2)=0
   eps1(3)=0.85
   eps1(4)=1
   eps1(5)=1
   eps2(1)=0
   eps2(2)=0
   eps2(3)=0.5
   eps2(4)=0.5
   eps2(5)=1
   eps2(6)=1.00
   eps1(6)=0

   call rdsolvff(n,'solvent')

   call valelff(vz)
   ! D3
   s8 = s8 / s6

   do i=1,max_elem
      do j=1,max_elem
         sr42ff(j,i)=3.0d0*s8*r2r4(i)*r2r4(j)
         r094ff(j,i)=a1*sqrt(3.0d0*r2r4(i)*r2r4(j))+a2
         zabff (j,i)=vz(i)*vz(j)
      enddo
   enddo
   call setr0ab(max_elem,0.52917726d0,r0abff)
   r0abff = 16.5 / r0abff**1.5

   allocate( c6ff(n,n), source = 0.0_wp )
   do i1=1,n
      iz1=at(i1)
      do i2=1,i1
         iz2=at(i2)
         call getc6(iz1,iz2,cn(i1),cn(i2),c6)
         c6ff(i2,i1)=c6 * s6 ! simulates solvent for s6 < 1
         c6ff(i1,i2)=c6 * s6
      enddo
   enddo


!  call ff_eg  (n,at,xyz,e,g)
!  call ff_nonb(n,at,xyz,e,g)
!  call ff_hb  (n,at,xyz,e,g)
!  write(*,*) 'QMDFF energy on input coordinates: ',e
!  write(*,*) 'QMDFF grad   on input coordinates: ',sqrt(sum(g**2))
!  stop

end subroutine ff_ini

!cccccccccccccccccccccccccccccccccccccccccccccccccccc
! read solvent coord and FF
!cccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine rdsolvff(nat,fname)
   implicit none
   character(len=*),intent(in) :: fname

   integer, intent(in) :: nat
   integer  :: at(nat)
   real(wp) :: xyz(3,nat)

   integer  :: i,j,n,nn,nterm,idum,imass,ich
   real(wp) :: qdum,r,xx(10),dum1,dum2
   character(len=80) atmp

!  include 'ffparam.f90'

   open(newunit=ich,file=fname)

   read (ich,'(a)') atmp
   call readl(atmp,xx,nn)
   n=idint(xx(1))
   if(n.ne.nat) then
      write(*,*) n,nat
      stop 'ff read error'
   endif
   qdum=xx(2)
   read (ich,'(a)') atmp
   write(output_unit,'(a)')
   write(output_unit,'(a)')'reading QMDFF intramolecular FF:'
   write(output_unit,'(a)')'method used for charges          '//trim(atmp)
   if(index(atmp,'morsethr=').ne.0) then
      call readl(atmp,xx,nn)
      morsethr=xx(nn)
      write(output_unit,'(a,f12.8)') 'taking morsethr from file ',morsethr
   else
      morsethr=99.
   endif

   allocate(qff(n), source = 0.0_wp)
   do i=1,n
      read (ich,*)at(i),xyz(1:3,i),qff(i),imass
   enddo
   read (ich,*) nbond,nangl,ntors,nhb,nnci,nbond12
!  if(nbond.gt.ndim.or.nangl.gt.ndim.or.ntors.gt.ndim) &
!  & call raise('E','too many FF terms',1)
   allocate ( bond(2,nbond), source = 0 )
   allocate ( vbond(3,nbond), source = 0.0_wp )
   do i=1,nbond
      read (ich,*)bond(1:2,i),vbond(1:3,i)
   enddo

   allocate ( angl(3,nangl), source = 0 )
   allocate ( vangl(2,nangl), source = 0.0_wp )
   do i=1,nangl
      read (ich,*) angl(1,i),angl(2,i),angl(3,i),vangl(1:2,i)
   enddo

   allocate ( tors(6,ntors), source = 0 )
   allocate ( vtors(3*ntterm+2,ntors), source = 0.0_wp )
   do i=1,ntors
      read (ich,*) tors(1:6,i),vtors(1:2,i), &
         &         (vtors(3*(j-1)+3,i), &
         &          vtors(3*(j-1)+4,i), &
         &          vtors(3*(j-1)+5,i),j=1,tors(5,i))
   enddo
   do i=1,ntors
      do j=1,tors(5,i)
         vtors(3*(j-1)+4,i)=vtors(3*(j-1)+4,i)*pi
      enddo
   enddo

   if(nhb.gt.0) then
      allocate( hb(3,nhb), source = 0 )
      allocate( vhb(2,nhb), source = 0.0_wp )
      read(ich,'(6(3i5,2x))')hb(1:3,1:nhb)
      do i=1,nhb
         if(at(hb(3,i)).eq.1)then
            call hbpara(10.0d0,5.0d0,qff(hb(1,i)),dum1)
            vhb(1,i)=dum1*scalehb(at(hb(1,i)))
            call hbpara(10.0d0,5.0d0,qff(hb(2,i)),dum1)
            vhb(2,i)=dum1*scalehb(at(hb(2,i)))
         else
            call hbpara(-6.5d0,1.0d0,qff(hb(3,i)),dum1)
            vhb(1,i)=scalexb(at(hb(3,i)))*dum1
         endif
      enddo
   endif

   allocate ( nci(3,nnci), source = 0 )
   read(ich,'(8(3i5,2x))')nci(1:3,1:nnci)

   close(ich)

   qff=qff*qdum

   write(output_unit,*)'nbond,nangl,ntors   :',nbond,nangl,ntors
   write(output_unit,*)'NCI terms added     :',nnci
   write(output_unit,*)'HB/XB  terms added  :',nhb
   return

   do i=1,nbond
      write(*,'(2i6,3F14.9)')bond(1,i),bond(2,i),vbond(1:3,i)
   enddo
   do i=1,nangl
      write(*,'(3i6,5F14.9)')angl(1,i),angl(2,i),angl(3,i), &
         &                        vangl(1:2,i)
   enddo
   do i=1,ntors
      write(*,'(5i6,6F12.8)')tors(1:5,i),vtors(1:2+ntterm,i)
   enddo

end subroutine rdsolvff

! evaluate the internal energy using the FF as
! specified in common (ffparam)

subroutine ff_eg(n,at,xyz,e,g)
   implicit none
   integer  :: n,at(n)
   real(wp) :: xyz(3,n),e,g(3,n)

!  include 'ffparam.f90'

   real(wp) :: vp(3),dt,deddt,rp,cosa,rab2,rcb2,rmul1,rmul2
   real(wp) :: cosna(4),sinna(4),v(4),sa(4),ca(4),phix(4),dphi(4)
   real(wp) :: va(3),vb(3),vc(3),vd(3),vba(3),vcb(3),vdc(3),vab(3)
   real(wp) :: vca(3),vdb(3),vt(3),vu(3),vtu(3),dedt(3),dedu(3)
   real(wp) :: deda(3),dedb(3),dedc(3),dedd(3),vtmp1(3),vtmp2(3)
   real(wp) :: rt2,ru2,rtu,rcb,rtru,damp,damp2,rac2,racut,omega
   real(wp) :: aai,aai2,r2,rjk,dampjk,damp2jk,dampij,damp2ij,dampjl
   real(wp) :: rkl,rjl,dampkl,damp2kl,damp2jl,term1(3),term2(3),term3(3)
   integer  :: ib,id
   real(wp) :: dda(3),ddb(3),ddc(3),ddd(3)

   real(wp) :: r,thab,thbc,ra(3),fac,eb,or,kij,rij,alpha,dij
   real(wp) :: dei(3),dej(3),dek(3),kijk,rb(3),ea,dedtheta,pi6,ban
   real(wp) :: c0,c1,c2,theta,sinth,costh,expo
   real(wp) :: dphidri(3),dphidrj(3),dphidrk(3),dphidrl(3),step
   real(wp) :: et,phi0,rn,kijkl,dedphi,phi,cosrph,cosrph0,er,el
   real(wp) :: dphi1,dphi2,x1cos,x2cos,x1sin,x2sin,phipi,e1,e2,ef
   integer  :: i,j,k,l,m,ic,ia,list(4),jj,mm,ii,kk,it,nt,it2


   e=0.0_wp
   g=0.0_wp

   ! strech
   do m=1,nbond
      i=bond(1,m)
      j=bond(2,m)
      r2=   ((xyz(1,i)-xyz(1,j))**2 &
      &     +(xyz(2,i)-xyz(2,j))**2 &
      &     +(xyz(3,i)-xyz(3,j))**2)
      r =sqrt(r2)
      do ic=1,3
         ra(ic)=xyz(ic,i)-xyz(ic,j)
      end do
      rij=vbond(1,m)
      kij=vbond(2,m)
      ! LJ
      if(bond(3,m).eq.0)then
         aai=vbond(3,m)
         aai2 =aai/2
         e=e+    kij*(1.+(rij/r)**aai - 2.*(rij/r)**aai2 )
         fac=aai*kij*(   -(rij/r)**aai  +  (rij/r)**aai2 )/r2
         ! Morse
      elseif(bond(3,m).eq.1)then
         aai=0.5d0*vbond(3,m)/rij
         e=e+    kij*(1.d0-exp(-aai*(r-rij)))**2
         fac=2.d0*aai*kij*exp(-aai*(r-rij))*(1.d0-exp(-aai*(r-rij)))/r
         ! Morse metal
!     elseif(bond(3,m).eq.1)then
!        aai=((7.d0*vbond(3,m)**4+ &
!    &        36.d0*vbond(3,m)**3 &
!    &       +44.d0*vbond(3,m)**2) &
!    &       /(192.d0*rij**4))**0.25d0
!        e=e+    kij*(1.d0-exp(-aai*(r-rij)))**4
!        fac=4.d0*aai*kij*exp(-aai*(r-rij))* &
!    &      (1.d0-exp(-aai*(r-rij)))**3/r
      endif

      do ic=1,3
         g(ic,i)=g(ic,i)+fac*ra(ic)
         g(ic,j)=g(ic,j)-fac*ra(ic)
      enddo
   enddo

   ! bend
   do m=1,nangl
      j = angl(1,m)
      i = angl(2,m)
      k = angl(3,m)
      c0  =vangl(1,m)
      kijk=vangl(2,m)
      va(1:3) = xyz(1:3,i)
      vb(1:3) = xyz(1:3,j)
      vc(1:3) = xyz(1:3,k)
      vab = va-vb
      vcb = vc-vb
      rab2 = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
      rcb2 = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
      call crprod(vcb,vab,vp)
      rp = norm2(vp)+1.d-14
      call impsc(vab,vcb,cosa)
      cosa = dble(min(1.0d0,max(-1.0d0,cosa)))
      theta= dacos(cosa)

      call abdamp(at(i),at(j),rab2,dampij,damp2ij)
      call abdamp(at(k),at(j),rcb2,dampjk,damp2jk)
      damp=dampij*dampjk

      if(pi-c0.lt.1.d-6)then
         dt  = theta - c0
         ea  = kijk * dt**2
         deddt = 2.d0 * kijk * dt
      else
         ea=kijk*(cosa-cos(c0))**2
         deddt=2.*kijk*sin(theta)*(cos(c0)-cosa)
      endif

      e = e + ea * damp
      call crprod(vab,vp,deda)
      rmul1 = -deddt / (rab2*rp)
      deda = rmul1*deda
      call crprod(vcb,vp,dedc)
      rmul2 =  deddt / (rcb2*rp)
      dedc = rmul2*dedc
      dedb = deda + dedc
      term1(1:3)=ea*damp2ij*dampjk*vab(1:3)
      term2(1:3)=ea*damp2jk*dampij*vcb(1:3)
      g(1:3,i) = g(1:3,i) + deda(1:3)*damp+term1(1:3)
      g(1:3,j) = g(1:3,j) - dedb(1:3)*damp-term1(1:3)-term2(1:3)
      g(1:3,k) = g(1:3,k) + dedc(1:3)*damp+term2(1:3)
   enddo

   ! torsion
   do m=1,ntors
      i=tors(1,m)
      j=tors(2,m)
      k=tors(3,m)
      l=tors(4,m)
      nt=tors(5,m)
      phi0 =vtors(1,m)

      if(tors(6,m).ne.2)then
!        ----------------------------------------------------
         vab(1:3) = xyz(1:3,i)-xyz(1:3,j)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,k)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rkl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         call abdamp(at(i),at(j),rij,dampij,damp2ij)
         call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
         call abdamp(at(k),at(l),rkl,dampkl,damp2kl)
         damp= dampjk*dampij*dampkl

         phi=valijklff(n,xyz,i,j,k,l)
         call dphidr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)

         et=0
         dij=0
         mm=3
         do it=1,nt
            rn=vtors(mm,m)
            dphi1=phi-phi0
            dphi2=phi+phi0-pi2
            c1=rn*dphi1+vtors(mm+1,m)
            c2=rn*dphi2+vtors(mm+1,m)
            x1cos=cos(c1)
            x2cos=cos(c2)
            x1sin=sin(c1)
            x2sin=sin(c2)
            phipi=phi-pi
            ef   =erf(phipi)
            e1   =vtors(mm+2,m)*(1.+x1cos)
            e2   =vtors(mm+2,m)*(1.+x2cos)
            et   =et+0.5*(1.-ef)*e1+(0.5+0.5*ef)*e2
            expo =exp(-phipi**2)/spi
            dij  =dij-expo*e1- 0.5*(1.-ef)*vtors(mm+2,m)*x1sin*rn+ &
               &      expo*e2-(0.5+0.5*ef)*vtors(mm+2,m)*x2sin*rn
            mm=mm+3
         enddo

         et=  et*vtors(2,m)
         dij=dij*vtors(2,m)*damp

         term1(1:3)=et*damp2ij*dampjk*dampkl*vab(1:3)
         term2(1:3)=et*damp2jk*dampij*dampkl*vcb(1:3)
         term3(1:3)=et*damp2kl*dampij*dampjk*vdc(1:3)

         g(1:3,i)=g(1:3,i)+dij*dda(1:3)+term1
         g(1:3,j)=g(1:3,j)+dij*ddb(1:3)-term1+term2
         g(1:3,k)=g(1:3,k)+dij*ddc(1:3)+term3-term2
         g(1:3,l)=g(1:3,l)+dij*ddd(1:3)-term3
         e=e+et*damp
!        ----------------------------------------------------
      else
!        ----------------------------------------------------
         vab(1:3) = xyz(1:3,j)-xyz(1:3,i)
         vcb(1:3) = xyz(1:3,j)-xyz(1:3,k)
         vdc(1:3) = xyz(1:3,j)-xyz(1:3,l)
         rij = vab(1)*vab(1) + vab(2)*vab(2) + vab(3)*vab(3)
         rjk = vcb(1)*vcb(1) + vcb(2)*vcb(2) + vcb(3)*vcb(3)
         rjl = vdc(1)*vdc(1) + vdc(2)*vdc(2) + vdc(3)*vdc(3)
         call abdamp(at(i),at(j),rij,dampij,damp2ij)
         call abdamp(at(k),at(j),rjk,dampjk,damp2jk)
         call abdamp(at(j),at(l),rjl,dampjl,damp2jl)
         damp= dampjk*dampij*dampjl

         phi=omega(n,xyz,i,j,k,l)
         call domegadr(n,xyz,i,j,k,l,phi,dda,ddb,ddc,ddd)

         rn=vtors(3,m)
         if(rn.gt.1.d-6)then
            dphi1=phi-phi0
            c1=dphi1+pi
            x1cos=cos(c1)
            x1sin=sin(c1)
            et   =(1.+x1cos)*vtors(2,m)
            dij  =-x1sin*vtors(2,m)*damp
         else
            et =   vtors(2,m)*(cos(phi) -cos(phi0))**2
            dij=2.*vtors(2,m)* sin(phi)*(cos(phi0)-cos(phi))*damp
         endif

         term1(1:3)=et*damp2ij*dampjk*dampjl*vab(1:3)
         term2(1:3)=et*damp2jk*dampij*dampjl*vcb(1:3)
         term3(1:3)=et*damp2jl*dampij*dampjk*vdc(1:3)

         g(1:3,i)=g(1:3,i)+dij*dda(1:3)-term1
         g(1:3,j)=g(1:3,j)+dij*ddb(1:3)+term1+term2+term3
         g(1:3,k)=g(1:3,k)+dij*ddc(1:3)-term2
         g(1:3,l)=g(1:3,l)+dij*ddd(1:3)-term3
         e=e+et*damp
      endif
!        ----------------------------------------------------

   enddo

end subroutine ff_eg


! nonbonded part

subroutine ff_nonb(n,at,xyz,enb,g)
   implicit none
   integer n,at(*)
   real(wp)xyz(3,n),enb,g(3,n)

!  include 'ffparam.f90'

   integer i1,i2,iz1,iz2,k,nk
   real(wp) r2,r,r4,r6,r06,R0,t6,t8,c6t6,c6t8,t27,drij,e
   real(wp) dx,dy,dz,c6,x,alpha,oner,aa,damp,damp2,e0

   e=0
   if(nnci.lt.1)return

   do k=1,nnci
      i1=nci(1,k)
      i2=nci(2,k)
      nk=nci(3,k)
      dx=xyz(1,i1)-xyz(1,i2)
      dy=xyz(2,i1)-xyz(2,i2)
      dz=xyz(3,i1)-xyz(3,i2)
      r2=dx*dx+dy*dy+dz*dz
      r =sqrt(r2)
      oner =1.0d0/r

      iz1=at(i1)
      iz2=at(i2)
      R0=r094ff(iz1,iz2)
      c6=c6ff(i2,i1)
      r4=r2*r2
      r6=r4*r2
      r06=R0**6
      t6=r6+r06
      t8=r6*r2+r06*R0*R0
      c6t6=c6/t6
      c6t8=c6/t8
      t27=sr42ff(iz1,iz2)*c6t8
      e0=c6t6+t27
      e=e-e0*eps2(nk)
      drij=eps2(nk)*(c6t6*6.0d0*r4/t6+8.0d0*t27*r6/t8)
      g(1,i1)=g(1,i1)+dx*drij
      g(2,i1)=g(2,i1)+dy*drij
      g(3,i1)=g(3,i1)+dz*drij
      g(1,i2)=g(1,i2)-dx*drij
      g(2,i2)=g(2,i2)-dy*drij
      g(3,i2)=g(3,i2)-dz*drij

      e0=qff(i1)*qff(i2)*oner*eps1(nk)
      e=e+e0
      drij=e0/r2
      g(1,i1)=g(1,i1)-dx*drij
      g(2,i1)=g(2,i1)-dy*drij
      g(3,i1)=g(3,i1)-dz*drij
      g(1,i2)=g(1,i2)+dx*drij
      g(2,i2)=g(2,i2)+dy*drij
      g(3,i2)=g(3,i2)+dz*drij

      if(r.lt.25)then
         x    =zabff(iz1,iz2)
         alpha=r0abff(iz1,iz2)
         t27  =x*dexp(-alpha*r)
         e0   =t27*oner
         e    =e + e0*eps2(nk)
         drij=eps2(nk)*t27*(alpha*r+1.0d0)*oner/r2
         g(1,i1)=g(1,i1)-dx*drij
         g(2,i1)=g(2,i1)-dy*drij
         g(3,i1)=g(3,i1)-dz*drij
         g(1,i2)=g(1,i2)+dx*drij
         g(2,i2)=g(2,i2)+dy*drij
         g(3,i2)=g(3,i2)+dz*drij
      endif

   enddo

   enb = enb + e

end subroutine ff_nonb

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! HB/XB energy and gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine ff_hb(n,at,xyz,eh,g)
   implicit none
   integer n,at(*)
   real(wp)xyz(3,n),g(3,n),eh
!  include 'ffparam.f90'

   integer i1,i2,i,k,j
   real(wp)c1,c2,r,er,el,step,edum,dum

   if(nhb.lt.1)return

   do k=1,nhb
      i1 =hb(1,k)
      i2 =hb(2,k)
!     r  =sqrt((xyz(1,i1)-xyz(1,i2))**2 &
!        &    +(xyz(2,i1)-xyz(2,i2))**2 &
!        &    +(xyz(3,i1)-xyz(3,i2))**2)
      r  = norm2(xyz(:,i1)-xyz(:,i2))
      if(r.gt.25.0d0)cycle
      i  =hb(3,k)
      if(at(i).eq.1)then
         c1 =vhb(1,k)
         c2 =vhb(2,k)
         call eabhag(n,i1,i2,i,xyz,c1,c2,eh,g)
      else
         c1 =vhb(1,k)
         call eabxag(n,i1,i2,i,xyz,c1,eh,g)
      endif
   enddo

end subroutine ff_hb

!cccccccccccccccccccccccccccccccccccccccccccccc
! damping of bend and torsion for long
! bond distances to allow proper dissociation
!cccccccccccccccccccccccccccccccccccccccccccccc

pure subroutine abdamp(ati,atj,r2,damp,ddamp)
   use aoparam, only : rad
   implicit none
   integer, intent(in)  :: ati,atj
   real(wp),intent(in)  :: r2
   real(wp),intent(out) :: damp,ddamp
   real(wp) :: rr,rcut

   ! cut-off at about double of R_cov
   rcut =3.0 * 3.5710642*(rad(ati)+rad(atj))**2
   rr   =(r2/rcut)**2
   damp = 1.0d0/(1.0d0+rr)
   ddamp=-2.d0*2*rr/(r2*(1.0d0+rr)**2)

end subroutine abdamp

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eabh0(n,A,B,H,rab,xyz,ca,cb,energy)
   implicit none
   integer A,B,H,n
   real(wp)xyz(3,n),energy,rab

   real(wp)ang,xy,cosabh,d2ij,d2ik,d2jk,term,temp,fa,fb
   real(wp)aterm,dampm,dampl,xm,ym,zm,rhm,rab2,ca,cb,da
   real(wp):: longcut=8.
   real(wp):: alp= 12
   real(wp):: alp3= 6

   ! AB distance
   rab2=rab*rab
   ! long  damping
   dampl=1./(1.+(rab/longcut)**alp)

   ! cos angle A-B and H (or B-H H) is term
!   D2IK = (XYZ(1,A)-XYZ(1,H))**2+ &
!      &   (XYZ(2,A)-XYZ(2,H))**2+ &
!      &   (XYZ(3,A)-XYZ(3,H))**2
!   D2JK = (XYZ(1,H)-XYZ(1,B))**2+ &
!      &   (XYZ(2,H)-XYZ(2,B))**2+ &
!      &   (XYZ(3,H)-XYZ(3,B))**2
   d2ik = sum((xyz(:,A)-xyz(:,H))**2)
   d2jk = sum((xyz(:,H)-xyz(:,B))**2)
   if(d2ik.gt.d2jk)then
      xy = sqrt(rab2*d2jk)
      term = 0.5d0 * (rab2+d2jk-d2ik) / xy
   else
      xy = sqrt(rab2*d2ik)
      term = 0.5d0 * (rab2+d2ik-d2jk) / xy
   endif
   ! angle damping term
   aterm = (0.5d0*(term+1.0d0))**alp3

   ! donor-acceptor term
   fa=d2jk/d2ik
   fb=d2ik/d2jk
   da=(ca*fb+cb*fa)/(fa+fb)

   energy=-da*dampl*aterm/(rab2*rab)

end subroutine eabh0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eabh(n,A,B,H,xyz,ca,cb,energy)
   implicit none
   integer A,B,H,n
   real(wp)xyz(3,n),energy

   real(wp)rab,ang,xy,cosabh,d2ij,d2ik,d2jk,term,temp,fa,fb
   real(wp)aterm,dampm,dampl,xm,ym,zm,rhm,rab2,ca,cb,da
   real(wp):: longcut=8.
   real(wp):: alp=12
   real(wp):: alp3=6

   ! AB distance
!   rab2=   ((xyz(1,A)-xyz(1,B))**2 &
!      &    +(xyz(2,A)-xyz(2,B))**2 &
!      &    +(xyz(3,A)-xyz(3,B))**2)
   rab2 = sum((xyz(:,A)-xyz(:,B))**2) ! same as above, but easier to read
   rab=sqrt(rab2) ! norm2 intrinsic would work here
   ! long  damping
   dampl=1./(1.+(rab/longcut)**alp)

   ! cos angle A-B and H (or B-H H) is term
!   D2IK = (XYZ(1,A)-XYZ(1,H))**2+ &
!      &   (XYZ(2,A)-XYZ(2,H))**2+ &
!      &   (XYZ(3,A)-XYZ(3,H))**2
!   D2JK = (XYZ(1,H)-XYZ(1,B))**2+ &
!      &   (XYZ(2,H)-XYZ(2,B))**2+ &
!      &   (XYZ(3,H)-XYZ(3,B))**2
   d2ik = sum((xyz(:,A)-xyz(:,H))**2)
   d2jk = sum((xyz(:,H)-xyz(:,B))**2)
   if(d2ik.gt.d2jk)then
      xy = sqrt(rab2*d2jk)
      term = 0.5d0 * (rab2+d2jk-d2ik) / xy
   else
      xy = sqrt(rab2*d2ik)
      term = 0.5d0 * (rab2+d2ik-d2jk) / xy
   endif
   ! angle damping term
   aterm = (0.5d0*(term+1.0d0))**alp3

   ! donor-acceptor term
   fa=d2jk/d2ik
   fb=d2ik/d2jk
   da=(ca*fb+cb*fa)/(fa+fb)

   ! r^3 only slightly better than r^4 (SG, 8/12)
   energy=-da*dampl*aterm/(rab2*rab)

end subroutine eabh

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! charge depdendent HB/XB scaling routine
! for HB a= 10, b=5
! for XB a=-10, b=1
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine hbpara(a,b,q,c)
   implicit none
   real(wp)q,c
   real(wp)a,b

   c=exp(-a*q)/(exp(-a*q)+b)

end subroutine hbpara

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! hydrogen bonding term with analytical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eabhag(n,A,B,H,xyz,ca,cb,energy,gdr)
   implicit none
   integer A,B,H,n
   real(wp)xyz(3,n),ca,cb,energy,gdr(3,n)

   real(wp)cosabh,aterm,rdampl,da
   real(wp)rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
   real(wp)drah(3),drbh(3),drab(3)
   real(wp)ga(3),gb(3),gh(3),dg(3),dg2(3)
   real(wp)gi,denom,ratio,apref,aprod
   real(wp)eabh

   real(wp):: longcut=8.d0
   real(wp):: alp= 12.d0
   real(wp):: alp3= 6.d0

   !     compute distances
   drah(1:3)=xyz(1:3,A)-xyz(1:3,H)
   drbh(1:3)=xyz(1:3,B)-xyz(1:3,H)
   drab(1:3)=xyz(1:3,A)-xyz(1:3,B)

   !     A-B distance
   rab2=sum(drab**2)
   rab =sqrt(rab2)

   !     A-H distance
   rah2=sum(drah**2)
   rah =sqrt(rah2)

   !     B-H distance
   rbh2=sum(drbh**2)
   rbh =sqrt(rbh2)

   !     long  damping
   ratio = (rab/longcut)**alp
   rdampl=1.d0/(1.d0+ratio)
   rdampl=rdampl/rab2/rab

   !     cos angle A-B and H (or B-H H) is term
   if(rah2.gt.rbh2)then
      aprod = 1.d0/rbh/rab
      cosabh = -sum(drbh*drab)*aprod
   else
      aprod = 1.d0/rah/rab
      cosabh = sum(drah*drab)*aprod
   endif

   !     angle damping term
   aterm = 0.5d0*(cosabh+1.d0)
   apref = aterm**(alp3-1)
   aterm = aterm*apref
   apref = alp3*0.5d0*apref

   !     donor-acceptor term
   rah4 = rah2*rah2
   rbh4 = rbh2*rbh2
   denom = 1.d0/(rah4+rbh4)
   da = (ca*rah4 + cb*rbh4)*denom

   !     r^3 only slightly better than r^4 (SG, 8/12)
   eabh = -da*rdampl*aterm

   if(eabh.gt.-1.d-8) return

   energy = energy + eabh

   !     gradient
   !     donor-acceptor part
   gi = 4.d0*(ca-cb)*rah2*rbh4*denom*denom
   gi = -gi*rdampl*aterm
   ga(1:3) = gi*drah(1:3)
   gi = 4.d0*(cb-ca)*rbh2*rah4*denom*denom
   gi = -gi*rdampl*aterm
   gb(1:3) = gi*drbh(1:3)
   gh(1:3) = -ga(1:3)-gb(1:3)

   !     long-range damping part
   gi = rdampl*rdampl*rab*(3.d0+(3.d0+alp)*ratio)
   gi = gi*da*aterm
   dg(1:3) = gi*drab(1:3)
   ga(1:3) = ga(1:3) + dg(1:3)
   gb(1:3) = gb(1:3) - dg(1:3)

   !     angle part
   if(rah2.gt.rbh2)then
      dg(1:3) = -aprod*drbh(1:3) - cosabh*drab(1:3)/rab2
      dg2(1:3) = aprod*drab(1:3) + cosabh*drbh(1:3)/rbh2
      gi = -da*rdampl*apref
      dg(1:3) = gi*dg(1:3)
      dg2(1:3) = gi*dg2(1:3)

      ga(1:3) = ga(1:3) + dg(1:3)
      gh(1:3) = gh(1:3) + dg2(1:3)
      gb(1:3) = gb(1:3) - dg(1:3) - dg2(1:3)

   else
      dg(1:3) = -aprod*drah(1:3) + cosabh*drab(1:3)/rab2
      dg2(1:3) = -aprod*drab(1:3) + cosabh*drah(1:3)/rah2
      gi = -da*rdampl*apref
      dg(1:3) = gi*dg(1:3)
      dg2(1:3) = gi*dg2(1:3)

      gb(1:3) = gb(1:3) + dg(1:3)
      gh(1:3) = gh(1:3) + dg2(1:3)
      ga(1:3) = ga(1:3) - dg(1:3) - dg2(1:3)

   endif

   !     move gradients into place

   gdr(1:3,A) = gdr(1:3,A) + ga(1:3)
   gdr(1:3,B) = gdr(1:3,B) + gb(1:3)
   gdr(1:3,H) = gdr(1:3,H) + gh(1:3)

   return
end subroutine eabhag

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! XB energy routine
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eabx(n,A,B,H,xyz,ca,energy)
   implicit none
   integer A,B,H,n
   real(wp)xyz(3,n),energy,rab

   real(wp)ang,xy,cosabh,d2ij,d2ik,d2jk,term,temp,fa,fb
   real(wp)aterm,dampm,dampl,xm,ym,zm,rhm,rab2,ca,cb,da
   real(wp):: longcut=120.
   real(wp):: alp= 6
   real(wp):: alp3=6

   ! AB distance
!   rab2 = (XYZ(1,A)-XYZ(1,B))**2+ &
!      &   (XYZ(2,A)-XYZ(2,B))**2+ &
!      &   (XYZ(3,A)-XYZ(3,B))**2
   rab2 = sum((xyz(:,A)-xyz(:,B))**2) ! same as above, but easier to read
   ! long  damping
   dampl=1./(1.+(rab2/longcut)**alp)

   ! cos angle A-B and H (or B-H H) is term
!   D2IK = (XYZ(1,A)-XYZ(1,H))**2+ &
!      &   (XYZ(2,A)-XYZ(2,H))**2+ &
!      &   (XYZ(3,A)-XYZ(3,H))**2
!   D2JK = (XYZ(1,H)-XYZ(1,B))**2+ &
!      &   (XYZ(2,H)-XYZ(2,B))**2+ &
!      &   (XYZ(3,H)-XYZ(3,B))**2
   d2ik = sum((xyz(:,A)-xyz(:,H))**2)
   d2jk = sum((xyz(:,H)-xyz(:,B))**2)
   if(d2ik.gt.d2jk)then
      xy = sqrt(rab2*d2jk)
      term = 0.5d0 * (rab2+d2jk-d2ik) / xy
   else
      xy = sqrt(rab2*d2ik)
      term = 0.5d0 * (rab2+d2ik-d2jk) / xy
   endif
   ! angle damping term
   aterm = (0.5d0*(term+1.0d0))**alp3

   energy=-ca*dampl*aterm/d2jk

end subroutine eabx

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! halogen bonding term with numerical gradient
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine eabxag(n,A,B,H,xyz,ca,energy,gdr)
   implicit none
   integer A,B,H,n
   real(wp)xyz(3,n),ca,energy,gdr(3,n)

   real(wp)cosabh,aterm,rdampl
   real(wp)rab,rah,rbh,rab2,rah2,rbh2,rah4,rbh4
   real(wp)drah(3),drbh(3),drab(3)
   real(wp)ga(3),gb(3),gh(3),dg(3),dg2(3)
   real(wp)gi,denom,ratio,apref,aprod
   real(wp)eabh

   real(wp):: longcut=120.d0
   real(wp):: alp= 6.d0
   real(wp):: alp3=6.d0

   real(wp)step,er,el,dum
   integer i,j,i1,i2

   i1=A
   i2=B
   i =H
   call eabx(n,i1,i2,i,xyz,ca,er)
   energy=energy+er

   step=1.d-6
   dum=1./(2.0d0*step)
   do j=1,3
      xyz(j,i1)=xyz(j,i1)+step
      call eabx(n,i1,i2,i,xyz,ca,er)
      xyz(j,i1)=xyz(j,i1)-step*2.0d0
      call eabx(n,i1,i2,i,xyz,ca,el)
      xyz(j,i1)=xyz(j,i1)+step
      gdr(j,i1)=gdr(j,i1)+(er-el)*dum
   enddo
   do j=1,3
      xyz(j,i2)=xyz(j,i2)+step
      call eabx(n,i1,i2,i,xyz,ca,er)
      xyz(j,i2)=xyz(j,i2)-step*2.0d0
      call eabx(n,i1,i2,i,xyz,ca,el)
      xyz(j,i2)=xyz(j,i2)+step
      gdr(j,i2)=gdr(j,i2)+(er-el)*dum
   enddo
   do j=1,3
      xyz(j,i )=xyz(j,i )+step
      call eabx(n,i1,i2,i,xyz,ca,er)
      xyz(j,i )=xyz(j,i )-step*2.0d0
      call eabx(n,i1,i2,i,xyz,ca,el)
      xyz(j,i )=xyz(j,i )+step
      gdr(j,i )=gdr(j,i )+(er-el)*dum
   enddo

   return
end subroutine eabxag


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! set all bonds to Morse function if De > thr
! ie weak bonds are still LJ
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

subroutine setmorse(n,at,thr,echo)
   implicit none
   integer n,at(n)
   logical echo
!  include 'ffparam.f90'
   real(wp)thr
   integer i,k,i1,i2

   k=0
   do i=1,nbond12
      i1=bond(1,i)
      i2=bond(2,i)
      if(vbond(2,i).gt.thr) then
         bond(3,i)=1
         k=k+1
      else
         bond(3,i)=0
      endif
   enddo

   if(echo)then
      write(*,'('' LJ->Morse switch De threshold'',f8.3)') thr
      write(*,'('' switching '',i4,'' 1,2 stretch potentials'')') k
   endif

end subroutine setmorse


real(wp)function valijklff(natoms,xyz,i,j,k,l)

!  .....................................................................

implicit none

external &
   &     vecnorm,valijk

integer &
   &     ic,i,j,k,l,natoms

real(wp)&
   &     xyz(3,natoms), &
   &     eps,ra(3),rb(3),rc(3),na(3),nb(3), &
   &     rab,rbc,thab,thbc,valijk, &
   &     vecnorm,nan,nbn,rcn,snanb,deter,pi

parameter (eps=1.0d-14)
data pi/3.1415926535897932384626433832795029d0/

! ... get torsion coordinate
do ic=1,3
   ra(ic)=xyz(ic,j)-xyz(ic,i)
   rb(ic)=xyz(ic,k)-xyz(ic,j)
   rc(ic)=xyz(ic,l)-xyz(ic,k)
end do

! ... determinante of rb,ra,rc
deter= ra(1)*(rb(2)*rc(3)-rb(3)*rc(2)) &
   &  -ra(2)*(rb(1)*rc(3)-rb(3)*rc(1)) &
   &  +ra(3)*(rb(1)*rc(2)-rb(2)*rc(1))

thab=valijk(natoms,xyz,i,k,j)
thbc=valijk(natoms,xyz,j,l,k)
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

valijklff=acos(snanb)

! the gradient dphir is only compatible with this subroutine
! if the statement below is commented out. If not, opt. and
! Hessian show large errors and imags. I don't understand
! this entirely but thats how it is.
! SG, Sat May 24 11:41:42 CEST 2014

!     if (deter.lt.0) then
!        valijkl=2.d0*pi-valijkl
!     end if

end function valijklff

subroutine valelff(z)
   real(wp)at,z(*)
   integer i

   z(1:max_elem)=0
   do i=1,54
      at=float(i)
      !   -----H
      z(i)=at
      if(i.le.10.and.i.gt.2)then
         !   ------Li-F
         z(i)=at-2
      endif
      if(i.le.13.and.i.gt.10)then
         !   ------Na-Al
         z(i)=at-10
      endif
      if(i.le.18.and.i.gt.13)then
         !   ------Si-Ar
         z(i)=at-10
      endif
      if((i.le.36.and.i.ge.30).or.i.eq.19.or. i.eq.20) then
         !   ------K,Ca,Zn-Kr
         !   4s4p
         z(i)=at-18
         if(i.ge.30)z(i)=at-28
      endif
      if(i.le.29.and.i.ge.21) then
         ! Sc-Cu
         z(i)=at-18
      endif
      if(i.le.47.and.i.ge.37) then
         ! Y-Ag
         z(i)=at-36
      endif
      if(i.gt.47.and.i.le.54) then
         ! In-Xe
         z(i)=at-46
      endif
   enddo

   ! not well define for these metals
   do i=57,80
      z(i)=4.0
   enddo

   do i=81,86
      z(i)=float(i)-78.0
   enddo

   ! modifications (fit)
   z(1:2)  =z(1:2)*  2.35
   z(3:10) =z(3:10)* 0.95
   z(11:18)=z(11:18)*0.75
   z(19:54)=z(19:54)*0.65
   ! just extrapolated
   z(55:max_elem)=z(55:max_elem)*0.60

   ! special for group 1 and 2, fitted to E_int and O-M distances
   ! at TPSS-D3/def2-TZVP level
   z(3) =1.7
   z(11)=2.5
   z(19)=3.0
   z(37)=3.0
   z(55)=3.0
   z(87)=3.0

   z( 4)=5.5
   z(12)=3.0
   z(20)=2.8
   z(38)=2.6
   z(56)=2.4
   z(88)=2.4

end subroutine valelff
end module qmdff
