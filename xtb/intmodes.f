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

! project cartesian on internal mode to determine str, bend, tors
      subroutine intmodestep(n,bmat,u,step,geo,na,nb,nc,coord)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: n
      integer, intent(in)    :: na(n),nb(n),nc(n)
      real(wp),intent(in)    :: bmat(3*n-6,3*n)
      real(wp),intent(in)    :: u(3*n)
      real(wp),intent(inout) :: geo(3,n)
      real(wp),intent(inout) :: coord(3,n)
      real(wp),intent(in)    :: step

      real(wp),allocatable :: dit(:)
      real(wp),allocatable :: geo2(:,:)
      real(wp) norm
      integer k,j,kl,l,lend,n36
      logical fail

      allocate(dit(3*n-6), geo2(3,n))

      dit = 0
      do k=1,3*n
         do j=1,3*n-6
            dit(j)=dit(j)+bmat(j,k)*u(k)
         enddo
      enddo

      kl=0
      geo2=0
      do k=2,n
         lend=3
         if(k.eq.2) lend=1
         if(k.eq.3) lend=2
         do l=1,lend
            kl=kl+1
            geo(l,k)=geo(l,k)+dit(kl)*step
         enddo
      enddo
!     new cartesians
      call gmetry(n,geo,coord,na,nb,nc,fail)

      end subroutine intmodestep

! project cartesian on internal mode to determine str, bend, tors
      pure subroutine modetyp(n,bmat,u,root,vtyp)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: n
      integer, intent(in) :: root
      real(wp),intent(in) :: bmat(3*n-6,3*n)
      real(wp),intent(in) :: u(3*n,3*n)
      real(wp),intent(inout) :: vtyp(3)

      real(wp) :: norm
      real(wp),allocatable :: dit(:)
      integer  :: k,j,kl,l,lend

      allocate(dit(3*n-6))

      dit = 0
      do k=1,3*n
         do j=1,3*n-6
            dit(j)=dit(j)+bmat(j,k)*u(k,root)
         enddo
      enddo

      vtyp = 0
      kl=0
      do k=2,n
         lend=3
         if(k.eq.2) lend=1
         if(k.eq.3) lend=2
         do l=1,lend
            kl=kl+1
            vtyp(l)=vtyp(l)+dit(kl)**2
         enddo
      enddo

      norm=sum(vtyp)
      vtyp=vtyp/norm

      end subroutine modetyp

! bmatrix dzmat/dxyz
      subroutine bzmat(n,at,xyzin,bmat)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: n,at(n)
      real(wp),intent(in) :: xyzin(3,n)
      real(wp),intent(inout) :: bmat(3*n-6,3*n)

      real(wp),allocatable :: xyz(:,:),geo(:,:),br(:),bl(:)
      real(wp) :: step,one
      parameter (one=1.0d0)
      integer, allocatable :: na(:),nb(:),nc(:)
      integer :: i,k,n3,ia,ii,ic,n36,j,kl
      integer lend,l,ij
      character(len=2),external :: asym

      allocate(xyz(3,n),geo(3,n),br(3*n-6),bl(3*n-6))
      allocate(na(n),nb(n),nc(n))

      n36=3*n-6
      n3=3*n
      xyz=xyzin*0.52916790d0
      step=1.d-5

! eq. zmat
      geo = 0
      call xyzint(xyz,n,na,nb,nc,one,geo)
!     do i=1,n
!        write(*,'(2x,a2,f14.8,f14.8,f14.8,3i5)')
!    .   asym(at(i)),geo(1,i),180.*geo(2:3,i)/3.14159,na(i),nb(i),nc(i)
!     enddo

! make bmat
      ij=0
      do i=1,n
         do j=1,3
            ij=ij+1
            xyz(j,i)=xyz(j,i)+step
            call xyzgeo(xyz,n,na,nb,nc,one,geo)
            kl=0
            do k=2,n
               lend=3
               if(k.eq.2) lend=1
               if(k.eq.3) lend=2
               do l=1,lend
                  kl=kl+1
                  br(kl)=geo(l,k)
               enddo
            enddo
            xyz(j,i)=xyz(j,i)-step*2.
            call xyzgeo(xyz,n,na,nb,nc,one,geo)
            xyz(j,i)=xyz(j,i)+step
            kl=0
            do k=2,n
               lend=3
               if(k.eq.2) lend=1
               if(k.eq.3) lend=2
               do l=1,lend
                  kl=kl+1
                  bl(kl)=geo(l,k)
               enddo
            enddo
            do k=1,n36
               bmat(k,ij)=0.5*(br(k)-bl(k))/step
            enddo
         enddo
      enddo

      end


!     *****************************************************************

      subroutine makenabc(xyzin,molvec,at,nat,n,n2,nmol,
     &                    fragind,na,nb,nc)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: n,n2
      integer, intent(inout) :: na(n2), nb(n2), nc(n2)
      integer, intent(in)    :: at(n), nmol, fragind(100)
      integer, intent(inout) :: nat(n2)
      integer, intent(in)    :: molvec(n)
      real(wp),intent(in)    :: xyzin(3,n)

      integer  :: i,j,k,m,nm,ntot,mm
      real(wp) :: r,sum3(3),one
      integer, allocatable :: ind(:,:),molv(:),idum(:)
      real(wp),allocatable :: xyz(:,:),rr(:),tmp(:,:)
      parameter (one=1.0d0)

      allocate(ind(n2,n2),molv(n2),idum(n))
      allocate(xyz(3,n2),rr(n2),tmp(3,n2))

      if(nmol.gt.1)then
      k=1
      ntot=0
      do nm=1,nmol  ! loop over fragments
         m=0
         do i=k,fragind(nm)
            m=m+1
            tmp(1:3,m)=xyzin(1:3,i)
            idum(m)=at(i)
         enddo
         call cma(m,idum,tmp,sum3)
         ntot=ntot+1
         nat(ntot)=107
         molv(ntot)=nm
         xyz(1:3,ntot)=sum3(1:3)
         do i=k,fragind(nm)
            ntot=ntot+1
            xyz(1:3,ntot)=xyzin(1:3,i)
            molv(ntot)=molvec(i)
            nat(ntot)=at(i)
         enddo
         k=k+fragind(nm)
      enddo
      else
      nat=at
      molv=molvec
      xyz=xyzin
      endif

      do i=1,n2
         do j=1,n2
            r=(xyz(1,i)-xyz(1,j))**2+
     1        (xyz(2,i)-xyz(2,j))**2+
     2        (xyz(3,i)-xyz(3,j))**2
            rr(j)=r+abs(molv(i)-molv(j))*100000.
            if(i.eq.j) rr(j)=1.d+42
            ind(j,i)=j
         enddo
         call qsort(rr,1,n,ind(1,i))
      enddo
      na(2)=1
      call bonded(3,na,nb,nc,ind,n2,molv)
      na(1)=0
      nb(1)=0
      nc(1)=0
      nb(2)=0
      nc(2)=0
      nc(3)=0

      write(*,*) 'z-matrix connectivity:'
      mm=1
      do i=1,n2
         na(i)=0
         if(i.gt.1)then
            if(molv(i).ne.molv(i-1)) then
               na(i)=mm
               mm=na(i)
            endif
         endif
         call bonded(i,na,nb,nc,ind,n2,molv)
!        write(*,*)'atom ',i,'  na,nb,nc :',na(i),nb(i),nc(i)
      enddo

      call xyzint(xyz,n2,na,nb,nc,1.d0,tmp)
      call xyzgeo(xyz,n2,na,nb,nc,one,tmp)
      call zmatpr(n2,nat,tmp,na,nb,nc,1)

      end subroutine makenabc

!     *****************************************************************

      subroutine cart2zmat(xyzin,molvec,at,nat,n,n2,nmol,
     &                    fragind,na,nb,nc,geo)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: n,n2
      integer, intent(in) :: na(n2), nb(n2), nc(n2)
      integer, intent(in) :: at(n), nmol, fragind(100)
      integer, intent(in) :: molvec(n),nat(n2)
      real(wp),intent(in) :: xyzin(3,n)
      real(wp),intent(in) :: geo(3,n2)

      integer  :: i,j,k,m,nm,ntot,mm
      real(wp) :: r,sum3(3),one
      integer, allocatable :: ind(:,:),molv(:),idum(:)
      real(wp),allocatable :: xyz(:,:),rr(:),tmp(:,:)
      parameter (one=1.0d0)
      allocate(ind(n2,n2),molv(n2),idum(n))
      allocate(xyz(3,n2),rr(n2),tmp(3,n))

      if(nmol.gt.1)then
      k=1
      ntot=0
      do nm=1,nmol  ! loop over fragments
         m=0
         do i=k,fragind(nm)
            m=m+1
            tmp(1:3,m)=xyzin(1:3,i)
            idum(m)=at(i)
         enddo
         call cma(m,idum,tmp,sum3)
         ntot=ntot+1
         xyz(1:3,ntot)=sum3(1:3)
         do i=k,fragind(nm)
            ntot=ntot+1
            xyz(1:3,ntot)=xyzin(1:3,i)
         enddo
         k=k+fragind(nm)
      enddo
      else
      xyz=xyzin
      endif

      call xyzgeo(xyz,n2,na,nb,nc,one,geo)

      end subroutine cart2zmat

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      pure subroutine bonded(i,na,nb,nc,ind,n,molvec)
      use iso_fortran_env, wp => real64
      implicit none
      integer,intent(in) :: i,n,ind(n,n)
      integer,intent(in) :: molvec(n)
      integer,intent(inout) :: na(n),nb(n),nc(n)
      integer :: j,is
      integer :: icon

      nb(i)=0
      nc(i)=0
      is=0
      if(na(i).ne.0) goto 10
      do j=1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)) then
            is=j
            na(i)=icon
            goto 10
         endif
      enddo
 10   continue
      if(na(i).eq.0)then
      do j=1,n
         icon=ind(j,i)
         if(icon.lt.i) then
            is=j
            na(i)=icon
            goto 20
         endif
      enddo
 20   continue
      endif
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)
     &   .and.icon.ne.na(i)) then
            is=j
            nb(i)=icon
            goto 30
         endif
      enddo
 30   continue
      if(nb(i).eq.0)then
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i
     &   .and.icon.ne.na(i)) then
            is=j
            nb(i)=icon
            goto 40
         endif
      enddo
      endif
 40   continue
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i.and.molvec(icon).eq.molvec(i)
     &   .and.icon.ne.na(i).and.icon.ne.nb(i)) then
            is=j
            nc(i)=icon
            goto 50
         endif
      enddo
 50   continue
      if(nc(i).eq.0)then
      do j=is+1,n
         icon=ind(j,i)
         if(icon.lt.i
     &   .and.icon.ne.na(i).and.icon.ne.nb(i)) then
            is=j
            nc(i)=icon
            return
         endif
      enddo
      endif

      end subroutine bonded

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine zmat2cart(n,n2,at,geo,xyz,na,nb,nc,fail)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: n,n2
      integer, intent(in) :: na(n2),nb(n2),nc(n2)
      integer, intent(in) :: at(n2)
      real(wp),intent(inout) :: xyz(3,n)
      real(wp),intent(in) :: geo(3,n2)
      logical, intent(out) :: fail

      integer k,m
      real(wp),allocatable :: tmp(:,:)
      allocate(tmp(3,n2))

      call gmetry(n2,geo,tmp,na,nb,nc,fail)

      m=0
      do k=1,n2
         if(at(k).lt.100) then
            m=m+1
            xyz(1:3,m)=tmp(1:3,k)
         endif
      enddo

      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!     subroutine int2cart(nall,s1,s2,xyzref,xyz,n,nmol,fragind)
!     use ls_rmsd
!     implicit none
!     integer n,nall,s1,s2,fragind(100),nmol
!     real(wp) xyzref(3,n,nall),xyz(3,n)

!     integer nm,i,j,k,m,mm
!     real(wp),allocatable::tmp(:,:),tmp2(:,:)
!     real(wp) gdum(3,3), u(3,3), x(3), y(3), dum
!
!     k=1
!     mm=0
!     do nm=1,nmol  ! loop over fragments
!        m=0
!        do i=k,fragind(nm)
!           m=m+1
!        enddo
!        mm=m
!        allocate(tmp(3,mm),tmp2(3,mm))
!        m=0
!        do i=k,fragind(nm)
!           m=m+1
!           tmp (1:3,m)=xyzref(1:3,i,1)
!           tmp2(1:3,m)=xyzref(1:3,i,s1)                 ! copy bad zmatrix
!        enddo
!        call rmsd(mm,tmp,tmp2,1,u1,x,y1,dum,.false.,gdum) ! det orient
!        m=0
!        do i=k,fragind(nm)
!           m=m+1
!           tmp (1:3,m)=xyzref(1:3,i,1)
!           tmp2(1:3,m)=xyzref(1:3,i,s2)                 ! copy bad zmatrix
!        enddo
!        call rmsd(mm,tmp,tmp2,1,u2,x,y2,dum,.false.,gdum) ! det orient
!        write(*,*) nm,dum
!        call prmat(6,u,3,3,'u')
!        write(*,*) x
!        write(*,*) y
!        do i=1,mm                                       ! shift
!           tmp2(1:3,i)=tmp2(1:3,i)-y(1:3)
!        enddo
!        tmp2=matmul(u,tmp2)                             ! rotate
!        m=0
!        do i=k,fragind(nm)
!           m=m+1
!           xyz(1:3,i)=tmp2(1:3,m)+x(1:3)                ! shift back
!        enddo
!        deallocate(tmp,tmp2)
!        k=k+fragind(nm)
!     enddo

!     end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine xyzgeo(xyz,numat,na,nb,nc,degree,geo)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: numat
      real(wp),intent(in)    :: xyz(3,numat)
      integer, intent(in)    :: na(numat), nb(numat), nc(numat)
      real(wp),intent(in)    :: degree
      real(wp),intent(inout) :: geo(3,numat)
!**********************************************************************
!
!   xyzgeo converts coordinates from cartesian to internal.
!
!     on input xyz  = array of cartesian coordinates
!              numat= number of atoms
!              na   = numbers of atom to which atoms are related
!                     by distance
!              nb   = numbers of atom to which atoms are related
!                     by angle
!              nc   = numbers of atom to which atoms are related
!                     by dihedral
!
!    on output geo  = internal coordinates in angstroms, radians,
!                     and radians
!
!**********************************************************************
      integer :: i,j,k,l
      integer :: ii
      do i=2,numat
         j=na(i)
         k=nb(i)
         l=nc(i)
         if(i.lt.3) cycle
         ii=i
         call bangle(xyz,ii,j,k,geo(2,i))
         geo(2,i)=geo(2,i)*degree
         if(i.lt.4) cycle
!
!   make sure dihedral is meaninglful
!
!        call bangle(xyz,j,k,l,angl)
!        tol=0.2617994d0
!        if(angl.gt.3.1415926d0-tol.or.angl.lt.tol)then
!
!  angle is unsatisfactory, let's search for another atom for
!  defining the dihedral.
! 88       sum=100.d0
!        do 74 i1=1,ii-1
!           r=(xyz(1,i1)-xyz(1,k))**2+
!    1          (xyz(2,i1)-xyz(2,k))**2+
!    2          (xyz(3,i1)-xyz(3,k))**2
!           if(r.lt.sum.and.i1.ne.j.and.i1.ne.k) then
!         call bangle(xyz,j,k,i1,angl)
!        if(angl.lt.3.1415926d0-tol.or.angl.gt.tol)then
!              sum=r
!              l=i1
!              nc(ii)=l
!           endif
!           endif
! 74  continue
!     if(sum.gt.99.d0.and.tol.gt.0.1d0)then
!
! anything within 5 degrees?
!
!     tol=0.087266d0
!     goto 88
!     endif
!           endif
         call dihed(xyz,ii,j,k,l,geo(3,i))
         geo(3,i)=geo(3,i)*degree
      geo(1,i)= sqrt((xyz(1,i)-xyz(1,j))**2+
     1               (xyz(2,i)-xyz(2,j))**2+
     2               (xyz(3,i)-xyz(3,j))**2)
      enddo
      geo(1,1)=0.d0
      geo(2,1)=0.d0
      geo(3,1)=0.d0
      geo(2,2)=0.d0
      geo(3,2)=0.d0
      geo(3,3)=0.d0
      end subroutine xyzgeo

!     *****************************************************************
      subroutine xyzint(xyz,numat,na,nb,nc,degree,geo)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: numat
      real(wp),intent(in)    :: xyz(3,numat)
      integer, intent(inout) :: na(numat), nb(numat), nc(numat)
      real(wp),intent(in)    :: degree
      real(wp),intent(inout) :: geo(3,numat)
!**********************************************************************
!
! xyzint works out the internal coordinates of a molecule.
!        the "rules" for the connectivity are as follows:
!        atom i is defined as being at a distance from the nearest
!        atom j, atom j already having been defined.
!        atom i makes an angle with atom j and the atom k, which has
!        already been defined, and is the nearest atom to j
!        atom i makes a dihedral angle with atoms j, k, and l. l having
!        been defined and is the nearest atom to k, and j, k and l
!        have a contained angle in the range 15 to 165 degrees,
!        if possible.
!
!        note that geo and xyz must not be the same in the call.
!
!   on input xyz    = cartesian array of numat atoms
!            degree = 1 if angles are to be in radians
!            degree = 57.29578 if angles are to be in radians
!
!**********************************************************************
      integer  :: nai1,nai2
      integer  :: i,j,k
      integer  :: im1
      real(wp) :: r,sum
      nai1=0
      nai2=0
      do i=1,numat
         na(i)=2
         nb(i)=3
         nc(i)=4
         im1=i-1
         if(im1.eq.0) cycle
         sum=100.d0
         do j=1,im1
            r=(xyz(1,i)-xyz(1,j))**2+
     1          (xyz(2,i)-xyz(2,j))**2+
     2          (xyz(3,i)-xyz(3,j))**2
            if(r.lt.sum.and.na(j).ne.j.and.nb(j).ne.j) then
               sum=r
               k=j
            endif
         enddo
!
!   atom i is nearest to atom k
!
         na(i)=k
         if(i.gt.2)nb(i)=na(k)
         if(i.gt.3)nc(i)=nb(k)
!
!   find any atom to relate to na(i)
!
      enddo
      na(1)=0
      nb(1)=0
      nc(1)=0
      nb(2)=0
      nc(2)=0
      nc(3)=0

      call xyzgeo(xyz,numat,na,nb,nc,degree,geo)
      end subroutine xyzint

      pure subroutine bangle(xyz,i,j,k,angle)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in) :: i,j,k
      real(wp),intent(in) :: xyz(3,*)
      real(wp),intent(out) :: angle
!********************************************************************
!
! bangle calculates the angle between atoms i,j, and k. the
!        cartesian coordinates are in xyz.
!
!********************************************************************
      real(wp) :: d2ij,d2jk,d2ik,xy,temp
      d2ij = (xyz(1,i)-xyz(1,j))**2+
     1       (xyz(2,i)-xyz(2,j))**2+
     2       (xyz(3,i)-xyz(3,j))**2
      d2jk = (xyz(1,j)-xyz(1,k))**2+
     1       (xyz(2,j)-xyz(2,k))**2+
     2       (xyz(3,j)-xyz(3,k))**2
      d2ik = (xyz(1,i)-xyz(1,k))**2+
     1       (xyz(2,i)-xyz(2,k))**2+
     2       (xyz(3,i)-xyz(3,k))**2
      xy = sqrt(d2ij*d2jk+1.d-14)
      temp = 0.5d0 * (d2ij+d2jk-d2ik) / xy
      if (temp .gt. 1.0d0) temp=1.0d0
      if (temp .lt. -1.0d0) temp=-1.0d0
      angle = acos( temp )
      end subroutine bangle
      subroutine dihed(xyz,i,j,k,l,angle)
      use iso_fortran_env, wp => real64
      implicit none
      integer, intent(in)    :: i,j,k,l
      real(wp),intent(inout) :: angle
      real(wp),intent(in)    :: xyz(3,*)
!********************************************************************
!
!      dihed calculates the dihedral angle between atoms i, j, k,
!            and l.  the cartesian coordinates of these atoms
!            are in array xyz.
!
!     dihed is a modified version of a subroutine of the same name
!           which was written by dr. w. theil in 1973.
!
!********************************************************************
      real(wp),parameter :: pi = 3.14159265358979d0
      real(wp) :: xi1,xj1,xl1,xi2,xj2,xl2
      real(wp) :: yi1,yj1,yl1,yi2,yj2,yl2,yi3,yl3
      real(wp) :: zi1,zj1,zl1
      real(wp) :: dist,yxdist
      real(wp) :: cosa,sinth,costh,cosph,sinph
      real(wp) :: ddd
      xi1=xyz(1,i)-xyz(1,k)
      xj1=xyz(1,j)-xyz(1,k)
      xl1=xyz(1,l)-xyz(1,k)
      yi1=xyz(2,i)-xyz(2,k)
      yj1=xyz(2,j)-xyz(2,k)
      yl1=xyz(2,l)-xyz(2,k)
      zi1=xyz(3,i)-xyz(3,k)
      zj1=xyz(3,j)-xyz(3,k)
      zl1=xyz(3,l)-xyz(3,k)
!      rotate around z axis to put kj along y axis
      dist= sqrt(xj1**2+yj1**2+zj1**2+1.d-14)
      cosa=zj1/dist
      if(cosa.gt.1.0d0) cosa=1.0d0
      if(cosa.lt.-1.0d0) cosa=-1.0d0
      ddd=1.0d0-cosa**2
      if(ddd.le.0.0) go to 10
      yxdist=dist* sqrt(ddd)
      if(yxdist.gt.1.0d-12) go to 20
   10 continue
      xi2=xi1
      xl2=xl1
      yi2=yi1
      yl2=yl1
      costh=cosa
      sinth=0.d0
      go to 30
   20 cosph=yj1/yxdist
      sinph=xj1/yxdist
      xi2=xi1*cosph-yi1*sinph
      xj2=xj1*cosph-yj1*sinph
      xl2=xl1*cosph-yl1*sinph
      yi2=xi1*sinph+yi1*cosph
      yj2=xj1*sinph+yj1*cosph
      yl2=xl1*sinph+yl1*cosph
!      rotate kj around the x axis so kj lies along the z axis
      costh=cosa
      sinth=yj2/dist
   30 continue
      yi3=yi2*costh-zi1*sinth
      yl3=yl2*costh-zl1*sinth
      call dang(xl2,yl3,xi2,yi3,angle)
!     if (angle .lt. 0.) angle=2.0d0*pi+angle
!sg   if (angle .ge. 2.0d0*pi    ) angle=0.d0
      end subroutine dihed
      pure subroutine dang(a1,a2,b1,b2,rcos)
      use iso_fortran_env, wp => real64
      implicit none
      real(wp),intent(inout) :: a1,a2,b1,b2
      real(wp),intent(out)   :: rcos
!*********************************************************************
!
!    dang  determines the angle between the points (a1,a2), (0,0),
!          and (b1,b2).  the result is put in rcos.
!
!*********************************************************************
      real(wp),parameter :: pi = 3.14159265358979d0
      real(wp) :: zero
      real(wp) :: anorm,bnorm
      real(wp) :: sinth,costh
      zero=1.0d-10
      if( abs(a1).lt.zero.and. abs(a2).lt.zero) go to 10
      if( abs(b1).lt.zero.and. abs(b2).lt.zero) go to 10
      anorm=1.0d0/ sqrt(a1**2+a2**2)
      bnorm=1.0d0/ sqrt(b1**2+b2**2)
      a1=a1*anorm
      a2=a2*anorm
      b1=b1*bnorm
      b2=b2*bnorm
      sinth=(a1*b2)-(a2*b1)
      costh=a1*b1+a2*b2
      if(costh.gt.1.0d0) costh=1.0d0
      if(costh.lt.-1.0d0) costh=-1.0d0
      rcos= acos(costh)
      if( abs(rcos).lt.zero) go to 10
      if(sinth.gt.0.d0) rcos=2.0d0*pi-rcos
      rcos=-rcos
      return
   10 rcos=0.0d0
      end subroutine dang

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      pure subroutine gmetry(natoms, geo, coord, na,nb,nc,fail)
      use iso_fortran_env, wp => real64
      implicit none

      integer, intent(in)    :: natoms
      real(wp),intent(inout) :: coord(3,natoms)
      real(wp),intent(in)    :: geo(3,natoms)
      integer, intent(in)    :: na(natoms), nb(natoms), nc(natoms)
      logical, intent(out)   :: fail
      real(wp) :: ccos
      real(wp) :: xa,ya,za
      real(wp) :: xb,yb,zb
      real(wp) :: xd,yd,zd
      real(wp) :: rbc,xpa,xpb,xyb
      real(wp) :: ypa,xqa,zqa,yza
      real(wp) :: ypd,zpd,xpd,zqd,xqd,yqd,xrd
      real(wp) :: costh,sinth,sinph,cosph
      real(wp) :: coskh,sinkh,sina,cosa,sind,cosd
      integer  :: i,k
      integer  :: ma,mb,mc

      fail=.false.

      coord(1,1)=0.0d00
      coord(2,1)=0.0d00
      coord(3,1)=0.0d00
      coord(1,2)=geo(1,2)
      coord(2,2)=0.0d00
      coord(3,2)=0.0d00
      if(natoms.eq.2) goto 110
      ccos=cos(geo(2,3))
      if(na(3).eq.1)then
         coord(1,3)=coord(1,1)+geo(1,3)*ccos
      else
         coord(1,3)=coord(1,2)-geo(1,3)*ccos
      endif
      coord(2,3)=geo(1,3)*sin(geo(2,3))
      coord(3,3)=0.0d00
      do 100 i=4,natoms
         cosa=cos(geo(2,i))
         mb=nb(i)
         mc=na(i)
         xb=coord(1,mb)-coord(1,mc)
         yb=coord(2,mb)-coord(2,mc)
         zb=coord(3,mb)-coord(3,mc)
         rbc=1.0d00/dsqrt(xb*xb+yb*yb+zb*zb)
         if (abs(cosa).lt.0.9999999999d0) go to 40
!
!     atoms mc, mb, and (i) are collinear
!
         rbc=geo(1,i)*rbc*cosa
         coord(1,i)=coord(1,mc)+xb*rbc
         coord(2,i)=coord(2,mc)+yb*rbc
         coord(3,i)=coord(3,mc)+zb*rbc
         go to 100
!
!     the atoms are not collinear
!
   40    ma=nc(i)
         xa=coord(1,ma)-coord(1,mc)
         ya=coord(2,ma)-coord(2,mc)
         za=coord(3,ma)-coord(3,mc)
!
!     rotate about the z-axis to make yb=0, and xb positive.  if xyb is
!     too small, first rotate the y-axis by 90 degrees.
!
         xyb=dsqrt(xb*xb+yb*yb)
         k=-1
         if (xyb.gt.0.1d00) go to 50
         xpa=za
         za=-xa
         xa=xpa
         xpb=zb
         zb=-xb
         xb=xpb
         xyb=dsqrt(xb*xb+yb*yb)
         k=1
!
!     rotate about the y-axis to make zb vanish
!
   50    costh=xb/xyb
         sinth=yb/xyb
         xpa=xa*costh+ya*sinth
         ypa=ya*costh-xa*sinth
         sinph=zb*rbc
         cosph=dsqrt(abs(1.d00-sinph*sinph))
         xqa=xpa*cosph+za*sinph
         zqa=za*cosph-xpa*sinph
!
!     rotate about the x-axis to make za=0, and ya positive.
!
         yza=dsqrt(ypa*ypa+zqa*zqa)
         if(yza.lt.1.d-4 )then
            if(yza.lt.1.d-4)goto 70
!           write(*,'(/9x,'' atoms'',i3,'','',i3,'', and'',i3,
!    1'' are within'',f7.4,'' angstroms of a straight line'')')
!    2mc,mb,ma,yza
            fail=.true.
            return
         endif
         coskh=ypa/yza
         sinkh=zqa/yza
         goto 80
   70    continue
!
!   angle too small to be important
!
         coskh=1.d0
         sinkh=0.d0
   80    continue
!
!     coordinates :-   a=(xqa,yza,0),   b=(rbc,0,0),  c=(0,0,0)
!     none are negative.
!     the coordinates of i are evaluated in the new frame.
!
         sina=sin(geo(2,i))
         sind=-sin(geo(3,i))
         cosd=cos(geo(3,i))
         xd=geo(1,i)*cosa
         yd=geo(1,i)*sina*cosd
         zd=geo(1,i)*sina*sind
!
!     transform the coordinates back to the original system.
!
         ypd=yd*coskh-zd*sinkh
         zpd=zd*coskh+yd*sinkh
         xpd=xd*cosph-zpd*sinph
         zqd=zpd*cosph+xd*sinph
         xqd=xpd*costh-ypd*sinth
         yqd=ypd*costh+xpd*sinth
         if (k.lt.1) go to 90
         xrd=-zqd
         zqd=xqd
         xqd=xrd
   90    coord(1,i)=xqd+coord(1,mc)
         coord(2,i)=yqd+coord(2,mc)
         coord(3,i)=zqd+coord(3,mc)
  100 continue
  110 continue

      end subroutine gmetry

      pure subroutine cma(nat,at,xyz,sum3)
      use iso_fortran_env, wp => real64
! atomic masses
      use splitparam
      implicit none
      integer, intent(in)  :: nat
      real(wp),intent(in)  :: xyz(3,nat)
      integer, intent(in)  :: at(nat)
      real(wp),intent(out) :: sum3(3)
      real(wp) :: sumw,sumwx,sumwy,sumwz
      integer  :: i

      sumw=0.0_wp
      sumwx=0.0_wp
      sumwy=0.0_wp
      sumwz=0.0_wp

      do i=1,nat
         sumw=sumw+atmass(i)
         sumwx=sumwx+atmass(i)*xyz(1,i)
         sumwy=sumwy+atmass(i)*xyz(2,i)
         sumwz=sumwz+atmass(i)*xyz(3,i)
      enddo

      sum3(1)=sumwx/sumw
      sum3(2)=sumwy/sumw
      sum3(3)=sumwz/sumw

      end subroutine cma
