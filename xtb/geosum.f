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

      subroutine geosum(n,ic,xyzin)
      use aoparam
      use splitparam
      implicit none
      integer n,ic(n)
      real*8 xyzin(3,n)
      character*2 asym                 

      real*8 ,allocatable :: x(:),y(:),z(:),rmin(:),rmax(:),dev(:)
      real*8 ,allocatable :: rlist(:,:),av(:) 
      integer,allocatable :: bondlist(:,:),nbond(:),nty(:)

      real*8 brad,r,xxx,yyy,zzz,avv
      real*8 xyz(3,3),ang,xyz4(3,4),ra(3),rb(3)
      integer i, j, nn, k, imin, jmin, iat1, iat2
      integer nm,ij,m,iat,jat, nsel
      integer nt,a,b,ia,ib,ja,jb,ii,jj,lin,linij
      character*31  lab,lab2,angsum1 (n*150)
      character*30  lab4,angsum2(n*150)
      character*50  lab3
      character*21  lab21
                                                        
c first nsel atoms are considered in angle/dihedreal printout
      nsel=min(10,n)

      nn=100*100                                                    
      allocate(x(n),y(n),z(n),rmin(nn),rmax(nn),dev(nn),av(nn),
     .         rlist(n,150),bondlist(n,150),nbond(n),nty(nn))
                                                        
      x(1:n)=xyzin(1,1:n)
      y(1:n)=xyzin(2,1:n)
      z(1:n)=xyzin(3,1:n)
      rmin=10000.
      k   =0
      do i=1,n
         nbond(i)=0
         do j=1,n
            if(i.eq.j)cycle
            k=k+1
            r=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
     .        *0.529177260d0
c the max bond dist in angst
            brad=rad(ic(i))+rad(ic(j))
            brad=brad+0.30*brad 
            if(r.le.brad) then
               nbond(i)=nbond(i)+1
               bondlist(i,nbond(i))=j
               rlist   (i,nbond(i))=r
            endif
         enddo
      enddo

      write(*,*)
      write(*,*)'Bond Distances (Angstroems)'
      write(*,*)'---------------------------'
      nt=0
      nty=0
      av=0
      dev=0
      rmin=1000
      rmax=-1
      do i=1,n
         do a=1,nbond(i)
            nt=nt+1
            write(lab21,'(a2,i3,''-'',a2,i3)')
     .      asym(ic(i)),i,asym(ic(bondlist(i,a))),bondlist(i,a) 
            call rmblank(lab21,lab4)
            write(angsum1(nt),'(a,''='',f6.4)')
     .      trim(lab4),rlist(i,a)
            iat1=ic(i)
            iat2=ic(bondlist(i,a))
            linij=lin(iat1,iat2)
            av(linij)=av(linij)+rlist(i,a)
            if(rlist(i,a).lt.rmin(linij))rmin(linij)=rlist(i,a)
            if(rlist(i,a).gt.rmax(linij))rmax(linij)=rlist(i,a)
            nty(linij)=nty(linij)+1
         enddo
      enddo
      
      do i=1,n
         do a=1,nbond(i)
            iat1=ic(i)
            iat2=ic(bondlist(i,a))
            linij=lin(iat1,iat2)
            avv=av(linij)/float(nty(linij))
            dev(linij)=dev(linij)+(avv-rlist(i,a))**2
         enddo
      enddo
      write(*,'(6a21)')(angsum1(i),i=1,nt)
      do i=1,100
         do j=1,i
            if(nty(lin(i,j)).gt.0)then
            write(*,'(2A3,'' Rav='',f6.4,'' sigma='',f6.4,
     .      ''  Rmin='',f6.4,''  Rmax='',f6.4,i6)')
     .      asym(i),asym(j),
     .      av(lin(i,j))/nty(lin(i,j)),
     .      sqrt(dev(lin(i,j))/float(nty(lin(i,j)))),
     .      rmin(lin(i,j)),rmax(lin(i,j)),int(0.5*nty(lin(i,j))+0.0001)
            endif
         enddo
      enddo

      write(*,*)
      write(*,*)'selected bond angles (degree)'
      write(*,*)'--------------------'
      nt=0
      do i=1,nsel
         m=nbond(i)
         if(m.gt.1)then
         xyz(1,2)=x(i)
         xyz(2,2)=y(i)
         xyz(3,2)=z(i)
         do a=1,m
            ia=bondlist(i,a)
            do b=1,a-1
               ib=bondlist(i,b)
               nt=nt+1
               write(lab,'(a2,i5,''-'',a2,i5,''-'',a2,i5)')
     .asym(ic(ia)),ia,asym(ic(i)),i,asym(ic(ib)),ib
               call rmblank(lab,lab2)
               xyz(1,1)=x(ia)
               xyz(2,1)=y(ia)
               xyz(3,1)=z(ia)
               xyz(1,3)=x(ib)
               xyz(2,3)=y(ib)
               xyz(3,3)=z(ib)
               call xbangle(xyz,ang)
               write(angsum1(nt),'(a,''='',f6.2)')trim(lab2),ang
            enddo
         enddo
         endif
      enddo
      write(*,'(4a)')(angsum1(i),i=1,nt)

      write(*,*)
      write(*,*)'selected dihedral angles (degree)'
      write(*,*)'---------------------------------'
      nt=0
      do i=1,nsel
         xyz4(1,2)=x(i)
         xyz4(2,2)=y(i)
         xyz4(3,2)=z(i)
         do j=1,i
            if(i.eq.j)cycle
            xyz4(1,3)=x(j)
            xyz4(2,3)=y(j)
            xyz4(3,3)=z(j)
            k=k+1
            r=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2+(z(i)-z(j))**2)
     .        *0.529177260d0
c the max bond dist in angst
            brad=rad(ic(i))+rad(ic(j))
            brad=brad+0.25*brad 
            if(r.le.brad) then
               do ii=1,nbond(i)
                  ib=bondlist(i,ii)
                  do jj=1,nbond(j)
                     jb=bondlist(j,jj)
                     if(ib.ne.jb.and.ib.ne.j.and.jb.ne.i)then
                        nt=nt+1
                        write(lab3,'(a2,i5,''-'',a2,i5,''-'',
     .                               a2,i5,''-'',a2,i5)')
     .                  asym(ic(ib)),ib,asym(ic(i)),i,
     .                  asym(ic(j)),j,  asym(ic(jb)),jb
                        call rmblank(lab3,lab4)
                        xyz4(1,1)=x(ib)
                        xyz4(2,1)=y(ib)
                        xyz4(3,1)=z(ib)
                        xyz4(1,4)=x(jb)
                        xyz4(2,4)=y(jb)
                        xyz4(3,4)=z(jb)
                        call xdihed(xyz4,ang)
                        write(angsum2(nt),'(a,''='',f6.2)')
     .                  trim(lab4),ang
                     endif
                  enddo
               enddo
            endif
         enddo
      enddo
      write(*,'(4a)')(angsum2(i),i=1,nt)

      if(iatf2.ne.0.and.iatf1.ne.0)then
      write(*,*)
      write(*,*)'CMA Distance (Angstroems)'
      write(*,*)'---------------------------'
      call cmafrag(n,ic,xyzin,ra,rb)
      write(*,'(''R(CMA):'',f8.4)')rcma*0.529177260d0
      endif

      end

C-------------------------------------------------

      SUBROUTINE RMBLANK(AS,RE)
      CHARACTER*(*) AS
      CHARACTER*(*) RE
  
      RE=' '
      K=0
      NSP=ICHAR(' ')
      DO 10 I=1,len(as)
         J=ICHAR(AS(I:I))
         IF(J.NE.NSP)THEN  
            K=K+1
            RE(K:K)=AS(I:I)
         ENDIF
  10  CONTINUE
      RETURN
      END

      SUBROUTINE XBANGLE(XYZ,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,3)
*********************************************************************
*
* BANGLE CALCULATES THE ANGLE BETWEEN ATOMS I,J, AND K. THE
*        CARTESIAN COORDINATES ARE IN XYZ.
*
*********************************************************************
      i=1
      j=2
      k=3
      D2IJ = (XYZ(1,I)-XYZ(1,J))**2+
     1       (XYZ(2,I)-XYZ(2,J))**2+
     2       (XYZ(3,I)-XYZ(3,J))**2
      D2JK = (XYZ(1,J)-XYZ(1,K))**2+
     1       (XYZ(2,J)-XYZ(2,K))**2+
     2       (XYZ(3,J)-XYZ(3,K))**2
      D2IK = (XYZ(1,I)-XYZ(1,K))**2+
     1       (XYZ(2,I)-XYZ(2,K))**2+
     2       (XYZ(3,I)-XYZ(3,K))**2
      XY = SQRT(D2IJ*D2JK)
      TEMP = 0.5D0 * (D2IJ+D2JK-D2IK) / XY
      IF (TEMP .GT. 1.0D0) TEMP=1.0D0
      IF (TEMP .LT. -1.0D0) TEMP=-1.0D0
      ANGLE = ACOS( TEMP )*180.0d0/3.14159265358979d0
      RETURN
      END
      SUBROUTINE XDIHED(XYZ,ANGLE)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION XYZ(3,4)
*********************************************************************
*
*      DIHED CALCULATES THE DIHEDRAL ANGLE BETWEEN ATOMS I, J, K,
*            AND L.  THE CARTESIAN COORDINATES OF THESE ATOMS
*            ARE IN ARRAY XYZ.
*
*     DIHED IS A MODIFIED VERSION OF A SUBROUTINE OF THE SAME NAME
*           WHICH WAS WRITTEN BY DR. W. THEIL IN 1973.
*
*********************************************************************
      DATA PI/3.14159265358979D0/
      i=1
      j=2
      k=3
      l=4
      XI1=XYZ(1,I)-XYZ(1,K)
      XJ1=XYZ(1,J)-XYZ(1,K)
      XL1=XYZ(1,L)-XYZ(1,K)
      YI1=XYZ(2,I)-XYZ(2,K)
      YJ1=XYZ(2,J)-XYZ(2,K)
      YL1=XYZ(2,L)-XYZ(2,K)
      ZI1=XYZ(3,I)-XYZ(3,K)
      ZJ1=XYZ(3,J)-XYZ(3,K)
      ZL1=XYZ(3,L)-XYZ(3,K)
C      ROTATE AROUND Z AXIS TO PUT KJ ALONG Y AXIS
      DIST= SQRT(XJ1**2+YJ1**2+ZJ1**2)
      COSA=ZJ1/DIST
      IF(COSA.GT.1.0D0) COSA=1.0D0
      IF(COSA.LT.-1.0D0) COSA=-1.0D0
      DDD=1.0D0-COSA**2
      IF(DDD.LE.0.0) GO TO 10
      YXDIST=DIST* SQRT(DDD)
      IF(YXDIST.GT.1.0D-9) GO TO 20
   10 CONTINUE
      XI2=XI1
      XL2=XL1
      YI2=YI1
      YL2=YL1
      COSTH=COSA
      SINTH=0.D0
      GO TO 30
   20 COSPH=YJ1/YXDIST
      SINPH=XJ1/YXDIST
      XI2=XI1*COSPH-YI1*SINPH
      XJ2=XJ1*COSPH-YJ1*SINPH
      XL2=XL1*COSPH-YL1*SINPH
      YI2=XI1*SINPH+YI1*COSPH
      YJ2=XJ1*SINPH+YJ1*COSPH
      YL2=XL1*SINPH+YL1*COSPH
C      ROTATE KJ AROUND THE X AXIS SO KJ LIES ALONG THE Z AXIS
      COSTH=COSA
      SINTH=YJ2/DIST
   30 CONTINUE
      YI3=YI2*COSTH-ZI1*SINTH
      YL3=YL2*COSTH-ZL1*SINTH
      CALL XDANG(XL2,YL3,XI2,YI3,ANGLE)
      IF (ANGLE .LT. 0.) ANGLE=2.0D0*PI+ANGLE
      IF (ANGLE .GE. 2.0d0*PI    ) ANGLE=0.D0
      ANGLE=ANGLE*180.0d0/pi
      RETURN
      END
      SUBROUTINE XDANG(A1,A2,B1,B2,RCOS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
**********************************************************************
*
*    DANG  DETERMINES THE ANGLE BETWEEN THE POINTS (A1,A2), (0,0),
*          AND (B1,B2).  THE RESULT IS PUT IN RCOS.
*
**********************************************************************
      DATA PI/3.14159265358979D0/
      ZERO=1.0D-6
      IF( ABS(A1).LT.ZERO.AND. ABS(A2).LT.ZERO) GO TO 10
      IF( ABS(B1).LT.ZERO.AND. ABS(B2).LT.ZERO) GO TO 10
      ANORM=1.0D0/ SQRT(A1**2+A2**2)
      BNORM=1.0D0/ SQRT(B1**2+B2**2)
      A1=A1*ANORM
      A2=A2*ANORM
      B1=B1*BNORM
      B2=B2*BNORM
      SINTH=(A1*B2)-(A2*B1)
      COSTH=A1*B1+A2*B2
      IF(COSTH.GT.1.0D0) COSTH=1.0D0
      IF(COSTH.LT.-1.0D0) COSTH=-1.0D0
      RCOS= ACOS(COSTH)
      IF( ABS(RCOS).LT.4.0D-4) GO TO 10
      IF(SINTH.GT.0.D0) RCOS=2.0D0*PI-RCOS
      RCOS=-RCOS
      RETURN
   10 RCOS=0.0D0
      RETURN
      END

