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
! You should have received a copy of the GNU Lesser General Public Licen
! along with xtb.  If not, see <https://www.gnu.org/licenses/>.

      subroutine preig(io,occ,f,e,istart,norbs)
      use xtb_mctc_accuracy, only : wp
      implicit none
      integer, intent(in) :: io
      integer, intent(in) :: istart
      real(wp),intent(in) :: f
      integer, intent(in) :: norbs
      real(wp),intent(in) :: occ(norbs)
      real(wp),intent(in) :: e(norbs)
      integer :: n,ntimes,nrest,j,n2,k,i

      write(io,'(/,10x,''eigenvalues'')')
      n=8
      ntimes=(norbs-istart+1)/n
      nrest=mod(norbs-istart+1,n)
      if(ntimes.eq.0) nrest=norbs-istart+1

      j=istart
      n2=n+istart-1

      do k=1,ntimes
      write(io,100)(i,i=j,n2)
      write(io,200)(occ(i),i=j,n2)
      write(io,300)(f*e(i),i=j,n2)
      j =j +n
      n2=n2+n
      enddo

      if(nrest.gt.0.or.ntimes.eq.0) then
      write(io,100)(i,i=j,j+nrest-1)
      write(io,200)(occ(i),i=j,j+nrest-1)
      write(io,300)(f*e(i),i=j,j+nrest-1)
      endif

      return

 100  format(' #    : ',2x,8(3x,i6,2x))
 200  format(' occ. : ',2x,8(4x,f6.3,1x))
 300  format(' eps  : ',2x,8f11.3)
      end subroutine preig

      SUBROUTINE PREIG2(IO,OCC,E,NORBS)
      use xtb_mctc_accuracy, only : wp
      implicit none
      integer, intent(in) :: io
      integer, intent(in) :: norbs
      real(wp),intent(in) :: occ(norbs)
      real(wp),intent(in) :: e(norbs)
      integer :: n,ntimes,nrest,j,n2,k,i

      n=6
      ntimes=norbs/n
      nrest=mod(norbs,n)
      if(ntimes.eq.0) nrest=norbs

      j=1
      n2=n

      do k=1,ntimes
      write(io,100)(i,i=j,n2)
      write(io,200)(occ(i),i=j,n2)
      write(io,300)(e(i),i=j,n2)
      j =j +n
      n2=n2+n
      enddo

      if(nrest.gt.0.or.ntimes.eq.0) then
      write(io,100)(i,i=j,j+nrest-1)
      write(io,200)(occ(i),i=j,j+nrest-1)
      write(io,300)(e(i),i=j,j+nrest-1)
      endif

      return

 100  format('#       :',2x,6(4x,i4,2x))
 200  format('# atoms :',2x,6(4x,f5.3,1x))
 300  format('shift ev:',2x,6f10.5)
      end subroutine preig2

      SUBROUTINE PREIG3(IO,E,NORBS)
      IMPLICIT INTEGER (I-N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION E(*)
      N=10
      NTIMES=NORBS/N
      NREST=MOD(NORBS,N)
      IF(NTIMES.EQ.0) NREST=NORBS
      J=1
      N2=N
      DO K=1,NTIMES
!     WRITE(IO,100)(I,I=J,N2)
      WRITE(IO,300)J,N2,(E(i),I=J,N2)
      J =J +N
      N2=N2+N
      ENDDO
      IF(NREST.GT.0.OR.NTIMES.EQ.0) THEN
      WRITE(IO,300)J,J+NREST-1,(E(i),I=J,J+NREST-1)
      ENDIF

 100  FORMAT('atoms : ',2X,6(4X,I4,2X))
 300  FORMAT(' value',i5,'-',i5,':',2X,12F6.2)
      END SUBROUTINE PREIG3

      SUBROUTINE PRMAT(IUOUT,R,N,M,HEAD)
      IMPLICIT INTEGER(I-N)
      REAL*8 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
!     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
!     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
!     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
!
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMAT4(IUOUT,R,N,M,HEAD)
      IMPLICIT INTEGER(I-N)
      REAL*4 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
!     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
!     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
!     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
!
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMATI(IUOUT,RR,N,M,HEAD)
      IMPLICIT INTEGER(I-N)
      CHARACTER*(*) HEAD
      integer RR(*)
      REAL*4 R
      DIMENSION R(n*n)
!     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
!     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
!     ((N+1)*N)/2 WHEN M IS ZERO

      R(1:n*n)=float(RR(1:n*n))
      WRITE(IUOUT,1001) HEAD
      NKPB=6
      IF(M)10,10,80
!
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',4X,6(3X,I4,3X),/)
 1003 FORMAT(' ',I4,6F10.5)
      END

      SUBROUTINE PRMATS(IUOUT,R,N,M,HEAD)
      IMPLICIT INTEGER(I-N)
      REAL*8 R
      CHARACTER*(*) HEAD
      DIMENSION R(*)
!     SUBROUTINE PRINTS MATRIX R,WHICH IS SUPPOSED
!     TO HAVE DIMENSION N,M  WHEN M IS NONZERO AND
!     ((N+1)*N)/2 WHEN M IS ZERO

      WRITE(IUOUT,1001) HEAD
      NKPB=12
      IF(M)10,10,80
!
   10 CONTINUE
      IBL=N/NKPB
      IR=N-IBL*NKPB
      J1=1
      K1S=1
      KD=0
      IF(IBL.EQ.0) GO TO 50
      J2=NKPB
      DO 40 I=1,IBL
      WRITE(IUOUT,1002)(J,J=J1,J2)
      K1=K1S
      K2=K1
      KK=0
      DO 20 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   20 K2=K1+KK
      J1=J1+NKPB
      IF(J1.GT.N) RETURN
      J2=J2+NKPB
      K2=K1-1
      K1=K2+1
      K2=K1+(NKPB-1)
      K1S=K2+1
      KK=KD+NKPB
      DO 30 J=J1,N
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KK
   30 K2=K2+KK
   40 KD=KD+NKPB
   50 IF(IR.EQ.0) GO TO 70
      K1=K1S
      J2=J1+IR-1
      KK=0
      K2=K1
      WRITE(IUOUT,1002)(J,J=J1,J2)
      WRITE(IUOUT,1003)
      DO 60 J=J1,J2
      WRITE(IUOUT,1003)J,(R(K),K=K1,K2)
      KK=KK+1
      K1=K1+KD+KK
   60 K2=K1+KK
   70 RETURN
   80 IBL=M/NKPB
      IR=M-IBL*NKPB
      I2=0
      K2=0
      IF(IBL.EQ.0) GO TO 100
      DO 90 I=1,IBL
      I1=(I-1)*N*NKPB+1
      I2=I1+(NKPB-1)*N
      K1=K2+1
      K2=K1+(NKPB-1)
      WRITE(IUOUT,1002)(K,K=K1,K2)
      DO 90 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
   90 I2=I1+(NKPB-1)*N
  100 IF(IR.EQ.0) GO TO 120
      I1=IBL*N*NKPB+1
      I2=I1+(IR-1)*N
      K1=K2+1
      K2=M
      WRITE(IUOUT,1002)(K,K=K1,K2)
      WRITE(IUOUT,1003)
      DO 110 J=1,N
      WRITE(IUOUT,1003)J,(R(IJ),IJ=I1,I2,N)
      I1=I1+1
      I2=I1+(IR-1)*N
  110 CONTINUE
  120 WRITE(IUOUT,1003)
      RETURN
 1001 FORMAT(/' MATRIX PRINTED:',2X,A)
 1002 FORMAT(/,' ',2X,12(2X,I3,1X),/)
 1003 FORMAT(' ',I4,12F6.2)
      END

      subroutine preigf(io,e,norbs)
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      dimension e(*)
      n=6
      ntimes=norbs/n
      nrest=mod(norbs,n)
      if(ntimes.eq.0) nrest=norbs
      j=1
      n2=n
      do k=1,ntimes
      write(io,300)     (e(i),i=j,n2)
      j =j +n
      n2=n2+n
      enddo
      if(nrest.gt.0.or.ntimes.eq.0) then
      write(io,300)            (e(i),i=j,j+nrest-1)
      endif

 100  format(' #     : ',2x,6(4x,i4,2x))
 300  format('eigval : '           ,2x,10f9.2)
      end subroutine preigf

      subroutine preigf0(io,e,norbs)
      implicit integer (i-n)
      implicit real*8 (a-h,o-z)
      dimension e(*)
      n=6
      ntimes=norbs/n
      nrest=mod(norbs,n)
      if(ntimes.eq.0) nrest=norbs
      j=1
      n2=n
      do k=1,ntimes
      write(io,300)     (e(i),i=j,n2)
      j =j +n
      n2=n2+n
      enddo
      if(nrest.gt.0.or.ntimes.eq.0) then
      write(io,300)            (e(i),i=j,j+nrest-1)
      endif

 100  format(' #     : ',2x,6(4x,i4,2x))
 300  format('eig.   : '           ,2x,10f9.2)
      end subroutine preigf0
