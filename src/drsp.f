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


c     real(wp) a(100),e(10),z(10,10),aux1(20),aux2(20)
c     a(1)=0.
c     a(2)=1.
c     a(3)=0.
c     a(4)=0.
c     a(5)=1.
c     a(6)=0.
c     n=3
c     m=3
c     call rsp (10,a,n,n,e,z,aux1,aux2)
c     do i=1,n
c        write(*,'(f10.5,5x,3F12.6)') e(i), (z(j,i)      ,j=1,n)
c     enddo
c     end

      SUBROUTINE RSP(A,N,MATZ,W,Z)
         use xtb_mctc_accuracy, only : wp
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(n*(n+1)/2),  W(n), Z(n,n)  
      DIMENSION :: FV1(2*n),FV2(2*n)
*******************************************************************
*
*   EISPACK DIAGONALIZATION ROUTINES: TO FIND THE EIGENVALUES AND
*           EIGENVECTORS (IF DESIRED) OF A REAL SYMMETRIC PACKED MATRIX.
* ON INPUT-      N  IS THE ORDER OF THE MATRIX  A,
*                A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
*                   PACKED MATRIX STORED ROW-WISE,
*             MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF ONLY
*                   EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO
*                   ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND
*                   EIGENVECTORS.
* ON OUTPUT-     W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,
*                Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,
*
*******************************************************************
* THIS SUBROUTINE WAS CHOSEN AS BEING THE MOST RELIABLE. (JJPS)
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      real(wp) esp,eta

      if(n.eq.1)then
         Z(1,1)=1.0d0 
         W(1)=A(1)
         return
      endif

c     CALL EPSETA(EPS,ETA)
c     write(*,*) eps,eta
      eps=1.000000000000000D-016 
      eta=1.000000000000000D-039
      
      NV=(N*(N+1))/2
      NM=N 
      CALL  TRED3(N,NV,A,W,FV1,FV2,EPS,EPS)
      IF (MATZ .NE. 0) GO TO 10
C     ********** FIND EIGENVALUES ONLY **********
      CALL  TQLRAT(N,W,FV2,IERR,EPS)
      GO TO 40
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********
   10 DO 30    I = 1, N
C
         DO 20    J = 1, N
            Z(J,I)=0.0D0
   20    CONTINUE
C
         Z(I,I)=1.0D0
   30 CONTINUE
C
      CALL  TQL2(NM,N,W,FV1,Z,IERR,EPS)
      IF (IERR .NE. 0) GO TO 40
      CALL  TRBAK3(NM,N,NV,A,N,Z)
C     ********** LAST CARD OF RSP **********
   40 RETURN
      END

      SUBROUTINE EPSETA(EPS,ETA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     COMPUTE AND RETURN ETA, THE SMALLEST REPRESENTABLE NUMBER,
C     AND EPS IS THE SMALLEST NUMBER FOR WHICH 1+EPS.NE.1.
C
C
      ETA = 1.D0
   10 IF((ETA/2.D0).EQ.0.D0) GOTO 20
      IF(ETA.LT.1.D-38) GOTO 20
      ETA = ETA / 2.D0
      GOTO 10
   20 EPS = 1.D0
   30 IF((1.D0+(EPS/2.D0)).EQ.1.D0) GOTO 40
      IF(EPS.LT.1.D-17) GOTO 40
      EPS = EPS / 2.D0
      GOTO 30
   40 RETURN
      END

      SUBROUTINE TQL2(NM,N,D,E,Z,IERR,EPS)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C               ===== PROCESSED BY AUGMENT, VERSION 4N =====
C     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
C               ----- LOCAL VARIABLES -----
C               ----- GLOBAL VARIABLES -----
      DIMENSION D(*), E(*), Z(NM,*)
C               ----- SUPPORTING PACKAGE FUNCTIONS -----
C               ===== TRANSLATED PROGRAM =====
C
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND
C     WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS
C     FULL MATRIX TO TRIDIAGONAL FORM.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN
C          THE IDENTITY MATRIX.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT
C          UNORDERED FOR INDICES 1,2,...,IERR-1,
C
C        E HAS BEEN DESTROYED,
C
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED
C          EIGENVALUES,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C
      IERR = 0
      IF (N .EQ. 1) GO TO 160
C
      DO 10   I = 2, N
   10 E(I-1) = E(I)
C
      F=0.0D0
      B=0.0D0
      E(N)=0.0D0
C
      DO 110   L = 1, N
         J = 0
         H=EPS*(ABS (D(L))+ABS (E(L)))
         IF (B .LT. H) B=H
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********
         DO 20   M = L, N
            IF (ABS (E(M)).LE.B)  GO TO 30
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
   20    CONTINUE
C
   30    IF (M .EQ. L) GO TO 100
   40    IF (J .EQ. 30) GO TO 150
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         G = D(L)
         P=(D(L1)-G)/(2.0D0*E(L))
         R=SQRT (P*P+1.0D0)
         D(L)=E(L)/(P+SIGN (R,P))
         H = G - D(L)
C
         DO 50   I = L1, N
   50    D(I) = D(I) - H
C
         F = F + H
C     ********** QL TRANSFORMATION **********
         P = D(M)
         C=1.0D0
         S=0.0D0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 90   II = 1, MML
            I = M - II
            G = C * E(I)
            H = C * P
            IF (ABS (P).LT.ABS (E(I)))  GO TO 60
            C = E(I) / P
            R=SQRT (C*C+1.0D0)
            E(I+1) = S * P * R
            S = C / R
            C=1.0D0/R
            GO TO 70
   60       C = P / E(I)
            R=SQRT (C*C+1.0D0)
            E(I+1) = S * E(I) * R
            S=1.0D0/R
            C = C * S
   70       P = C * D(I) - S * G
            D(I+1) = H + S * (C * G + S * D(I))
C     ********** FORM VECTOR **********
            DO 80   K = 1, N
               H = Z(K,I+1)
               Z(K,I+1) = S * Z(K,I) + C * H
               Z(K,I) = C * Z(K,I) - S * H
   80       CONTINUE
C
   90    CONTINUE
C
         E(L) = S * P
         D(L) = C * P
         IF (ABS (E(L)).GT.B)  GO TO 40
  100    D(L) = D(L) + F
  110 CONTINUE
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********
      DO 140   II = 2, N
         I = II - 1
         K = I
         P = D(I)
C
         DO 120   J = II, N
            IF (D(J) .GE. P) GO TO 120
            K = J
            P = D(J)
  120    CONTINUE
C
         IF (K .EQ. I) GO TO 140
         D(K) = D(I)
         D(I) = P
C
         DO 130   J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
  130    CONTINUE
C
  140 CONTINUE
C
      GO TO 160
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  150 IERR = L
  160 RETURN
C     ********** LAST CARD OF TQL2 **********
      END
      SUBROUTINE TQLRAT(N,D,E2,IERR,EPS)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C               ===== PROCESSED BY AUGMENT, VERSION 4N =====
C     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
C               ----- LOCAL VARIABLES -----
C               ----- GLOBAL VARIABLES -----
      DIMENSION D(*), E2(*)
C               ----- SUPPORTING PACKAGE FUNCTIONS -----
C               ===== TRANSLATED PROGRAM =====
C
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,
C
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY.
C
C      ON OUTPUT-
C
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE
C          THE SMALLEST EIGENVALUES,
C
C        E2 HAS BEEN DESTROYED,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN
C                     DETERMINED AFTER 30 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C
      IERR = 0
      IF (N .EQ. 1) GO TO 140
C
      DO 10   I = 2, N
   10 E2(I-1) = E2(I)
C
      F=0.0D0
      B=0.0D0
      E2(N)=0.0D0
C
      DO 120   L = 1, N
         J = 0
         H=EPS*(ABS (D(L))+SQRT (E2(L)))
         IF (B .GT. H) GO TO 20
         B = H
         C = B * B
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
   20    DO 30   M = L, N
            IF (E2(M) .LE. C) GO TO 40
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
C                THROUGH THE BOTTOM OF THE LOOP **********
   30    CONTINUE
C
   40    IF (M .EQ. L) GO TO 80
   50    IF (J .EQ. 30) GO TO 130
         J = J + 1
C     ********** FORM SHIFT **********
         L1 = L + 1
         S=SQRT (E2(L))
         G = D(L)
         P=(D(L1)-G)/(2.0D0*S)
         R=SQRT (P*P+1.0D0)
         D(L)=S/(P+SIGN (R,P))
         H = G - D(L)
C
         DO 60   I = L1, N
   60    D(I) = D(I) - H
C
         F = F + H
C     ********** RATIONAL QL TRANSFORMATION **********
         G = D(M)
         IF (G.EQ.0.0D0) G=B
         H = G
         S=0.0D0
         MML = M - L
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
         DO 70   II = 1, MML
            I = M - II
            P = G * H
            R = P + E2(I)
            E2(I+1) = S * R
            S = E2(I) / R
            D(I+1) = H + S * (H + D(I))
            G = D(I) - E2(I) / G
            IF (G.EQ.0.0D0) G=B
            H = G * P / R
   70    CONTINUE
C
         E2(L) = S * G
         D(L) = H
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST **********
         IF (H.EQ.0.0D0)  GO TO 80
         IF (ABS (E2(L)).LE.ABS (C/H))  GO TO 80
         E2(L) = H * E2(L)
         IF (E2(L).NE.0.0D0)  GO TO 50
   80    P = D(L) + F
C     ********** ORDER EIGENVALUES **********
         IF (L .EQ. 1) GO TO 100
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
         DO 90   II = 2, L
            I = L + 2 - II
            IF (P .GE. D(I-1)) GO TO 110
            D(I) = D(I-1)
   90    CONTINUE
C
  100    I = 1
  110    D(I) = P
  120 CONTINUE
C
      GO TO 140
C     ********** SET ERROR -- NO CONVERGENCE TO AN
C                EIGENVALUE AFTER 30 ITERATIONS **********
  130 IERR = L
  140 RETURN
C     ********** LAST CARD OF TQLRAT **********
      END
      SUBROUTINE TRBAK3(NM,N,NV,A,M,Z)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C               ===== PROCESSED BY AUGMENT, VERSION 4N =====
C     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
C               ----- LOCAL VARIABLES -----
C               ----- GLOBAL VARIABLES -----
      DIMENSION A(*), Z(NM,*)
C               ===== TRANSLATED PROGRAM =====
C
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRBAK3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL SYMMETRIC
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     SYMMETRIC TRIDIAGONAL MATRIX DETERMINED BY  TRED3.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL TRANSFORMATIONS
C          USED IN THE REDUCTION BY  TRED3  IN ITS FIRST
C          N*(N+1)/2 POSITIONS,
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED,
C
C        Z CONTAINS THE EIGENVECTORS TO BE BACK TRANSFORMED
C          IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT-
C
C        Z CONTAINS THE TRANSFORMED EIGENVECTORS
C          IN ITS FIRST M COLUMNS.
C
C     NOTE THAT TRBAK3 PRESERVES VECTOR EUCLIDEAN NORMS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 50
      IF (N .EQ. 1) GO TO 50
C
      DO 40   I = 2, N
         L = I - 1
         IZ = (I * L) / 2
         IK = IZ + I
         H = A(IK)
         IF (H.EQ.0.0D0)  GO TO 40
C
         DO 30   J = 1, M
            S=0.0D0
            IK = IZ
C
            DO 10   K = 1, L
               IK = IK + 1
               S = S + A(IK) * Z(K,J)
   10       CONTINUE
C     ********** DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW **********
            S = (S / H) / H
            IK = IZ
C
            DO 20   K = 1, L
               IK = IK + 1
               Z(K,J) = Z(K,J) - S * A(IK)
   20       CONTINUE
C
   30    CONTINUE
C
   40 CONTINUE
C
   50 RETURN
C     ********** LAST CARD OF TRBAK3 **********
      END
      SUBROUTINE TRED3(N,NV,A,D,E,E2,EPS,ETA)
      IMPLICIT INTEGER (I-N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C               ===== PROCESSED BY AUGMENT, VERSION 4N =====
C     APPROVED FOR VAX 11/780 ON MAY 6,1980.  J.D.NEECE
C               ----- LOCAL VARIABLES -----
C               ----- GLOBAL VARIABLES -----
      DIMENSION A(*), D(*), E(*), E2(*)
C               ----- SUPPORTING PACKAGE FUNCTIONS -----
C               ===== TRANSLATED PROGRAM =====
C
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        N IS THE ORDER OF THE MATRIX,
C
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,
C
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.
C
C     ON OUTPUT-
C
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL
C          TRANSFORMATIONS USED IN THE REDUCTION,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 100   II = 1, N
         I = N + 1 - II
         L = I - 1
         IZ = ( I * L ) / 2
         H=0.0D0
         SCALE=0.0D0
         DO 10   K = 1, L
            IZ = IZ + 1
            D(K) = A(IZ)
            SCALE=SCALE+ABS( D(K) )
   10    CONTINUE
C
         IF ( SCALE.NE.0.D0 ) GO TO 20
         E(I)=0.0D0
         E2(I)=0.0D0
         GO TO 90
C
   20    DO 30   K = 1, L
            D(K) = D(K) / SCALE
            H = H + D(K) * D(K)
   30    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         F = D(L)
         G=-SIGN (SQRT (H),F)
         E(I) = SCALE * G
         H = H - F * G
         D(L) = F - G
         A(IZ) = SCALE * D(L)
         IF (L .EQ. 1) GO TO 90
         F=0.0D0
C
         DO 70   J = 1, L
            G=0.0D0
            JK = (J * (J-1)) / 2
C     ********** FORM ELEMENT OF A*U **********
            K = 0
   40       K = K + 1
            JK = JK + 1
            G = G + A(JK) * D(K)
            IF ( K .LT. J ) GO TO 40
            IF ( K .EQ. L ) GO TO 60
   50       JK = JK + K
            K = K + 1
            G = G + A(JK) * D(K)
            IF ( K .LT. L ) GO TO 50
C     ********** FORM ELEMENT OF P **********
   60       CONTINUE
            E(J) = G / H
            F = F + E(J) * D(J)
   70    CONTINUE
C
         HH = F / (H + H)
         JK = 0
C     ********** FORM REDUCED A **********
         DO 80   J = 1, L
            F = D(J)
            G = E(J) - HH * F
            E(J) = G
C
            DO 80   K = 1, J
               JK = JK + 1
               A(JK) = A(JK) - F * E(K) - G * D(K)
   80    CONTINUE
C
   90    D(I) = A(IZ+1)
         A(IZ+1)=SCALE*SQRT (H)
  100 CONTINUE
C
      RETURN
C     ********** LAST CARD OF TRED3 **********
      END
