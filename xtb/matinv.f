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

      SUBROUTINE SMATINV (A,LDM,N,D)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LDM, N
      REAL*8, INTENT(OUT)   :: D
      REAL*4, INTENT(INOUT) :: A(LDM,*)
      INTEGER             :: I, J, K, L(N), M(N)
      REAL*8                :: BIGA, TEMP
      REAL*8, PARAMETER     :: TOL = 1.0D-12
 
      D = 1.0
 
      DO K = 1, N
        L(K) = K
        M(K) = K
        BIGA = A(K,K)
        DO J = K, N
          DO I = K, N
            IF ( ABS(BIGA).LT.ABS(A(J,I)) ) THEN
              BIGA = A(J,I)
              L(K) = I
              M(K) = J
            END IF
          END DO
        END DO
        J = L(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = -A(I,K)
            A(I,K) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        I = M(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = -A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        IF ( ABS(BIGA).LT.TOL ) THEN
          D = 0.0
          RETURN
        END IF
        DO I = 1, N
          IF ( I.NE.K ) A(K,I) = A(K,I)/(-BIGA)
        END DO
        DO I = 1, N
          DO J = 1, N
            IF ( I.NE.K ) THEN
              IF ( J.NE.K ) A(J,I) = A(K,I)*A(J,K) + A(J,I)
            END IF
          END DO
        END DO
        DO J = 1, N
          IF ( J.NE.K ) A(J,K) = A(J,K)/BIGA
        END DO
        D = MAX(-1.0D25,MIN(1.0D25,D))
        D = D*BIGA
        A(K,K) = 1.0/BIGA
      END DO
!
      K = N
      DO
!
        K = K - 1
        IF ( K.LE.0 ) EXIT
        I = L(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = A(K,J)
            A(K,J) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        J = M(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = A(I,K)
            A(I,K) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
      END DO
!
      END

      SUBROUTINE DMATINV (A,LDM,N,D)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: LDM, N
      REAL*8, INTENT(OUT)   :: D
      REAL*8, INTENT(INOUT) :: A(LDM,*)
      INTEGER             :: I, J, K, L(N), M(N)
      REAL*8                :: BIGA, TEMP
      REAL*8, PARAMETER     :: TOL = 1.0D-12
!
      D = 1.0
!
      DO K = 1, N
        L(K) = K
        M(K) = K
        BIGA = A(K,K)
        DO J = K, N
          DO I = K, N
            IF ( ABS(BIGA).LT.ABS(A(J,I)) ) THEN
              BIGA = A(J,I)
              L(K) = I
              M(K) = J
            END IF
          END DO
        END DO
        J = L(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = -A(I,K)
            A(I,K) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        I = M(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = -A(K,J)
            A(K,J) = A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        IF ( ABS(BIGA).LT.TOL ) THEN
          D = 0.0
          RETURN
        END IF
        DO I = 1, N
          IF ( I.NE.K ) A(K,I) = A(K,I)/(-BIGA)
        END DO
        DO I = 1, N
          DO J = 1, N
            IF ( I.NE.K ) THEN
              IF ( J.NE.K ) A(J,I) = A(K,I)*A(J,K) + A(J,I)
            END IF
          END DO
        END DO
        DO J = 1, N
          IF ( J.NE.K ) A(J,K) = A(J,K)/BIGA
        END DO
        D = MAX(-1.0D25,MIN(1.0D25,D))
        D = D*BIGA
        A(K,K) = 1.0/BIGA
      END DO
!
      K = N
      DO
!
        K = K - 1
        IF ( K.LE.0 ) EXIT
        I = L(K)
        IF ( I.GT.K ) THEN
          DO J = 1, N
            TEMP = A(K,J)
            A(K,J) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
        J = M(K)
        IF ( J.GT.K ) THEN
          DO I = 1, N
            TEMP = A(I,K)
            A(I,K) = -A(I,J)
            A(I,J) = TEMP
          END DO
        END IF
      END DO
!
      END
