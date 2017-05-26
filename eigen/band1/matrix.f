C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----------------------<<<      matrix.f      >>>-----------------------------
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       SUBROUTINE MATRIX(MMAX,MM,A,B,F,FLAG,TPLUS,TMINUS)
       INCLUDE 'indat.h' 
       INTEGER  MMAX
       DOUBLE PRECISION A(MMAX+1,MMAX+1), B(MMAX+1,MMAX+1)
       DOUBLE PRECISION TPLUS(2*MMAX+2), TMINUS(2*MMAX+2)
       DOUBLE PRECISION F(MMAX+1,*)
       DOUBLE PRECISION TEMP, ANORM, BNORM
       INTEGER  I, J, MM, FLAG(*)
C   Initialize coefficient arrays
C
C ---- Define Matrices A and B; eigenvalue problem is solved
C      in the form: 
C  ******      A * X = LAMBDA * B * X     ******
C      where A is the Tridiagonal matrix and B the upper
C ---- Call routine defining matrices.
C ---- Definition of actual arrays.
       IF (DEBUG .EQ. 0) THEN
          CALL DENDRITE(MMAX,MM,A,B,F,FLAG,TPLUS,TMINUS)
       END IF
       IF (DEBUG .EQ. 1) THEN
C ---- Use quantum harmonic oscillator as a test case.
          CALL SHO(MMAX,MM,A,B,F)
       END IF
       IF (DEBUG .EQ. 2) THEN
C ---- Matrix for convolution by Chebyshev polynomial.
          CALL CHEBY(MMAX,MM,A,B,F)
       END IF
       IF (TRANS .EQ. 1) THEN 
C Transpose matrices
          DO 10 J = 1, MM
             DO 20 I = J+1, MM+1
C ----    Transpose A
                TEMP   = A(J,I)
                A(J,I) = A(I,J)
                A(I,J) = TEMP
C ----    Transpose B
                TEMP   = B(J,I)
                B(J,I) = B(I,J)
                B(I,J) = TEMP
 20          CONTINUE
 10       CONTINUE
       ENDIF
       IF (REVERS .EQ. 1) THEN
C Reverse roles between A and B: solve 1/l A X = B X
          DO 30 J = 1, MM+1
             DO 40 I = 1, MM+1
                TEMP   = A(I,J)
                A(I,J) = B(I,J)
                B(I,J) = TEMP
 40          CONTINUE
 30       CONTINUE
       ENDIF
C-----------------Print matrices
c      WRITE(*,*) '--------- A ----------'
c       DO 777 I = 1, MM+1
c          WRITE(*,*) 'I = ', I, (A(I,J), J = 1, MM+1)
c777     CONTINUE
c       WRITE(*,*) '--------- B ----------'
c        DO 778 I = 1, MM+1
c           WRITE(*,*) 'I = ', I, (B(I,J), J = 1, MM+1)
c778     CONTINUE
C__________________________________________________
C   Compute Frobenius norms of matrices A and B
C__________________________________________________
       CALL NORM(MMAX,MM,A,ANORM)
       CALL NORM(MMAX,MM,B,BNORM)
       WRITE(*,*) ' ANORM = ', ANORM, ' BNORM = ', BNORM
       END
C------------------------------------------------------------------
