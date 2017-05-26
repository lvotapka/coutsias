c        program costest
cimplicit none
c     INCLUDE 'param.h'
c     INCLUDE 'cost.h'
c     DOUBLE PRECISION PI
c     PARAMETER(PI=3.1415926535897932384626D0)
cdouble precision DATA(M_MAX+2,N_MAX+2), TRANSFORM(M_MAX+2)
cdouble precision DAT2(M_MAX+2,N_MAX+2)
cinteger M, N, I, J, II
c       N = 1
c       READ(*,*) M
c       do i = 0,M
c          ii = i+1
c          do j = 1,N+1
c             DATA(II,J) = 0.d0
c             DAT2(II,J) = 0.d0
c          end do
c          DATA(II,1) = 1.d0
c          DATA(II,1) = cos(PI*I/M)
c          DATA(II,1) = 2.d0*cos(PI*I/M)**2-1.d0
c          DAT2(II,1) = DATA(II,1)
c       end do
c       CALL COSFTINI(M)
cCALL COST(DATA,M,TRANSFORM,1)
cCALL COST(DATA,M,TRANSFORM,-1)
c       write(*,*) (DATA(I,1)-DAT2(I,1), I=1,M+1)
c       write(*,*) (DATA(I,1), I=1,M+1)
c       END
C------------------------------------------------------
C                ---- COSFOU.FOR ----
C------------------------------------------------------
      SUBROUTINE COSFTINI(M)
      INTEGER M, I,J
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932384626D0)
      INCLUDE 'limits.h'
      INCLUDE 'cost.h'
      DO I=0,M
         DO J=0,M
            CHEBT(I,J) = cos(PI*I*J/M)
         END DO
      END DO
      END
C--------------------------------------------------------
      SUBROUTINE COST (DATA, M, TRANSFORM, ISIGN)
C--------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'limits.h'
      INCLUDE 'cost.h'
      INTEGER M, ISIGN, I, J, II, JJ
      DOUBLE PRECISION DATA(M_MAX+2), TRANSFORM(M_MAX+2)
      DO I = 0, M
         II = I+1
         IF(ISIGN .EQ. 1) THEN
            TRANSFORM(II) = .5D0*DATA(1)
         ELSE
            TRANSFORM(II) =      DATA(1)
         END IF
         DO J = 1, M-1
            JJ = J+1
            TRANSFORM(II) = TRANSFORM(II)+ CHEBT(I,J)*DATA(JJ)
         END DO
         IF(ISIGN .EQ. 1) THEN
        TRANSFORM(II) = (2.d0/M)*(TRANSFORM(II)+.5D0*(-1)**I *DATA(M+1))
         ELSE
        TRANSFORM(II) =           TRANSFORM(II)+     (-1)**I *DATA(M+1)
         END IF
      END DO
      IF (ISIGN .EQ. 1) THEN
         TRANSFORM(  1) = .5d0*TRANSFORM(  1)
         TRANSFORM(M+1) = .5d0*TRANSFORM(M+1)
      ENDIF

      END
C--------------------------------------------------------------------
