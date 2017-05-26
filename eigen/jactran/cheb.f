c        program costest
c implicit none
c     integer     m_max, n_max
c       parameter(m_max = 512, n_max=0)
c     DOUBLE PRECISION CHEBT(0:M_MAX,0:M_MAX)
c     DOUBLE PRECISION PI
c     PARAMETER(PI=3.1415926535897932384626D0)
c double precision DATA(0:M_MAX,0:N_MAX+1), TRANSFORM(0:M_MAX)
c double precision DAT2(0:M_MAX)
c integer M, N, I, J
c       N = 0
c       READ(*,*) M
c       do i = 0,M
c          do j = 0,N+1
c          DATA(I,J) = 0.d0
c          end do
c          DATA(I,0) = 1.d0
c          DATA(I,0) = cos(PI*I/M)
c          DATA(I,0) = 2.d0*cos(PI*I/M)**2-1.d0
c          DAT2(I) = DATA(I,0)
c       end do
c       CALL COSFTINI(M_MAX,M,CHEBT)
c CALL COST(M_MAX,M,DATA,TRANSFORM,1,CHEBT)
c       write(*,*) (TRANSFORM(I), I=0,M)
c CALL COST(M_MAX,M,DATA,TRANSFORM,-1,CHEBT)
c       write(*,*) (DATA(I,0)-DAT2(I), I=0,M)
c       write(*,*) (DATA(I,0), I=0,M)
c       END
C------------------------------------------------------
C                ---- COSFOU.FOR ----
C------------------------------------------------------
      SUBROUTINE COSFTINI(M_MAX,M,CHEBT)
      INTEGER M_MAX,M, I,J
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932384626D0)
      DOUBLE PRECISION CHEBT(0:M_MAX,0:*)
      DO I=0,M
         DO J=0,M
            CHEBT(I,J) = cos(PI*I*J/M)
         END DO
      END DO
      END
C--------------------------------------------------------
      SUBROUTINE COST (M_MAX,M,DATA,TRANSFORM, ISIGN, CHEBT)
C--------------------------------------------------------
      IMPLICIT NONE
      integer     m_max
      DOUBLE PRECISION CHEBT(0:M_MAX,0:*)
      INTEGER M, ISIGN, I, J
      DOUBLE PRECISION DATA(0:*), TRANSFORM(0:*)
      DO I = 0, M
         IF(ISIGN .EQ. 1) THEN
            TRANSFORM(I) = .5D0*DATA(0)
         ELSE
            TRANSFORM(I) =      DATA(0)
         END IF
         DO J = 1, M-1
            TRANSFORM(I) = TRANSFORM(I)+ CHEBT(I,J)*DATA(J)
         END DO
         IF(ISIGN .EQ. 1) THEN
        TRANSFORM(I) = (2.d0/M)*(TRANSFORM(I)+.5D0*(-1)**I *DATA(M))
         ELSE
        TRANSFORM(I) =           TRANSFORM(I)+     (-1)**I *DATA(M)
         END IF
      END DO
      IF (ISIGN .EQ. 1) THEN
         TRANSFORM(0) = .5d0*TRANSFORM(0)
         TRANSFORM(M) = .5d0*TRANSFORM(M)
      ENDIF

      END
C--------------------------------------------------------------------
