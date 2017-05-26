C--------------------------------------------------------
      SUBROUTINE JACT(M_MAX,M,DATA,TRANSFORM, ISIGN,JACM,W,G,b)
C--------------------------------------------------------
      IMPLICIT NONE
      integer     m_max
      DOUBLE PRECISION JACM(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION W(0:M_MAX), G(0:M_MAX),b(0:M_MAX)
      INTEGER M, ISIGN, I, J
      DOUBLE PRECISION DATA(0:*), TRANSFORM(0:*)
      IF(ISIGN .EQ. 1) THEN
        DO J = 0,M
         b(j) = data(j)*w(j)
        END DO
        DO I = 0, M
          transform(i) = 0.d0
          DO J = 0, M
            TRANSFORM(i) = TRANSFORM(i)+JACM(i,j)*b(j)
          END DO
          transform(i) = transform(i)*g(i)
        END DO   
      ELSE
        DO I = 0, M
          data(i) = 0.d0
          DO J = 0, M
            data(i) = data(i)+JACM(j,i)*TRANSFORM(j)
          END DO
        END DO   
      ENDIF
      RETURN
      END
C--------------------------------------------------------------------
