C_______________________________________________________________________
C                ---- JACINI.F ----
C_______________________________________________________________________
C Performs all the initializations necessary to work with arbitrary
C orthogonal polynomial family. 
C On input
C 	  M_MAX
C	  M
C 	  AL
C 	  BE
C 	  Kind: 1,2,3,4,5,6,7:Legendre,Chebyshev,Chebyshev2,Hermite
C                             Jacobi,Laguerre,Gegenbauer
C 	  kpts
C 	  endpts
C On output
C	  JACM
C	  T
C  	  W
C 	  G
C	  X
C	  B
C_______________________________________________________________________
      SUBROUTINE JACTINI(M_MAX,M,AL,BE,kind,kpts,endpts,JACM,T,W,G,b)
      implicit none
      INTEGER M_MAX,M, I,J,k
      double precision al, be, endpts(0:1)
      integer kind, kpts
      DOUBLE PRECISION JACM(0:M_MAX,0:M_MAX)
      DOUBLE PRECISION AA(0:512,0:512)
      DOUBLE PRECISION T(0:M_MAX), W(0:M_MAX), G(0:M_MAX), b(0:M_MAX)
C--- Generate nodes and weights
c       write(*,*) 'kind = ', kind,al,be, kpts, endpts(0), endpts(1)
      call gaussq(kind, m+1, al, be, kpts, endpts(0), b(0), t(0), w(0))
C--- Now T contains quadrature nodes and W the weights
C    Evaluate the polynomials at the nodes to get the 
C    matrix Q(j,k) = P_k(t_j)
      if (kind .eq. 6 .or. kind .eq. 4) goto 99
C-- Reverse arrays; 0 entry corresponds to +1 endpoint, M to -1.
      DO J=0,M
        b(j) = T(M-j)
      END DO
      DO J=0,M
        T(j) = b(j)
        b(j) = W(M-j)
      END DO
      DO J=0,M
        W(j) = b(j)
      END DO
 99   continue
c     write(*,*) '    j        T         W     '
c     DO J=0,M
c       write(*,*) j, T(j), W(j)
c     END DO
C--- Calculate polynomial values at each of the nodes
      DO J=0,M
        CALL JAC_Q(kind,T(j),al,be,0,M,JACM(0,j))
      END DO
C--- Compute the normalization constants
      DO K = 0,M
          G(k) = 0.d0
        DO J = 0,M
          G(k) = G(k)+W(j)*JACM(k,j)*JACM(k,j)
        END DO
          G(k) = 1/G(k)
      END DO
C--- Test orthonormality
c     Do i = 0,m
c       Do j = 0,m
c           AA(i,j) = 0.d0
c         Do k= 0,m
c           AA(i,j) = AA(i,j)+JACM(i,k)*G(j)*JACM(j,k)*W(k)
c         end do
c       end do
c     end do
c     call mug(AA,M_MAX,0,8,0,8)
      return
      END
C--------------------------------------------------------
