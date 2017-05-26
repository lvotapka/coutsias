C------------------------------------------------------------------
      SUBROUTINE D1MAT(ML, A1)
       IMPLICIT NONE
       INTEGER ML, I, J
       DOUBLE PRECISION  A1(0:ML, 0:ML)
C---- Fill array with zeroes
        DO I = 0, ML
           DO J = 0, ML
              A1(I,J) = 0.D0
           END DO
        END DO
C---- Insert non-zero elements
        DO I = 0, ML-1
           DO J = I+1, ML, 2
              A1(I,J) = 2.d0*dble(J)
           END DO
        END DO
C---- Special care for first row:
        DO J = 0, ML
           A1( 0,J) = A1(0,J)/2.D0
           A1(ML,J) = 0.d0
        END DO
       RETURN
       END
C------------------------------------------------------------------
      SUBROUTINE R1D1MAT(ML, A1, aspect)
       IMPLICIT NONE
       INTEGER ML, I, J
       DOUBLE PRECISION  A1(0:ML, 0:ML)
       DOUBLE PRECISION  aspect
C---- Fill array with zeroes
        DO I = 0, ML
           DO J = 0, ML
              A1(I,J) = 0.D0
           END DO
        END DO  
C---- Insert non-zero elements
        DO I = 0, ML-1
           DO J = I, ML
              A1(I,J) = 2.d0*dble(J)
           END DO
           DO J = I+1, ML, 2
              A1(I,J) = A1(I,J)*aspect
           END DO
        END DO  
C---- Special care for first row:
        DO J = 0, ML
           A1( 0,J) = A1(0,J)/2.D0
           A1(ML,J) = 0.d0
           A1( J,J) = A1(J,J)/2.D0
        END DO
       RETURN 
       END
C------------------------------------------------------------------
      SUBROUTINE R2LAPNMAT(ML,MODE, R2LAP,ASPECT)
C------------------------------------------------------------------
       IMPLICIT NONE
       INTEGER ML, I, J, MODE
       DOUBLE PRECISION  R2LAP(0:ML, 0:ML), ASPECT, A1, A2
       A1 = 2.d0*ASPECT
       A2 = ASPECT*ASPECT+1.d0
       DO I = 0,ML
         DO J = 0,ML
           R2LAP(I,J) = 0.d0
         END DO
       END DO
       DO I = 0,ML-1,2
         DO J = I+1,ML,2
           R2LAP(I  ,J  ) = A1*DBLE( J   *( J   * J   - I   * I   ))
           R2LAP(I  ,J+1) = A2*DBLE((J+1)*((J+1)*(J+1)- I   * I   ))
           R2LAP(I+1,J  ) = A2*DBLE( J   *( J   * J   -(I+1)*(I+1)))
           R2LAP(I+1,J+1) = A1*DBLE((J+1)*((J+1)*(J+1)-(I+1)*(I+1)))
         END DO
       END DO
       DO I = 0, ML
          R2LAP(0,I) = R2LAP(0,I)/2.d0
       END DO
       DO I = 0, ML
          R2LAP(I,I) = DBLE(I)**2 - DBLE(MODE)**2
       END DO
       RETURN
      END
C------------------------------------------------------------------
      SUBROUTINE R2LAPNMATMUL(ML,MODE, RDR,R2LAP)
       IMPLICIT NONE
       INTEGER ML, I, J, MODE
       DOUBLE PRECISION  RDR(0:ML, 0:ML)
       DOUBLE PRECISION  R2LAP(0:ML, 0:ML)
       DO I = 0,ML
         DO J = 0,ML
           R2LAP(I,J) = 0.d0
         END DO
       END DO
       CALL FFMUL(ML+1, ML+1, ML+1, ML+1, RDR, RDR, R2LAP)
       DO I = 0, ML
          R2LAP(I,I) = R2LAP(I,I) - DBLE(MODE)**2
       END DO
       RETURN
      END
C------------------------------------------------------------------
C-----------------------------------------------------------------------
C               -----------<    deriv.f   >------------
C       compute the derivatives of a Fourier (theta)-Chebyshev (x)
C       expansion
C-----------------------------------------------------------------------
	subroutine d1th(mmax,m,n,A,B,period)
C-----------------------------------------------------------------------
	implicit none
	integer mmax,m,n,period,i,j
	double precision A(0:mmax,0:1,0:*)
	double precision B(0:mmax,0:1,0:*)
        double precision tmp
C-----------------------------------------------------------------------
        do i = 0,m
          do j =0,n
            tmp      = dble(j)*A(i,1,j)
            B(i,1,j) =-dble(j)*A(i,0,j)
            B(i,0,j) = tmp
          end do
        end do
	return
	end
C-----------------------------------------------------------------------
        subroutine d1r(mmax,m,n,A,B)
C-----------------------------------------------------------------------
	implicit none
	integer mmax,m,n,i,j
	double precision A(0:mmax,0:*)
	double precision B(0:mmax,0:*)
C-----------------------------------------------------------------------
        do j = 0,n
            B(m  ,j) = 0.d0
            B(m-1,j) = dble(2*m)*A(m,j)
          do i = m-2,1,-1
            B(i,j) = B(i+2,j)+dble(2*(i+1))*A(i+1,j)
          end do
            B(0,j) = B(2,j)/2.d0+A(1,j)
        end do
	return
	end
C-----------------------------------------------------------------------
        subroutine r1d1r(mmax,m,n,A,B,aspect)
C-----------------------------------------------------------------------
        implicit none
        integer mmax,m,n,i,j
        double precision A(0:mmax,0:*)
        double precision B(0:mmax,0:*)
        double precision aspect
C-----------------------------------------------------------------------
        do j = 0,n
            B(m  ,j) = dble(m)*A(m,j)
            B(m-1,j) = dble(m-1)*A(m-1,j)+aspect*dble(2*m)*A(m,j)
          do i = m-2,2,-2
            B(i  ,j) = dble(i  )*A(i  ,j)
     _                +dble(i+1)*A(i+1,j)*2.d0*aspect
     _                +dble(i+2)*A(i+2,j)
     _                +          B(i+2,j)
            B(i-1,j) = dble(i-1)*A(i-1,j)
     _                +dble(i  )*A(i  ,j)*2.d0*aspect
     _                +dble(i+1)*A(i+1,j)
     _                +          B(i+1,j)
          end do
            B(0,j) = (2.d0*A(2,j)+B(2,j))/2.d0+aspect*A(1,j)
        end do
        return
        end
C-----------------------------------------------------------------------
      SUBROUTINE CHEBINT(mmax, m, n, A,V)
C--- Accepts an array whose columns are the Chebyshev coeffs. of a set
C    of functions. Produces the integrals of these functions over [-1,1].
       IMPLICIT NONE
       INTEGER mmax, m, n, I, J
       DOUBLE PRECISION  A(0:mmax, 0:*), V(*), II
       DO J = 0, n
             V(J) = 2.D0*A(0,J)
          DO I = 2, m, 2
             II = DBLE(I)
             V(J) = V(J)-2.D0*A(I,J)/(II**2-1.D0)
          END DO
       END DO
       RETURN
      END
C------------------------------------------------------------------

C-----------------------------------------------------------------------
	subroutine convol(mmax,m,vec,mat)
C-----------------------------------------------------------------------
        implicit none
        integer mmax, m, i, k
        double precision vec(0:2*mmax), mat(0:mmax,0:*)
        do i = 0,mmax
          do k = 0,mmax
            mat(i,k) = 0.d0
          end do
        end do
        do k = 0,m
          do i = k,m
            mat(i,i-k) = mat(i,i-k)+.5d0*vec(k)
            mat(i-k,i) = mat(i-k,i)+.5d0*vec(k)
          end do
          do i = 1, k
            mat(i,k-i) = mat(i,k-i)+.5d0*vec(k)
          end do
        end do
        do k = m+1, 2*m
          do i = k-m, m
            mat(i,k-i) = mat(i,k-i)+.5d0*vec(k) 
          end do
        end do
        return
        end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C                 -------<  c o n d 1 . f  >-------
C-----------------------------------------------------------------------
		subroutine cond1_mat(maxm,pp,kk,OP,aspect,debug,buff)
C-----------------------------------------------------------------------
C     initialize various banded preconditioners
C     PHASE 1: Chebyshev family
C     
C     PHASE 2: Jacobi of arbitrary type.
C       will characterize by giving recurrence functions for X and B.
C-----------------------------------------------------------------------
C-- Build type 1 (integration) preconditioner matrix B_[2]^pp (x+a)^kk
C   with bands -(pp+kk):(pp+kk) 
C   All bands are introduced, even those that may contain zeroes.
C   The computation is carried over enough additional modes so that 
C   truncation at maxm is exact.
C-----------------------------------------------------------------------
      implicit none
      integer maxm,pp,kk,debug,i,j,k,mmax,buff
      double precision aspect
      double precision OP(0:maxm+buff,-pp-kk:pp+kk)
      double precision a1(-16:16)
      double precision a2(-16:16)
C-----------------------------------------------------------------------
C-- Zero-out arrays
      do i = 0,maxm+buff
        do j = -pp-kk, pp+kk
          OP(i,j) = 0.d0
        end do
      end do
      do i = -16,16
        a1(i) = 0.d0
        a2(i) = 0.d0
      end do
C-----------------------------------------------------------------------
C- the default for op is I_[2], the identity with first two rows set to zero
C  but this is left to the programmer; this routine produces arrays correctly
C  with no special entries.
      do i = 0,maxm+buff
        OP(i,0) = 1.d0
      end do
C-------------------
C - extend the computation adequately to ensure that truncation at maxm
C   produces exact bands. Thus, treat last rows same as the rest;
C   must extend if higher powers/order are to be used.
c     mmax = maxm+buff
      mmax = maxm
C-------------------
C-- First build (x+a)^kk; bypass if kk = 0
      if (kk .gt. 0) then
        do k = 1, kk
c- zero help arrays
          do j = -k,k
            a1(j) = 0.d0
            a2(j) = 0.d0
          end do
c-   0-row
          do j = 0, k
            a2(j) = .5d0*op(1,j-1)+aspect*op(0,j)
          end do
C-  1-row
          do j = -1, k-1
            a1(j) = .5d0*op(2,j-1)+op(0,j+1)+aspect*op(1,j)
          end do
            a1(k  ) = .5d0*op(2,k-1)
          do j = 0,k
            op(0,j) = a2(j)
          end do
          do j = -1,k
            a2(j) = a1(j)
          end do
C-  all other rows
C- first compute row entries and store; old contents still needed for next row
          do i = 2, mmax+buff-k
            a1(-k  ) = .5d0*op(i-1,-k+1)
            do j = -k+1, k-1
              a1(j) = .5d0*(op(i-1,j+1)+op(i+1,j-1))+aspect*op(i,j)
            end do        
            a1(k  ) = .5d0*op(i+1,k-1)
C- Now, previous row no longer needed; replace with stored contents
C  from a2; then put a1 into a2, freeing it for next cycle
            do j = -k, k
              op(i-1,j) = a2(j)
              a2(j)     = a1(j)
            end do
          end do        
c-  last row
            do j = -k, k
              op(mmax+buff-k  ,j) = a2(j)
            end do
        end do        
      end if
C-----------------------------------------------------------------------
C-- Now multiply pp times by B
      do k = 1, pp
c- zero help arrays
        do j = -kk-k,kk+k
          a1(j) = 0.d0
          a2(j) = 0.d0
        end do
        i = 1
            a2(-kk-k) =  op(i-1,-kk-k+1)
          do j = -kk-k+1,kk+k-1
            a2( j  ) = (op(i-1, j+1)-op(i+1,j-1)/2.d0)/dble(i)
          end do
            a2( kk+k) =           -op(i+1,kk+k-1)/(2.d0*dble(i))
        do i = 2, mmax+buff-kk-k
            a1(-kk-k) = op(i-1,-kk-k+1)          /(2.d0*dble(i))
          do j = -kk-k+1,kk+k-1
            a1( j  ) = (op(i-1, j+1)-op(i+1,j-1))/(2.d0*dble(i))
          end do
            a1( kk+k) =             -op(i+1,kk+k-1)/(2.d0*dble(i))
          do j = -kk-k,kk+k
            op(i-1,j) = a2(j) 
            a2(j)     = a1(j)
          end do
        end do
          do j = -kk-k,kk+k
            op(mmax+buff-kk-k,j) = a2(j) 
          end do
      end do
C- Zero out first pp rows, regardless of content
C-- Characteristic of B_[2]^2 preconditioning. 
C-- here it is done to avoid any special entries in the B
      do j = -kk-pp,kk+pp
        do i = 0,pp-1
          op(i,j) = 0.d0
        end do
      end do

c       do j = 1,4*maxm+4
c         do i = 1,4*maxm+4
c           write(*,*) i,j, op(i,j)
c         end do
c       end do
c       STOP

      return
      end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
          subroutine rdiv(mmax,m,n,U,U0,a)
C Experimental division routine. Intended for the Fourier-Chebyshev
C Divides U by x+a, puts result in U0
C-----------------------------------------------------------------------
          implicit none
          integer mmax,m,n,i,j
          double precision a,ap,am,bp,bm,br
          double precision U(0:mmax,0:*)
          double precision U0(0:mmax,0:*)
          double precision z
C-----------------------------------------------------------------------
C if a is too small to resolve to double precision, give a-1 in one of the
C parameters, and consider it here
          if (a .lt. 1.d0) then
            am = a
            ap = 2.d0
            a  = 1.d0
          else
            ap = a+1.d0
            am = a-1.d0
          end if
          z  = dsqrt(am*ap)
          bp = a+z
          bm = a-z
          br = (bp**2+1.d0)/(bp**2-1.d0)
          write(*,*)  '   bp = ', bp
          write(*,*)  '   bm = ', bm
          write(*,*)  '   br = ', br
          do j = 0,n
              U0(m,j) = U(m,j)
            do i = m-1,0,-1
              U0(i,j) = U(i,j)-U(i+1,j)*bm
            end do
              U0(0,j) = U0(0,j)*bm
            do i = 1,m
              U0(i,j) = (U0(i,j)-U0(i-1,j))*bm
            end do
            z      = -U0(0,j)*br*bm
            U0(0,j) = 2.d0*bp**2*U0(0,j)/(bp**2-1.d0)
            do i = 1, m
              U0(i,j) = (U0(i,j) + z)*2.d0
              z = -z*bm
            end do
          end do
          return
          end
C-----------------------------------------------------------------------
c koper.f
c========================================================================
              subroutine koper(ORD,MAXM,MC,DOP1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine that builds a square Chebyshev operator matrix
c  ASSUMING MC IS AN EVEN NUMBER
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer maxm,mc
      integer i, j, ORD
      double precision dop1(MAXM+1,mc)
      double precision c1,c2,i1,i2,j2,j4
c------------------------------------------------------------------------
c....Initializations
c------------------------------------------------------------------------
      do 20 i=1,maxm+1
        do 15 j=1,mc
          dop1(i,j  ) = 0.0d0
15      continue
20    continue
c....build operator matrices
c------------------------------------------------------------------------
c....explicitly finding fourth derivative (assuming mc even) 
c------------------------------------------------------------------------
      if (ORD .eq. 4) then
            c1 = 1.d0/48.d0
            c2 = 1.d0/24.d0
        do 25 j=4,mc-2,2  
            j2 = j**2  
            dop1(1,j+1) = c1*j*(j2*(j2 - 4.d0)**2)
25      continue
        do 35 i=2,mc-4
          do 30 j=i+3,mc-1,2
            i2 = ((i*1.d0) - 1.d0)**2
            j2 = j**2
            dop1(i,j+1) = c2*j*(j2*(j2 - 4.d0)**2 - 3.d0*i2*(j2**2)
     _                           + 3.d0*(i2**2)*j2 - i2*(i2 - 4.d0)**2)
30        continue        
35      continue          
      end if
c------------------------------------------------------------------------
c....explicitly finding third derivative (assuming mc even)
c------------------------------------------------------------------------
      if (ORD .eq. 3) then
        do 40 j=3,mc-1,2
            j2 = j**2
            j4 = j2**2
            dop1(1,j+1) = .125d0*j*(j4 - 2.0d0*j2 + 1.d0)
40      continue
        do 50 i=2,mc-3
            i2 = ((i*1.d0) - 1.d0)**2
            i1 = i2 + 1.d0
            i2 = (i2 - 1.d0)**2
          do 45 j=i+2,mc-1,2
            j2 = j**2
            j4 = j2**2
            dop1(i,j+1) = .25d0*j*(j4 - 2.d0*j2*i1 + i2)
45        continue
50      continue
      end if
c------------------------------------------------------------------------
c....explicitly finding second derivative
c------------------------------------------------------------------------
      if (ORD .eq. 2) then
        do 65 j=2, mc-2,2
          dop1(1,j+1) =  .5d0*(j**3)
65      continue
        do 75 i=2, mc-2
          i1 = (1.d0*(i-1))**2
          do 70 j=i+1,mc-1,2
            j2 = (1.d0*j)**2
            dop1(i  ,j+1) = j*(j2-i1)
70        continue
75      continue
      end if
c------------------------------------------------------------------------
c....explicitly finding first derivative
c------------------------------------------------------------------------
      if (ORD .eq. 1) then
        do 90 j=1, mc-1,2
          dop1(1,j+1) = j*1.0d0
90      continue
        do 100 i=2, mc-1
          do 95 j=i,mc-1,2
            dop1(i  ,j+1) = j*2.0d0
95        continue
100     continue
      end if

      end
c=======================================================================
c mdivmul.f
c==============================================================

c--------------------------------------------------------------------
      subroutine rsqmul(rho,rsqrho,maxm,n,m,a)
c--------------------------------------------------------------------
      integer          n, m, i, j,maxm
      double precision  rho(maxm+1,n+1), rsqrho(maxm+1,n+1)
      double precision a, a1, a2

      a1= 4.0d0*a
      a2= 4.0d0*a**2+2.0D0

      do 10 j=1, n+1
        rsqrho(1,j) = (a2*rho(1,j)+a1*rho(2,j)+rho(3,j))/4.0d0
        rsqrho(2,j) = (2.0d0*a1*rho(1,j)+(a2+1.0d0)*rho(2,j)
     -                 +a1*rho(3,j)+rho(4,j))/4.0d0
        rsqrho(3,j) = (2.0d0*rho(1,j)+a1*rho(2,j)+a2*rho(3,j)
     -                 +a1*rho(4,j)+rho(5,j))/4.0d0
10    continue

      do 30 i=4, m-1
        do 20 j=1, n+1
          rsqrho(i,j)=(rho(i-2,j)+a1*rho(i-1,j)+a2*rho(i,j)
     -                 +a1*rho(i+1,j)+rho(i+2,j))/4.0d0
20      continue
30    continue

      do 40 j=1, n+1
        rsqrho(m,j)=(rho(m-2,j)+a1*rho(m-1,j)+a2*rho(m,j)
     -                 +a1*rho(m+1,j))/4.0d0
        rsqrho(m+1,j)= (rho(m-1,j)+a1*rho(m,j)+a2*rho(m+1,j))/4.0d0
40    continue

      end    
c====================================================================

c--------------------------------------------------------------------
      subroutine rmul(rho,rrho,maxm,n,m,a)
c--------------------------------------------------------------------
      integer          n, m, i, j,maxm
      double precision  rho(maxm+1,n+1), rrho(maxm+1,n+1)
      double precision a

      do 10 j=1, n+1
        rrho(1,j) = a*rho(1,j)+0.5d0*rho(2,j)
        rrho(2,j) = (rho(1,j)+a*rho(2,j))+0.5d0*rho(3,j)
10    continue

      do 30 i=3, m
        do 20 j=1, n+1
          rrho(i,j)=((rho(i-1,j)+rho(i+1,j))/2.0d0)+a*rho(i,j)
20      continue
30    continue

      do 40 j=1, n+1
        rrho(m+1,j)=a*rho(m+1,j)+0.5d0*rho(m,j)
40    continue

      end
c=======================================================================
