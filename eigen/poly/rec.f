C-----------------------------------------------------------------------
C              ---------< r e c X . f >--------
C  returns the recurrence coefficients (standard scaling) for the selected
C  polynomial family
C  All cases: k not 0 means the computation is intended for BC and result
C  will be shifted, hence unsuitable for Fourier matrix calculation
C  currently, only k=1 can be used with Chebyshev. 
C  Kind=3;  works only for BC
C  Kind=7 works here, but fourier matrix calculation not im plemented.
C   jacobi
C             a1*Q_(n+1) = (a2+a3*x)*Q_n - a4*Q_(n-1)
C     a1 = 2(n+1)(n+a+b+1)(2n+a+b)
C     a2 = (2n+a+b+1)(a^2-b^2)
C     a3 = (2n+a+b)_3
C     a4 = 2(n+a)(n+b)(2n+a+b+2)
C     Q_0 = 1, Q_1 = (a-b)/2 + (a+b+2)x/2
C   chebyshev
C             Q_(n+1) = 2*x*Q_n - Q_(n-1)
C     Q_0 = 1, Q_1 = x
C   gegenbauer
C             a1*Q_(n+1) = a3*x*Q_n - a4*Q_(n-1)
C     a1 = (n+1)
C     a3 = 2(n+v)
C     a4 = (n+2v-1)
C     Q_0 = 1, Q_1 = 2vx
C   For the derivatives, 
C        D Q_k^(al) = 2*al*Q_(k-1)^(al+1)
C-----------------------------------------------------------------------
	subroutine rec_X(m,A,B,kind,al,be)
	implicit none
        integer m, i, kind, j
	double precision al,be, a0,aa,ab
	double precision A(0:m,-1:1), B(0:m,-1:1)
C-----------------------------------------------------------------------
        do i = 0, m
          do j = -1,1
            A(i,j) = 0.d0
            B(i,j) = 0.d0
          end do
        end do
C-----------------------------------------------------------------------
        go to (10, 20, 30, 40, 50, 60, 10) kind
C--- Legendre or Gegenbauer
 10       A(0, 1) = al/(al+1.d0)
        do i = 1, m
               aa = .5d0/(al+dble(i)-1.d0)
               ab = .5d0/(al+dble(i)+1.d0)
          A(i,-1) = aa* dble(i)
          A(i, 1) = ab*(dble(i)+2.d0*al)
          B(i,-1) = aa
          B(i, 1) =-ab
        end do
        return
C-----------------------------------------------------------------------
C--- Chebyshev
 20       A(0, 1) = 0.5d0
        do i = 1, m
              aa = .5d0/dble(i)
          A(i,-1) = .5d0
          A(i, 1) = .5d0
          B(i,-1) = aa
          B(i, 1) =-aa
        end do
          A(1,-1) = 1.d0
          B(1,-1) = 1.d0
        return
C for the chebyshev polynomials of the second kind, use scaled gegenbauers.
 30     write(*,*) 'ANOTHER HOLE TO PATCH YOU MORON!'
        stop
        return
C-----------------------------------------------------------------------
C---- Hermite polynomials
 40       A(0, 1) = dble(i)+1.d0
        do i = 1, m
          A(i,-1) = .5d0
          A(i, 1) = dble(i)+1.d0
          B(i,-1) = .5d0/dble(i)
        end do
        return
C-----------------------------------------------------------------------
C---- Jacobi polynomials
 50     aa   = al+be
        ab   = al*al-be*be
        A(0, 0) = -ab/(aa*(aa+2))
        A(0, 1) = 2.d0*(al+1)*(be+1)/((aa+2)*(aa+3)) 
        do i = 1, m
          a0 = 2*dble(i)+aa
          A(i,-1) = 2*dble(i)*(dble(i)+aa)/(a0*(a0-1.d0))
          A(i, 0) = -ab/(a0*(a0+2))
          A(i, 1) = 2*(dble(i)+al+1)*(dble(i)+be+1)/((a0+2.d0)*(a0+3))
          B(i,-1) =  A(i,-1)/dble(i)
          B(i, 1) = -A(i, 1)/(dble(i)+aa+1.d0)
          B(i, 0) = -2.d0*A(i,0)/aa
        end do
        return
C-----------------------------------------------------------------------
C---- Laguerre polynomials
 60            aa = al+1.d0
          A(0, 1) =-aa
          A(0, 0) = aa
        do i = 1, m
          A(i,-1) = -dble(i)
          A(i, 1) = -(aa+dble(i))
          A(i, 0) = 2.d0*dble(i)+aa
          B(i,-1) = -1.d0
          B(i, 0) =  1.d0
        end do
        return
        end
C-----------------------------------------------------------------------
