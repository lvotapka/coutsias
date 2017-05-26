C-----------------------------------------------------------------------
C              ---------< p o l y . f >--------
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
	subroutine jac_q(kind,x,al,be,k,n,q)
	implicit none
	double precision q(0:*), x,al,be, a1,a2,a3,a4,aa,ab
	integer n,k,i,ii, kind
C-----------------------------------------------------------------------
        go to (10, 20, 30, 40, 50, 60, 70) kind
 20     if (k .eq. 0) then
        q(0) = 1.d0
        q(1) = x
        do i = 1, n-1
          q(i+1) = 2*x*q(i)-q(i-1)
        end do
        else
C for the chebyshev polynomials of the second kind, use scaled gegenbauers.
C The present version is not suitable for computations, as it shifts the 
C orders. Rather, it is useful only for BC, as are all forms in this routine.
C Should include a switch, that causes all things to be shifted by k places
C if they are to be used for point evaluation functionals of boundary
C conditions, or not shifted if they are to be used to construct transform
C matrices. ANOTHER THING UNDER CONSTRUCTION.
 30       al = 1.d0
          k = 0
          go to 70
 35       continue
          do i = n,1,-1
            q(i) = i*q(i-1)
          end do
            q(0) = 0.d0
        endif
        return
C-----------------------------------------------------------------------
 10     continue
        al   = .5d0
 70     al   = al+k
        q(0) = 1.d0
        q(1) = 2*al*x
        do i = 1, n-1
          a1 = (i+1)
          a3 = 2.d0*(al+i)
          a4 = (i+2*al-1)
          q(i+1) = (a3*x*q(i)-a4*q(i-1))/a1
        end do
          aa = 1.d0
        do i = 1, k
          aa = 2*aa*(al-i)
        end do
        do i = n,k,-1
          q(i) = aa*q(i-k)
        end do
        do i = 0, k-1
          q(i) = 0.d0
        end do
        al = al-k
        if (kind .eq. 3) goto 35
        return
C-----------------------------------------------------------------------
 50     al   = al+k
        be   = be+k
        q(0) = 1.d0
        q(1) = (al-be)/2+x*(al+be+2)/2
        ab   = al*al-be*be
        do i = 1, n-1
          aa = 2*i+al+be
          a1 = 2.d0*(i+1)*(i+al+be+1)*aa
          a2 = (aa+1)*ab
          a3 = aa*(aa+1)*(aa+2)
          a4 = 2*(i+al)*(i+be)*(aa+2)
          q(i+1) = ((a2+a3*x)*q(i)-a4*q(i-1))/a1
        end do
C-- compute derivative factors and shift up.
          al = al-k
          be = be-k
          do ii = n,k,-1
            a1 = (ii+al+be)/2.d0
            a2 = 1.d0
            do i = 1,k
              a2 = a2*(a1+i/2.d0)
            end do
            q(ii) = a2*q(ii-k)
          end do
          do i = 0,k-1
            q(i) = 0.d0
          end do
        return
C---- Hermite 
 40     q(0) = 1.d0
        q(1) = 2*x
        do i = 1, n-1
          q(i+1) = (x*q(i) - .5d0*q(i-1))/(i+1.d0)
        end do
        return
C---- Laguerre 
 60     q(0) = 1.d0
        q(1) = 1.d0+al-x
        do i = 1, n-1
          q(i+1) = ((2*i+al+1.d0-x)*q(i)-(al+i)*q(i-1))/(i+1.d0)
        end do
        return
        end
C-----------------------------------------------------------------------
