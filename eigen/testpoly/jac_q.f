C-----------------------------------------------------------------------
C              ---------< j a c _ q . f >--------
C   Compute the jacobi polynomial point-evaluation operator Q_x^n
C   i.e. (q_0(x),q_1(x),...,q_n(x) by recurrence relations. Results
C   remain accurate to ~10^-13 for n=10,000.
C   The recurrence relation used:
C             a1*Q_(n+1) = (a2+a3*x)*Q_n - a4*Q_(n-1)
C     a1 = 2(n+1)(n+a+b+1)(2n+a+b)
C     a2 = (2n+a+b+1)(a^2-b^2)
C     a3 = (2n+a+b)_3
C     a4 = 2(n+a)(n+b)(2n+a+b+2)
C     Q_0 = 1, Q_1 = (a-b)/2 + (a+b+2)x/2
C-----------------------------------------------------------------------
	subroutine jac_q(x,al,be,k,n,q)
C-----------------------------------------------------------------------
	implicit none
	double precision q(0:*), x,al,be, a1,a2,a3,a4,aa,ab
	integer n,k,i,ii
C-----------------------------------------------------------------------
        al   = al+k
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
        end
C-----------------------------------------------------------------------
