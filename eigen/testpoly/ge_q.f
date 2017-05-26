C-----------------------------------------------------------------------
C              ---------< g e _ q . f >--------
C   Compute the gegenbauer polynomial point-evaluation operator Q_x^(k)
C   i.e. (q_0^(k)(x),q_1^(k)(x),...,q_n^(k)(x)) by recurrence relations. 
C   Results remain accurate to ~10^-13 for n=10,000.
C   The recurrence relation used:
C             a1*Q_(n+1) = a3*x*Q_n - a4*Q_(n-1)
C     a1 = (n+1)
C     a3 = 2(n+v)
C     a4 = (n+2v-1)
C     Q_0 = 1, Q_1 = 2vx
C   For the derivatives, 
C        D Q_k^(nu) = nu*Q_(k-1)^(nu+1)
C-----------------------------------------------------------------------
	subroutine ge_q(x,nu,k,n,q)
C-----------------------------------------------------------------------
	implicit none
	double precision q(0:*), x,nu, a1,a3,a4,aa,ab
	integer n, k, i
C-----------------------------------------------------------------------
        nu   = nu+k
        q(0) = 1.d0
        q(1) = 2*nu*x
        do i = 1, n-1
          a1 = (i+1)
          a3 = 2.d0*(nu+i)
          a4 = (i+2*nu-1)
          q(i+1) = (a3*x*q(i)-a4*q(i-1))/a1
        end do
          aa = 1.d0
        do i = 1, k
          aa = 2*aa*(nu-i)
        end do
        do i = n,k,-1
          q(i) = aa*q(i-k)
        end do
        do i = 0
          q(i) = 0.d0
        end do
        nu = nu-k
        return
        end
C-----------------------------------------------------------------------
