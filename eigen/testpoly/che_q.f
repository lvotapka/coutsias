C-----------------------------------------------------------------------
C              ---------< c h e _ q . f >--------
C   Compute the chebyshev polynomial point-evaluation operator Q_x^n
C   i.e. (q_0(x),q_1(x),...,q_n(x) by recurrence relations. Results
C   remain accurate to ~10^-13 for n=10,000.
C   The recurrence relation used:
C             Q_(n+1) = 2*x*Q_n - Q_(n-1)
C     Q_0 = 1, Q_1 = x
C-----------------------------------------------------------------------
	subroutine che_q(x,k,n,q)
C-----------------------------------------------------------------------
	implicit none
	double precision q(0:*), x
	integer n, k, i
C-----------------------------------------------------------------------
        if (k .eq. 0) then
        q(0) = 1.d0
        q(1) = x
        do i = 1, n-1
          q(i+1) = 2*x*q(i)-q(i-1)
        end do
        else
          call ge_q(x,1.d0,k-1,n,q)
          do i = n,1,-1
            q(i) = i*q(i-1)
          end do
            q(0) = 0.d0
        endif
        return
        end
C-----------------------------------------------------------------------
