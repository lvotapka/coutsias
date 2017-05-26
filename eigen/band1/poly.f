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
C        D Q_k^(nu) = 2*nu*Q_(k-1)^(nu+1)
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
        do i = 0, k-1
          q(i) = 0.d0
        end do
        nu = nu-k
        return
        end
C-----------------------------------------------------------------------
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
        end
C-----------------------------------------------------------------------
