C-----------------------------------------------------------------------
C             -----------< j a c _ p o l . f >----------
C Horner evaluation of jacobi polynomials, using table on p.789,
C of Abramowitz and Stegun. No better than 10^-13 up to 10,000 order
C-----------------------------------------------------------------------
	subroutine jac_pol(x,al,be,n,jac)
C-----------------------------------------------------------------------
	implicit none
	double precision x,al,be,jac
	integer n
	double precision a,b,c,d
        integer i
C-----------------------------------------------------------------------
          a = 1.d0
	do i = n,1,-1
          b = (n-i+1.d0)*(al+be+n+i)
          c = 2.d0*i*(al+i)
          a = 1.d0-b*(1.d0-x)*a/c
        end do
          d = 1.d0
        do i = 1,n
          d = d*(i+al)/i
        end do
        jac = d*a
        return
        end
C-----------------------------------------------------------------------
