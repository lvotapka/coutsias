C-----------------------------------------------------------------------
C             -----------< c h e _ p o l . f >----------
C Horner evaluation of chebyshev polynomials, using table on p.789,
C of Abramowitz and Stegun. No better than 10^-13 up to 10,000 order
C at x=1; however everywhere else this method is UNSTABLE!
C-----------------------------------------------------------------------
	subroutine che_pol(x,n,che)
C-----------------------------------------------------------------------
	implicit none
	double precision x,che
	integer n
	double precision a,b,c,d
        integer i,nn, ni
C-----------------------------------------------------------------------
          a = 1.d0
          ni = n/2
          nn = n-ni*2
        if (nn .eq. 0) then
          do i = ni,1,-1
            b = 2*(ni-i+1.d0)*(ni+i-1)
            c = i*(2*i-1)
            a = 1.d0-b*x*x*a/c
          end do
            d = 1.d0
        else
          do i = ni,1,-1
            b = 2*(ni-i+1.d0)*(ni+i)
            c = i*(2*i+1)
            a = 1.d0-b*x*x*a/c
          end do
            d = x*(2*ni+1)
        endif
          do i = 0,ni-1,1
            d = -d
          end do
          ge = d*a
        return
        end
C-----------------------------------------------------------------------
