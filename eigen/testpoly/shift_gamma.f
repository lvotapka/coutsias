	program Gamma
	implicit none
        double precision pi, e, ex, sq, x, ln, p, Ga
        integer n, i, k
        integer*4 kk
	pi=3.1415926535897932384626d0
        write(*,*) ' sqrt(pi) = ', dsqrt(pi)
        e = dexp(1.d0)**2
c       x = .8d0
          write(*,*) ' x ='
          read(*,*) x
           ex = dble(int(x))-1
           kk = 1
           do i = 1, ex
             kk = kk*i
           end do
           write(*,*) ' int(x)! = ', kk
           n = 2
        do k = 1, 20
           sq= dsqrt(2.d0*pi/(n+x))
           ln= dlog(n+x)
           ex=dexp(x*(ln-1))
           p = 0.d0
          do i = 0,n-1
            p = p+dlog(x+i)
c           write(*,*) i, p
          end do
c           write(*,*) p+n, n*ln
           p = p+n*(1-ln)
c         write(*,*) ' p = ', p
c         write(*,*) 'sq*ex ', sq*ex
           p = dexp(-p)
           Ga = sq*ex*p
          write(*,*) 'n = ',n, ' Gamma(x) = ', Ga
          n = n*2
        end do
         stop
        end
