C-----------------------------------------------------------------------
C -- driver program for testing various routines that generate
c    polynomial values for Jacobi-type polynomials. The recurrence
c    evaluation appears stable (proof?), but the Horner evaluation
c    is quite unstable, and can only be used without stabilization
c    modifications to generate the first few polynomial values, with
c    a few exceptions. Since it is as computationally intensive as the
c    recursive evaluation, it is not worth pursuing. It is however a
c    possible source of fun problems for numerical analysis course.
C-----------------------------------------------------------------------
                         program gam
C-----------------------------------------------------------------------
        implicit none
        double precision qj(0:10000), qg(0:10000), qc(0:10000)
        double precision qa(0:10000), qd(0:10000)
        double precision x, al, be, nu, jac, ge, che
        double precision dgamma, aa, a1, a2
        integer n, i, ii, case, kk
C-----------------------------------------------------------------------
c       open(unit=13,file='in',status='old',form='formatted')
c         read(13,*) x
c         read(13,*) n
c         read(13,*) al
c         read(13,*) be
c         read(13,*) nu
c         read(13,*) case
c         read(13,*) kk
c       close(13)
c       write(*,*) 'case = ', case,nu,al,be, ' deriv = ', kk
C-----------------------------------------------------------------------
 1111   write(*,*) 'x =  ? '
        read(*,*) x
        aa = dgamma(x)
        write(*,*) x, aa
        if (x .ne. 111.d0) goto 1111
C-----------------------------------------------------------------------
        end
C-----------------------------------------------------------------------
