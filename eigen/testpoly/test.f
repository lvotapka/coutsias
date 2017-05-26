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
                         program test
C-----------------------------------------------------------------------
        implicit none
        double precision qj(0:10000), qg(0:10000), qc(0:10000)
        double precision qa(0:10000), qd(0:10000)
        double precision endpts(2)
        double precision x, al, be, nu, jac, ge, che
        double precision aa, a1, a2
        integer n, i, ii, case, kk
C-----------------------------------------------------------------------
        DOUBLE PRECISION PI, PIN, CHEB, XCHEB
        INTEGER IN, NP
        PI = ACOS(-1.0D0)
        open(unit=13,file='in',status='old',form='formatted')
          read(13,*) x
          read(13,*) n
          read(13,*) al
          read(13,*) be
          read(13,*) nu
          read(13,*) case
          read(13,*) kk
        close(13)
        write(*,*) 'case = ', case,nu,al,be, ' deriv = ', kk
C-----------------------------------------------------------------------
      endpts(1) =-1.d0
      endpts(2) = 1.d0
        DO IN = 2, 6
          N  = 2**IN
          WRITE(6,*) N, AL, BE
          NP = N + 1
          PIN= PI / N
      call gaussq(case, np, al, be, 2, endpts, qa(0), qd(0), qj(0))
          CALL JACOBL(N, AL, BE, qc(0))
          DO i = 0, N
            XCHEB = COS(I*PIN)
            WRITE(6,*) I, qc(I), qd(i), XCHEB
          END DO
        END DO
        end
C-----------------------------------------------------------------------
