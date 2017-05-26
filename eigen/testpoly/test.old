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
        double precision x, al, be, nu, jac, ge, che
        double precision aa, a1, a2
        integer n, i, ii, case, kk
C-----------------------------------------------------------------------
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
        if (case .eq. 1) then
        call jac_q(x,al,be,kk,n,qj)
        do i = kk, n
            a1 = (i+al+be)/2.d0
            a2 = 1.d0
            do ii = 1,kk
              a2 = a2*(a1+ii/2.d0)
            end do     
          call jac_pol(x,al+kk,be+kk,i-kk,jac)
            jac = a2*jac
          write(*,*) i, jac, qj(i)
        end do
        endif
C-----------------------------------------------------------------------
        if (case .eq. 2) then
        call ge_q(x,nu,kk,n,qg)
            aa = 1.d0
          do i = 0,kk-1
            aa = 2*aa*(nu+i)
          end do
        do i = kk, n
          call ge_pol(x,nu+kk,i-kk,ge)
          write(*,*) i, aa*ge, qg(i)
        end do
        endif
C-----------------------------------------------------------------------
        if (case .eq. 3) then
        call che_q(x,kk,n,qc)
            aa = 1.d0
          do i = 1,kk-1
            aa = aa*i
          end do
        do i = kk, n
        if (kk .eq. 0) then
          call che_pol(x,i,che)
        else
          call ge_pol(x,dble(kk),i-kk,che)
          che = i*aa*che
        end if
          write(*,*) i, che, qc(i)
        end do
        endif
        end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C             -----------< g e _ p o l . f >----------
C Horner evaluation of gegenbauer polynomials, using table on p.789,
C of Abramowitz and Stegun. No better than 10^-13 up to 10,000 order
C-----------------------------------------------------------------------
	subroutine ge_pol(x,nu,n,ge)
C-----------------------------------------------------------------------
	implicit none
	double precision x,nu,ge
	integer n
	double precision a,b,c,d
        integer i,nn, ni
C-----------------------------------------------------------------------
          a = 1.d0
          ni = n/2
          nn = n-ni*2
        if (nn .eq. 0) then
          do i = ni,1,-1
            b = 2*(ni-i+1.d0)*(nu+ni+i-1)
            c = i*(2*i-1)
            a = 1.d0-b*x*x*a/c
          end do
            d = 1.d0
        else
          do i = ni,1,-1
            b = 2*(ni-i+1.d0)*(nu+ni+i)
            c = i*(2*i+1)
            a = 1.d0-b*x*x*a/c
          end do
            d = 2*x*(nu+ni)
        endif
          do i = 0,ni-1,1
            d = -d*(i+nu)/(i+1)
          end do
          ge = d*a
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
          che = d*a
        return
        end
C-----------------------------------------------------------------------
C
C  Tests the routine to compute Jacobi roots
C
        SUBROUTINE TEST_JACOL(N,ALPHA,BETA)
	IMPLICIT NONE
	DOUBLE PRECISION PI, PIN, CHEB, XCHEB
	DOUBLE PRECISION ALPHA, BETA
	DOUBLE PRECISION XJAC(65)
	INTEGER IN, N, NP, I
	PI = ACOS(-1.0D0)
	ALPHA = -.5D0
	BETA  = -.5D0
	DO IN = 2, 6
	  N  = 2**IN
	  NP = N + 1
	  PIN= PI / N
	  CALL JACOBL(N, ALPHA, BETA, XJAC)
	  WRITE(6,100) N, ALPHA, BETA
	  DO I = 1, NP
	    XCHEB = COS((I-1)*PIN)
	    WRITE(6,101) I, XJAC(I), XCHEB
	  END DO
	END DO
	STOP
 100	FORMAT(//5X,'N =',I6,5X,'ALPHA =',D10.3,5X,'BETA =',D10.3/)
 101	FORMAT(I5,5X,D25.15,5X,D25.15)
	END
C-----------------------------------------------------------------------
	SUBROUTINE JACOBL(N,ALPHA,BETA,XJAC)
C
C  Computes the gauss-Lobato collocation points for Jacobi polynomials
C
C  N:		Degree of approximation
C  ALPHA:	parameter in Jacobi weight
C  BETA:	parameter in Jacobi weight
C
C  XJAC:	output array with the Gauss-Lobatto roots
C		(they are ordered from largest (+1.) to smallest (-1.)
C
	IMPLICIT NONE
	DOUBLE PRECISION XJAC(65)
	DOUBLE PRECISION ALP, BET, RV
	DOUBLE PRECISION ALPHA, BETA, EPS
	COMMON /JACPAR/ALP, BET, RV
	INTEGER N, KSTOP, NP, NH, J, K, I, JM, I, NPP
	DOUBLE PRECISION PNP1P,PDNP1P,PNP,PDNP,PNM1P,PDNM1
	DOUBLE PRECISION PNP1M,PDNP1M,PNM,PDNM,PNM1M
	DOUBLE PRECISION PNP1, PDNP1, PN, PDN, PNM1
	DOUBLE PRECISION POLY, PDER
	DOUBLE PRECISION RECSUM, DELX
	DOUBLE PRECISION DET, RP, RM, A, B, DTH, CD, CS, SD, SS
	DOUBLE PRECISION X, CSSAVE
	DATA KSTOP/10/
	DATA EPS/1.0D-12/
	ALP = ALPHA
	BET = BETA
	RV = 1.d0 + ALP
	NP = N + 1
C
C  Compute the parameters in the polynomial whose roots are desired
C
	CALL JACOBF(NP,PNP1P,PDNP1P,PNP,PDNP,PNM1P,PDNM1, 1.0D0)
	CALL JACOBF(NP,PNP1M,PDNP1M,PNM,PDNM,PNM1M,PDNM1,-1.0D0)
	DET = PNP*PNM1M-PNM*PNM1P
	RP  =-PNP1P
	RM  =-PNP1M
	A   =(RP*PNM1M-RM*PNM1P)/DET 
	B   =(RM*PNP  -RP*PNM  )/DET 
C
	XJAC(1) = 1.D0
	NH = (N+1)/2
C
C  Set-up recursion relation for initial guess for the roots
C
	DTH = 3.1415926535897932384626D0/(2*N+1)
	CD = COS(2.D0*DTH)
	SD = SIN(2.D0*DTH)
	CS = COS(DTH)
	SS = SIN(DTH)
C
C  Compute the first half of the roots by polynomial deflation
C
	DO J = 2, NH
	  X = CS
	  DO K = 1, KSTOP
	    CALL JACOBF(NP,PNP1,PDNP1,PN,PDN,PNM1,PDNM1,X)
	    POLY = PNP1 +A*PN +B*PNM1
	    PDER = PDNP1+A*PDN+B*PDNM1
	    RECSUM = 0.D0
	    JM = J - 1
	    DO I = 1, JM
	      RECSUM = RECSUM + 1.D0/(X - XJAC(I))
            END DO
	    DELX = -POLY/(PDER-RECSUM*POLY)
            X = X + DELX
            IF(ABS(DELX) .LT. EPS) GO TO 30
	  END DO
 30	  CONTINUE
	  XJAC(J) = X
	  CSSAVE = CS*CD-SS*SD
	   SS    = CS*SD+SS*CD
	   CS    = CSSAVE
	END DO
	XJAC(NP) = -1.D0
	NPP      = N+2
C
C  use symmetry for second half of the roots
C
	DO I = 2, NH
	  XJAC(NPP-I) = -XJAC(I)
	END DO
	IF (N .NE. 2*(N/2)) RETURN
	XJAC(NH+1) = 0.D0
	RETURN
	END
C-----------------------------------------------------------------------
	SUBROUTINE JACOBF(N,POLY,PDER,POLYM1,PDERM1,POLYM2,PDERM2,X)
C
C  computes the Jacobi polynomial (POLY) and its derivative
C  (PDER) of degree N at X
C
	IMPLICIT NONE
	DOUBLE PRECISION ALP, BET, RV
	COMMON /JACPAR/ALP,BET,RV
	DOUBLE PRECISION POLY, PDER, POLYM1, PDERM1, POLYM2, PDERM2
	DOUBLE PRECISION X, APB, POLYLST, PDERLST
	DOUBLE PRECISION A1, A2, A3, B3, A4
	DOUBLE PRECISION PSAVE, PDSAVE, POLYN, PDERN
	INTEGER N, K
	APB = ALP + BET
	POLY = 1.D0
	PDER = 0.D0
	IF (N .EQ. 0) RETURN   
	POLYLST = POLY
	PDERLST = PDER
	POLY = RV*X
	PDER = RV
	IF (N .EQ. 1) RETURN
	DO K = 2, N
	  A1 = 2.D0*K*(K+APB)*(2.D0*K+APB-2.D0)
	  A2 = (2.D0*K+APB-1.D0)*(ALP**2-BET**2)
	  B3 = (2.D0*K+APB-2.D0)
	  A3 = B3*(B3+1.D0)*(B3+2.D0)
	  A4 = 2.D0*(K+ALP-1.D0)*(K+BET-1.D0)*(2.D0*K+APB)
	  POLYN = ((A2+A3*X)*POLY-A4*POLYLST)/A1
	  PDERN = ((A2+A3*X)*PDER-A4*PDERLST+A3*POLY)/A1
	  PSAVE  = POLYLST
	  PDSAVE = PDERLST
	  POLYLST= POLY
	  POLY   = POLYN
	  PDERLST= PDER
	  PDER   = PDERN
	END DO
	POLYM1 = POLYLST
	PDERM1 = PDERLST
	POLYM2 = PSAVE
	PDERM2 = PDSAVE
	RETURN
	END
C-----------------------------------------------------------------------
c        program costest
cimplicit none
c     INCLUDE 'param.h'
c     INCLUDE 'cost.h'
c     DOUBLE PRECISION PI
c     PARAMETER(PI=3.1415926535897932384626D0)
cdouble precision DATA(M_MAX+2,N_MAX+2), TRANSFORM(M_MAX+2)
cdouble precision DAT2(M_MAX+2,N_MAX+2)
cinteger M, N, I, J, II
c       N = 1
c       READ(*,*) M
c       do i = 0,M
c          ii = i+1
c          do j = 1,N+1
c             DATA(II,J) = 0.d0
c             DAT2(II,J) = 0.d0
c          end do
c          DATA(II,1) = 1.d0
c          DATA(II,1) = cos(PI*I/M)
c          DATA(II,1) = 2.d0*cos(PI*I/M)**2-1.d0
c          DAT2(II,1) = DATA(II,1)
c       end do
c       CALL COSFTINI(M)
cCALL COST(DATA,M,TRANSFORM,1)
cCALL COST(DATA,M,TRANSFORM,-1)
c       write(*,*) (DATA(I,1)-DAT2(I,1), I=1,M+1)
c       write(*,*) (DATA(I,1), I=1,M+1)
c       END
C------------------------------------------------------
C                ---- COSFOU.FOR ----
C------------------------------------------------------
      SUBROUTINE COSFTINI(M)
      INTEGER M, I,J
      DOUBLE PRECISION PI
      PARAMETER(PI=3.1415926535897932384626D0)
      INCLUDE 'limits.h'
      INCLUDE 'cost.h'
      DO I=0,M
         DO J=0,M
            CHEBT(I,J) = cos(PI*I*J/M)
         END DO
      END DO
      END
C--------------------------------------------------------
      SUBROUTINE COST (DATA, M, TRANSFORM, ISIGN)
C--------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'limits.h'
      INCLUDE 'cost.h'
      INTEGER M, ISIGN, I, J, II, JJ
      DOUBLE PRECISION DATA(M_MAX+2), TRANSFORM(M_MAX+2)
      DO I = 0, M
         II = I+1
         IF(ISIGN .EQ. 1) THEN
            TRANSFORM(II) = .5D0*DATA(1)
         ELSE
            TRANSFORM(II) =      DATA(1)
         END IF
         DO J = 1, M-1
            JJ = J+1
            TRANSFORM(II) = TRANSFORM(II)+ CHEBT(I,J)*DATA(JJ)
         END DO
         IF(ISIGN .EQ. 1) THEN
        TRANSFORM(II) = (2.d0/M)*(TRANSFORM(II)+.5D0*(-1)**I *DATA(M+1))
         ELSE
        TRANSFORM(II) =           TRANSFORM(II)+     (-1)**I *DATA(M+1)
         END IF
      END DO
      IF (ISIGN .EQ. 1) THEN
         TRANSFORM(  1) = .5d0*TRANSFORM(  1)
         TRANSFORM(M+1) = .5d0*TRANSFORM(M+1)
      ENDIF

      END
C--------------------------------------------------------------------
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
