C
C  Tests the routine to compute Jacobi roots
C
cIMPLICIT NONE
cDOUBLE PRECISION PI, PIN, CHEB, XCHEB
cDOUBLE PRECISION ALPHA, BETA
cDOUBLE PRECISION XJAC(65)
cINTEGER IN, N, NP, I
cPI = ACOS(-1.0D0)
cALPHA = -.5D0
cBETA  = -.5D0
cDO IN = 2, 6
c  N  = 2**IN
c  NP = N + 1
c  PIN= PI / N
c  CALL JACOBL(N, ALPHA, BETA, XJAC)
c  WRITE(6,100) N, ALPHA, BETA
c  DO I = 1, NP
c    XCHEB = COS((I-1)*PIN)
c    WRITE(6,101) I, XJAC(I), XCHEB
c  END DO
cEND DO
cSTOP
c100	FORMAT(//5X,'N =',I6,5X,'ALPHA =',D10.3,5X,'BETA =',D10.3/)
c101	FORMAT(I5,5X,D25.15,5X,D25.15)
cEND
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
	CALL JACOBF(NP,ALP,BET,RV,PNP1P,PDNP1P,PNP,PDNP,PNM1P,PDNM1, 1.0D0)
	CALL JACOBF(NP,ALP,BET,RV,PNP1M,PDNP1M,PNM,PDNM,PNM1M,PDNM1,-1.0D0)
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
	    CALL JACOBF(NP,ALP,BET,RV,PNP1,PDNP1,PN,PDN,PNM1,PDNM1,X)
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
C  Hence: this routine only valid for symmetric Jacobis, with alpha=beta.
C  Perhaps eliminating this restriction might allow more generality,
C  But the half computed properly does not seem accurate enough
	DO I = 2, NH
	  XJAC(NPP-I) = -XJAC(I)
	END DO
	IF (N .NE. 2*(N/2)) RETURN
	XJAC(NH+1) = 0.D0
	RETURN
	END
C-----------------------------------------------------------------------
       SUBROUTINE JACOBF(N,ALP,BET,RV,POLY,PDER,POLYM1,PDERM1,POLYM2
     _                                                       ,PDERM2,X)
C
C  computes the Jacobi polynomial (POLY) and its derivative
C  (PDER) of degree N at X
C
	IMPLICIT NONE
	DOUBLE PRECISION ALP, BET, RV
cCOMMON /JACPAR/ALP,BET,RV
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
