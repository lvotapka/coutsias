C-----------------------------------------------------------------------
                PROGRAM EIGEN
C-----------------------------------------------------------------------
C                                                              Version I
	implicit none
	include 'limits.h'
	include 'indata.h'
C--- Include banded arrays for various solvers
	include '3579.h'
	include 'ops.h'
c-----------------------------------------------------------------------
c....Variables for Evals
c-----------------------------------------------------------------------
      DOUBLE PRECISION   AT(0:m_max,0:m_max),BT(0:m_max,0:m_max)
      DOUBLE PRECISION    Z(0:m_max,0:m_max)
      DOUBLE PRECISION ALFR(0:m_max),        ALFI(0:m_max)
      DOUBLE PRECISION BETA(0:m_max),        DIST(0:m_max)
      DOUBLE PRECISION                    EVAL_BETA(        2)
      DOUBLE PRECISION PREVR(0:m_max),      PREVI(0:m_max)
      INTEGER            TAG(0:m_max),       TYPE(0:m_max)
      INTEGER          EVLFIL,   EVCFIL
             PARAMETER(EVLFIL=23,EVCFIL=25)
c-----------------------------------------------------------------------
c....Variables for iteration 
c-----------------------------------------------------------------------
      INTEGER COUNT, CONV
C-----------------------------------------------------------------------
C various work arrays
	double precision    A(0:m_max,0:m_max)
	double precision    B(0:m_max,0:m_max)
        double precision   Y1(0:2*m_max,0:1)
        double precision   Y2(0:2*m_max,0:1)
C--------------------------
C-- Initializations
        integer         lu_in
              parameter(lu_in = 17)
	call lin_indata(lu_in)
C-----------------------------------------------------------------------
C PART C      Main Loop
c-----------------------------------------------------------------------
      if (oper_case .eq. 1 .or. oper_case .eq. 11) call oper1(M_F,A,B)
      if (oper_case .eq. 2) call oper2(M_F,A,B)
      if (oper_case .eq. 3) call oper3(M_F,A,B)
      if (oper_case .eq. 4) call oper4(M_F,A,B)
      CALL COMP_EVAL(m_max+1,M_0,M_F,A,B,STEP,WINDOW,MATZ,PERTURB,
     _               AT,BT,TRANS,REVERS,DDEBUG,EVL_WRITE,EVLNAM,EVCNAM,
     _               EVLFIL,EVCFIL,COUNT,COUNT_MAX,DIST,
     _               Z,ALFR,ALFI,BETA,PREVR,PREVI,EVAL_BETA,TAG,TYPE)
      write(*,*) eval_beta(1), eval_beta(2)
C-----------------------------------------------------------------------
C PART D
C-----------------------------------------------------------------------
        WRITE(*,*) 'Real part of eval =',EVAL_BETA(1)
        WRITE(*,*) 'Imag part of eval =',EVAL_BETA(2)
        end
C-----------------------------------------------------------------------
