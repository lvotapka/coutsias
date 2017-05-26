C-----------------------------------------------------------------------
C           --------<   precond.h   >--------
C-----------------------------------------------------------------------
        double precision     B1(0:m_max+buff,-1:1)
        double precision   B1R1(0:m_max+buff,-2:2)
        double precision   B1R2(0:m_max+buff,-3:3)
        double precision   B1R3(0:m_max+buff,-4:4)
        double precision   B1R4(0:m_max+buff,-5:5)
        double precision   B1R5(0:m_max+buff,-6:6)
        double precision     B2(0:m_max+buff,-2:2)
        double precision   B2R1(0:m_max+buff,-3:3)
        double precision   B2R2(0:m_max+buff,-4:4)
        double precision   B2R3(0:m_max+buff,-5:5)
        double precision   B2R4(0:m_max+buff,-6:6)
        double precision     R1(0:m_max+buff,-1:1)
        double precision     R2(0:m_max+buff,-2:2)
        double precision     R3(0:m_max+buff,-3:3)
        double precision     R4(0:m_max+buff,-4:4)
        double precision     R5(0:m_max+buff,-5:5)
        double precision     R6(0:m_max+buff,-6:6)
        double precision     D1(0:m_max+buff,0:m_max+buff)
        double precision   R1D1(0:m_max+buff,0:m_max+buff)
        double precision R2LAPN(0:m_max+buff,0:m_max+buff)
C-----------------------------------------------------------------------
        common/prec1/B1,B1R1,B1R2,B1R3,B1R4,B1R5,B2,B2R1,B2R2,B2R3,B2R4
     _              ,R1,R2,R3,R4,R5,R6,D1,R1D1,R2LAPN
C-----------------------------------------------------------------------
