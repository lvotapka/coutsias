C-----------------------------------------------------------------------
C            ------------<   3 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator P could be arbitrary second order operator, of
C  bandwidth <= l_max
C-----------------------------------------------------------------------
        double precision     H3(0:m_max     ,  0:1       , 0:n_max/2)
        double precision     O3(0:m_max+buff,  1:3       , 0:n_max/2)
        double precision     C3(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V3(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH3(1:2         ,  0:1       , 0:n_max/2)
C-----------------------------------------------------------------------
	common/tri/O3,H3,C3,V3,CH3
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C            ------------<   5 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator P could be arbitrary second order operator, of
C  bandwidth <= l_max
C-----------------------------------------------------------------------
        double precision     H5(0:m_max     ,  0:1       , 0:n_max/2)
        double precision     O5(0:m_max+buff,  1:5       , 0:n_max/2)
        double precision     C5(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V5(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH5(1:2         ,  0:1       , 0:n_max/2)
        double precision  alpha(0:2         ,  0:2                  )
C-----------------------------------------------------------------------
	common/penta/O5,H5,C5,V5,CH5,alpha
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C            ------------<   7 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator H could be arbitrary second order operator, of
C  bandwidth <= l_max = 7
C-----------------------------------------------------------------------
        double precision     H7(0:m_max     ,  0:1       , 0:n_max/2)
        double precision     O7(0:m_max+buff,  1:7       , 0:n_max/2)
        double precision     C7(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V7(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH7(1:2         ,  0:1       , 0:n_max/2)
C-----------------------------------------------------------------------
	common/hepta/O7,H7,C7,V7,CH7
C-----------------------------------------------------------------------
C            ------------<   9 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator H could be arbitrary second order operator, of
C  bandwidth <= l_max = 9
C-----------------------------------------------------------------------
        double precision     H9(0:m_max     ,  0:1       , 0:n_max/2)
        double precision     O9(0:m_max+buff,  1:9       , 0:n_max/2)
        double precision     C9(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V9(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH9(1:2         ,  0:1       , 0:n_max/2)
C-----------------------------------------------------------------------
	common/ennea/O9,H9,C9,V9,CH9
C-----------------------------------------------------------------------
C            ------------<  bc . h   >------------
C-----------------------------------------------------------------------
        double precision     CD(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     CN(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     NS(0:m_max     ,  1:2       , 0:n_max/2)
C-----------------------------------------------------------------------
	common/constraints/CD,CN,NS
C-----------------------------------------------------------------------
