C-----------------------------------------------------------------------
C            ------------<   5 . h   >------------
C-----------------------------------------------------------------------
C  Arrays used for the solution of annular poisson problem
C  written so operator P could be arbitrary second order operator, of
C  bandwidth <= l_max
C-----------------------------------------------------------------------
        double precision     H5(0:m_max+1   ,  0:1       , 0:n_max/2)
        double precision     O5(0:m_max+buff,  1:5       , 0:n_max/2)
        double precision     C5(0:m_max     ,  1:2       , 0:n_max/2)
        double precision     V5(1:2         ,  0:1       , 0:n_max/2)
        double precision    CH5(1:2         ,  0:1       , 0:n_max/2)
        double precision  alpha(0:2         ,  0:2                  )
C-----------------------------------------------------------------------
	common/penta/O5,H5,C5,V5,CH5,alpha
C-----------------------------------------------------------------------
