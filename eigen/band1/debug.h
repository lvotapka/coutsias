C-----------------------------------------------------------------------
C            ------------<  d e b u g . h   >-------------
C-----------------------------------------------------------------------
C  Arrays used for debugging purposes
C-----------------------------------------------------------------------
        double precision  PH0(  0:m_max     , 0:1         , 0:n_max/2  )
        double precision   P0(  0:m_max+buff, 1:2*q_max+1 , 0:n_max/2  )
        double precision CHP0(  1:2         , 0:1         , 0:n_max/2  )
        double precision   X0(  0:m_max     , 0:1         , 0:n_max/2+1)
        double precision   bb0(0:m_max)
        double precision   cc0(0:m_max)
C-----------------------------------------------------------------------
        common/trash/P0,PH0,CHP0,X0,bb0,cc0
C-----------------------------------------------------------------------
