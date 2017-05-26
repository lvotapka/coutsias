      integer         m_max   , n_max
            parameter(m_max=256, n_max=2)
      integer         p_max
            parameter(p_max = max(m_max,n_max))
      integer         q_max
            parameter(q_max = 4)
      integer         buff
            parameter(buff = 36)

      integer         m_max_act, n_max_act
            parameter(m_max_act  = (m_max/3)*2 + 2)
            parameter(n_max_act  = (n_max/3)*2 + 2)


