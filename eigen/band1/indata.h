c---------------------------------------------------
c           ------<<<  indata.h   >>>-----
c---------------------------------------------------
c... variables read in 'in'-file
c---------------------------------------------------
c......physical variables
c---------------------------------------------------
      double precision pi, nu, al, be
      integer          oper_case,bc_case,mode
c---------------------------------------------------
c......control variables
c---------------------------------------------------
      double precision window,perturb,alpha_tol
      integer          ddebug,matz,evl_write,m_0,m_f
      integer          step,revers,trans
      integer          debug, count_max
      character*36     evlnam, evcnam
c---------------------------------------------------
c......common definition
c---------------------------------------------------
      common /din_dat/ pi,nu,al,be,window,perturb,alpha_tol
      common /iin_dat/ count_max,oper_case,bc_case,mode,
     _               ddebug,matz,evl_write,m_0,m_f,
     _               step,revers,trans,debug 
      common/ cin_dat/ evlnam, evcnam
