ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine oper4(MM,A,B)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c               A * x = (lambda) * B * x                                c
c                                                                       c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'limits.h'
      include 'indata.h'
      integer mm,mc
      integer i, j
      double precision A(0:m_max,0:m_max),   B(0:m_max,0:m_max)
        double precision band1(0:m_max+buff, 2*q_max+1)
        double precision band2(0:m_max+buff, 2*q_max+1)
        double precision band3(0:m_max+buff, 2*q_max+1)
        double precision     H(0:m_max,0:1)
        double precision    CD(0:m_max     ,  1:2)
        double precision    CN(0:m_max     ,  1:2)
        double precision    X1(0:m_max)
      double precision tau,j2,j3,j4,j5,n2,pmode,ii,mi,a1,a0
c------------------------------------------------------------------------
c....Initializations
c------------------------------------------------------------------------
       pmode      = dble(mode)
       n2         = pmode**2
C-----------------------------------------------------------------------
C--- define Dirichlet & Neumann constraints and grid for GET
C-------Still under development;do not use!
        a1 = 1.d0
          do i = 0,2*m_max
               X1(i) = cos(pi*i/(2.d0*m_max))
                 ii = dble(i*i)
                 mi = -a1*ii
            CD(i,1) = 1.d0
            CD(i,2) = a1
            CN(i,1) = ii
            CN(i,2) = mi
            a1        =-a1
          end do
C-----------------------------------------------------------------------
c....build operator matrices
c-----------------------------------------------------------------------
        call cond1_ge(m_max, 0, 1,nu, band1, 0.d0, debug, buff)
        call bfcon(m_max,mm,-1,1,band1,A,buff)
        call cond1_ja(m_max, 0, 1,al,be, band2, 1.d0, debug, buff)
        call bfcon(m_max,mm,-1,1,band2,B,buff)
c       call cond1_ja(m_max, 0, 0,al,be, band1, 1.d0, debug, buff)
c       call bfcon(m_max,mm,0,0,band1,A,buff)
        A(0,0) = 0.d0
        A(1,1) = 0.d0
        write(*,*) 'A'
        call mug(A,m_max,0,6,0,6)
        call mug(A,m_max,mm-10,mm,mm-10,mm)
c       call cond1_mat(m_max, 2, 0, al,be, band2, 1.d0, debug, buff)
c       call bfcon(m_max,mm,-2,2,band2,B,buff)
        B(0,0) = 1.d0
        B(1,1) = 1.d0
c       do i = 0,mm
c         B(0,i) = CD(i,1)
c         B(1,i) = CD(i,2)
c       end do
        write(*,*) 'B'
        call mug(B,m_max,0,6,0,6)
        call mug(B,m_max,mm-10,mm,mm-10,mm)
        stop
c       call mug(B,m_max,15,20,15,20)
        a1 = pi*pi
        a0 = a1/4.d0
        do i = 0,mm/2
         ii = i*i
         j2 = (2*i+1)*(2*i+1)
         write(*,*) a0*j2, a1*ii
        end do
c-----------------------------------------------------------------------
      return
      end
c======================================================================
