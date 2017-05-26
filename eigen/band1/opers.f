ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine oper1(MM,A,B)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                       c
c               A * x = (lambda) * B * x                                c
c     Chebyshev Laplace operator, Integrator version                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      include 'limits.h'
      include 'indata.h'
      integer mm,mc
      integer i, j, k
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
C--- define Dirichlet & Neumann constraints and grid for FFT
          do i = 0,m_max
               X1(i) = cos(pi*i/(2.d0*m_max))
          end do
          call che_q( 1.d0, 0, m_max, CD(0,1))
          call che_q(-1.d0, 0, m_max, CD(0,2))
          call che_q( 1.d0, 1, m_max, CN(0,1))
          call che_q(-1.d0, 1, m_max, CN(0,2))
C-----------------------------------------------------------------------
c....build operator matrices
c-----------------------------------------------------------------------
        call cond1_mat(m_max, 0, 0, band1, 0.d0, debug, buff)
        call cond1_mat(m_max, 2, 0, band2, 0.d0, debug, buff)
        if (oper_case .eq. 11) then
          do i = 2,mm
              band1(i,1) = band1(i,1)*dble(i)
            do k = 1,5
              band2(i,k) = band2(i,k)*dble(i)
            end do
c           write(*,*) i, (band2(i,k), k=1,5)
          end do
        end if
        call bfcon(m_max,mm,-2,2,band2,B,buff)
        call bfcon(m_max,mm, 0,0,band1,A,buff)
        do i = 0,mm
          B(0,i) = 0.d0
          B(1,i) = 0.d0
          A(0,0) = 0.d0
          A(1,1) = 0.d0
          if (bc_case .eq. 1) then
            B(0,i) = CD(i,1)
            B(1,i) = CD(i,2)
          end if
          if (bc_case .eq. -1) then
            A(0,i) = CD(i,1)
            A(1,i) = CD(i,2)
          end if
          if (bc_case .eq. 2) then
            B(0,i) = CN(i,1)
            B(1,i) = CN(i,2)
          end if
          if (bc_case .eq. -2) then
            A(0,i) = CN(i,1)
            A(1,i) = CN(i,2)
          end if
        end do
        if (bc_case .eq. 0) then
          B(0,0) = 1.d0
          B(1,1) = 1.d0
        end if
        if (bc_case .eq. -100) then
          A(0,0) = 1.d0
          A(1,1) = 1.d0
        end if
        write(*,*) 'A'
        call mug(A,m_max,0,6,0,6)
c       call mug(A,m_max,mm-10,mm,mm-10,mm)
        write(*,*) 'B'
        call mug(B,m_max,0,10,0,10)
c       call mug(B,m_max,mm-10,mm,mm-10,mm)
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
c matout.f
c======================================================================
      subroutine matmult1(a,b,h,maxm,mm,ma,nb,namb)

      double precision a(maxm,mm), b(maxm,mm),
     _  h(maxm,mm)
      integer maxm,mm,ma, nb, namb

      do 410 i=1, ma
        do 405 j=1, nb
          h(i,j)=0.0d0
          do 400 k=1, namb
            h(i,j)=h(i,j)+a(i,k)*b(k,j) 
400       continue
405     continue
410   continue

      end
c========================================================================
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine oper2(MM,A,B)
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
      double precision C(0:m_max,0:m_max)
        double precision band1(0:m_max+buff, 2*q_max+1)
        double precision band2(0:m_max+buff, 2*q_max+1)
        double precision band3(0:m_max+buff, 2*q_max+1)
        double precision     H(0:m_max,0:1)
        double precision    CD(0:m_max     ,  1:2)
        double precision    CN(0:m_max     ,  1:2)
        double precision    X1(0:m_max)
      double precision bb(0:m_max)
      double precision cc(0:m_max)
      double precision tau,j2,j3,j4,j5,n2,pmode,ii,mi,a1,a0
c------------------------------------------------------------------------
c....Initializations
c------------------------------------------------------------------------
       pmode      = dble(mode)
       n2         = pmode**2
C-----------------------------------------------------------------------
C--- define Dirichlet & Neumann constraints and grid for FFT
        a1 = 1.d0
          do i = 0,m_max
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
        call koper(2,m_max, mm+1, A)
        write(*,*) 'A'
        call mug(A,m_max,0,10,0,10)
        call mug(A,m_max,mm-10,mm,mm-10,mm)
        call cond1_mat(m_max, 2, 0, band2, 1.d0, debug, buff)
        call bfcon(m_max,mm,-2,2,band2,B,buff)
        write(*,*) 'B'
        call mug(B,m_max,0,10,0,10)
        call mug(B,m_max,mm-10,mm,mm-10,mm)
        call ffmul(m_max+1,MM+1,MM+1,MM+1,A,B,C)
        write(*,*) 'C'
        call mug(C,m_max,0,10,0,10)
        call mug(C,m_max,mm-10,mm,mm-10,mm)
        call cond1_mat(m_max, 0, 0, band2, 1.d0, debug, buff)
        call bfcon(m_max,mm,0,0,band2,B,buff)
        do i = mm-2,0,-1
          do j = 0, mm
             A(i+2,j) = A(i,j)
             B(i+2,j) = B(i,j)
          end do
        end do
        do i = 0,mm                   
            B(0,i) = 0.d0               
            B(1,i) = 0.d0
            A(0,i) = 0.d0               
            A(1,i) = 0.d0
          if (bc_case .eq. 1) then
            A(0,i) = CD(i,1)
            A(1,i) = CD(i,2)
          end if    
          if (bc_case .eq. 2) then
            A(0,i) = CN(i,1)
            A(1,i) = CN(i,2)
          end if
        end do           
        if (bc_case .eq. 0) then
          B(0,0) = 1.d0
          B(1,1) = 1.d0
        end if
        write(*,*) 'A'
        call mug(A,m_max,0,10,0,10)
        call mug(A,m_max,mm-10,mm,mm-10,mm)
        write(*,*) 'B'
        call mug(B,m_max,0,10,0,10)
        call mug(B,m_max,mm-10,mm,mm-10,mm)
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
              subroutine oper3(MM,A,B)
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
          do i = 0,m_max
               X1(i) = cos(pi*i/(2.d0*m_max))
          end do
          call ge_q( 1.d0,nu,0,m_max,CD(0,1))
          call ge_q(-1.d0,nu,0,m_max,CD(0,2))
          call ge_q( 1.d0,nu,1,m_max,CN(0,1))
          call ge_q(-1.d0,nu,1,m_max,CN(0,2))
C-----------------------------------------------------------------------
c....build operator matrices
c-----------------------------------------------------------------------
        call cond1_ge(m_max, 0, 0,nu, band1, 1.d0, debug, buff)
        call cond1_ge(m_max, 2, 0, nu, band2, 1.d0, debug, buff)
        call bfcon(m_max,mm,0,0,band1,A,buff)
        A(0,0) = 0.d0
        A(1,1) = 0.d0
        write(*,*) 'A'
        call mug(A,m_max,0,6,0,6)
c       call mug(A,m_max,mm-10,mm,mm-10,mm)
        call bfcon(m_max,mm,-2,2,band2,B,buff)
        do i = 0,mm                   
          B(0,i) = 0.d0               
          B(1,i) = 0.d0
          if (bc_case .eq. 1) then
          B(0,i) = CD(i,1)
          B(1,i) = CD(i,2)
          end if    
          if (bc_case .eq. 2) then
          B(0,i) = CN(i,1)
          B(1,i) = CN(i,2)
          end if
        end do           
        if (bc_case .eq. 0) then
          B(0,0) = 1.d0
          B(1,1) = 1.d0
        end if
        write(*,*) 'B'
        call mug(B,m_max,0,6,0,6)
c       call mug(B,m_max,mm-10,mm,mm-10,mm)
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
          do i = 0,m_max
               X1(i) = cos(pi*i/(2.d0*m_max))
          end do
          call jac_q( 1.d0, al, be, 0, m_max, CD(0,1))
          call jac_q(-1.d0, al, be, 0, m_max, CD(0,2))
          call jac_q( 1.d0, al, be, 1, m_max, CN(0,1))
          call jac_q(-1.d0, al, be, 1, m_max, CN(0,2))
C-----------------------------------------------------------------------
c....build operator matrices
c-----------------------------------------------------------------------
C--Test: compare Jacobi(0,0) to Gegenbauer(1/2) (i.e. Legendre)
c       call cond1_ge(m_max, 2, 0,nu, band1, 0.d0, debug, buff)
c       call bfcon(m_max,mm,-2,2,band1,A,buff)
c       call cond1_ja(m_max, 1, 0,al,be, band2, 0.d0, debug, buff)
c       call bfcon(m_max,mm,-1,1,band2,B,buff)
        call cond1_ja(m_max, 0, 0,al,be, band1, 0.d0, debug, buff)
        call bfcon(m_max,mm,0,0,band1,A,buff)
        A(0,0) = 0.d0
        A(1,1) = 0.d0
        write(*,*) 'A'
        call mug(A,m_max,0,6,0,6)
c       call mug(A,m_max,mm-10,mm,mm-10,mm)
        call cond1_ja(m_max, 2, 0, al,be, band2, 1.d0, debug, buff)
        call bfcon(m_max,mm,-2,2,band2,B,buff)
        do i = 0,mm                   
            B(0,i) = 0.d0               
            B(1,i) = 0.d0
          if (bc_case .eq. 1) then
            B(0,i) = CD(i,1)
            B(1,i) = CD(i,2)
          end if    
          if (bc_case .eq. 2) then
            B(0,i) = CN(i,1)
            B(1,i) = CN(i,2)
          end if
        end do           
        if (bc_case .eq. 0) then
          B(0,0) = 1.d0
          B(1,1) = 1.d0
        end if
        write(*,*) 'B'
        call mug(B,m_max,0,6,0,6)
c       call mug(B,m_max,mm-10,mm,mm-10,mm)
c       a1 = pi*pi
c       a0 = a1/4.d0
c       do i = 0,mm/2
c        ii = i*i
c        j2 = (2*i+1)*(2*i+1)
c        write(*,*) a0*j2, a1*ii
c       end do
c-----------------------------------------------------------------------
      return
      end
c======================================================================
