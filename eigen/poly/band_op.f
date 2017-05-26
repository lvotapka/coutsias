C-----------------------------------------------------------------------
C INITIALIZE ARRAYS FOR 2*p+1 BANDED SOLVER: UL factorization,
C                                            Constraints,
C   1st or 2nd order diff. operator          Homogeneous solutions, 
C       B^2_[2](x+a)^2 preconditioning       adjustment operators.
C
C      ***    P^(p)_[2]+B_[2] P^(p-1)+B^2_[2] P^(p-2)   ***
C
C  P^(p-k) = alpha(k,p-k)T_(p-k)+alpha(k,p-k-1)T_(p-k-1)+...+alpha(k,0)T_0
C
C  Caller must precalculate the polynomial expansions of the P^(k) (whose
C  index gives the max. allowable degree for each) and provide them in the
C  argument array alpha(0:pn,0:p).
C  In this geometry, r=(X+aspect) and for a>1 (annulus) or a=1 (disk)
C  all expressions are derived in the same way. However for the magic
C  combinations (x+1)^2 D^2 and (x+1) D an entirely different preconditioner 
C  is available. (This may be related to some connection between even or odd
C  subclass and the interval (0,1)?)
C   This is not in the pattern of the present banded solver, which
C  assumes B^2_[2] preconditioning.
C-----------------------------------------------------------------------
      subroutine band_op(mmax,m,q,n,P,OP,pn,CP,VP,HP,CHP,debug,buff1,
     _                   U,bb,cc)
C-----------------------------------------------------------------------
C Arguments
C           mmax        : leading column dimension of arrays (0:m_max)
C                
C           m,n         : effective  dimensions are 0:m,*,0:n
C           q           : (symmetric) bandwidth of operator 
C           P           : The (eventually factored) operator in band form
C           OP          : The original operator stored in band form
C           pn          : number of top null-rows for operator (order)
C           CP          : The two constraints, assumed identical for sine
C                           and cosine modes (1: outer, 2: inner)
C                         This aspect, peculiar to Chebyshev expansions,
C                         must be dealt with for the general case.
C                                                   0 : 1
C           VP          : The "BC" values, (out:in,cos:sin,0:n/2)
C           HP          : The homos: (0:m,hom_0:hom_1,0:n/2)
C           CHP         : CHP(i,j) = <hom_i , CP_j> (inverted)
C  ****     alpha(k,p-k): the (supplied) array of expansion coeffs. for the
C     (not implemented)     polynomial coefficients P^(p-k),k=0,1,2
C-----------------------------------------------------------------------
	implicit none
        integer mmax,m,q,n,pn,debug,sing
        integer buff1
        double precision     HP(  0:mmax  ,  0:1       , 0:*)
        double precision      U(  0:mmax  ,  0:1       , 0:*)
        double precision      P(  0:mmax+buff1 , -q:  q     , 0:*)
        double precision     OP(  0:mmax+buff1 ,  1:2*q+1   , 0:*)
        double precision     CP(  0:mmax  ,  1:2       , 0:*)
        double precision     VP(  1:2     ,  0:1       , 0:*)
        double precision    CHP(  1:2     ,  0:1       , 0:*)
        double precision     bb(  0:*), cc(  0:*)
	integer i,j,k,kk,l,cond, case
        double precision det, u0
C-----------------------------------------------------------------------
C -- Define the RHS for the homogeneous solution: the problem is
C    nonsingular and Q_0^k (k=pn-1) has a rank pn component in Null space. 
C    Solve in Q_pn^M by assuming, in turn, eg. HP = (1,0,..,0) + f E Q_2^M
        do j = 0,n
          do k = -q,q
            do i = 0,pn-1
C-- zero out top PN rows of operator (safety measure!)
              P(i,k,j) = 0.d0
            end do
            do i = 0,m
              OP(i,k+q+1,j) = P(i,k,j)
            end do
          end do
          do k = 0,pn-1
c$        do k = 0,1
            do i = 0, m
              HP(i,k,j) = 0.d0
            end do
            do i = pn,q+k
c$          do i = 2,q+k
              HP(i,k,j) = -P(i,-i+k,j)
            end do
          end do
        end do
        if (debug .gt. 2) then
          write(*,*) ' P = '
          call mug(P,mmax+buff1,0,debug,0,2*q)
        end if
C- The UL factorization;
c                           -----  U : P(i, 0.. q,j), unit diagonal  ---
c                           -----  L : P(i,-1..-q,j)                 ---
	do j = 0, n
	  do i = m, pn+1, -1
c$ do i = m, 3, -1
            kk = min(q,i-pn)
c$          kk = min(q,i-2)
            do k = 1,kk
c--U-factors
                P(i-k,k,j) = P(i-k,k,j)/P(i,0,j)
              do l = 1, min(q,i)
c-- new L-factors
                P(i-k,k-l,j) = P(i-k,k-l,j)-P(i-k,k,j)*P(i,-l,j)
              end do
            end do
          end do
	end do
      if (debug .gt. 2) then
        do j = 0,n
          call bbmul(mmax,m,1,-q,-q,q,0,q,P(0,1,j),P(0,-q,j)
     _                ,OP(0,1,j),buff1)
          call bbadd(mmax,m,-q,-q,-q,q,0,q,OP(0,1,j),P(0,-q,j)
     _                ,OP(0,1,j),1.d0,1.d0,buff1)
          write(*,*) ' test UL factorization: U*L =? P, mode = ',j
          call mug(OP,mmax+buff1,0,debug,0,2*q)
        end do
          write(*,*) ' P in band_op after U-L (mode 0)'
          call mug(P,mmax+buff1,0,debug,0,2*q)
C-----------------------------------------------------------------------
C-- The homogeneous solutions
C-- Compute homogeneous solutions in Q_2^M
          write(*,*) ' HP before homos: first ',q,' entries of P '
          call mug(HP,mmax,0,debug,0,2*n+1)
      end if
      call ul_sol(mmax,m,q,n,P,OP,pn,CP,VP,HP,CHP,HP,cond,0,0,debug,
     _            buff1,U,bb,cc)
c-- adjust 0-modes
        do j = 0,n
          do i = 0,pn-1
            do k = 0,pn-1
              HP(i,k,j)=0.d0
            end do
              HP(i,i,j)=1.d0
          end do
        end do
      if (debug .gt. 2) then
          write(*,*) ' HP adjusted '
          call mug(HP,mmax,0,debug,0,2*n+1)
          write(*,*) 'after adjustment Op*HP =? 0 '
      call opmul(mmax,m,n,-q,q,OP,HP,U,bb,cc,buff1)
          call mug(U,mmax,0,debug,0,2*n+1)
      end if
C-- Compute inner products of constraints with homo. solutions
        do j=0, n
          do k=1, pn
            do l=0, pn-1
              CHP(k,l,j) = CP(0,k,j)*HP(0,l,j) 
            end do
          end do
          do i = 1, m
            do k=1, pn
              do l=0, pn-1
                CHP(k,l,j) = CHP(k,l,j)+CP(i,k,j)*HP(i,l,j)
              end do
            end do
          end do
            if (debug .gt. 2) then
              write(*,*) 'BC matrix: '
              write(*,*) '1,0 = ',CHP(1,0,j)
          if (pn .eq. 2) then
              write(*,*) '2,0 = ',CHP(2,0,j)
              write(*,*) '1,1 = ',CHP(1,1,j)
              write(*,*) '2,1 = ',CHP(2,1,j)
          endif
            end if
C-- and invert to get 2X2Xj constraint operator
          if (pn .eq. 1) then
            CHP(1,0,j) =  1.d0/CHP(1,0,j)
          else
            det = CHP(1,0,j)*CHP(2,1,j)-CHP(2,0,j)*CHP(1,1,j)
            if (abs(det) .ge. .1d-20 ) then
              CHP(2,0,j) = -CHP(2,0,j)/det
              CHP(1,1,j) = -CHP(1,1,j)/det
              u0         =  CHP(1,0,j)/det
              CHP(1,0,j) =  CHP(2,1,j)/det
              CHP(2,1,j) = u0
              if (debug .gt. 2) then
                write(*,*) ' j = ',j,'  det  =  ',det
                write(*,*) 'Inverse: '
                write(*,*) '1,0 = ',CHP(1,0,j)
                write(*,*) '2,0 = ',CHP(2,0,j)
                write(*,*) '1,1 = ',CHP(1,1,j)
                write(*,*) '2,1 = ',CHP(2,1,j)
              end if
            else
              u0 = CHP(1,1,j)**2+CHP(2,1,j)**2
              CHP(2,0,j) =  CHP(1,1,j)/u0
              CHP(2,1,j) =  CHP(2,1,j)/u0
              CHP(1,1,j) =  0.d0
              CHP(1,0,j) =  0.d0
              if (debug .gt. 2) then
                write(*,*) ' j = ',j,'  det  =  ',det
                write(*,*) 'Pseudo-Inverse: '
                write(*,*) '1,0 = ',CHP(1,0,j)
                write(*,*) '2,0 = ',CHP(2,0,j)
                write(*,*) '1,1 = ',CHP(1,1,j)
                write(*,*) '2,1 = ',CHP(2,1,j)
              end if
            end if
          end if
        end do
        return
        end
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C SOLVE BANDED SYSTEM in U*L factored form: 
c                           -----  U : P(i, 1.. q,j), unit diagonal  ---
c                           -----  L : P(i, 0..-q,j)                 ---
C-----------------------------------------------------------------------
      subroutine ul_sol(mmax,m1,q,n1,P,OP,pn,CP,VP,HP,CHP,U,cond,sing,
     _                  case,debug,buff1,U0,bb,cc)
C-----------------------------------------------------------------------
C        m_max: max column dimension
C           m,n         : effective  dimensions are 0:m,*,0:n
C           q           : bandwidth of operator 
C           U           : input-output array
C           P           : The operator stored in band form
C           pn          : order of operator and number of constraints
C-----------------------------------------------------------------------
        implicit none
C-----------------------------------------------------------------------
        integer buff1, q, pn
        integer m1, n1, k, case, i, j, l
        double precision    U(  0:mmax  ,  0:1  , 0:*)
        double precision   U0(  0:mmax  ,  0:1  , 0:*)
        double precision   HP(  0:mmax  ,  0:1  , 0:*)
        double precision    P(  0:mmax+buff1  , -q:  q  , 0:*)
        double precision  OP(  0:mmax+buff1  ,  1:2*q+1, 0:*)
        double precision   CP(  0:mmax  ,  1:2  , 0:*)
        double precision   VP(  1:2     ,  0:1  , 0:*)
        double precision  CHP(  1:2     ,  0:1  , 0:*)
        double precision   bb(0:mmax), cc(0:mmax)
        double precision  u1, u2, u4, u3, tau
	integer           debug, cond,sing,mmax,ll
C-----------------------------------------------------------------------
C NOTE: postconditioned triangular factors must be adjusted! Still to do...
C-----------------------------------------------------------------------
C--                                   U z =  f
C--- Inefficient programming, with compares in do loop arguments. FIX!!!
        if (debug .gt. 2) then
          write(*,*) ' The BC matrix C',2*q+1
          call mug(CP,mmax,0,debug,0,2*n1+1)
          write(*,*) ' f before Uz = f '
          call mug(U,mmax,0,debug,0,2*n1+1)
        end if
        do j = 0, n1
          do k = 0,pn-1
c$        do k = 0,1
            do i = m1,pn,-1
c$          do i = m1,2,-1
              ll = min(q,m1-i)
              do l = 1, ll
                U(i,k,j) = U(i,k,j)-P(i,l,j)*U(i+l ,k,j)
              end do
            end do
              do i = 0,pn-1
                U(i,k,j) = 0.d0
c$              U(0,k,j) = 0.d0
c$              U(1,k,j) = 0.d0
              end do
          end do
        end do
        if (debug .gt. 2) then
      do j = 0, n1
      call bfmul(mmax,m1,0,1,q,P(0,1,j),U(0,0,j),U0(0,0,j)
     _                                          ,bb,cc,buff1)
      end do
      call addm(1.d0,U0,1.d0,U,U0,mmax,m1,2*n1+1)
          write(*,*) 'check correctness of Uz =? f, mode = ',j
          call mug(U0,mmax,0,debug,0,2*n1+1)
          write(*,*) ' z after Uz = f '
          call mug(U,mmax,0,debug,0,2*n1+1)
        end if
C--                                   L x = z
        do j = 0, n1
          do k = 0,pn-1
C$        do k = 0,1
            do i = pn,m1
C$          do i = 2,m1
              ll = min(q,i)
              do l = 1,ll
                U(i,k,j)=U(i,k,j)-P(i,-l,j)*U(i-l,k,j)
              end do
                U(i,k,j)=U(i,k,j)/P(i,0,j)
            end do
          end do
        end do
        if (debug .gt. 2) then
      do j = 0, n1
      call bfmul(mmax,m1,0,-q,0,P(0,-q,j),U(0,0,j),U0(0,0,j)
     _                                          ,bb,cc,buff1)
      end do
          write(*,*) 'check correctness of Lx =? z'
          call mug(U0,mmax,0,debug,0,2*n1+1)
          write(*,*) ' x after Px = z '
          call mug(U,mmax,0,debug,0,2*n1+1)
        if (debug .gt. 2) then
          call mug(U,mmax,m1-2,m1+1,0,2*n1+1)
        end if
      call opmul(mmax,m1,n1,-q,q,OP,U,U0,bb,cc,buff1)
          write(*,*) 'check correctness of ULx =? f'
          call mug(U0,mmax,0,debug,0,2*n1+1)
        end if
C-----------------------------------------------------------------------
C-- Adjustment for constraints: sing =0, no singularity
C-- must not do when computing homo. solutions
      if (case .eq. 0) then
        return
      end if
      if (sing .eq. 0) then
        do j = 0, n1
          do k = 0,pn-1
c$        do k = 0,1
              u1 = VP(1,k,j)
            if (pn .eq. 2) then
              u2 = VP(2,k,j)
            end if
            do i = 0,m1
              u1 = u1-CP(i,1,j)*U(i,k,j)
              u2 = 0.d0
            if (pn .eq. 2) then
              u2 = u2-CP(i,2,j)*U(i,k,j)
            end if
            end do
              u3 = CHP(1,0,j)*u1+CHP(1,1,j)*u2
              u4 = 0.d0
            if (pn .eq. 2) then
              u4 = CHP(2,0,j)*u1+CHP(2,1,j)*u2
            end if
c           write(*,*) u3, u4
            do i = 0,m1
              U(i,k,j)=U(i,k,j)+u3*HP(i,0,j)+u4*HP(i,1,j)
            end do
          end do
        end do
        if (debug .gt. 2) then
          write(*,*) 'check correctness of final solution incl. BC.'
          call opmul(mmax,m1,n1,-q,q,OP,U,U0,bb,cc,buff1)
          call mug(U0,mmax,0,debug,0,2*n1+1)
          write(*,*) 'final solution incl. BC.'
          call mug(U,mmax,0,debug,0,2*n1+1)
          write(*,*) ' check conformity with BC'
          u1 = 0.d0
          u2 = 0.d0
          u3 = 0.d0
          u4 = 0.d0
          tau = 1.d0
          do i = 0,m1
           u1 = u1 + CP(i,1,0)*U(i,0,0) 
           u2 = u2 + U(i,0,0) 
           u3 = u3 + CP(i,2,0)*U(i,0,0) 
           u4 = u4 + tau*U(i,0,0) 
           tau = -tau
          end do
           write(*,*) 'outer BC ', VP(1,0,0),u1
           write(*,*) 'inner BC ', VP(2,0,0),u3
           write(*,*) ' Boundary values (out-in) ', u2, u4
        end if
      end if
      return
      end
C----------------------------------------------------------------------
