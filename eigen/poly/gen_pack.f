C-----------------------------------------------------------------------
C                 -------<  c o n d 1_P O L . f  >-------
C-----------------------------------------------------------------------
      subroutine cond1(maxm,pp,kk,OP,aspect,debug,buff,X,B)
C-----------------------------------------------------------------------
C  initialize various banded preconditioners for polynomial family
C       will characterize by giving recurrence functions for X and B.
C       these are computed in POLINI, the routine that initializes
C       arrays for the polynomial transform.
C-----------------------------------------------------------------------
C-- Build type 1 (integration) preconditioner matrix B_[2]^pp (X+a)^kk
C   with bands -(pp+kk):(pp+kk); independent of polynomial family 
C   All bands are introduced, even those that may contain zeroes.
C   The computation is carried over enough additional modes so that 
C   truncation at maxm is exact.
C-----------------------------------------------------------------------
      implicit none
      integer maxm,pp,kk,debug,i,j,k,mmax,buff
      double precision aspect
      double precision OP(0:maxm+buff,-pp-kk:pp+kk)
      double precision  X(0:maxm+buff,    -1:1    )
      double precision  B(0:maxm+buff,    -1:1    )
      double precision a1(-16:16)
      double precision a2(-16:16)
C-----------------------------------------------------------------------
C-- Zero-out arrays
      do i = 0,maxm+buff
        do j = -pp-kk, pp+kk
          OP(i,j) = 0.d0
        end do
      end do
      do i = -16,16
        a1(i) = 0.d0
        a2(i) = 0.d0
      end do
C-----------------------------------------------------------------------
C- the default for op is I_[2], the identity with first two rows set to zero
C  but this is left to the programmer; this routine produces arrays correctly
C  with no special entries.
      do i = 0,maxm+buff
        OP(i,0) = 1.d0
      end do
C-------------------
C - extend the computation adequately to ensure that truncation at maxm
C   produces exact bands. Thus, treat last rows same as the rest;
C   must extend if higher powers/order are to be used.
c     mmax = maxm+buff
      mmax = maxm
C-------------------
C-- First build (x+a)^kk; bypass if kk = 0
      if (kk .gt. 0) then
        do k = 1, kk
c- zero help arrays
          do j = -k,k
            a1(j) = 0.d0
            a2(j) = 0.d0
          end do
c-   0-row
          i = 0
          do j = 0, k
            a2(j) = (aspect+X(i, 0))*op(i  ,j  )
     _             +        X(i, 1) *op(i+1,j-1)
          end do
C-  1-row
          i = 1
          do j = -1, k-1
            a1(j) = (aspect+X(i, 0))*op(i  ,j  )
     _             +        X(i, 1) *op(i+1,j-1)
     _             +        X(i,-1) *op(i-1,j+1)
          end do
            a1(k) =         X(i, 1) *op(i+1,k-1)
          do j = 0,k
            op(0,j) = a2(j)
          end do
          do j = -1,k
            a2(j) = a1(j)
          end do
C-  all other rows
C- first compute row entries and store; old contents still needed for next row
          do i = 2, mmax+buff-k
            a1(-k  ) =      X(i,-1) *op(i-1,-k+1)
            do j = -k+1, k-1
            a1(j) = (aspect+X(i, 0))*op(i  ,j  )
     _             +        X(i, 1) *op(i+1,j-1)
     _             +        X(i,-1) *op(i-1,j+1)
            end do        
            a1(k) =         X(i, 1) *op(i+1,k-1)
C- Now, previous row no longer needed; replace with stored contents
C  from a2; then put a1 into a2, freeing it for next cycle
            do j = -k, k
              op(i-1,j) = a2(j)
              a2(j)     = a1(j)
            end do
          end do        
c-  last row
            do j = -k, k
              op(mmax+buff-k  ,j) = a2(j)
            end do
        end do        
      end if
C-----------------------------------------------------------------------
C-- Now multiply pp times by B
      do k = kk+1, kk+pp
c- zero help arrays
        do j = -k,k
          a2(j) = 0.d0
        end do
        do i = 1, mmax+buff-k
            a1(-k) = op(i-1,-k+1)*B(i,-1)
          do j = -k+1,k-1
            a1( j) = op(i-1, j+1)*B(i,-1)
     _              +op(i  , j  )*B(i, 0)
     _              +op(i+1, j-1)*B(i, 1)
          end do
            a1( k) = op(i+1, k-1)*B(i, 1)
          do j = -k,k
            op(i-1,j) = a2(j) 
            a2(j)     = a1(j)
          end do
        end do
          do j = -k,k
            op(mmax+buff-k,j) = a2(j) 
          end do
      end do
C- Zero out first pp rows, regardless of content
C-- Characteristic of B_[2]^2 preconditioning. 
C-- here it is done to avoid any special entries in the B
      do j = -kk-pp,kk+pp
        do i = 0,pp-1
          op(i,j) = 0.d0
        end do
      end do
      return
      end
C-----------------------------------------------------------------------
