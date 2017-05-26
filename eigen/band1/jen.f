C-----------------------------------------------------------------------
C                 -------<  c o n d_P O L . f  >-------
C-----------------------------------------------------------------------
	subroutine cond_pol(maxm,pp,kk,kind,al,be,OP,aspect,debug,buff)
C-----------------------------------------------------------------------
C  initialize various banded pre/post-conditioners for polynomial family
C     KIND 1: Legendre                                  (-: post)
C     KIND 2: Chebyshev
C     KIND 3: Chebyshev 2nd kind
C     KIND 4: Hermite 
C     KIND 5: Jacobi 
C     KIND 6: Laguerre
C     KIND 7: Gegenbauer
C       will characterize by giving recurrence functions for X and B.
C       these are computed in POLINI, the routine that initializes
C       arrays for the polynomial transform.
C-----------------------------------------------------------------------
C-- Build type 1 (integration) preconditioner matrix G_[2]^pp (x+a)^kk
C   with bands -(pp+kk):(pp+kk) 
C   All bands are introduced, even those that may contain zeroes.
C   The computation is carried over enough additional modes so that 
C   truncation at maxm is exact.
C-----------------------------------------------------------------------
      implicit none
      integer maxm,pp,kk,debug,i,j,k,mmax,buff, kind
      double precision al,be, aspect, y0,y1,y2,y3,y4
      double precision nu
      double precision OP(0:maxm+buff,-pp-kk:pp+kk)
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
               y0 = al*al-be*be
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
               y1 = 2.d0*i+al+be
               y2 =      i+al+be
               y3 =      i+al
               y4 =      i   +be
            if (y1 .eq. 0.d0) then
            a2(j) =  aspect                             *op(i  ,j  )
     _             +(2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,j-1)
            else
            a2(j) = (aspect-y0/(y1*(y1+2)))             *op(i  ,j  )
     _             +(2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,j-1)
            endif
          end do
C-  1-row
          i = 1
               y1 = 2.d0*i+al+be
               y2 =      i+al+be
               y3 =      i+al
               y4 =      i   +be
          do j = -1, k-1
            a1(j) = (aspect-y0/(y1*(y1+2)))             *op(i  ,j  )
     _             +(2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,j-1)
     _             +(2.d0* i    * y2   /((y1-1)* y1   ))*op(i-1,j+1)
          end do
            a1(k) = (2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,k-1)
          do j = 0,k
            op(0,j) = a2(j)
          end do
          do j = -1,k
            a2(j) = a1(j)
          end do
C-  all other rows
C- first compute row entries and store; old contents still needed for next row
          do i = 2, mmax+buff-k
               y1 = 2.d0*i+al+be
               y2 =      i+al+be
               y3 =      i+al
               y4 =      i   +be
            a1(-k  ) = (2.d0* i    * y2   /((y1-1)* y1   ))*op(i-1,-k+1)
            do j = -k+1, k-1
               a1(j) = (aspect-y0/(y1*(y1+2)))             *op(i  ,j  )
     _                +(2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,j-1)
     _                +(2.d0* i    * y2   /((y1-1)* y1   ))*op(i-1,j+1)
            end do        
            a1(k) = (2.d0*(y3+1)*(y4+1)/((y1+2)*(y1+3)))*op(i+1,k-1)
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
               y0 = al-be
      do k = 1, pp
c- zero help arrays
        do j = -kk-k,kk+k
          a1(j) = 0.d0
          a2(j) = 0.d0
        end do
        i = 1
               y1 = 2.d0*i+al+be
               y2 =      i+al+be
               y3 =      i+al   +1
               y4 =      i   +be+1
            a2(-kk-k) = op(i-1,-kk-k+1)*(2*y2/(y1*(y1-1.d0)))
          do j = -kk-k+1,kk+k-1
            a2( j   ) = op(i-1, j   +1)*(2*y2/(y1*(y1-1.d0)))
     _                 +op(i  , j     )*(2*y0/((y1+2)*y1))
     _                 -op(i+1, j   -1)*(2*y3*y4/((y2+1)*(y1+2)*(y1+3)))
          end do
            a2( kk+k) =-op(i+1, kk+k-1)*(2*y3*y4/((y2+1)*(y1+2)*(y1+3)))
        do i = 2, mmax+buff-kk-k
               y1 = 2.d0*i+al+be
               y2 =      i+al+be
               y3 =      i+al   +1
               y4 =      i   +be+1
            a1(-kk-k) = op(i-1,-kk-k+1)*(2*y2/(y1*(y1-1.d0)))
          do j = -kk-k+1,kk+k-1
            a1( j   ) = op(i-1, j   +1)*(2*y2/(y1*(y1-1.d0)))
     _                 +op(i  , j     )*(2*y0/((y1+2)*y1))
     _                 -op(i+1, j   -1)*(2*y3*y4/((y2+1)*(y1+2)*(y1+3)))
          end do
            a1( kk+k) =-op(i+1, kk+k-1)*(2*y3*y4/((y2+1)*(y1+2)*(y1+3)))
          do j = -kk-k,kk+k
            op(i-1,j) = a2(j) 
            a2(j)     = a1(j)
          end do
        end do
          do j = -kk-k,kk+k
            op(mmax+buff-kk-k,j) = a2(j) 
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

c       do i = 0,maxm
c           write(*,*) i, (op(i,j), j=-kk-pp,kk+pp)
c       end do

      return
      end
C-----------------------------------------------------------------------
