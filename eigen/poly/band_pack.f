C-----------------------------------------------------------------------
C           ----------<   B  B  A  D  D  .  f    >----------
C-----------------------------------------------------------------------
        subroutine bbadd(mmax,m,al,bl,cl,au,bu,cu,A,B,C,alpha,beta,buff)
C-----------------------------------------------------------------------
C Add two banded matrices A,B and put result in C. Compute to make
C                 C = alpha*A + beta*B
C sure that cu/l >= max(au/l,bu/l) else C has wrong bandwidth
C Normally caller is responsible for supplying correct bandwidth specs.
C This operation does not deteriorate correctness in given truncation.
C A can be overwritten by C
      implicit none
      integer m,mmax,al,bl,cl,au,bu,cu
      integer mm,i,k,buff, ml, mu
      double precision alpha,beta
      double precision A(0:mmax+buff,al:au)
      double precision B(0:mmax+buff,bl:bu)
      double precision C(0:mmax+buff,cl:cu)
C-----------------------------------------------------------------------
        ml = min(al,bl) 
        mu = max(au,bu) 
        mm = m+buff
      if (cl .le. ml .and. cu .ge. mu) then
        do i = 0, mm
C-- zero-out C if it is different from A
          if (cl .ne. al .or. cu .ne. au) then
            do k = cl,cu
              C(i,k) = 0.d0
            end do
          end if
C-- accumulate A (or, sibmply, multiply by a number alpha)
          do k = al,au
            C(i,k) = alpha*A(i,k)
          end do
C-- add in B's multiple by beta
          do k = bl,bu
            C(i,k) = C(i,k)+beta*B(i,k)
          end do
        end do
      else
        write(*,*) 'in BBADD: attempt to add into inadequate bands'
        STOP
      end if
C-----------------------------------------------------------------------
	return
	end
C-----------------------------------------------------------------------
C           ----------<   B  B  M  U  L  .  f    >----------
C-----------------------------------------------------------------------
          subroutine bbmul(mmax,m,al,bl,cl,au,bu,cu,A,B,C,buff)
C-----------------------------------------------------------------------
C Multiply two banded matrices A,B and put result in C. Compute to make
C sure that au/l+bu/l=cu/l else C has wrong bandwidth
C ***>>> If all are saved by diagonal (i.e. the first mm rows are correct),
C        then after multiplication the last au rows are contaminated: hence
C        don't compute them.
C        upper au rows return with zeroes; must pad banded matrices with enough
C        extra correct entries so product is correct to desired truncation,
C        which means the padding depends on how many iteration levels a certain
C        matrix is going to seed.
C buff   : a global variable characterizing the code. Set as parameter,
C          since it determines true dimension of all banded matrices in code.
C-----------------------------------------------------------------------
      implicit none
      integer m,mmax,al,bl,cl,au,bu,cu
      integer i,j,k,mm,buff
      double precision A(0:mmax+buff,al:au)
      double precision B(0:mmax+buff,bl:bu)
      double precision C(0:mmax+buff,cl:cu)
      mm = m+buff
c- zero-out C
      do j = cl,cu
        do i = 0, mm
          C(i,j) = 0.d0
        end do
      end do
      do j = al, au
        do k = bl,bu
          do i = max(0,-j,-j-k),min(mm,mm-au)
C-as matrices C(i,i+j+k)=C(i,i+j+k)+A(i,i+j)*B(i+j,i+j+k); transl. to banded:
            C(i,j+k) = C(i,j+k)+A(i,j)*B(i+j,k)
          end do
        end do
      end do
C-----------------------------------------------------------------------
	return
	end
C-----------------------------------------------------------------------
C           ----------<   B  F  C  O  N  .  f    >----------
C-----------------------------------------------------------------------
          subroutine bfcon(mmax,m,al,au,A,B,buff)
C-----------------------------------------------------------------------
C al, au: lower and upper bands of banded matrix A
C only one page of an operator is copied per call 
      implicit none  
      integer buff,mmax,m,al,au
      integer i,k,l,ll,lu
      double precision A(0:mmax+buff,al:au)
      double precision B(0:mmax,0:mmax)
C-----------------------------------------------------------------------
      do k = 0,m   
        do i = 0,m
          B(i,k) = 0.d0
        end do
      end do
      do i = 0,m
        ll = max(-i,al)
        lu = min(m-i,au)
        do l = ll,lu
          B(i,i+l) = A(i,l)
        end do
      end do
      return              
      end
C-----------------------------------------------------------------------
C           ----------<   B  F  M  U  L  .  f    >----------
C-----------------------------------------------------------------------
          subroutine bfmul(mmax,m,n,al,au,A,B,C,bb,cc,buff)
C-----------------------------------------------------------------------
C Multiply banded matrix A times full matrix B and put result in C. 
C Can have B = C so that second input may or may not be overwritten
C Use, e.g. for various preconditioning operations.
C al, au: lower and upper bands of banded matrix A
C Note that only one page of an operator is done per call.
      implicit none
      integer buff,mmax,m,n,al,au
      integer i,j,k,kk
      double precision A(0:mmax+buff,al:au)
      double precision B(0:mmax,0:1,0:*)
      double precision C(0:mmax,0:1,0:*)
      double precision bb(0:mmax)
      double precision cc(0:mmax)
C-----------------------------------------------------------------------
      do j = 0,n
        do k = 0,1
          do i = 0,m
            bb(i) = B(i,k,j)
            cc(i) = 0.d0
          end do
          if (al .le. 0 .and. au .ge. 0) then
            do kk = 0,au
              do i = 0, m-kk
                cc(i) = cc(i)+A(i,kk)*bb(i+kk)
              end do
            end do
            do kk = al,-1
              do i = -kk, m
                cc(i) = cc(i)+A(i,kk)*bb(i+kk)
              end do
            end do
          else
            if (al .ge. 0) then
            do kk = al,au
              do i = 0, m-kk
                cc(i) = cc(i)+A(i,kk)*bb(i+kk)
              end do
            end do
            end if
            if (au .le. 0) then
            do kk = al,au
              do i = -kk, m
                cc(i) = cc(i)+A(i,kk)*bb(i+kk)
              end do
            end do
            end if
          end if
          do i = 0,m
            C(i,k,j) = cc(i)
          end do
        end do
      end do
C-----------------------------------------------------------------------
	return
	end
C-----------------------------------------------------------------------
C           ----------<   O  P  M  U  L  .  f    >----------
C-----------------------------------------------------------------------
          subroutine opmul(mmax,m,n,al,au,A,B,C,bb,cc,buff)
C-----------------------------------------------------------------------
C Multiply banded operator A times full matrix B and put result in C. 
C Can have B = C so that second input may or may not be overwritten
C Use, e.g. for various preconditioning operations.
C al, au: lower and upper bands of banded matrix A
C Note that only one page of an operator is done per call.
      implicit none
      integer buff,mmax,m,n,al,au
      integer i,j,k,kk
      double precision A(0:mmax+buff,al:au,0:*)
      double precision B(0:mmax,0:1,0:*)
      double precision C(0:mmax,0:1,0:*)
      double precision bb(0:mmax)
      double precision cc(0:mmax)
C-----------------------------------------------------------------------
      do j = 0,n
        do k = 0,1
          do i = 0,m
            bb(i) = B(i,k,j)
            cc(i) = 0.d0
          end do
          if (al .le. 0 .and. au .ge. 0) then
            do kk = 0,au
              do i = 0, m-kk
                cc(i) = cc(i)+A(i,kk,j)*bb(i+kk)
              end do
            end do
            do kk = al,-1
              do i = -kk, m
                cc(i) = cc(i)+A(i,kk,j)*bb(i+kk)
              end do
            end do
          else
            if (al .ge. 0) then
            do kk = al,au
              do i = 0, m-kk
                cc(i) = cc(i)+A(i,kk,j)*bb(i+kk)
              end do
            end do
            end if
            if (au .le. 0) then
            do kk = al,au
              do i = -kk, m
                cc(i) = cc(i)+A(i,kk,j)*bb(i+kk)
              end do
            end do
            end if
          end if
          do i = 0,m
            C(i,k,j) = cc(i)
          end do
        end do
      end do
C-----------------------------------------------------------------------
	return
	end
C-----------------------------------------------------------------------
