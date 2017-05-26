      subroutine cosfou (m_max,n_max,p_max,data, m, n, isign
     -             ,umc, vmc, wmc, vnc, wnc, copyc, ibitmc, ibitnc)

c....  this subroutine calculates the cosine-fourier transform of
c....  a two-dimensional array of mxn real numbers, data(i,j),
c....  i=1,...,m+1,  j=1,...,n.
      implicit none
      integer m_max,n_max, p_max, m, n, isign, i, j
      integer          ibitmc(m_max)
      integer          ibitnc(n_max)
c     double precision    umc(m_max)
      double precision    umc(m_max-2)
      double precision    vmc(m_max/2+2)
      double precision    wmc(m_max/2)
      double precision    vnc(n_max/2+2)
      double precision    wnc(n_max/2)
      double precision  copyc(p_max+2)
      double precision   data(m_max+1,*)

      if (n .eq. 0) then
        do i=1,m+1
          copyc(i)=data(i,1)
        end do
        call cosft(m,copyc,umc,vmc,wmc,ibitmc,isign)
        do i=1,m+1
          data(i,1)=copyc(i)
        end do
        return
      endif
      if (n .ge. 2) then
      do j=1,n+1
        do i=1,m+1
          copyc(i)=data(i,j)
        end do
        call cosft(m,copyc,umc,vmc,wmc,ibitmc,isign)
        do i=1,m+1
          data(i,j)=copyc(i)
        end do
      end do
      endif
      if (m .gt. 0) then
      do i=1,m+1
        do j=1,n+2
          copyc(j)=data(i,j)
        end do
        call realftc(n,copyc,vnc,wnc,ibitnc,isign)
        do j=1,n+2
          data(i,j)=copyc(j)
        end do
      end do
      end if
      return
      end

c *************************************************

      subroutine cosft(n,y,u,v,w,ibit,isign)

c....   this subroutine calculates the cosine transform of
c....   an array of n+1 real numbers, f_j (j = 0,1,...,n) :
c....
c....                 n            pi
c....        f_j = sum    a_k*cos( -- *j*k)
c....                 k=0          n
c....
c....
c....   variables: (i=input, o=output)
c....   ==========
c....   n     integer          (i)   n+1=number of real numbers
c....   y     real             (i,o) dimension: nmax+2
c....                                array containing the n+1 real numbers
c....   u     real array       (i)   dimension: nmax-2
c....                                trigonometric factors used in cosft
c....   v     real array       (i)   dimension: (nmax/2)+2
c....                                trigonometric factors used in realft
c....   w     real array       (i)   dimension: nmax/2
c....                                trigonometric factors used in cfft
c....   ibit  integer array    (i)   dimension: nmax
c....                                bit-reversal array used in cfft
c....   isign integer          (i)   specifies direction of transformation (see
c....                                below)
c....
c....   isign=1:  forward transformation
c....   ========
c....      input:  y(i)=f_j    ... n+1 real function values
c....      output: y(i)=a_k    ... n+1 cosine coefficients,
c....
c....   isign=-1: invers transformation
c....   =========
c....      input:  y(i)=a_k    ... n+1 cosine coefficients
c....      output: y(i)=f_j    ... n+1 real function values
c....
c....               -------------------------------------------
c....
c....   note : a_0 and a_(n) do not contain any weight factor (see
c....            the above def. of f_j)
c....
c.........................................................................

      implicit none

      integer n, ibit(*), m, j, isign
      real*8   wr,wi,sum,y1,y2,c1
      real*8   y(*), u(*), v(*), w(*)

      if (isign.eq.-1) then
        y(1)=y(1)*2.d0
        y(n+1)=y(n+1)*2.d0
      end if
      sum =0.5d0*(y(1)-y(n+1))
      y(1)=0.5d0*(y(1)+y(n+1))
      m=n/2
      do 11 j=1,m-1
        wr=u(2*(j-1)+1)
        wi=u(2*(j-1)+2)
        y1=0.5d0*(y(j+1)+y(n-j+1))
        y2=(y(j+1)-y(n-j+1))
        y(j+1)  =y1-wi*y2
        y(n-j+1)=y1+wi*y2
        sum=sum+wr*y2
11    continue
      call realftc(n,y,v,w,ibit,1)
      sum=sum*2.d0/n
      y(2)=sum
      do 12 j=4,n,2
        sum=sum+y(j)
        y(j)=sum
12    continue
      if (isign.eq.-1) then
        y(1)=n*y(1)
        y(n+1)=n*y(n+1)
        c1=n/2.d0
        do j=2,n
          y(j)=c1*y(j)
        end do
      end if
      return
      end

c ******************************************************

      subroutine realftc(n,data,v,w,ibit,isign)

c....   this subroutine calculates the fourier transform of
c....   an array of n real numbers, f_j (j = 0,1,...,n-1) :
c....
c....        f_j = a_0 + a_(n/2)*cos(pi*j) +
c....
c....                 (n/2)-1           2*pi                   2*pi
c....              sum       { a_k*cos( ----- *j*k) + b_k*sin( ---- *j*k) }
c....                  k=1                n                      n
c....
c....
c....   variables: (i=input, o=output)
c....   ==========
c....   n     integer          (i)   number of real numbers
c....   data  real             (i,o) dimension: nmax+2
c....                                array containing the n real numbers
c....   v     real array       (i)   dimension: (nmax/2)+2
c....                                trigonometric factors used in realft
c....   w     real array       (i)   dimension: nmax/2
c....                                trigonometric factors used in cfft
c....   ibit  integer array    (i)   dimension: nmax
c....                                bit-reversal array used in cfft
c....   isign integer          (i)   specifies direction of transformation (see
c....                                below)
c....
c....   isign=1:  forward transformation
c....   ========
c....      input:  data(i)=f_j           ... n real function values
c....      output: data(i)=a_k and b_k   ... n real fourier coefficients,
c....                                        organized as:
c....                                         data(1)    = a_0
c....                                         data(2)    = 0
c....                                         data(2k+1) = a_k, k=1 ... (n/2)-1
c....                                         data(2k+2) = b_k, k=1 ... (n/2)-1
c....                                         data(n+1)  = a_(n/2)
c....                                         data(n+2)  = 0
c....
c....   isign=-1: invers transformation
c....   =========
c....      input:  data(i)=a_k and b_k    ... n real fourier coefficients
c....      output: data(i)=f_j            ... n real function values
c....
c....               -------------------------------------------
c....
c....   note : - definition of "forwards" and "inverse"
c....          - the full value of n is passed in the call to realft
c....          - the array data must be dimensioned to nmax+2
c....          - cosine and sine formulation rather than exponential (which
c....              means a change in sign of the sin terms)
c....          - a_0 and a_(n/2) do not contain any weight factor (see
c....              the above def. of f_j)
c....
c.........................................................................

      implicit none

      integer n, isign, ibit(*), np3, i, i1, i2, i3, i4
      real*8   data(*), v(*), w(*)
      real*8   c1, c2, h1r, h1i, h2r, h2i, wr, wi

      c1=0.5d0
      if (isign.eq.1) then
        c2=-0.5d0
        call cfftc(n/2,data,w,ibit,+1)
        data(n+1)=data(1)
        data(n+2)=data(2)
      else
        c2=0.5d0
        data(1)  =data(1)*2.d0
        data(n+1)=data(n+1)*2.d0
      endif
      np3=n+3
      do 11 i=1,n/4+1
        i1=2*i-1
        i2=i1+1
        i3=np3-i2
        i4=i3+1
        h1r= c1*(data(i1)+data(i3))
        h1i= c1*(data(i2)-data(i4))*isign
        h2r=-c2*(data(i2)+data(i4))*isign
        h2i= c2*(data(i1)-data(i3))
        wr=v(2*(i-1)+1)
        wi=v(2*(i-1)+2)*isign
        data(i1)= h1r+wr*h2r+wi*h2i
        data(i2)=isign*(-h1i-wr*h2i+wi*h2r)
        data(i3)= h1r-wr*h2r-wi*h2i
        data(i4)=isign*(h1i-wr*h2i+wi*h2r)
11    continue
      if (isign.eq.1) then
        data(2)=  0.d0
        data(n+2)=0.d0
        data(1)  =data(1)/2.d0
        data(n+1)=data(n+1)/2.d0
      else
        call cfftc(n/2,data,w,ibit,-1)
        data(n+1)=0.d0
        data(n+2)=0.d0
      endif
      return
      end

 
c ******************************************************

	subroutine cfftc(n,d,w,ibit,isign)

c....   this subroutine calculates the complex fourier transform of
c....   an array of n complex numbers, f_j (j = 0,1,...,n-1) :
c....
c....                 n-1               2*pi
c....        f_j = sum     c_k * exp(i* ----- * j*k)
c....                 k=0                 n
c....
c....
c....   variables: (i=input, o=output)
c....   ==========
c....   n     integer          (i)   number of complex numbers
c....   d     complex or real  (i,o) dimension: nmax(=number of real numbers)
c....                                array containing n complex numbers,
c....                                or equivalently 2n real numbers
c....                                designating alternately the
c....                                real and imaginary parts
c....   w     real array       (i)   dimension: nmax/2
c....                                trigonometric factors
c....   ibit  integer array    (i)   dimension: nmax
c....                                bit-reversal array
c....   isign integer          (i)   specifies direction of transformation (see
c....                                below)
c....
c....   isign=1:  forward transformation
c....   ========
c....      input:  d(i)=f_j   ... n complex function values
c....      output: d(i)=c_k   ... n complex fourier coefficients
c....
c....   isign=-1: invers transformation
c....   =========
c....      input:  d(i)=c_k   ... n complex fourier coefficients
c....      output: d(i)=f_j   ... n complex function values
c....
c....               -------------------------------------------
c....
c....   note !!!: the definition of forward and invers transformations
c....               is opposite of the definition in "numerical recipies" !!
c....             the values of c_k have been normalized by division with n,
c....               as they should.
c....
c.........................................................................

        implicit none

	integer n,m,i,j,mmax,k,kd,istep,isign,ibit(*),nn
	real*8   d(*),tempr,tempi,tr,ti,w(*)

        nn=2*n
	do 10 i=1,nn,2
                j=ibit(i)
		if(j.gt.i) then
			tempr = d(j)
			tempi = d(j+1)
			d(j)  = d(i)
			d(j+1)= d(i+1)
			d(i)  = tempr
			d(i+1)= tempi
		endif
   10	continue
	mmax = 2
	kd= n/2
  11	if(nn.gt.mmax) then
		istep = 2*mmax
		do 13 m=1,mmax,2
			k = 1+(m-1)*kd
			do 12 i=m,nn,istep
                                j=i+mmax
				tr=d(j)*w(k)+d(j+1)*w(k+1)*isign
				ti=-d(j)*w(k+1)*isign+d(j+1)*w(k)
				d(j)=d(i)-tr
				d(j+1)=d(i+1)-ti
				d(i)=d(i)+tr
				d(i+1)=d(i+1)+ti
   12			continue
   13		continue
		kd = kd/2
		mmax = istep
		goto 11
	endif

        if (isign.eq.1) then
          do 20 i=1,nn
            d(i) = d(i)/n
20        continue
        end if

	return
	end

 
c *************************************************

      subroutine fftinic(m_max,n_max,p_max,m,n,umc, vmc, wmc, vnc, wnc, 
     _                   copyc, ibitmc, ibitnc)

c....   this subroutine initializes the bit-reversal
c....   and trigonometric arrays used in
c....   cosft, realft, and cfft (and therefore also in cosfou).
c....   the index m refers to cosft, and n refers to realft, which
c....   are the two different transformations in cosfou.

      implicit none

      integer m_max, n_max, p_max, m, n
      integer i, j, mm
      real*8   pi, factor
      parameter(pi=3.1415926535897932384626d0)

      integer          ibitmc(m_max)
      integer          ibitnc(n_max)
c     double precision    umc(m_max)
      double precision    umc(m_max-2)
      double precision    vmc(m_max/2+2)
      double precision    wmc(m_max/2)
      double precision    vnc(n_max/2+2)
      double precision    wnc(n_max/2)
      double precision  copyc(p_max+2)

      j=1
      do i=1,m,2
        ibitmc(i)=j
        mm=m/2
 1      if (mm.ge.2 .and. j.gt.mm) then
          j=j-mm
          mm=mm/2
          go to 1
        end if
        j=j+mm
      end do

      j=1
      do i=1,n,2
        ibitnc(i)=j
        mm=n/2
 2      if (mm.ge.2 .and. j.gt.mm) then
          j=j-mm
          mm=mm/2
          go to 2
        end if
        j=j+mm
      end do

      factor = 2*pi/n
      do i=1, n/4+1
        vnc(2*(i-1)+1) = cos(factor*(i-1))
        vnc(2*(i-1)+2) = sin(factor*(i-1))
      end do

      factor= 2*pi/(n/2)
      do i=1, n/4
        wnc(2*(i-1)+1) = cos(factor*(i-1))
        wnc(2*(i-1)+2) = sin(factor*(i-1))
      end do

      factor = pi/m
      do i=1, m/2-1
        umc(2*(i-1)+1) = cos(factor*i)
        umc(2*(i-1)+2) = sin(factor*i)
      end do

      factor = 2*pi/m
      do i=1, m/4+1
        vmc(2*(i-1)+1) = cos(factor*(i-1))
        vmc(2*(i-1)+2) = sin(factor*(i-1))
      end do

      factor = 2*pi/(m/2)
      do i=1, m/4
        wmc(2*(i-1)+1) = cos(factor*(i-1))
        wmc(2*(i-1)+2) = sin(factor*(i-1))
      end do

      end

