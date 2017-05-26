      subroutine d9gaml (xmin, xmax)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c calculate the minimum and maximum legal bounds for x in gamma(x).
c xmin and xmax are not the only bounds, but they are the only non-
c trivial ones to calculate.
c
c             output arguments --
c xmin   dble prec minimum legal value of x in gamma(x).  any smaller
c        value of x might result in underflow.
c xmax   dble prec maximum legal value of x in gamma(x).  any larger
c        value of x might cause overflow.
c
      double precision xmin, xmax, alnbig, alnsml, xln, xold, d1mach,
     1  dlog
      external d1mach, dlog
c
      alnsml = dlog(d1mach(1))
      xmin = -alnsml
      do 10 i=1,10
        xold = xmin
        xln = dlog(xmin)
        xmin = xmin - xmin*((xmin+0.5d0)*xln - xmin - 0.2258d0 + alnsml)
     1    / (xmin*xln+0.5d0)
        if (dabs(xmin-xold).lt.0.005d0) go to 20
 10   continue
      call seteru (27hd9gaml  unable to find xmin, 27, 1, 2)
c
 20   xmin = -xmin + 0.01d0
c
      alnbig = dlog (d1mach(2))
      xmax = alnbig
      do 30 i=1,10
        xold = xmax
        xln = dlog(xmax)
        xmax = xmax - xmax*((xmax-0.5d0)*xln - xmax + 0.9189d0 - alnbig)
     1    / (xmax*xln-0.5d0)
        if (dabs(xmax-xold).lt.0.005d0) go to 40
 30   continue
      call seteru (27hd9gaml  unable to find xmax, 27, 2, 2)
c
 40   xmax = xmax - 0.01d0
      xmin = dmax1 (xmin, -xmax+1.d0)
c
      return
      end
