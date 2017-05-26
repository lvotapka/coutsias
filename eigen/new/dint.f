      double precision function dint (x)
c december 1983 edition. w. fullerton, c3, los alamos scientific lab.
c
c dint is the double precision equivalent of aint.  this portable
c version is quite efficient when the argument is reasonably small (a
c common case), and so no faster machine-dependent version is needed.
c
      double precision x, xscl, scale, xbig, xmax, part, d1mach,
     1  dlog
      external d1mach, dlog, i1mach, r1mach
      data npart, scale, xbig, xmax / 0, 3*0.0d0 /
c
      if (npart.ne.0) go to 10
      ibase = i1mach(10)
      xmax = 1.0d0/d1mach (4)
      xbig = amin1 (float (i1mach(9)), 1.0/r1mach(4))
      scale = ibase**int(dlog(xbig)/dlog(dble(float(ibase)))-0.5d0)
      npart = dlog(xmax)/dlog(scale) + 1.0d0
c
 10   if (x.lt.(-xbig) .or. x.gt.xbig) go to 20
c
      dint = int(sngl(x))
      return
c
 20   xscl = dabs(x)
      if (xscl.gt.xmax) go to 50
c
      do 30 i=1,npart
        xscl = xscl/scale
 30   continue
c
      dint = 0.0d0
      do 40 i=1,npart
        xscl = xscl*scale
        ipart = xscl
        part = ipart
        xscl = xscl - part
        dint = dint*scale + part
 40   continue
c
      if (x.lt.0.0d0) dint = -dint
      return
c
 50   call seteru (
     1 57hdint dabs(x) may be too big to be repd a an exact integer, 57, 1, 1)
      dint = x
      return
c
      end
