      double precision function d9lgmc (x)
c august 1977 edition.  w. fullerton, c3, los alamos scientific lab.
c
c compute the log gamma correction factor for x .ge. 10. so that
c dlog (dgamma(x)) = dlog(dsqrt(2*pi)) + (x-.5)*dlog(x) - x + d9lgmc(x)
c
      double precision x, algmcs(15), xbig, xmax, dcsevl, d1mach,
     1  dexp, dlog, dsqrt
      external d1mach, dcsevl, dexp, dlog, dsqrt, initds
c
c series for algm       on the interval  0.          to  1.00000e-02
c                                        with weighted error   1.28e-31
c                                         log weighted error  30.89
c                               significant figures required  29.81
c                                    decimal places required  31.48
c
      data algmcs(  1) / +.1666389480 4518632472 0572965082 2 d+0      /
      data algmcs(  2) / -.1384948176 0675638407 3298605913 5 d-4      /
      data algmcs(  3) / +.9810825646 9247294261 5717154748 7 d-8      /
      data algmcs(  4) / -.1809129475 5724941942 6330626671 9 d-10     /
      data algmcs(  5) / +.6221098041 8926052271 2601554341 6 d-13     /
      data algmcs(  6) / -.3399615005 4177219443 0333059966 6 d-15     /
      data algmcs(  7) / +.2683181998 4826987489 5753884666 6 d-17     /
      data algmcs(  8) / -.2868042435 3346432841 4462239999 9 d-19     /
      data algmcs(  9) / +.3962837061 0464348036 7930666666 6 d-21     /
      data algmcs( 10) / -.6831888753 9857668701 1199999999 9 d-23     /
      data algmcs( 11) / +.1429227355 9424981475 7333333333 3 d-24     /
      data algmcs( 12) / -.3547598158 1010705471 9999999999 9 d-26     /
      data algmcs( 13) / +.1025680058 0104709120 0000000000 0 d-27     /
      data algmcs( 14) / -.3401102254 3167487999 9999999999 9 d-29     /
      data algmcs( 15) / +.1276642195 6300629333 3333333333 3 d-30     /
c
      data nalgm, xbig, xmax / 0, 2*0.d0 /
c
      if (nalgm.ne.0) go to 10
      nalgm = initds (algmcs, 15, sngl(d1mach(3)) )
      xbig = 1.0d0/dsqrt(d1mach(3))
      xmax = dexp (dmin1(dlog(d1mach(2)/12.d0), -dlog(12.d0*d1mach(1))))
c
 10   if (x.lt.10.d0) call seteru (23hd9lgmc  x must be ge 10, 23, 1, 2)
      if (x.ge.xmax) go to 20
c
      d9lgmc = 1.d0/(12.d0*x)
      if (x.lt.xbig) d9lgmc = dcsevl (2.0d0*(10.d0/x)**2-1.d0, algmcs,
     1  nalgm) / x
      return
c
 20   d9lgmc = 0.d0
      call seteru (34hd9lgmc  x so big d9lgmc underflows, 34, 2, 0)
      return
c
      end
