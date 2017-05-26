      double precision function dexp (x)
c may 1980 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, expcs(14), twon16(17), aln216, f, xint, xmax,
     1  xmin, y,  d1mach, dint, d9pak, dcsevl, dlog
      external d1mach, d9pak, dcsevl, dint, dlog, initds
c
c series for exp        on the interval -1.00000e+00 to  1.00000e+00
c                                        with weighted error   2.30e-34
c                                         log weighted error  33.64
c                               significant figures required  32.28
c                                    decimal places required  34.21
c
      data exp cs(  1) / +.8665694933 1498571273 3404647266 231 d-1    /
      data exp cs(  2) / +.9384948692 9983956189 6336579701 203 d-3    /
      data exp cs(  3) / +.6776039709 9816826407 4353014653 601 d-5    /
      data exp cs(  4) / +.3669312003 9380592780 1891250687 610 d-7    /
      data exp cs(  5) / +.1589590536 1746184464 1928517821 508 d-9    /
      data exp cs(  6) / +.5738598786 3020660125 2990815262 106 d-12   /
      data exp cs(  7) / +.1775744485 9142151180 2306980226 000 d-14   /
      data exp cs(  8) / +.4807991668 4237242267 5950244533 333 d-17   /
      data exp cs(  9) / +.1157163768 8182857280 9260000000 000 d-19   /
      data exp cs( 10) / +.2506506102 5549771993 2458666666 666 d-22   /
      data exp cs( 11) / +.4935717081 4049582848 0000000000 000 d-25   /
      data exp cs( 12) / +.8909295727 4063424000 0000000000 000 d-28   /
      data exp cs( 13) / +.1484480629 0799786666 6666666666 666 d-30   /
      data exp cs( 14) / +.2296789166 3018666666 6666666666 666 d-33   /
c
c twon16(i) is 2.0**((i-1)/16) - 1.0
      data twon16(  1) / 0.0d0                                         /
      data twon16(  2) / +.4427378242 7413840321 9664787399 29 d-1     /
      data twon16(  3) / +.9050773266 5257659207 0106557607 07 d-1     /
      data twon16(  4) / +.1387886347 5669165370 3830283841 51 d+0     /
      data twon16(  5) / +.1892071150 0272106671 7499970560 47 d+0     /
      data twon16(  6) / +.2418578120 7348404859 3677468726 59 d+0     /
      data twon16(  7) / +.2968395546 5100966593 3754117792 45 d+0     /
      data twon16(  8) / +.3542555469 3689272829 8014740140 70 d+0     /
      data twon16(  9) / +.4142135623 7309504880 1688724209 69 d+0     /
      data twon16( 10) / +.4768261459 3949931138 6907480374 04 d+0     /
      data twon16( 11) / +.5422108254 0794082361 2291862090 73 d+0     /
      data twon16( 12) / +.6104903319 4925430817 9520667357 40 d+0     /
      data twon16( 13) / +.6817928305 0742908606 2250952466 42 d+0     /
      data twon16( 14) / +.7562521603 7329948311 2160619375 31 d+0     /
      data twon16( 15) / +.8340080864 0934246348 7083189588 28 d+0     /
      data twon16( 16) / +.9152065613 9714729387 2611270295 83 d+0     /
      data twon16( 17) / 1.d0                                          /
c
c aln216 is 16.0/alog(2.) - 23.0
      data aln216 / +.8312065422 3414517758 7948960302 74 d-1     /
      data nterms, xmin, xmax /0, 2*0.d0 /
c
      if (nterms.ne.0) go to 10
      nterms = initds (expcs, 14, 0.1*sngl(d1mach(3)))
      xmin = dlog (d1mach(1)) + .01d0
      xmax = dlog (d1mach(2)) - 0.001d0
c
 10   if (x.lt.xmin) go to 20
      if (x.gt.xmax) call seteru (
     1  31hdexp    x so big dexp overflows, 31, 2, 2)
c
      xint = dint (x)
      y = x - xint
c
      y = 23.d0*y + x*aln216
      n = y
      f = y - dble(float(n))
      n = 23.d0*xint + dble(float(n))
      n16 = n/16
      if (n.lt.0) n16 = n16 - 1
      ndx = n - 16*n16 + 1
c
      dexp = 1.0d0 + (twon16(ndx) + f*(1.0d0 + twon16(ndx)) *
     1  dcsevl (f, expcs, nterms) )
c
      dexp = d9pak (dexp, n16)
      return
c
 20   call seteru (34hdexp    x so small dexp underflows, 34, 1, 0)
      dexp = 0.d0
      return
c
      end
