      function alog (x)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      dimension alncs(6), center(4), alncen(5)
      external csevl, inits, r1mach
c
c series for aln        on the interval  0.          to  3.46021d-03
c                                        with weighted error   1.50e-16
c                                         log weighted error  15.82
c                               significant figures required  15.65
c                                    decimal places required  16.21
c
      data aln cs( 1) /   1.3347199877 973882e0 /
      data aln cs( 2) /    .0006937562 83284112e0 /
      data aln cs( 3) /    .0000004293 40390204e0 /
      data aln cs( 4) /    .0000000002 89338477e0 /
      data aln cs( 5) /    .0000000000 00205125e0 /
      data aln cs( 6) /    .0000000000 00000150e0 /
c
      data center(1) / 1.0 /
      data center(2) / 1.25 /
      data center(3) / 1.50 /
      data center(4) / 1.75 /
c
      data alncen(  1) / 0.0e0                                         /
      data alncen(  2) / +.2231435513 14209755 e+0                     /
      data alncen(  3) / +.4054651081 08164381 e+0                     /
      data alncen(  4) / +.5596157879 35422686 e+0                     /
      data alncen(  5) / +.6931471805 59945309 e+0                     /
c
c aln2 = alog(2.0) - 0.625
      data aln2 / 0.0681471805 59945309e0 /
      data nterms / 0 /
c
      if (nterms.eq.0) nterms = inits (alncs, 6, 28.9*r1mach(3))
c
      if (x.le.0.) call seteru (
     1  29halog    x is zero or negative, 29, 1, 2)
c
      call r9upak (x, y, n)
c
      xn = n - 1
      y = 2.0*y
      ntrval = 4.0*y - 2.5
      if (ntrval.eq.5) t = ((y-1.0)-1.0) / (y+2.0)
      if (ntrval.lt.5) t = (y-center(ntrval))/(y+center(ntrval))
      t2 = t*t
c
      alog = 0.625*xn + (aln2*xn + alncen(ntrval) + 2.0*t +
     1  t*t2*csevl(578.0*t2-1.0, alncs, nterms) )
c
      return
      end
