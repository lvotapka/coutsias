      double precision function dlog (x)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, alncs(11), center(4), alncen(5), aln2, y, t,
     1  t2, xn,  dcsevl, d1mach
      external d1mach, dcsevl, initds
c
c series for aln        on the interval  0.          to  3.46021e-03
c                                        with weighted error   4.15e-32
c                                         log weighted error  31.38
c                               significant figures required  31.21
c                                    decimal places required  31.90
c
      data aln cs(  1) / +.1334719987 7973881561 6893860471 87 d+1     /
      data aln cs(  2) / +.6937562832 8411286281 3724383542 25 d-3     /
      data aln cs(  3) / +.4293403902 0450834506 5592108036 62 d-6     /
      data aln cs(  4) / +.2893384779 5432594580 4664403875 87 d-9     /
      data aln cs(  5) / +.2051251753 0340580901 7418134477 26 d-12    /
      data aln cs(  6) / +.1503971705 5497386574 6151533199 99 d-15    /
      data aln cs(  7) / +.1129454069 5636464284 5216133333 33 d-18    /
      data aln cs(  8) / +.8635578867 1171868881 9466666666 66 d-22    /
      data aln cs(  9) / +.6695299053 4350370613 3333333333 33 d-25    /
      data aln cs( 10) / +.5249155744 8151466666 6666666666 66 d-28    /
      data aln cs( 11) / +.4153054068 0362666666 6666666666 66 d-31    /
c
      data center(1) / 1.0d0 /
      data center(2) / 1.25d0 /
      data center(3) / 1.50d0 /
      data center(4) / 1.75d0 /
c
      data alncen(  1) / 0.0d0                                         /
      data alncen(  2) / +.2231435513 1420975576 6295090309 83 d+0     /
      data alncen(  3) / +.4054651081 0816438197 8013115464 34 d+0     /
      data alncen(  4) / +.5596157879 3542268627 0888500526 82 d+0     /
      data alncen(  5) / +.6931471805 5994530941 7232121458 17 d+0     /
c
c aln2 = alog(2.0) - 0.625
      data aln2 / 0.0681471805 5994530941 7232121458 18d0 /
      data nterms / 0 /
c
      if (nterms.eq.0) nterms = initds (alncs, 11, 28.9*sngl(d1mach(3)))
c
      if (x.le.0.d0) call seteru (
     1  29hdlog    x is zero or negative, 29, 1, 2)
c
      call d9upak (x, y, n)
c
      xn = n - 1
      y = 2.0d0*y
      ntrval = 4.0d0*y - 2.5d0
c
      if (ntrval.eq.5) t = ((y-1.0d0)-1.0d0) / (y+2.0d0)
      if (ntrval.lt.5) t = (y-center(ntrval)) / (y+center(ntrval))
      t2 = t*t
      dlog = 0.625d0*xn + (aln2*xn + alncen(ntrval) + 2.0d0*t +
     1  t*t2*dcsevl(578.d0*t2-1.0d0, alncs, nterms) )
c
      return
      end
