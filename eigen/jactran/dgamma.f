C-----------------------------------------------------------------------
      double precision function dgamma (x)
C-----------------------------------------------------------------------
c jan 1984 edition.  w. fullerton, c3, los alamos scientific lab.
c jan 1994 wpp@ips.id.ethz.ch, ehg@research.att.com   declare xsml
      double precision x, gamcs(42), dxrel, pi, sinpiy, sq2pil, xmax,
     1  xmin, y, d9lgmc, dcsevl, d1mach, dexp, dint, dlog,
     2  dsin, dsqrt, xsml
      external d1mach, d9lgmc, dcsevl, dexp, dint, dlog, dsin, dsqrt,
     1  initds
c
c series for gam        on the interval  0.          to  1.00000e+00
c                                        with weighted error   5.79e-32
c                                         log weighted error  31.24
c                               significant figures required  30.00
c                                    decimal places required  32.05
c
      data gam cs(  1) / +.8571195590 9893314219 2006239994 2 d-2      /
      data gam cs(  2) / +.4415381324 8410067571 9131577165 2 d-2      /
      data gam cs(  3) / +.5685043681 5993633786 3266458878 9 d-1      /
      data gam cs(  4) / -.4219835396 4185605010 1250018662 4 d-2      /
      data gam cs(  5) / +.1326808181 2124602205 8400679635 2 d-2      /
      data gam cs(  6) / -.1893024529 7988804325 2394702388 6 d-3      /
      data gam cs(  7) / +.3606925327 4412452565 7808221722 5 d-4      /
      data gam cs(  8) / -.6056761904 4608642184 8554829036 5 d-5      /
      data gam cs(  9) / +.1055829546 3022833447 3182350909 3 d-5      /
      data gam cs( 10) / -.1811967365 5423840482 9185589116 6 d-6      /
      data gam cs( 11) / +.3117724964 7153222777 9025459316 9 d-7      /
      data gam cs( 12) / -.5354219639 0196871408 7408102434 7 d-8      /
      data gam cs( 13) / +.9193275519 8595889468 8778682594 0 d-9      /
      data gam cs( 14) / -.1577941280 2883397617 6742327395 3 d-9      /
      data gam cs( 15) / +.2707980622 9349545432 6654043308 9 d-10     /
      data gam cs( 16) / -.4646818653 8257301440 8166105893 3 d-11     /
      data gam cs( 17) / +.7973350192 0074196564 6076717535 9 d-12     /
      data gam cs( 18) / -.1368078209 8309160257 9949917230 9 d-12     /
      data gam cs( 19) / +.2347319486 5638006572 3347177168 8 d-13     /
      data gam cs( 20) / -.4027432614 9490669327 6657053469 9 d-14     /
      data gam cs( 21) / +.6910051747 3721009121 3833697525 7 d-15     /
      data gam cs( 22) / -.1185584500 2219929070 5238712619 2 d-15     /
      data gam cs( 23) / +.2034148542 4963739552 0102605193 2 d-16     /
      data gam cs( 24) / -.3490054341 7174058492 7401294910 8 d-17     /
      data gam cs( 25) / +.5987993856 4853055671 3505106602 6 d-18     /
      data gam cs( 26) / -.1027378057 8722280744 9006977843 1 d-18     /
      data gam cs( 27) / +.1762702816 0605298249 4275966074 8 d-19     /
      data gam cs( 28) / -.3024320653 7353062609 5877211204 2 d-20     /
      data gam cs( 29) / +.5188914660 2183978397 1783355050 6 d-21     /
      data gam cs( 30) / -.8902770842 4565766924 4925160106 6 d-22     /
      data gam cs( 31) / +.1527474068 4933426022 7459689130 6 d-22     /
      data gam cs( 32) / -.2620731256 1873629002 5732833279 9 d-23     /
      data gam cs( 33) / +.4496464047 8305386703 3104657066 6 d-24     /
      data gam cs( 34) / -.7714712731 3368779117 0390152533 3 d-25     /
      data gam cs( 35) / +.1323635453 1260440364 8657271466 6 d-25     /
      data gam cs( 36) / -.2270999412 9429288167 0231381333 3 d-26     /
      data gam cs( 37) / +.3896418998 0039914493 2081663999 9 d-27     /
      data gam cs( 38) / -.6685198115 1259533277 9212799999 9 d-28     /
      data gam cs( 39) / +.1146998663 1400243843 4761386666 6 d-28     /
      data gam cs( 40) / -.1967938586 3451346772 9510399999 9 d-29     /
      data gam cs( 41) / +.3376448816 5853380903 3489066666 6 d-30     /
      data gam cs( 42) / -.5793070335 7821357846 2549333333 3 d-31     /
c
      data pi / 3.1415926535 8979323846 2643383279 50 d0 /
c sq2pil is 0.5*alog(2*pi) = alog(sqrt(2*pi))
      data sq2pil / 0.9189385332 0467274178 0329736405 62 d0 /
      data ngam, xmin, xmax, xsml, dxrel / 0, 4*0.d0 /
c
      if (ngam.ne.0) go to 10
      ngam = initds (gamcs, 42, 0.1*sngl(d1mach(3)) )
c
      call d9gaml (xmin, xmax)
      xsml = dexp (dmax1 (dlog(d1mach(1)), -dlog(d1mach(2)))+0.01d0)
      dxrel = dsqrt (d1mach(4))
c
 10   y = dabs(x)
      if (y.gt.10.d0) go to 50
c
c compute gamma(x) for -xbnd .le. x .le. xbnd.  reduce interval and find
c gamma(1+y) for 0.0 .le. y .lt. 1.0 first of all.
c
      n = x
      if (x.lt.0.d0) n = n - 1
      y = x - dble(float(n))
      n = n - 1
      dgamma = 0.9375d0 + dcsevl (2.d0*y-1.d0, gamcs, ngam)
      if (n.eq.0) return
c
      if (n.gt.0) go to 30
c
c compute gamma(x) for x .lt. 1.0
c
      n = -n
      if (x.eq.0.d0) call seteru (14hdgamma  x is 0, 14, 4, 2)
      if (x.lt.0.0d0 .and. x+dble(float(n-2)).eq.0.d0) call seteru (
     1  31hdgamma  x is a negative integer, 31, 4, 2)
      if (x.lt.(-0.5d0) .and. dabs((x-dint(x-0.5d0))/x).lt.dxrel) call
     1  seteru (
     2 58hdgamma  answer lt half prec because x too near neg integer
     3  , 58, 1, 1)
      if (y.lt.xsml) call seteru (
     1  54hdgamma  x is so close to 0.0 that the result overflows,
     2  54, 5, 2)
c
      do 20 i=1,n
        dgamma = dgamma/(x+dble(float(i-1)) )
 20   continue
      return
c
c gamma(x) for x .ge. 2.0 and x .le. 10.0
c
 30   do 40 i=1,n
        dgamma = (y+dble(float(i))) * dgamma
 40   continue
      return
c
c gamma(x) for dabs(x) .gt. 10.0.  recall y = dabs(x).
c
 50   if (x.gt.xmax) call seteru (32hdgamma  x so big gamma overflows,
     1  32, 3, 2)
c
      dgamma = 0.d0
      if (x.lt.xmin) call seteru (35hdgamma  x so small gamma underflows
     1  , 35, 2, 0)
      if (x.lt.xmin) return
c
      dgamma = dexp ((y-0.5d0)*dlog(y) - y + sq2pil + d9lgmc(y) )
      if (x.gt.0.d0) return
c
      if (dabs((x-dint(x-0.5d0))/x).lt.dxrel) call seteru (
     1  61hdgamma  answer lt half precision, x too near negative integer
     2  , 61, 1, 1)
c
      sinpiy = dsin (pi*y)
      if (sinpiy.eq.0.d0) call seteru (
     1  31hdgamma  x is a negative integer, 31, 4, 2)
c
      dgamma = -pi/(y*sinpiy*dgamma)
c
      return
      end
C-----------------------------------------------------------------------
      function alog (x)
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
      function csevl (x, cs, n)
C-----------------------------------------------------------------------
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c evaluate the n-term chebyshev series cs at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).  also see fox
c and parker, chebyshev polys in numerical analysis, oxford press, p.56.
c
c             input arguments --
c x      value at which the series is to be evaluated.
c cs     array of n terms of a chebyshev series.  in eval-
c        uating cs, only half the first coef is summed.
c n      number of terms in array cs.
c
      dimension cs(1)
c
      if (n.lt.1) call seteru (28hcsevl   number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hcsevl   number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1) .or. x.gt.1.1) call seteru (
     1  25hcsevl   x outside (-1,+1), 25, 1, 1)
c
      b1 = 0.
      b0 = 0.
      twox = 2.*x
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n + 1 - i
        b0 = twox*b1 - b2 + cs(ni)
 10   continue
c
      csevl = 0.5 * (b0-b2)
c
      return
      end
C-----------------------------------------------------------------------
      DOUBLE PRECISION FUNCTION D1MACH(I)
C-----------------------------------------------------------------------
C
C  DOUBLE-PRECISION MACHINE CONSTANTS
C
C  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  D1MACH( 5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
C  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
C  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.  IF YOU DO NOT
C  KNOW WHICH SET TO USE, TRY BOTH AND SEE WHICH GIVES PLAUSIBLE
C  VALUES.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR D1MACH.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER SC
C
      DOUBLE PRECISION DMACH(5)
C
      EQUIVALENCE (DMACH(1),SMALL(1))
      EQUIVALENCE (DMACH(2),LARGE(1))
      EQUIVALENCE (DMACH(3),RIGHT(1))
      EQUIVALENCE (DMACH(4),DIVER(1))
      EQUIVALENCE (DMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR BIG-ENDIAN IEEE ARITHMETIC (BINARY FORMAT)
C     MACHINES IN WHICH THE MOST SIGNIFICANT BYTE IS STORED FIRST,
C     SUCH AS THE AT&T 3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G.
C     SUN 3), AND MACHINES THAT USE SPARC, HP, OR IBM RISC CHIPS.
C
       DATA SMALL(1),SMALL(2) /    1048576,          0 /
       DATA LARGE(1),LARGE(2) / 2146435071,         -1 /
       DATA RIGHT(1),RIGHT(2) / 1017118720,          0 /
       DATA DIVER(1),DIVER(2) / 1018167296,          0 /
       DATA LOG10(1),LOG10(2) / 1070810131, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR LITTLE-ENDIAN (BINARY) IEEE ARITHMETIC
C     MACHINES IN WHICH THE LEAST SIGNIFICANT BYTE IS STORED FIRST,
C     E.G. IBM PCS AND OTHER MACHINES THAT USE INTEL 80X87 OR DEC
C     ALPHA CHIPS.
C
C      DATA SMALL(1),SMALL(2) /          0,    1048576 /
C      DATA LARGE(1),LARGE(2) /         -1, 2146435071 /
C      DATA RIGHT(1),RIGHT(2) /          0, 1017118720 /
C      DATA DIVER(1),DIVER(2) /          0, 1018167296 /
C      DATA LOG10(1),LOG10(2) / 1352628735, 1070810131 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1),SMALL(2) /    1048576,          0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,         -1 /
C      DATA RIGHT(1),RIGHT(2) /  856686592,          0 /
C      DATA DIVER(1),DIVER(2) /  873463808,          0 /
C      DATA LOG10(1),LOG10(2) / 1091781651, 1352628735 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA SMALL(1) / ZC00800000 /
C      DATA SMALL(2) / Z000000000 /
C
C      DATA LARGE(1) / ZDFFFFFFFF /
C      DATA LARGE(2) / ZFFFFFFFFF /
C
C      DATA RIGHT(1) / ZCC5800000 /
C      DATA RIGHT(2) / Z000000000 /
C
C      DATA DIVER(1) / ZCC6800000 /
C      DATA DIVER(2) / Z000000000 /
C
C      DATA LOG10(1) / ZD00E730E7 /
C      DATA LOG10(2) / ZC77800DC0 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O0000000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O0007777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA SMALL(1) / O1771000000000000 /
C      DATA SMALL(2) / O7770000000000000 /
C
C      DATA LARGE(1) / O0777777777777777 /
C      DATA LARGE(2) / O7777777777777777 /
C
C      DATA RIGHT(1) / O1461000000000000 /
C      DATA RIGHT(2) / O0000000000000000 /
C
C      DATA DIVER(1) / O1451000000000000 /
C      DATA DIVER(2) / O0000000000000000 /
C
C      DATA LOG10(1) / O1157163034761674 /
C      DATA LOG10(2) / O0006677466732724 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / 00564000000000000000B /
C      DATA SMALL(2) / 00000000000000000000B /
C
C      DATA LARGE(1) / 37757777777777777777B /
C      DATA LARGE(2) / 37157777777777777774B /
C
C      DATA RIGHT(1) / 15624000000000000000B /
C      DATA RIGHT(2) / 00000000000000000000B /
C
C      DATA DIVER(1) / 15634000000000000000B /
C      DATA DIVER(2) / 00000000000000000000B /
C
C      DATA LOG10(1) / 17164642023241175717B /
C      DATA LOG10(2) / 16367571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA SMALL(1) / O"00564000000000000000" /
C      DATA SMALL(2) / O"00000000000000000000" /
C
C      DATA LARGE(1) / O"37757777777777777777" /
C      DATA LARGE(2) / O"37157777777777777774" /
C
C      DATA RIGHT(1) / O"15624000000000000000" /
C      DATA RIGHT(2) / O"00000000000000000000" /
C
C      DATA DIVER(1) / O"15634000000000000000" /
C      DATA DIVER(2) / O"00000000000000000000" /
C
C      DATA LOG10(1) / O"17164642023241175717" /
C      DATA LOG10(2) / O"16367571421742254654" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1
C
C      DATA SMALL(1),SMALL(2) / '00100000'X, '00000000'X /
C      DATA LARGE(1),LARGE(2) / '7FFFFFFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) / '3CC00000'X, '00000000'X /
C      DATA DIVER(1),DIVER(2) / '3CD00000'X, '00000000'X /
C      DATA LOG10(1),LOG10(2) / '3FF34413'X, '509F79FF'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA SMALL(1) / 201354000000000000000B /
C      DATA SMALL(2) / 000000000000000000000B /
C
C      DATA LARGE(1) / 577767777777777777777B /
C      DATA LARGE(2) / 000007777777777777776B /
C
C      DATA RIGHT(1) / 376434000000000000000B /
C      DATA RIGHT(2) / 000000000000000000000B /
C
C      DATA DIVER(1) / 376444000000000000000B /
C      DATA DIVER(2) / 000000000000000000000B /
C
C      DATA LOG10(1) / 377774642023241175717B /
C      DATA LOG10(2) / 000007571421742254654B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C
C     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED
C     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC DMACH(5)
C
C      DATA SMALL/20K,3*0/,LARGE/77777K,3*177777K/
C      DATA RIGHT/31420K,3*0/,DIVER/32020K,3*0/
C      DATA LOG10/40423K,42023K,50237K,74776K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '37777577 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000333 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000334 /
C      DATA LOG10(1),LOG10(2) / '23210115, '10237777 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O402400000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O376777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O604400000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O606400000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O776464202324, O117571775714 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA SMALL(1),SMALL(2) / Z00100000, Z00000000 /
C      DATA LARGE(1),LARGE(2) / Z7FFFFFFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z33100000, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z34100000, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z41134413, Z509F79FF /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA SMALL(1),SMALL(2) / Z'00100000', Z'00000000' /
C      DATA LARGE(1),LARGE(2) / Z'7EFFFFFF', Z'FFFFFFFF' /
C      DATA RIGHT(1),RIGHT(2) / Z'33100000', Z'00000000' /
C      DATA DIVER(1),DIVER(2) / Z'34100000', Z'00000000' /
C      DATA LOG10(1),LOG10(2) / Z'41134413', Z'509F79FF' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "033400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "344777777777 /
C      DATA RIGHT(1),RIGHT(2) / "113400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "114400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "144117571776 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA SMALL(1),SMALL(2) / "000400000000, "000000000000 /
C      DATA LARGE(1),LARGE(2) / "377777777777, "377777777777 /
C      DATA RIGHT(1),RIGHT(2) / "103400000000, "000000000000 /
C      DATA DIVER(1),DIVER(2) / "104400000000, "000000000000 /
C      DATA LOG10(1),LOG10(2) / "177464202324, "047674776746 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /    8388608,           0 /
C      DATA LARGE(1),LARGE(2) / 2147483647,          -1 /
C      DATA RIGHT(1),RIGHT(2) /  612368384,           0 /
C      DATA DIVER(1),DIVER(2) /  620756992,           0 /
C      DATA LOG10(1),LOG10(2) / 1067065498, -2063872008 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O00040000000, O00000000000 /
C      DATA LARGE(1),LARGE(2) / O17777777777, O37777777777 /
C      DATA RIGHT(1),RIGHT(2) / O04440000000, O00000000000 /
C      DATA DIVER(1),DIVER(2) / O04500000000, O00000000000 /
C      DATA LOG10(1),LOG10(2) / O07746420232, O20476747770 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     SMALL, LARGE, RIGHT, DIVER, LOG10 SHOULD BE DECLARED
C     INTEGER SMALL(4), LARGE(4), RIGHT(4), DIVER(4), LOG10(4)
C
C      DATA SMALL(1),SMALL(2) /    128,      0 /
C      DATA SMALL(3),SMALL(4) /      0,      0 /
C
C      DATA LARGE(1),LARGE(2) /  32767,     -1 /
C      DATA LARGE(3),LARGE(4) /     -1,     -1 /
C
C      DATA RIGHT(1),RIGHT(2) /   9344,      0 /
C      DATA RIGHT(3),RIGHT(4) /      0,      0 /
C
C      DATA DIVER(1),DIVER(2) /   9472,      0 /
C      DATA DIVER(3),DIVER(4) /      0,      0 /
C
C      DATA LOG10(1),LOG10(2) /  16282,   8346 /
C      DATA LOG10(3),LOG10(4) / -31493, -12296 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA SMALL(3),SMALL(4) / O000000, O000000 /
C
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA LARGE(3),LARGE(4) / O177777, O177777 /
C
C      DATA RIGHT(1),RIGHT(2) / O022200, O000000 /
C      DATA RIGHT(3),RIGHT(4) / O000000, O000000 /
C
C      DATA DIVER(1),DIVER(2) / O022400, O000000 /
C      DATA DIVER(3),DIVER(4) / O000000, O000000 /
C
C      DATA LOG10(1),LOG10(2) / O037632, O020232 /
C      DATA LOG10(3),LOG10(4) / O102373, O147770 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WITH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA SMALL(1),SMALL(2) / :10000000000, :00000100001 /
C      DATA LARGE(1),LARGE(2) / :17777777777, :37777677775 /
C      DATA RIGHT(1),RIGHT(2) / :10000000000, :00000000122 /
C      DATA DIVER(1),DIVER(2) / :10000000000, :00000000123 /
C      DATA LOG10(1),LOG10(2) / :11504046501, :07674600177 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000
C
C      DATA SMALL(1),SMALL(2) / $00000000,  $00100000 /
C      DATA LARGE(1),LARGE(2) / $FFFFFFFF,  $7FEFFFFF /
C      DATA RIGHT(1),RIGHT(2) / $00000000,  $3CA00000 /
C      DATA DIVER(1),DIVER(2) / $00000000,  $3CB00000 /
C      DATA LOG10(1),LOG10(2) / $509F79FF,  $3FD34413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA SMALL(1),SMALL(2) / O000040000000, O000000000000 /
C      DATA LARGE(1),LARGE(2) / O377777777777, O777777777777 /
C      DATA RIGHT(1),RIGHT(2) / O170540000000, O000000000000 /
C      DATA DIVER(1),DIVER(2) / O170640000000, O000000000000 /
C      DATA LOG10(1),LOG10(2) / O177746420232, O411757177572 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER
C
C      DATA SMALL(1),SMALL(2) /        128,           0 /
C      DATA LARGE(1),LARGE(2) /     -32769,          -1 /
C      DATA RIGHT(1),RIGHT(2) /       9344,           0 /
C      DATA DIVER(1),DIVER(2) /       9472,           0 /
C      DATA LOG10(1),LOG10(2) /  546979738,  -805796613 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER
C
C      DATA SMALL(1),SMALL(2) / Z00000080, Z00000000 /
C      DATA LARGE(1),LARGE(2) / ZFFFF7FFF, ZFFFFFFFF /
C      DATA RIGHT(1),RIGHT(2) / Z00002480, Z00000000 /
C      DATA DIVER(1),DIVER(2) / Z00002500, Z00000000 /
C      DATA LOG10(1),LOG10(2) / Z209A3F9A, ZCFF884FB /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2
C
C      DATA SMALL(1),SMALL(2) /       '80'X,        '0'X /
C      DATA LARGE(1),LARGE(2) / 'FFFF7FFF'X, 'FFFFFFFF'X /
C      DATA RIGHT(1),RIGHT(2) /     '2480'X,        '0'X /
C      DATA DIVER(1),DIVER(2) /     '2500'X,        '0'X /
C      DATA LOG10(1),LOG10(2) / '209A3F9A'X, 'CFF884FB'X /, SC/987/
C
C  ***  ISSUE STOP 779 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) STOP 779
C  ***  ISSUE STOP 778 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (DMACH(4) .GE. 1.0D0) STOP 778
      IF (I .LT. 1  .OR.  I .GT. 5) GOTO 999
      D1MACH = DMACH(I)
      RETURN
  999 WRITE(*,1999) I
 1999 FORMAT(' D1MACH - I OUT OF BOUNDS',I10)
      STOP
      END
C-----------------------------------------------------------------------
      subroutine d9gaml (xmin, xmax)
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
      double precision function d9lgmc (x)
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
      double precision function d9pak (y, n)
C-----------------------------------------------------------------------
c december 1979 edition. w. fullerton, c3, los alamos scientific lab.
c
c pack a base 2 exponent into floating point number x.  this routine is
c almost the inverse of d9upak.  it is not exactly the inverse, because
c dabs(x) need not be between 0.5 and 1.0.  if both d9pak and 2.d0**n
c were known to be in range we could compute
c                d9pak = x * 2.0d0**n
c
      double precision y, aln2b, aln210, d1mach
      external d1mach, i1mach
      data nmin, nmax / 2*0 /
      data aln210 / 3.321928094 8873623478 7031942948 9 d0 /
c
      if (nmin.ne.0) go to 10
      aln2b = 1.0d0
      if (i1mach(10).ne.2) aln2b = d1mach(5)*aln210
      nmin = aln2b*dble(float(i1mach(15)))
      nmax = aln2b*dble(float(i1mach(16)))
c
 10   call d9upak (y, d9pak, ny)
c
      nsum = n + ny
      if (nsum.lt.nmin) go to 40
      if (nsum.gt.nmax) call seteru (
     1  31hd9pak   packed number overflows, 31, 1, 2)
c
      if (nsum.eq.0) return
      if (nsum.gt.0) go to 30
c
 20   d9pak = 0.5d0*d9pak
      nsum = nsum + 1
      if (nsum.ne.0) go to 20
      return
c
 30   d9pak = 2.0d0*d9pak
      nsum = nsum - 1
      if (nsum.ne.0) go to 30
      return
c
 40   call seteru (32hd9pak   packed number underflows, 32, 1, 0)
      d9pak = 0.0d0
      return
c
      end
C-----------------------------------------------------------------------
      subroutine d9upak (x, y, n)
C-----------------------------------------------------------------------
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      double precision x, y, absx
c
      absx = dabs(x)
      n = 0
      y = 0.0d0
      if (x.eq.0.0d0) return
c
 10   if (absx.ge.0.5d0) go to 20
      n = n - 1
      absx = absx*2.0d0
      go to 10
c
 20   if (absx.lt.1.0d0) go to 30
      n = n + 1
      absx = absx*0.5d0
      go to 20
c
 30   y = dsign (absx, x)
      return
c
      end
C-----------------------------------------------------------------------
      double precision function dcsevl (x, a, n)
C-----------------------------------------------------------------------
c evaluate the n-term chebyshev series a at x.  adapted from
c r. broucke, algorithm 446, c.a.c.m., 16, 254 (1973).
c
c             input arguments --
c x      dble prec value at which the series is to be evaluated.
c a      dble prec array of n terms of a chebyshev series.  in eval-
c        uating a, only half the first coef is summed.
c n      number of terms in array a.
c
      double precision a(n), x, twox, b0, b1, b2
c
      if (n.lt.1) call seteru (28hdcsevl  number of terms le 0, 28, 2,2)
      if (n.gt.1000) call seteru (31hdcsevl  number of terms gt 1000,
     1  31, 3, 2)
      if (x.lt.(-1.1d0) .or. x.gt.1.1d0) call seteru (
     1  25hdcsevl  x outside (-1,+1), 25, 1, 1)
c
      twox = 2.0d0*x
      b1 = 0.d0
      b0 = 0.d0
      do 10 i=1,n
        b2 = b1
        b1 = b0
        ni = n - i + 1
        b0 = twox*b1 - b2 + a(ni)
 10   continue
c
      dcsevl = 0.5d0 * (b0-b2)
c
      return
      end
C-----------------------------------------------------------------------
      double precision function dexp (x)
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
      double precision function dint (x)
C-----------------------------------------------------------------------
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
     1 57hdint dabs(x) may be too big to be repd a an exact integer
     2    , 57, 1, 1)
      dint = x
      return
c
      end
C-----------------------------------------------------------------------
      double precision function dlog (x)
C-----------------------------------------------------------------------
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
C-----------------------------------------------------------------------
      double precision function dsin (x)
C-----------------------------------------------------------------------
c august 1980 edition.  w. fullerton, los alamos scientific lab.
c
c this routine is based on the algorithm of cody and waite in
c argonne tm-321, software manual working note number 1
c
      double precision x, sincs(15), pihi, pilo, pirec, pi2rec, xsml,
     1  xwarn, xmax, y, xn, sgn, f, dint, dcsevl, d1mach, dsqrt
      external d1mach, dcsevl, dint, dsqrt, initds
c
c series for sin    on the interval  0.00000e+00 to  2.46740e+00
c                                        with weighted error   2.56e-34
c                                         log weighted error  33.59
c                               significant figures required  33.01
c                                    decimal places required  34.18
c
      data sin cs(  1) / -0.3749911549 5587317583 9919279977 323464d0/
      data sin cs(  2) / -0.1816031552 3725020186 3830316158 004754d0/
      data sin cs(  3) /  0.0058047092 7459863355 9427341722 857921d0/
      data sin cs(  4) / -0.0000869543 1177934075 7113212316 353178d0/
      data sin cs(  5) /  0.0000007543 7014808885 1481006839 927030d0/
      data sin cs(  6) / -0.0000000042 6712966505 5961107126 829906d0/
      data sin cs(  7) /  0.0000000000 1698042294 5488168181 824792d0/
      data sin cs(  8) / -0.0000000000 0005012057 8889961870 929524d0/
      data sin cs(  9) /  0.0000000000 0000011410 1026680010 675628d0/
      data sin cs( 10) / -0.0000000000 0000000020 6437504424 783134d0/
      data sin cs( 11) /  0.0000000000 0000000000 0303969595 918706d0/
      data sin cs( 12) / -0.0000000000 0000000000 0000371357 734157d0/
      data sin cs( 13) /  0.0000000000 0000000000 0000000382 486123d0/
      data sin cs( 14) / -0.0000000000 0000000000 0000000000 336623d0/
      data sin cs( 15) /  0.0000000000 0000000000 0000000000 000256d0/
c
c pihi + pilo = pi.  pihi is exactly representable on all machines
c with at least 8 bits of precision.  whether it is exactly
c represented depends on the compiler.  this routine is more
c accurate if it is exactly represented.
      data pihi / 3.140625d0 /
      data pilo / 9.676535897 9323846264 3383279502 88d-4/
      data pirec / 0.3183098861 8379067153 7767526745 03d0 /
      data pi2rec / 0.6366197723 6758134307 5535053490 06d0 /
      data ntsn, xsml, xwarn, xmax / 0, 3*0.0d0 /
c
      if (ntsn.ne.0) go to 10
      ntsn = initds (sincs, 15, 0.1*sngl(d1mach(3)))
c
      xsml = dsqrt (2.0d0*d1mach(3))
      xmax = 1.0d0/d1mach(4)
      xwarn = dsqrt (xmax)
c
 10   y = dabs (x)
      if (y.gt.xmax) call seteru (
     1  42hdsin    no precision because abs(x) is big, 42, 2, 2)
      if (y.gt.xwarn) call seteru (
     1  54hdsin    answer lt half precision because abs(x) is big,
     2  54, 1, 1)
c
      dsin = x
      if (y.lt.xsml) return
c
      xn = dint (y*pirec+0.5d0)
      n2 = dmod (xn, 2.0d0) + 0.5d0
      sgn = x
      if (n2.ne.0) sgn = -sgn
      f = (y-xn*pihi) - xn*pilo
c
      dsin = f + f*dcsevl(2.0d0*(f*pi2rec)**2-1.0d0, sincs, ntsn)
      if (sgn.lt.0.0d0) dsin = -dsin
      if (dabs(dsin).gt.1.0d0) dsin = dsign (1.0d0, dsin)
c
      return
      end
C-----------------------------------------------------------------------
      double precision function dsqrt (x)
C-----------------------------------------------------------------------
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
      double precision x, sqrt2(3), y,  d9pak, d1mach
      external alog, d1mach, d9pak
      data sqrt2(1) / 0.7071067811 8654752440 0844362104 85 d0 /
      data sqrt2(2) / 1.0 d0 /
      data sqrt2(3) / 1.4142135623 7309504880 1688724209 70 d0 /
c
      data niter / 0 /
c
      if (niter.eq.0) niter = 1.443*alog(-0.104*alog(0.1*sngl(d1mach(3))
     1  )) + 1.0
c
      if (x.le.0.d0) go to 20
c
      call d9upak (x, y, n)
      ixpnt = n/2
      irem = n - 2*ixpnt + 2
c
c the approximation below has accuracy of 4.16 digits.
      z = y
      dsqrt = .261599e0 + z*(1.114292e0 + z*(-.516888e0 + z*.141067e0))
c
      do 10 iter=1,niter
        dsqrt = dsqrt + 0.5d0*(y - dsqrt*dsqrt) / dsqrt
 10   continue
c
      dsqrt = d9pak (sqrt2(irem)*dsqrt, ixpnt)
      return
c
 20   if (x.lt.0.d0) call seteru (21hdsqrt   x is negative, 21, 1, 1)
      dsqrt = 0.0d0
      return
c
      end
C-----------------------------------------------------------------------
      subroutine e9rint(messg,nw,nerr,save)
C-----------------------------------------------------------------------
c
c  this routine stores the current error message or prints the old one,
c  if any, depending on whether or not save = .true. .
c
      integer messg(nw)
      logical save
      external i1mach, i8save
c
c  messgp stores at least the first 72 characters of the previous
c  message. its length is machine dependent and must be at least
c
c       1 + 71/(the number of characters stored per integer word).
c
      integer messgp(36),fmt(14),ccplus
c
c  start with no previous message.
c
      data messgp(1)/1h1/, nwp/0/, nerrp/0/
c
c  set up the format for printing the error message.
c  the format is simply (a1,14x,72axx) where xx=i1mach(6) is the
c  number of characters stored per integer word.
c
      data ccplus  / 1h+ /
c
      data fmt( 1) / 1h( /
      data fmt( 2) / 1ha /
      data fmt( 3) / 1h1 /
      data fmt( 4) / 1h, /
      data fmt( 5) / 1h1 /
      data fmt( 6) / 1h4 /
      data fmt( 7) / 1hx /
      data fmt( 8) / 1h, /
      data fmt( 9) / 1h7 /
      data fmt(10) / 1h2 /
      data fmt(11) / 1ha /
      data fmt(12) / 1hx /
      data fmt(13) / 1hx /
      data fmt(14) / 1h) /
c
      if (.not.save) go to 20
c
c  save the message.
c
        nwp=nw
        nerrp=nerr
        do 10 i=1,nw
 10     messgp(i)=messg(i)
c
        go to 30
c
 20   if (i8save(1,0,.false.).eq.0) go to 30
c
c  print the message.
c
        iwunit=i1mach(4)
        write(iwunit,9000) nerrp
 9000   format(7h error ,i4,4h in )
c
        call s88fmt(2,i1mach(6),fmt(12))
        write(iwunit,fmt) ccplus,(messgp(i),i=1,nwp)
c
 30   return
c
      end
C-----------------------------------------------------------------------
      subroutine eprint
C-----------------------------------------------------------------------
c
c  this subroutine prints the last error message, if any.
c
      integer messg(1)
c
      call e9rint(messg,1,1,.false.)
      return
c
      end
      INTEGER FUNCTION I1MACH(I)
      INTEGER I
C
C  I/O UNIT NUMBERS.
C
C    I1MACH( 1) = THE STANDARD INPUT UNIT.
C
C    I1MACH( 2) = THE STANDARD OUTPUT UNIT.
C
C    I1MACH( 3) = THE STANDARD PUNCH UNIT.
C
C    I1MACH( 4) = THE STANDARD ERROR MESSAGE UNIT.
C
C  WORDS.
C
C    I1MACH( 5) = THE NUMBER OF BITS PER INTEGER STORAGE UNIT.
C
C    I1MACH( 6) = THE NUMBER OF CHARACTERS PER CHARACTER STORAGE UNIT.
C                 FOR FORTRAN 77, THIS IS ALWAYS 1.  FOR FORTRAN 66,
C                 CHARACTER STORAGE UNIT = INTEGER STORAGE UNIT.
C
C  INTEGERS.
C
C    ASSUME INTEGERS ARE REPRESENTED IN THE S-DIGIT, BASE-A FORM
C
C               SIGN ( X(S-1)*A**(S-1) + ... + X(1)*A + X(0) )
C
C               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,S-1.
C
C    I1MACH( 7) = A, THE BASE.
C
C    I1MACH( 8) = S, THE NUMBER OF BASE-A DIGITS.
C
C    I1MACH( 9) = A**S - 1, THE LARGEST MAGNITUDE.
C
C  FLOATING-POINT NUMBERS.
C
C    ASSUME FLOATING-POINT NUMBERS ARE REPRESENTED IN THE T-DIGIT,
C    BASE-B FORM
C
C               SIGN (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C               WHERE 0 .LE. X(I) .LT. B FOR I=1,...,T,
C               0 .LT. X(1), AND EMIN .LE. E .LE. EMAX.
C
C    I1MACH(10) = B, THE BASE.
C
C  SINGLE-PRECISION
C
C    I1MACH(11) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(12) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(13) = EMAX, THE LARGEST EXPONENT E.
C
C  DOUBLE-PRECISION
C
C    I1MACH(14) = T, THE NUMBER OF BASE-B DIGITS.
C
C    I1MACH(15) = EMIN, THE SMALLEST EXPONENT E.
C
C    I1MACH(16) = EMAX, THE LARGEST EXPONENT E.
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.  ALSO, THE VALUES OF
C  I1MACH(1) - I1MACH(4) SHOULD BE CHECKED FOR CONSISTENCY
C  WITH THE LOCAL OPERATING SYSTEM.  FOR FORTRAN 77, YOU MAY WISH
C  TO ADJUST THE DATA STATEMENT SO IMACH(6) IS SET TO 1, AND
C  THEN TO COMMENT OUT THE EXECUTABLE TEST ON I .EQ. 6 BELOW.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE, EXCEPT PERHAPS
C  FOR IMACH(1) - IMACH(4).
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR I1MACH.
C
      INTEGER CRAY1, IMACH(16), OUTPUT, SANITY, SMALL(2)
      COMMON /D8MACH/ CRAY1
C/6S
C/7S
      SAVE IMACH, SANITY
C/
      REAL RMACH
C
      EQUIVALENCE (IMACH(4),OUTPUT), (RMACH,SMALL(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   32 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
       DATA IMACH(11) /   24 /
       DATA IMACH(12) / -125 /
       DATA IMACH(13) /  128 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
       DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA IMACH( 1) /    7 /
C      DATA IMACH( 2) /    2 /
C      DATA IMACH( 3) /    2 /
C      DATA IMACH( 4) /    2 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   33 /
C      DATA IMACH( 9) / Z1FFFFFFFF /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -256 /
C      DATA IMACH(13) /  255 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) / -256 /
C      DATA IMACH(16) /  255 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -50 /
C      DATA IMACH(16) /  76 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  48 /
C      DATA IMACH( 6) /   6 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  39 /
C      DATA IMACH( 9) / O0007777777777777 /
C      DATA IMACH(10) /   8 /
C      DATA IMACH(11) /  13 /
C      DATA IMACH(12) / -50 /
C      DATA IMACH(13) /  76 /
C      DATA IMACH(14) /  26 /
C      DATA IMACH(15) / -32754 /
C      DATA IMACH(16) /  32780 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / 00007777777777777777B /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   47 /
C      DATA IMACH(12) / -929 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   94 /
C      DATA IMACH(15) / -929 /
C      DATA IMACH(16) / 1069 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   60 /
C      DATA IMACH( 6) /   10 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   48 /
C      DATA IMACH( 9) / O"00007777777777777777" /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   47 /
C      DATA IMACH(12) / -929 /
C      DATA IMACH(13) / 1070 /
C      DATA IMACH(14) /   94 /
C      DATA IMACH(15) / -929 /
C      DATA IMACH(16) / 1069 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   53 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA IMACH( 1) /     5 /
C      DATA IMACH( 2) /     6 /
C      DATA IMACH( 3) /   102 /
C      DATA IMACH( 4) /     6 /
C      DATA IMACH( 5) /    64 /
C      DATA IMACH( 6) /     8 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    63 /
C      DATA IMACH( 9) /  777777777777777777777B /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    47 /
C      DATA IMACH(12) / -8189 /
C      DATA IMACH(13) /  8190 /
C      DATA IMACH(14) /    94 /
C      DATA IMACH(15) / -8099 /
C      DATA IMACH(16) /  8190 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C      DATA IMACH( 1) /   11 /
C      DATA IMACH( 2) /   12 /
C      DATA IMACH( 3) /    8 /
C      DATA IMACH( 4) /   10 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) /32767 /
C      DATA IMACH(10) /   16 /
C      DATA IMACH(11) /    6 /
C      DATA IMACH(12) /  -64 /
C      DATA IMACH(13) /   63 /
C      DATA IMACH(14) /   14 /
C      DATA IMACH(15) /  -64 /
C      DATA IMACH(16) /   63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA IMACH( 1) /       5 /
C      DATA IMACH( 2) /       6 /
C      DATA IMACH( 3) /       0 /
C      DATA IMACH( 4) /       6 /
C      DATA IMACH( 5) /      24 /
C      DATA IMACH( 6) /       3 /
C      DATA IMACH( 7) /       2 /
C      DATA IMACH( 8) /      23 /
C      DATA IMACH( 9) / 8388607 /
C      DATA IMACH(10) /       2 /
C      DATA IMACH(11) /      23 /
C      DATA IMACH(12) /    -127 /
C      DATA IMACH(13) /     127 /
C      DATA IMACH(14) /      38 /
C      DATA IMACH(15) /    -127 /
C      DATA IMACH(16) /     127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /   43 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   63 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   7 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z7FFFFFFF /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  63 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  63 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA IMACH( 1) /   5 /
C      DATA IMACH( 2) /   6 /
C      DATA IMACH( 3) /   6 /
C      DATA IMACH( 4) /   6 /
C      DATA IMACH( 5) /  32 /
C      DATA IMACH( 6) /   4 /
C      DATA IMACH( 7) /   2 /
C      DATA IMACH( 8) /  31 /
C      DATA IMACH( 9) / Z'7FFFFFFF' /
C      DATA IMACH(10) /  16 /
C      DATA IMACH(11) /   6 /
C      DATA IMACH(12) / -64 /
C      DATA IMACH(13) /  62 /
C      DATA IMACH(14) /  14 /
C      DATA IMACH(15) / -64 /
C      DATA IMACH(16) /  62 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   54 /
C      DATA IMACH(15) / -101 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    5 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / "377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   62 /
C      DATA IMACH(15) / -128 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGER ARITHMETIC.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   16 /
C      DATA IMACH( 6) /    2 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   15 /
C      DATA IMACH( 9) / 32767 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE PRIME 50 SERIES SYSTEMS
C     WTIH 32-BIT INTEGERS AND 64V MODE INSTRUCTIONS,
C     SUPPLIED BY IGOR BRAY.
C
C      DATA IMACH( 1) /            1 /
C      DATA IMACH( 2) /            1 /
C      DATA IMACH( 3) /            2 /
C      DATA IMACH( 4) /            1 /
C      DATA IMACH( 5) /           32 /
C      DATA IMACH( 6) /            4 /
C      DATA IMACH( 7) /            2 /
C      DATA IMACH( 8) /           31 /
C      DATA IMACH( 9) / :17777777777 /
C      DATA IMACH(10) /            2 /
C      DATA IMACH(11) /           23 /
C      DATA IMACH(12) /         -127 /
C      DATA IMACH(13) /         +127 /
C      DATA IMACH(14) /           47 /
C      DATA IMACH(15) /       -32895 /
C      DATA IMACH(16) /       +32637 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA IMACH( 1) /     0 /
C      DATA IMACH( 2) /     0 /
C      DATA IMACH( 3) /     7 /
C      DATA IMACH( 4) /     0 /
C      DATA IMACH( 5) /    32 /
C      DATA IMACH( 6) /     1 /
C      DATA IMACH( 7) /     2 /
C      DATA IMACH( 8) /    31 /
C      DATA IMACH( 9) /  2147483647 /
C      DATA IMACH(10) /     2 /
C      DATA IMACH(11) /    24 /
C      DATA IMACH(12) /  -125 /
C      DATA IMACH(13) /   128 /
C      DATA IMACH(14) /    53 /
C      DATA IMACH(15) / -1021 /
C      DATA IMACH(16) /  1024 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C     NOTE THAT THE PUNCH UNIT, I1MACH(3), HAS BEEN SET TO 7
C     WHICH IS APPROPRIATE FOR THE UNIVAC-FOR SYSTEM.
C     IF YOU HAVE THE UNIVAC-FTN SYSTEM, SET IT TO 1.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   36 /
C      DATA IMACH( 6) /    6 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   35 /
C      DATA IMACH( 9) / O377777777777 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   27 /
C      DATA IMACH(12) / -128 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   60 /
C      DATA IMACH(15) /-1024 /
C      DATA IMACH(16) / 1023 /, SANITY/987/
C
C     MACHINE CONSTANTS FOR VAX.
C
C      DATA IMACH( 1) /    5 /
C      DATA IMACH( 2) /    6 /
C      DATA IMACH( 3) /    7 /
C      DATA IMACH( 4) /    6 /
C      DATA IMACH( 5) /   32 /
C      DATA IMACH( 6) /    4 /
C      DATA IMACH( 7) /    2 /
C      DATA IMACH( 8) /   31 /
C      DATA IMACH( 9) / 2147483647 /
C      DATA IMACH(10) /    2 /
C      DATA IMACH(11) /   24 /
C      DATA IMACH(12) / -127 /
C      DATA IMACH(13) /  127 /
C      DATA IMACH(14) /   56 /
C      DATA IMACH(15) / -127 /
C      DATA IMACH(16) /  127 /, SANITY/987/
C
C  ***  ISSUE STOP 775 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SANITY .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      (SMALL(1) .EQ. 1117925532
     *           .AND. SMALL(2) .EQ. -448790528)
     *       .OR.     (SMALL(2) .EQ. 1117925532
     *           .AND. SMALL(1) .EQ. -448790528)) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*               *** VAX WITH D_FLOATING ***
               IMACH(10) = 2
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
            IMACH(11) = IMACH(14)
            IMACH(12) = IMACH(15)
            IMACH(13) = IMACH(16)
         ELSE
            RMACH = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*               *** IEEE ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -125
               IMACH(13) = 128
               IMACH(14) = 53
               IMACH(15) = -1021
               IMACH(16) = 1024
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*               *** VAX ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -127
               IMACH(13) = 127
               IMACH(14) = 56
               IMACH(15) = -127
               IMACH(16) = 127
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*               *** IBM MAINFRAME ***
               IMACH(10) = 16
               IMACH(11) = 6
               IMACH(12) = -64
               IMACH(13) = 63
               IMACH(14) = 14
               IMACH(15) = -64
               IMACH(16) = 63
               SANITY = 987
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               IMACH(10) = 2
               IMACH(11) = 24
               IMACH(12) = -128
               IMACH(13) = 127
               IMACH(14) = 53
               IMACH(15) = -1024
               IMACH(16) = 1023
               SANITY = 987
            ELSE
*              CRAY1 = 4617762693716115456
               CRAY1 = 4617762
               CRAY1 = 1000000*CRAY1 + 693716
               CRAY1 = 1000000*CRAY1 + 115456
               IF (SMALL(1) .NE. CRAY1) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
               IMACH(1) = 5
               IMACH(2) = 6
               IMACH(3) = 102
               IMACH(4) = 6
               IMACH(5) = 64
               IMACH(6) = 8
               IMACH(7) = 2
               IMACH(8) = 63
*              IMACH(9) = 9223372036854775807
               IMACH(9) = 9223372
               IMACH(9) = 1000000*IMACH(9) + 36854
               IMACH(9) = 1000000*IMACH(9) + 775807
               IMACH(10) = 2
               IMACH(11) = 47
               IMACH(12) = -8189
               IMACH(13) = 8190
               IMACH(14) = 94
               IMACH(15) = -8099
               IMACH(16) = 8190
               SANITY = 987
               GO TO 10
               END IF
            END IF
         IMACH( 1) = 5
         IMACH( 2) = 6
         IMACH( 3) = 7
         IMACH( 4) = 6
         IMACH( 5) = 32
         IMACH( 6) = 4
         IMACH( 7) = 2
         IMACH( 8) = 31
         IMACH( 9) = 2147483647
         SANITY = 987
         END IF
C/6S
C9010 FORMAT(/47H Adjust autodoubled I1MACH by uncommenting data/
C    * 52H statements appropriate for your machine and setting/
C    * 46H IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.)
C9020 FORMAT(/46H Adjust I1MACH by uncommenting data statements/
C    * 30H appropriate for your machine.)
C/7S
 9010 FORMAT(/' Adjust autodoubled I1MACH by uncommenting data'/
     * ' statements appropriate for your machine and setting'/
     * ' IMACH(I) = IMACH(I+3) for I = 11, 12, and 13.')
 9020 FORMAT(/' Adjust I1MACH by uncommenting data statements'/
     * ' appropriate for your machine.')
C/
 10   IF (I .LT. 1  .OR.  I .GT. 16) GO TO 30
C
      I1MACH = IMACH(I)
C/6S
C/7S
      IF (I .EQ. 6) I1MACH = 1
C/
      RETURN
C
 30   WRITE(*,*) 'I1MACH(I): I =',I,' is out of bounds.'
C
*     CALL FDUMP
C
      STOP
C
* /* C source for I1MACH -- remove the * in column 1 */
* /* Note that some values may need changing -- see the comments below. */
*#include <stdio.h>
*#include <float.h>
*#include <limits.h>
*#include <math.h>
*
*long i1mach_(long *i)
*{
*	switch(*i){
*	  case 1:  return 5;	/* standard input  unit -- may need changing */
*	  case 2:  return 6;	/* standard output unit -- may need changing */
*	  case 3:  return 7;	/* standard punch  unit -- may need changing */
*	  case 4:  return 0;	/* standard error  unit -- may need changing */
*	  case 5:  return 32;	/* bits per integer -- may need changing */
*	  case 6:  return 1;	/* Fortran 77 value: 1 character */
*	  			/*    per character storage unit */
*	  case 7:  return 2;	/* base for integers -- may need changing */
*	  case 8:  return 31;	/* digits of integer base -- may need changing */
*	  case 9:  return LONG_MAX;
*	  case 10: return FLT_RADIX;
*	  case 11: return FLT_MANT_DIG;
*	  case 12: return FLT_MIN_EXP;
*	  case 13: return FLT_MAX_EXP;
*	  case 14: return DBL_MANT_DIG;
*	  case 15: return DBL_MIN_EXP;
*	  case 16: return DBL_MAX_EXP;
*	  }
*
*	fprintf(stderr, "invalid argument: i1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END
      integer function i8save(isw,ivalue,set)
c
c  if (isw = 1) i8save returns the current error number and
c               sets it to ivalue if set = .true. .
c
c  if (isw = 2) i8save returns the current recovery switch and
c               sets it to ivalue if set = .true. .
c
      logical set
c
      integer iparam(2)
c  iparam(1) is the error number and iparam(2) is the recovery switch.
c
c  start execution error free and with recovery turned off.
c
      data iparam(1) /0/,  iparam(2) /2/
c
      i8save=iparam(isw)
      if (set) iparam(isw)=ivalue
c
      return
c
      end
      function initds (dos, nos, eta)
c june 1977 edition.   w. fullerton, c3, los alamos scientific lab.
c
c initialize the double precision orthogonal series dos so that initds
c is the number of terms needed to insure the error is no larger than
c eta.  ordinarily eta will be chosen to be one-tenth machine precision.
c
c             input arguments --
c dos    dble prec array of nos coefficients in an orthogonal series.
c nos    number of coefficients in dos.
c eta    requested accuracy of series.
c
      double precision dos(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinitds  number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(sngl(dos(i)))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinitds  eta may be too small, 28,
     1  1, 2)
      initds = i
c
      return
      end
C-----------------------------------------------------------------------
      function inits (os, nos, eta)
C-----------------------------------------------------------------------
c april 1977 version.  w. fullerton, c3, los alamos scientific lab.
c
c initialize the orthogonal series so that inits is the number of terms
c needed to insure the error is no larger than eta.  ordinarily, eta
c will be chosen to be one-tenth machine precision.
c
c             input arguments --
c os     array of nos coefficients in an orthogonal series.
c nos    number of coefficients in os.
c eta    requested accuracy of series.
c
      dimension os(nos)
c
      if (nos.lt.1) call seteru (
     1  35hinits   number of coefficients lt 1, 35, 2, 2)
c
      err = 0.
      do 10 ii=1,nos
        i = nos + 1 - ii
        err = err + abs(os(i))
        if (err.gt.eta) go to 20
 10   continue
c
 20   if (i.eq.nos) call seteru (28hinits   eta may be too small, 28,
     1  1, 2)
      inits = i
c
      return
      end
C-----------------------------------------------------------------------
      REAL FUNCTION R1MACH(I)
C-----------------------------------------------------------------------
      INTEGER I
C
C  SINGLE-PRECISION MACHINE CONSTANTS
C
C  R1MACH(1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
C
C  R1MACH(2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
C
C  R1MACH(3) = B**(-T), THE SMALLEST RELATIVE SPACING.
C
C  R1MACH(4) = B**(1-T), THE LARGEST RELATIVE SPACING.
C
C  R1MACH(5) = LOG10(B)
C
C  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
C  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
C  REMOVING THE C FROM COLUMN 1.
C
C  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), THE FIRST
C  SET OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
C
C  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
C  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
C  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
C  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
C  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
C  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
C
C  COMMENTS JUST BEFORE THE END STATEMENT (LINES STARTING WITH *)
C  GIVE C SOURCE FOR R1MACH.
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
      INTEGER CRAY1, SC
      COMMON /D8MACH/ CRAY1
C/6S
C/7S
      SAVE SMALL, LARGE, RIGHT, DIVER, LOG10, SC
C/
      REAL RMACH(5)
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
C     3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
C     PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).
C
       DATA SMALL(1) /     8388608 /
       DATA LARGE(1) /  2139095039 /
       DATA RIGHT(1) /   864026624 /
       DATA DIVER(1) /   872415232 /
       DATA LOG10(1) /  1050288283 /, SC/987/
C
C     MACHINE CONSTANTS FOR AMDAHL MACHINES.
C
C      DATA SMALL(1) /    1048576 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  990904320 /
C      DATA DIVER(1) / 1007681536 /
C      DATA LOG10(1) / 1091781651 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
C
C      DATA RMACH(1) / Z400800000 /
C      DATA RMACH(2) / Z5FFFFFFFF /
C      DATA RMACH(3) / Z4E9800000 /
C      DATA RMACH(4) / Z4EA800000 /
C      DATA RMACH(5) / Z500E730E8 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS.
C
C      DATA RMACH(1) / O1771000000000000 /
C      DATA RMACH(2) / O0777777777777777 /
C      DATA RMACH(3) / O1311000000000000 /
C      DATA RMACH(4) / O1301000000000000 /
C      DATA RMACH(5) / O1157163034761675 /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN4 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / 00564000000000000000B /
C      DATA RMACH(2) / 37767777777777777776B /
C      DATA RMACH(3) / 16414000000000000000B /
C      DATA RMACH(4) / 16424000000000000000B /
C      DATA RMACH(5) / 17164642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR FTN5 ON THE CDC 6000/7000 SERIES.
C
C      DATA RMACH(1) / O"00564000000000000000" /
C      DATA RMACH(2) / O"37767777777777777776" /
C      DATA RMACH(3) / O"16414000000000000000" /
C      DATA RMACH(4) / O"16424000000000000000" /
C      DATA RMACH(5) / O"17164642023241175720" /, SC/987/
C
C     MACHINE CONSTANTS FOR CONVEX C-1.
C
C      DATA RMACH(1) / '00800000'X /
C      DATA RMACH(2) / '7FFFFFFF'X /
C      DATA RMACH(3) / '34800000'X /
C      DATA RMACH(4) / '35000000'X /
C      DATA RMACH(5) / '3F9A209B'X /, SC/987/
C
C     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
C
C      DATA RMACH(1) / 200034000000000000000B /
C      DATA RMACH(2) / 577767777777777777776B /
C      DATA RMACH(3) / 377224000000000000000B /
C      DATA RMACH(4) / 377234000000000000000B /
C      DATA RMACH(5) / 377774642023241175720B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
C
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING LINE -
C     STATIC RMACH(5)
C
C      DATA SMALL/20K,0/,LARGE/77777K,177777K/
C      DATA RIGHT/35420K,0/,DIVER/36020K,0/
C      DATA LOG10/40423K,42023K/, SC/987/
C
C     MACHINE CONSTANTS FOR THE HARRIS SLASH 6 AND SLASH 7.
C
C      DATA SMALL(1),SMALL(2) / '20000000, '00000201 /
C      DATA LARGE(1),LARGE(2) / '37777777, '00000177 /
C      DATA RIGHT(1),RIGHT(2) / '20000000, '00000352 /
C      DATA DIVER(1),DIVER(2) / '20000000, '00000353 /
C      DATA LOG10(1),LOG10(2) / '23210115, '00000377 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE HONEYWELL DPS 8/70 SERIES.
C
C      DATA RMACH(1) / O402400000000 /
C      DATA RMACH(2) / O376777777777 /
C      DATA RMACH(3) / O714400000000 /
C      DATA RMACH(4) / O716400000000 /
C      DATA RMACH(5) / O776464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9 AND THE SEL SYSTEMS 85/86.
C
C      DATA RMACH(1) / Z00100000 /
C      DATA RMACH(2) / Z7FFFFFFF /
C      DATA RMACH(3) / Z3B100000 /
C      DATA RMACH(4) / Z3C100000 /
C      DATA RMACH(5) / Z41134413 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE INTERDATA 8/32
C     WITH THE UNIX SYSTEM FORTRAN 77 COMPILER.
C
C     FOR THE INTERDATA FORTRAN VII COMPILER REPLACE
C     THE Z'S SPECIFYING HEX CONSTANTS WITH Y'S.
C
C      DATA RMACH(1) / Z'00100000' /
C      DATA RMACH(2) / Z'7EFFFFFF' /
C      DATA RMACH(3) / Z'3B100000' /
C      DATA RMACH(4) / Z'3C100000' /
C      DATA RMACH(5) / Z'41134413' /, SC/987/
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR).
C
C      DATA RMACH(1) / "000400000000 /
C      DATA RMACH(2) / "377777777777 /
C      DATA RMACH(3) / "146400000000 /
C      DATA RMACH(4) / "147400000000 /
C      DATA RMACH(5) / "177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1) /    8388608 /
C      DATA LARGE(1) / 2147483647 /
C      DATA RIGHT(1) /  880803840 /
C      DATA DIVER(1) /  889192448 /
C      DATA LOG10(1) / 1067065499 /, SC/987/
C
C      DATA RMACH(1) / O00040000000 /
C      DATA RMACH(2) / O17777777777 /
C      DATA RMACH(3) / O06440000000 /
C      DATA RMACH(4) / O06500000000 /
C      DATA RMACH(5) / O07746420233 /, SC/987/
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRANS SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C      DATA SMALL(1),SMALL(2) /   128,     0 /
C      DATA LARGE(1),LARGE(2) / 32767,    -1 /
C      DATA RIGHT(1),RIGHT(2) / 13440,     0 /
C      DATA DIVER(1),DIVER(2) / 13568,     0 /
C      DATA LOG10(1),LOG10(2) / 16282,  8347 /, SC/987/
C
C      DATA SMALL(1),SMALL(2) / O000200, O000000 /
C      DATA LARGE(1),LARGE(2) / O077777, O177777 /
C      DATA RIGHT(1),RIGHT(2) / O032200, O000000 /
C      DATA DIVER(1),DIVER(2) / O032400, O000000 /
C      DATA LOG10(1),LOG10(2) / O037632, O020233 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
C
C      DATA SMALL(1) / $00800000 /
C      DATA LARGE(1) / $7F7FFFFF /
C      DATA RIGHT(1) / $33800000 /
C      DATA DIVER(1) / $34000000 /
C      DATA LOG10(1) / $3E9A209B /, SC/987/
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES.
C
C      DATA RMACH(1) / O000400000000 /
C      DATA RMACH(2) / O377777777777 /
C      DATA RMACH(3) / O146400000000 /
C      DATA RMACH(4) / O147400000000 /
C      DATA RMACH(5) / O177464202324 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX UNIX F77 COMPILER.
C
C      DATA SMALL(1) /       128 /
C      DATA LARGE(1) /    -32769 /
C      DATA RIGHT(1) /     13440 /
C      DATA DIVER(1) /     13568 /
C      DATA LOG10(1) / 547045274 /, SC/987/
C
C     MACHINE CONSTANTS FOR THE VAX-11 WITH
C     FORTRAN IV-PLUS COMPILER.
C
C      DATA RMACH(1) / Z00000080 /
C      DATA RMACH(2) / ZFFFF7FFF /
C      DATA RMACH(3) / Z00003480 /
C      DATA RMACH(4) / Z00003500 /
C      DATA RMACH(5) / Z209B3F9A /, SC/987/
C
C     MACHINE CONSTANTS FOR VAX/VMS VERSION 2.2.
C
C      DATA RMACH(1) /       '80'X /
C      DATA RMACH(2) / 'FFFF7FFF'X /
C      DATA RMACH(3) /     '3480'X /
C      DATA RMACH(4) /     '3500'X /
C      DATA RMACH(5) / '209B3F9A'X /, SC/987/
C
C  ***  ISSUE STOP 777 IF ALL DATA STATEMENTS ARE COMMENTED...
      IF (SC .NE. 987) THEN
*        *** CHECK FOR AUTODOUBLE ***
         SMALL(2) = 0
         RMACH(1) = 1E13
         IF (SMALL(2) .NE. 0) THEN
*           *** AUTODOUBLED ***
            IF (      SMALL(1) .EQ. 1117925532
     *          .AND. SMALL(2) .EQ. -448790528) THEN
*              *** IEEE BIG ENDIAN ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2146435071
               LARGE(2) = -1
               RIGHT(1) = 1017118720
               RIGHT(2) = 0
               DIVER(1) = 1018167296
               DIVER(2) = 0
               LOG10(1) = 1070810131
               LOG10(2) = 1352628735
            ELSE IF ( SMALL(2) .EQ. 1117925532
     *          .AND. SMALL(1) .EQ. -448790528) THEN
*              *** IEEE LITTLE ENDIAN ***
               SMALL(2) = 1048576
               SMALL(1) = 0
               LARGE(2) = 2146435071
               LARGE(1) = -1
               RIGHT(2) = 1017118720
               RIGHT(1) = 0
               DIVER(2) = 1018167296
               DIVER(1) = 0
               LOG10(2) = 1070810131
               LOG10(1) = 1352628735
            ELSE IF ( SMALL(1) .EQ. -2065213935
     *          .AND. SMALL(2) .EQ. 10752) THEN
*              *** VAX WITH D_FLOATING ***
               SMALL(1) = 128
               SMALL(2) = 0
               LARGE(1) = -32769
               LARGE(2) = -1
               RIGHT(1) = 9344
               RIGHT(2) = 0
               DIVER(1) = 9472
               DIVER(2) = 0
               LOG10(1) = 546979738
               LOG10(2) = -805796613
            ELSE IF ( SMALL(1) .EQ. 1267827943
     *          .AND. SMALL(2) .EQ. 704643072) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               SMALL(2) = 0
               LARGE(1) = 2147483647
               LARGE(2) = -1
               RIGHT(1) = 856686592
               RIGHT(2) = 0
               DIVER(1) = 873463808
               DIVER(2) = 0
               LOG10(1) = 1091781651
               LOG10(2) = 1352628735
            ELSE
               WRITE(*,9010)
               STOP 777
               END IF
         ELSE
            RMACH(1) = 1234567.
            IF (SMALL(1) .EQ. 1234613304) THEN
*              *** IEEE ***
               SMALL(1) = 8388608
               LARGE(1) = 2139095039
               RIGHT(1) = 864026624
               DIVER(1) = 872415232
               LOG10(1) = 1050288283
            ELSE IF (SMALL(1) .EQ. -1271379306) THEN
*              *** VAX ***
               SMALL(1) = 128
               LARGE(1) = -32769
               RIGHT(1) = 13440
               DIVER(1) = 13568
               LOG10(1) = 547045274
            ELSE IF (SMALL(1) .EQ. 1175639687) THEN
*              *** IBM MAINFRAME ***
               SMALL(1) = 1048576
               LARGE(1) = 2147483647
               RIGHT(1) = 990904320
               DIVER(1) = 1007681536
               LOG10(1) = 1091781651
            ELSE IF (SMALL(1) .EQ. 1251390520) THEN
*              *** CONVEX C-1 ***
               SMALL(1) = 8388608
               LARGE(1) = 2147483647
               RIGHT(1) = 880803840
               DIVER(1) = 889192448
               LOG10(1) = 1067065499
            ELSE
*              CRAY1 = 4617762693716115456
               CRAY1 = 4617762
               CRAY1 = 1000000*CRAY1 + 693716
               CRAY1 = 1000000*CRAY1 + 115456
               IF (SMALL(1) .NE. CRAY1) THEN
                  WRITE(*,9020)
                  STOP 777
                  END IF
*              *** CRAY 1, XMP, 2, AND 3 ***
*              SMALL(1) = 2306828171632181248
               SMALL(1) = 2306828
               SMALL(1) = 1000000*SMALL(1) + 171632
               SMALL(1) = 1000000*SMALL(1) + 181248
*              LARGE(1) = 6917247552664371198
               LARGE(1) = 6917247
               LARGE(1) = 1000000*LARGE(1) + 552664
               LARGE(1) = 1000000*LARGE(1) + 371198
*              RIGHT(1) = 4598878906987053056
               RIGHT(1) = 4598878
               RIGHT(1) = 1000000*RIGHT(1) + 906987
               RIGHT(1) = 1000000*RIGHT(1) + 053056
*              DIVER(1) = 4599160381963763712
               DIVER(1) = 4599160
               DIVER(1) = 1000000*DIVER(1) + 381963
               DIVER(1) = 1000000*DIVER(1) + 763712
*              LOG10(1) = 4611574008272714704
               LOG10(1) = 4611574
               LOG10(1) = 1000000*LOG10(1) + 008272
               LOG10(1) = 1000000*LOG10(1) + 714704
               END IF
            END IF
         SC = 987
         END IF
C
C  ***  ISSUE STOP 776 IF ALL DATA STATEMENTS ARE OBVIOUSLY WRONG...
      IF (RMACH(4) .GE. 1.0) STOP 776
*C/6S
*C     IF (I .LT. 1  .OR.  I .GT. 5)
*C    1   CALL SETERR(24HR1MACH - I OUT OF BOUNDS,24,1,2)
*C/7S
*      IF (I .LT. 1  .OR.  I .GT. 5)
*     1   CALL SETERR('R1MACH - I OUT OF BOUNDS',24,1,2)
*C/
C
      IF (I .LT. 1 .OR. I .GT. 5) THEN
         WRITE(*,*) 'R1MACH(I): I =',I,' is out of bounds.'
         STOP
         END IF
      R1MACH = RMACH(I)
      RETURN
C/6S
C9010 FORMAT(/42H Adjust autodoubled R1MACH by getting data/
C    *42H appropriate for your machine from D1MACH.)
C9020 FORMAT(/46H Adjust R1MACH by uncommenting data statements/
C    *30H appropriate for your machine.)
C/7S
 9010 FORMAT(/' Adjust autodoubled R1MACH by getting data'/
     *' appropriate for your machine from D1MACH.')
 9020 FORMAT(/' Adjust R1MACH by uncommenting data statements'/
     *' appropriate for your machine.')
C/
C
* /* C source for R1MACH -- remove the * in column 1 */
*#include <stdio.h>
*#include <float.h>
*#include <math.h>
*
*float r1mach_(long *i)
*{
*	switch(*i){
*	  case 1: return FLT_MIN;
*	  case 2: return FLT_MAX;
*	  case 3: return FLT_EPSILON/FLT_RADIX;
*	  case 4: return FLT_EPSILON;
*	  case 5: return log10(FLT_RADIX);
*	  }
*
*	fprintf(stderr, "invalid argument: r1mach(%ld)\n", *i);
*	exit(1);
*	return 0; /* for compilers that complain of missing return values */
*	}
      END
C-----------------------------------------------------------------------
      subroutine r9upak (x, y, n)
C-----------------------------------------------------------------------
c august 1980 portable edition.  w. fullerton, los alamos scientific lab
c
c unpack floating point number x so that x = y * 2.0**n, where
c 0.5 .le. abs(y) .lt. 1.0 .
c
      absx = abs(x)
      n = 0
      y = 0.0
      if (x.eq.0.0) return
c
 10   if (absx.ge.0.5) go to 20
      n = n - 1
      absx = absx*2.0
      go to 10
c
 20   if (absx.lt.1.0) go to 30
      n = n + 1
      absx = absx*0.5
      go to 20
c
 30   y = sign (absx, x)
      return
c
      end
C-----------------------------------------------------------------------
      subroutine s88fmt( n, w, ifmt )
C-----------------------------------------------------------------------
c
c  s88fmt  replaces ifmt(1), ... , ifmt(n) with
c  the characters corresponding to the n least significant
c  digits of w.
c
      integer n,w,ifmt(n)
c
      integer nt,wt,digits(10)
c
      data digits( 1) / 1h0 /
      data digits( 2) / 1h1 /
      data digits( 3) / 1h2 /
      data digits( 4) / 1h3 /
      data digits( 5) / 1h4 /
      data digits( 6) / 1h5 /
      data digits( 7) / 1h6 /
      data digits( 8) / 1h7 /
      data digits( 9) / 1h8 /
      data digits(10) / 1h9 /
c
      nt = n
      wt = w
c
 10   if (nt .le. 0) return
        idigit = mod( wt, 10 )
        ifmt(nt) = digits(idigit+1)
        wt = wt/10
        nt = nt - 1
        go to 10
c
      end
      subroutine seterr (messg, nmessg, nerr, iopt)
c
c  this version modified by w. fullerton to dump if iopt = 1 and
c  not recovering.
c  seterr sets lerror = nerr, optionally prints the message and dumps
c  according to the following rules...
c
c    if iopt = 1 and recovering      - just remember the error.
c    if iopt = 1 and not recovering  - print, dump and stop.
c    if iopt = 2                     - print, dump and stop.
c
c  input
c
c    messg  - the error message.
c    nmessg - the length of the message, in characters.
c    nerr   - the error number. must have nerr non-zero.
c    iopt   - the option. must have iopt=1 or 2.
c
c  error states -
c
c    1 - message length not positive.
c    2 - cannot have nerr=0.
c    3 - an unrecovered error followed by another error.
c    4 - bad value for iopt.
c
c  only the first 72 characters of the message are printed.
c
c  the error handler calls a subroutine named fdump to produce a
c  symbolic dump. to complete the package, a dummy version of fdump
c  is supplied, but it should be replaced by a locally written version
c  which at least gives a trace-back.
c
      integer messg(1)
      external i1mach, i8save
c
c  the unit for error messages.
c
      iwunit=i1mach(4)
c
      if (nmessg.ge.1) go to 10
c
c  a message of non-positive length is fatal.
c
        write(iwunit,9000)
 9000   format(52h1error    1 in seterr - message length not positive.)
        go to 60
c
c  nw is the number of words the message occupies.
c
 10   nw=(min0(nmessg,72)-1)/i1mach(6)+1
c
      if (nerr.ne.0) go to 20
c
c  cannot turn the error state off using seterr.
c
        write(iwunit,9001)
 9001   format(42h1error    2 in seterr - cannot have nerr=0//
     1         34h the current error message follows///)
        call e9rint(messg,nw,nerr,.true.)
        itemp=i8save(1,1,.true.)
        go to 50
c
c  set lerror and test for a previous unrecovered error.
c
 20   if (i8save(1,nerr,.true.).eq.0) go to 30
c
        write(iwunit,9002)
 9002   format(23h1error    3 in seterr -,
     1         48h an unrecovered error followed by another error.//
     2         48h the previous and current error messages follow.///)
        call eprint
        call e9rint(messg,nw,nerr,.true.)
        go to 50
c
c  save this message in case it is not recovered from properly.
c
 30   call e9rint(messg,nw,nerr,.true.)
c
      if (iopt.eq.1 .or. iopt.eq.2) go to 40
c
c  must have iopt = 1 or 2.
c
        write(iwunit,9003)
 9003   format(42h1error    4 in seterr - bad value for iopt//
     1         34h the current error message follows///)
        go to 50
c
c  test for recovery.
c
 40   if (iopt.eq.2) go to 50
c
      if (i8save(2,0,.false.).eq.1) return
c
c     call eprint
c     stop
c
 50   call eprint
c60   call fdump
 60   continue
      stop
c
      end
C-----------------------------------------------------------------------
      subroutine seteru (messg, nmessg, nerr, iopt)
C-----------------------------------------------------------------------
      common /cseter/ iunflo
      integer messg(1)
      data iunflo / 0 /
c
      if (iopt.ne.0) call seterr (messg, nmessg, nerr, iopt)
      if (iopt.ne.0) return
c
      if (iunflo.le.0) return
      call seterr (messg, nmessg, nerr, 1)
c
      return
      end
C-----------------------------------------------------------------------
