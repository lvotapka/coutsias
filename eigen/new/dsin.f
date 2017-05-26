      double precision function dsin (x)
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
