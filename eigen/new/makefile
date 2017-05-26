PROG   = gam
OBJS   = gam.o csevl.o d9gaml.o d9lgmc.o d9pak.o d9upak.o dcsevl.o dexp.o \
         dgamma.o dint.o dlog.o dsin.o dsqrt.o e9rint.o eprint.o i8save.o \
         initds.o inits.o r9upak.o s88fmt.o seterr.o seteru.o alog.o d1mach.o\
         i1mach.o r1mach.o 

FFLAGS = 
DFLAG  = 
COMP   =  f77 -c -e -fast
LINK   = f77 -O -cg89



.SUFFIXES:
.SUFFIXES: .f .o
.f.o:
	$(COMP) $(FFLAGS) $(DFLAG) $< 

$(PROG): $(OBJS) 
	$(LINK) $(DFLAG) $(OBJS) -o $@ 
	chmod 771 $@

