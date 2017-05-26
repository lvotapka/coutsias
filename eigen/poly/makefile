PROG   = test
BIN    = ./Bin
AUX    = .
OBJS   = jac.o poly.o gauss.o dgamma.o cheb.o jacini.o test.o full_pack.o\
         gen_pack.o rec.o band_op.o band_pack.o cheb_pack.o
LINK   = f77 
COMP   = f77 -c +es +e +O3
FFLAGS = 

# .SUFFIXES:
# .SUFFIXES: .f .o
# .f.o:
# 	$(COMP) $(FFLAGS) $(DFLAG) $< 

$(PROG): $(OBJS) 
	$(LINK) $(DFLAG) $(OBJS) -o $@
	chmod 771 $@

