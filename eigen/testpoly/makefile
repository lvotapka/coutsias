PROG   = test
BIN    = ./Bin
AUX    = .
OBJS   = test.o jacobi.o gauss.o dgamma.o
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

