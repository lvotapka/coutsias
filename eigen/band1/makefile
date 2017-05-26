PROG   = eigen
BIN    = ./Bin
AUX    = .
OBJS   = $(BIN)/eigen.o $(BIN)/indata.o $(BIN)/cosfou.o    \
         $(BIN)/band_pack.o $(BIN)/band_op.o $(BIN)/full_pack.o \
	 $(BIN)/opers.o  $(BIN)/cheb_pack.o  $(BIN)/ge_pack.o \
	 $(BIN)/comp_eval.o  $(BIN)/eval.o \
	 $(BIN)/poly.o  $(BIN)/jac_pack.o $(BIN)/eval_matrix.o 
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

$(BIN)/eigen.o: eigen.f limits.h  indata.h ops.h 3579.h
	$(COMP) -c eigen.f
	mv *.o $(BIN)

$(BIN)/opers.o: opers.f limits.h  indata.h 
	$(COMP) -c opers.f
	mv *.o $(BIN)

$(BIN)/cosfou.o: $(AUX)/cosfou.f
	$(COMP) -c $(AUX)/cosfou.f
	mv *.o $(BIN)

$(BIN)/indata.o: indata.f limits.h indata.h
	$(COMP) -c indata.f
	mv *.o $(BIN)

$(BIN)/band_pack.o: $(AUX)/band_pack.f
	$(COMP) -c $(AUX)/band_pack.f
	mv *.o $(BIN)

$(BIN)/band_op.o: $(AUX)/band_op.f
	$(COMP) -c $(AUX)/band_op.f
	mv *.o $(BIN)

$(BIN)/full_pack.o: $(AUX)/full_pack.f limits.h
	$(COMP) -c $(AUX)/full_pack.f
	mv *.o $(BIN)

$(BIN)/cheb_pack.o: $(AUX)/cheb_pack.f 
	$(COMP) -c $(AUX)/cheb_pack.f
	mv *.o $(BIN)

$(BIN)/poly.o: $(AUX)/poly.f 
	$(COMP) -c $(AUX)/poly.f
	mv *.o $(BIN)

$(BIN)/ge_pack.o: $(AUX)/ge_pack.f 
	$(COMP) -c $(AUX)/ge_pack.f
	mv *.o $(BIN)

$(BIN)/jac_pack.o: $(AUX)/jac_pack.f 
	$(COMP) -c $(AUX)/jac_pack.f
	mv *.o $(BIN)

$(BIN)/eval.o: $(AUX)/eval.f 
	$(COMP) -c $(AUX)/eval.f
	mv *.o $(BIN)

$(BIN)/comp_eval.o: $(AUX)/comp_eval.f 
	$(COMP) -c $(AUX)/comp_eval.f
	mv *.o $(BIN)

$(BIN)/eval_matrix.o: $(AUX)/eval_matrix.f 
	$(COMP) -c $(AUX)/eval_matrix.f
	mv *.o $(BIN)
