SUNDIALS_ROOT = ${PWD}/../..
include ../../../Makefile.head

all: $(SUNDIALS_LIB_DIR)/libsundials.so $(SUNDIALS_LIB_DIR)/libsundials.a


OBJ_FILES = sundials_band.o sundials_dense.o sundials_direct.o sundials_iterative.o sundials_math.o sundials_nvector.o sundials_pcg.o sundials_sparse.o sundials_spbcgs.o sundials_spfgmr.o sundials_spgmr.o sundials_sptfqmr.o
     
$(SUNDIALS_LIB_DIR)/libsundials.so: $(OBJ_FILES)
				 gcc -shared -o $(SUNDIALS_LIB_DIR)/libsundials.so $(OBJ_FILES)

$(SUNDIALS_LIB_DIR)/libsundials.a: $(OBJ_FILES)
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_band.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_dense.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_direct.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_iterative.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_math.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_nvector.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_pcg.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_sparse.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_spbcgs.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_spfgmr.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_spgmr.o
			  $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials.a sundials_sptfqmr.o


sundials_band.o: sundials_band.c $(SUNDIALS_INC_DIR)/sundials/sundials_band.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_band.c

sundials_dense.o: sundials_dense.c $(SUNDIALS_INC_DIR)/sundials/sundials_dense.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_dense.c

sundials_direct.o: sundials_direct.c $(SUNDIALS_INC_DIR)/sundials/sundials_direct.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_direct.c

sundials_iterative.o: sundials_iterative.c $(SUNDIALS_INC_DIR)/sundials/sundials_iterative.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_iterative.c

sundials_math.o: sundials_math.c $(SUNDIALS_INC_DIR)/sundials/sundials_math.h 
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_math.c

sundials_nvector.o: sundials_nvector.c $(SUNDIALS_INC_DIR)/sundials/sundials_nvector.h 
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_nvector.c

sundials_pcg.o: sundials_pcg.c $(SUNDIALS_INC_DIR)/sundials/sundials_pcg.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_pcg.c

sundials_sparse.o: sundials_sparse.c $(SUNDIALS_INC_DIR)/sundials/sundials_sparse.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_sparse.c

sundials_spbcgs.o: sundials_spbcgs.c $(SUNDIALS_INC_DIR)/sundials/sundials_spbcgs.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_spbcgs.c

sundials_spfgmr.o: sundials_spfgmr.c $(SUNDIALS_INC_DIR)/sundials/sundials_spfgmr.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_spfgmr.c

sundials_spgmr.o: sundials_spgmr.c $(SUNDIALS_INC_DIR)/sundials/sundials_spgmr.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_spgmr.c

sundials_sptfqmr.o: sundials_sptfqmr.c $(SUNDIALS_INC_DIR)/sundials/sundials_sptfqmr.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. sundials_sptfqmr.c

clean:
	/bin/rm -f *.o *.log $(SUNDIALS_LIB_DIR)/libsundials.so $(SUNDIALS_LIB_DIR)/libsundials.a 