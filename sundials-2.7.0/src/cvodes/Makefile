SUNDIALS_ROOT = ${PWD}/../..
include ../../../Makefile.head

#all: $(SUNDIALS_LIB_DIR)/libsundials_cvodes.so $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a
all: $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a



#OBJ_FILES = cvodea.o cvodea_io.o cvodes.o cvodes_band.o cvodes_bandpre.o cvodes_bbdpre.o cvodes_dense.o cvodes_diag.o cvodes_direct.o cvodes_io.o cvodes_klu.o cvodes_lapack.o cvodes_sparse.o cvodes_spbcgs.o cvodes_spgmr.o cvodes_spils.o cvodes_sptfqmr.o cvodes_superlumt.o

OBJ_FILES = cvodea.o cvodea_io.o cvodes.o cvodes_band.o cvodes_bandpre.o cvodes_bbdpre.o cvodes_dense.o cvodes_diag.o cvodes_direct.o cvodes_io.o cvodes_lapack.o cvodes_sparse.o cvodes_spbcgs.o cvodes_spgmr.o cvodes_spils.o cvodes_sptfqmr.o 
     
#$(SUNDIALS_LIB_DIR)/libsundials_cvodes.so: $(OBJ_FILES)
#				 gcc -shared -o $(SUNDIALS_LIB_DIR)/libsundials_cvodes.so $(OBJ_FILES)

$(SUNDIALS_LIB_DIR)/libsundials_cvodes.a: $(OBJ_FILES)
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodea.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodea_io.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_band.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_bandpre.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_bbdpre.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_dense.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_diag.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_direct.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_io.o
#				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_klu.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_lapack.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_sparse.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_spbcgs.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_spgmr.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_spils.o
				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_sptfqmr.o
#				$(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a cvodes_superlumt.o

cvodea.o: cvodea.c cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_types.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodea.c

cvodea_io.o: cvodea_io.c cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_types.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodea_io.c

cvodes.o: cvodes.c cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_types.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes.c

cvodes_band.o: cvodes_band.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_band.h cvodes_direct_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_band.c
	       
cvodes_bandpre.o: cvodes_bandpre.c cvodes_impl.h cvodes_bandpre_impl.h cvodes_spils_impl.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_sptfqmr.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_spbcgs.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_spgmr.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_bandpre.c

cvodes_bbdpre.o: cvodes_bbdpre.c cvodes_impl.h cvodes_bbdpre_impl.h cvodes_spils_impl.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_sptfqmr.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_spbcgs.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_spgmr.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_bbdpre.c

cvodes_dense.o: cvodes_dense.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_dense.h cvodes_direct_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_dense.c

cvodes_diag.o: cvodes_diag.c cvodes_diag_impl.h cvodes_impl.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_diag.c

cvodes_direct.o: cvodes_direct.c cvodes_direct_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_direct.c

cvodes_io.o: cvodes_io.c cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_types.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_io.c

#cvodes_klu.o: cvodes_klu.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_klu.h cvodes_impl.h cvodes_sparse_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_klu_impl.h
#	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_klu.c

cvodes_lapack.o: cvodes_lapack.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_lapack.h cvodes_direct_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_lapack.c

cvodes_sparse.o: cvodes_sparse.c cvodes_sparse_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_sparse.c

cvodes_spbcgs.o: cvodes_spbcgs.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_spbcgs.h cvodes_spils_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_spbcgs.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_spbcgs.c

cvodes_spgmr.o: cvodes_spgmr.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_spgmr.h cvodes_spils_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_spgmr.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_spgmr.c

cvodes_spils.o: cvodes_spils.c cvodes_impl.h cvodes_spils_impl.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_spils.c

cvodes_sptfqmr.o: cvodes_sptfqmr.c $(SUNDIALS_INC_DIR)/cvodes/cvodes_sptfqmr.h cvodes_spils_impl.h cvodes_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h $(SUNDIALS_INC_DIR)/sundials/sundials_sptfqmr.h
	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_sptfqmr.c

#cvodes_superlumt.o: cvodes_superlumt.c cvodes_impl.h cvodes_sparse_impl.h $(SUNDIALS_INC_DIR)/cvodes/cvodes_superlumt.h $(SUNDIALS_INC_DIR)/sundials/sundials_superlumt_impl.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
#	  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) cvodes_superlumt.c

clean:
	/bin/rm -f *.o ${SUNDIALS_LIB_DIR}/libsundials_cvodes.so ${SUNDIALS_LIB_DIR}/libsundials_cvodes.a
