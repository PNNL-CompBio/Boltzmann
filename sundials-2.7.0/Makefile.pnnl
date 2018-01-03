SUNDIALS_ROOT = ${PWD}
include ../Makefile.head

#all: $(SUNDIALS_LIB_DIR)/libsundials_cvodes.so $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so $(SUNDIALS_LIB_DIR)/libsundials.so $(SUNDIALS_LIB_DIR)/libsundials.a
all: $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a $(SUNDIALS_LIB_DIR)/libsundials.a

#$(SUNDIALS_LIB_DIR)/libsundials_cvodes.so $(SUNDIALS_LIB_DIR)/libsundials_cvodes.a: src/cvodes/cvodea.c src/cvodes/cvodea_io.c src/cvodes/cvodes.c src/cvodes/cvodes_band.c src/cvodes/cvodes_bandpre.c src/cvodes/cvodes_bbdpre.c src/cvodes/cvodes_dense.c src/cvodes/cvodes_diag.c src/cvodes/cvodes_direct.c src/cvodes/cvodes_io.c src/cvodes/cvodes_lapack.c src/cvodes/cvodes_sparse.c src/cvodes/cvodes_spbcgs.c src/cvodes/cvodes_spgmr.c src/cvodes/cvodes_spils.c src/cvodes/cvodes_sptfqmr.c
$(SUNDIALS_LIB_DIR)/libsundials_cvodes.a: src/cvodes/cvodea.c src/cvodes/cvodea_io.c src/cvodes/cvodes.c src/cvodes/cvodes_band.c src/cvodes/cvodes_bandpre.c src/cvodes/cvodes_bbdpre.c src/cvodes/cvodes_dense.c src/cvodes/cvodes_diag.c src/cvodes/cvodes_direct.c src/cvodes/cvodes_io.c src/cvodes/cvodes_lapack.c src/cvodes/cvodes_sparse.c src/cvodes/cvodes_spbcgs.c src/cvodes/cvodes_spgmr.c src/cvodes/cvodes_spils.c src/cvodes/cvodes_sptfqmr.c
	cd src/cvodes && $(MAKE)

#$(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so: src/nvec_ser/fnvector_serial.c src/nvec_ser/nvector_serial.c
$(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a: src/nvec_ser/fnvector_serial.c src/nvec_ser/nvector_serial.c
	cd src/nvec_ser && $(MAKE)

#$(SUNDIALS_LIB_DIR)/libsundials.so $(SUNDIALS_LIB_DIR)/libsundials.a: src/sundials/sundials_band.c src/sundials/sundials_dense.c src/sundials/sundials_direct.c src/sundials/sundials_iterative.c src/sundials/sundials_math.c src/sundials/sundials_nvector.c src/sundials/sundials_pcg.c src/sundials/sundials_sparse.c src/sundials/sundials_spbcgs.c src/sundials/sundials_spfgmr.c src/sundials/sundials_spgmr.c src/sundials/sundials_sptfqmr.c
$(SUNDIALS_LIB_DIR)/libsundials.a: src/sundials/sundials_band.c src/sundials/sundials_dense.c src/sundials/sundials_direct.c src/sundials/sundials_iterative.c src/sundials/sundials_math.c src/sundials/sundials_nvector.c src/sundials/sundials_pcg.c src/sundials/sundials_sparse.c src/sundials/sundials_spbcgs.c src/sundials/sundials_spfgmr.c src/sundials/sundials_spgmr.c src/sundials/sundials_sptfqmr.c
	cd src/sundials && $(MAKE)

clean:
	/bin/rm -f  lib/*
	/bin/rm  -f src/cvodes/*.o
	/bin/rm  -f src/nvec_ser/*.o
	/bin/rm  -f src/sundials/*.o

