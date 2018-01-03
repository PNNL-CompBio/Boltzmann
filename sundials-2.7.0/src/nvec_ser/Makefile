SUNDIALS_ROOT = ${PWD}/../..
include ../../../Makefile.head

all: $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so

NVEC_SER_OBJS = fnvector_serial.o nvector_serial.o


$(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a: $(NVEC_SER_OBJS)
			       $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a fnvector_serial.o
			       $(AR) $(ARFLAGS) $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a nvector_serial.o

$(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so: $(NVEC_SER_OBJS)
	$(LD) -shared -o $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so $(NVEC_SER_OBJS)

fnvector_serial.o: fnvector_serial.c fnvector_serial.h
		   $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. fnvector_serial.c

nvector_serial.o: nvector_serial.c $(SUNDIALS_INC_DIR)/nvector/nvector_serial.h $(SUNDIALS_INC_DIR)/sundials/sundials_math.h
		  $(CC) -c $(OPT_FLAGS) -I$(SUNDIALS_INC_DIR) -I. nvector_serial.c

clean: 
	/bin/rm -f *.o *.log $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.a $(SUNDIALS_LIB_DIR)/libsundials_nvec_ser.so