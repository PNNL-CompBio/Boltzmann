#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_cvdls_errors.h"
int boltzmann_check_cvdls_errors(int flag,
				 void *cvode_mem,
				 struct state_struct *state,
				 const char *routine) {
  /*
    Called: by boltzmann_set_cvodes_linear_solver
    Calls:  fprintf,fflush
  */
  FILE *lfp;
  FILE *efp;
  int success;
  int padi;
  success = 1;
  if (flag < 0) {
    success = 0;
    lfp = state->lfp;
    if (lfp) {
      switch (flag) {
      case CVDLS_MEM_NULL:
	fprintf(lfp,"boltzmann_cvodes: %s CV_MEM_NULL error\n",routine);
	fflush(lfp);
	break;
      case CVDLS_MEM_FAIL:
	fprintf(lfp,"boltzmann_cvodes: %s memory allocate request failed\n",
		routine);
	fflush(lfp);
	break;
      case CVDLS_ILL_INPUT:
	fprintf(lfp,"boltzmann_cvodes: %s, illegal input error\n",routine);
	if (strcmp(routine,"CVDense") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s not compatable with the current Nvector module\n",routine);
	  fflush(lfp);
	} else if (strcmp(routine,"CVLapackDense") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s not compatable with the current Nvector module\n",routine);
	  fflush(lfp);
	} else if (strcmp(routine,"CVBand") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s not compatable with the current Nvector module or at least one of mlower or mupper is outsizde of [0..N-1]\n",routine);
	  fflush(lfp);
	} else if (strcmp(routine,"CVLapackBand") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s not compatable with the current Nvector module or at least one of mlower or mupper is outsizde of [0..N-1]\n",routine);
	  fflush(lfp);
	} 
	break;
	/*
      case CVDIAG_MEM_NULL:
	fprintf(lfp,"boltzmann_cvodes: %s CV_MEM_NULL error\n",routine);
	fflush(lfp);
	break;
      case CVDIAG_MEM_FAIL:
	fprintf(lfp,"boltzmann_cvodes: %s memory allocate request failed\n",
		routine);
	fflush(lfp);
	break;
      case CVDIAG_ILL_INPUT:
	fprintf(lfp,"boltzmann_cvodes: %s, illegal input error\n");
	fprintf(lfp,"boltzmann_cvodes: CVDiag is not compatable with the current Nvector module\n");
	fflush(lfp);
	break;
	*/
      } /* end switch (flag) */
    } /* end if lfp */
  } /* ennd if (flag < 0) */
  return(success);
}
