#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_cvspils_errors.h"
int boltzmann_check_cvspils_errors(int flag,
				   void *cvode_mem,
				   struct state_struct *state,
				   const char*routine) {
  /*
    Check for error returns from CVSpgmr, CVSpcbcg, or CVSptfqmr
  */
  FILE *lfp;
  FILE *efp;
  int success;
  int padi;
  success = 1;
  lfp     = state->lfp;
  if (flag != CVSPILS_SUCCESS) {
    success = 0;
    if (lfp) {
      switch (flag) {
      case CVSPILS_MEM_NULL:
	fprintf(lfp,"boltzmann_cvodes %s: Error cvode_mem not initialized\n",
		routine);
	break;
      case CVSPILS_MEM_FAIL:
	fprintf(lfp,"boltzmann_cvodes %s: Error memory allocation failed\n",
		routine);
	break;
      case CVSPILS_LMEM_NULL:
	fprintf(lfp,"boltmzann_cvodes %s: not initialized.\n",routine);
	break;
      case CVSPILS_ILL_INPUT:
	if (strcmp(routine,"CVSpgmr") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s invalid preconditioer type\n",
		routine);
	} else if (strcmp(routine,"CVSpbcg") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s invalid preconditioer type\n",
		routine);
	} else if (strcmp(routine,"CVSptfqmr") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s invalid preconditioer type\n",
		routine);
	} else if (strcmp(routine,"SetGSType") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s invalid gstype\n",routine);
	} else if (strcmp(routine,"SetEpsLin") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: %s eplifac is negative\n",routine);
	} 
	break;
      default:
	fprintf(lfp,"boltzmann_cvodes: Unknown error from %s was %d\n",
		routine,flag);
      } /* end switch(flag) */
    } /* end if (lfp) */
  } /* end if (flag != CV_SUCCESS) */
  return (success);
}
