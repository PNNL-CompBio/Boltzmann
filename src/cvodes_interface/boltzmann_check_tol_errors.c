#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_tol_errors.h"
int boltzmann_check_tol_errors(int flag, void *cvode_mem,
			       struct state_struct *state) {
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
      case CV_MEM_NULL:
	fprintf(lfp,"boltzmann_cvodes: Error: cvode_mem not initialized\n");
	break;
      case CV_NO_MALLOC:
	fprintf(lfp,"boltzmann_cvodes: Error: allocation error CVodeInit not called\n");
	break;
      case CV_ILL_INPUT:
	fprintf(lfp,"boltzmann_cvodes: negative tolerance for CVodeSStolerances\n");
	break;
      default:
	fprintf(lfp,"boltzmann_cvodes: Unknown error from CVodeSStolerances was %d\n",
		flag);
      } /* end switch(flag) */
      fflush(lfp);
    } /* end if (lfp) */
  }
  return(success);
}
