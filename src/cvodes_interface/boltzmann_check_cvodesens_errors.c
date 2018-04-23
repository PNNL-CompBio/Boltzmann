#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_cvodesens_errors.h"
int boltzmann_check_cvodesens_errors(int flag,
				     void *cvode_mem,
				     struct state_struct *state,
				     const char *routine) {
  /*
    Check error retrun codes for CVodeSens* routines.
    Called by: boltzmann_cvodes_init, boltzmann_print_sensitivities
    Calls:     fprint, fflush
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
      case CV_MEM_NULL:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: CV_MEM_NULL error\n",routine);
	fflush(lfp);
	break;
      case CV_MEM_FAIL:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: CV_MEM_FAIL error\n",routine);
	fflush(lfp);
	break;
      case CV_ILL_INPUT:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: CV_ILL_INPUT error\n",routine);
	fflush(lfp);
	break;
      case CV_NO_SENS:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: CvodeSensInit not called previoisly\n",routine);
	fflush(lfp);
	break;
      case CV_BAD_DKY:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: ys is NULL\n",routine);
	fflush(lfp);
	break;
      case CV_BAD_K:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: invalid k argument \n",routine);
	fflush(lfp);
	break;
      case CV_BAD_T:
	fprintf(lfp,"boltzmann_cvodes: CVode%s: invalid t argument \n",routine);
	fflush(lfp);
	break;
      }
    }
  }
  return(success);
}
