#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_cvodeset_errors.h"
int boltzmann_check_cvodeset_errors(int flag,
				    void *cvode_mem,
				    struct state_struct *state,
				    const char *routine) {
  /*
    Print error messages for cvodeset calls.
    Called by: boltzmann_cvodes_init
    Calls:     strcmp,fprintf,fflush
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
	fprintf(lfp,"boltzmann_cvodes: CVodeSet%s CV_MEM_NULL error\n",routine);
	break;
      case CV_ILL_INPUT:
	fprintf(lfp,"boltzmann_cvodes: CVodeSet%s CV_ILL_INPUT error\n",routine);
	if (strcmp(routine,"MaxOrder") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: max_ord must be positive\n");
	} else if (strcmp(routine,"StabLimDet") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: stability limit detection only good for lmm = CV_BDF\n");
	} else if (strcmp(routine,"MinStep") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: hmin was < 0 or > hmax\n");
	} else if (strcmp(routine,"MaxStep") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: hmax was <= 0 or < hmin\n");
	} else if (strcmp(routine,"StopTime") == 0) {
	  fprintf(lfp,"boltzmann_cvodes: tstop not beyond current time\n");
	}
      } /* end switch(flag) */
      fflush(lfp);
    } /* end if (lfp) */
  }
  return(success);
}
