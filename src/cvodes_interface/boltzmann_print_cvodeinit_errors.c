#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_print_cvodeinit_errors.h"
void boltzmann_print_cvodeinit_errors(int flag, 
				      void *cvode_mem, 
				      struct state_struct * state) {
  /*
    Print error return messages from CVodeInit call.
    Called by: boltzmann_cvodes
    Calls:     fprintf,fflush
  */
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  if (lfp) {
    switch (flag) {
    case CV_MEM_NULL:
      fprintf(lfp,"boltzmann_cvodes: Error: cvode_mem not initialized\n");
      break;
    case CV_MEM_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error: allocation error\n");
      break;
    case CV_ILL_INPUT:
      fprintf(lfp,"boltzmann_cvodes: Invalid argument to CVodeInit\n");
      break;
    default:
      fprintf(lfp,"boltzmann_cvodes: Unknown error from CVodeInit was %d\n",
	      flag);
    } /* end switch(flag) */
    fflush(lfp);
  } /* end if (lfp) */
}
