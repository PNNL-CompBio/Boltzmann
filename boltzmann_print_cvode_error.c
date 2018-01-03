#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_print_cvode_error.h"
void boltzmann_print_cvode_error(int flag,void *cvode_mem,  struct state_struct *state) {
  /*
    Print error return messages from CVode call.
    Called by: boltzmann_cvodes
    Calls:     fprintf,fflush
  */
  struct cvodes_params_struct * cvodes_params;
  int mxsteps;
  int mxnef;
  int mxncf;
  int padi;
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  cvodes_params = (struct cvodes_params_struct *)state->cvodes_params;
  mxsteps = cvodes_params->mxsteps;
  mxnef  = cvodes_params->mxnef;
  mxncf  = cvodes_params->mxncf;
  if (lfp) {
    switch(flag) {
    case CV_MEM_NULL:
      fprintf(lfp,"boltzmann_cvodes: Error cvode_mem was null in Cvode call\n");
      break;
    case CV_NO_MALLOC:
      fprintf(lfp,"boltzmann_cvodes: Error cvode_mem Not allocate by call to CVodeInit\n");
      break;
    case CV_ILL_INPUT:
      fprintf(lfp,"boltzmann_cvodes: Error from Cvode call was CV_ILL_INPUT\n");
      break;
    case CV_TOO_CLOSE:
      fprintf(lfp,"boltzmann_cvodes: Error from Cvode call T0 and TOUT too close\n");
      break;
    case CV_TOO_MUCH_WORK:
      fprintf(lfp,"boltzmann_cvodes: Error CV_TOO_MUCH_WORK from Cvode call, mxsteps = %d\n",mxsteps);
      break;
    case CV_TOO_MUCH_ACC:
      fprintf(lfp,"boltzmann_cvodes: Error CV_TOO_MUCH_ACC (too much accuracy requeste) from Cvode call.\n");
      break;
    case CV_ERR_FAILURE:
      /*
	Probably also want to print out h and hmin
      */
      fprintf(lfp,"boltzmann_cvodes: Error CV_ERR_FAILURE from Cvode call. MXNEF = %d\n",mxnef);
      break;
    case CV_CONV_FAILURE:
      /*
	Probably also want to print out h and hmin
      */
      fprintf(lfp,"boltzmann_cvodes: Error CV_CONV_FAILURE from Cvode call. MXNCF = %d\n",mxncf);
      break;
    case CV_LINIT_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error CV_LINIT_FAIL from Cvode call.\n");
      break;
    case CV_LSETUP_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error CV_LSETUP_FAIL from Cvode call.\n");
      break;
    case CV_LSOLVE_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error CV_LSOLVE_FAIL from Cvode call.\n");
      break;
    case CV_RHSFUNC_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error CV_RHSFUNC_FAIL from Cvode call.\n");
      break;
    case CV_FIRST_RHSFUNC_ERR:
      fprintf(lfp,"boltzmann_cvodes: Error CV_FIRST_RHSFUNC_FAIL from Cvode call.\n");
      break;
    case CV_REPTD_RHSFUNC_ERR:
      fprintf(lfp,"boltzmann_cvodes: Error CV_REPTD_RHSFUNC_ERR from Cvode call.\n");
      break;
    case CV_UNREC_RHSFUNC_ERR:
      fprintf(lfp,"boltzmann_cvodes: Error CV_UNREC_RHSFUNC_ERR from Cvode call.\n");
      break;
    case CV_RTFUNC_FAIL:
      fprintf(lfp,"boltzmann_cvodes: Error CV_RTFUNC_FAIL from Cvode call.\n");
      break;
    } /* end switch(flag) */
    fflush(lfp);
  } /* end if (lfp) */
}
