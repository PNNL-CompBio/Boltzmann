#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_check_cvodesens_errors.h"
#include "boltzmann_print_sensitivities.h"
void boltzmann_print_sensitivities(struct state_struct *state) {
  struct cvodes_params_struct *cvodes_params;
  double *vdata;
  N_Vector *ys0;
  N_Vector *dys;
  void *cvode_mem;
  double tret;
  int  ns;
  int  ny;
  int  flag;
  int  success;
  int  i;
  int  j;
  char *ode_sens_file;
  char *ode_dsens_file;
  FILE *sens_fp;
  FILE *dsens_fp;
  FILE *lfp;
  FILE *efp;
  success        = 1;
  ny             = state->nunique_molecules;
  ns             = state->number_reactions;
  cvodes_params  = state->cvodes_params;
  ode_sens_file  = state->ode_sens_file;
  ode_dsens_file = state->ode_dsens_file;
  lfp            = state->lfp;
  ys0            = cvodes_params->ys0;
  dys            = cvodes_params->dys;
  cvode_mem      = cvodes_params->cvode_mem;
  sens_fp        = fopen(ode_sens_file,"w");
  if (sens_fp) {
    /*
      Retrieve the sensitivites matrix.
    */
    flag = CVodeGetSens(cvode_mem,&tret,ys0);
    success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,"GetSens");
    if (success) {
      fprintf(sens_fp,"prop\\mlcl");
      for (j=0;j<ny;j++) {
	fprintf(sens_fp,"\t%d",j);
      }
      fprintf(sens_fp,"\n");
      for (i=0;i<ns;i++) {
	fprintf(sens_fp,"%d",i);
	vdata = NV_DATA_S(ys0[i]);
	for (j=0;j<ny;j++) {
	  fprintf(sens_fp,"\t%le",vdata[j]);
	}
	fprintf(sens_fp,"\n");
      }
    }
    fclose(sens_fp);
  } else {
    if (lfp) {
      fprintf(lfp,"Boltmzann_print_sensitivities: Error unable to open file %s for writing.\n",ode_sens_file);
      fflush(lfp);
    }
  }
  dsens_fp       = fopen(ode_dsens_file,"w");
  if (dsens_fp) {
    /*
      Retrieve the derivatives of the sensitivites matrix.
    */
    i = 1;
    flag = CVodeGetSensDky(cvode_mem,tret,i,dys);
    success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,"GetSensDky");
    if (success) {
      fprintf(dsens_fp,"prop\\mlcl");
      for (j=0;j<ny;j++) {
	fprintf(dsens_fp,"\t%d",j);
      }
      fprintf(dsens_fp,"\n");
      for (i=0;i<ns;i++) {
	fprintf(dsens_fp,"%d",i);
	vdata = NV_DATA_S(dys[i]);
	for (j=0;j<ny;j++) {
	  fprintf(dsens_fp,"\t%le",vdata[j]);
	}
	fprintf(dsens_fp,"\n");
      }
    }
    fclose(dsens_fp);
  } else {
    if (lfp) {
      fprintf(lfp,"Boltmzann_print_sensitivities: Error unable to open file %s for writing.\n",ode_dsens_file);
      fflush(lfp);
    }
  }
}
