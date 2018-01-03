#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "boltzmann_check_cvdls_errors.h"
/*
#include "boltzmann_check_cvsls_errors.h"
*/
#include "boltzmann_check_cvspils_errors.h"
#include "boltzmann_set_cvodes_linear_solver.h"
int boltzmann_set_cvodes_linear_solver(struct state_struct *state,
				       void *cvode_mem) {

  /*
    Called by: boltzmann_cvodes_init
    Calls:     CVDense,CVLapackDense,CVBand,CVLapackBand,CVDiag,
               CVKLU,CVSuperLUMT,CVSpgmr,CVSpbcg,CVSptfqmr,
	       boltzmann_check_cvdls_errors,
	       boltzmann_check_cvsls_errors,
	       boltzmann_check_cvspils_errors
  */
  struct cvodes_params_struct *cvodes_params;
  int linear_solver_method;
  int flag;

  int ny;
  int mupper;

  int mlower;
  int pretype;

  /*
  int nnz;
  int sparsetype;
  */

  int maxl;
  int success;

  success = 1;
  ny = state->nunique_molecules;
  cvodes_params = state->cvodes_params;
  linear_solver_method = cvodes_params->linear_solver_method;
  switch(linear_solver_method) {
  case 0: 
    flag = CVDense(cvode_mem,ny);
    success = boltzmann_check_cvdls_errors(flag,cvode_mem,state,"CVDense");
    break;
    /*
  case 1:
    flag = CVLapackDense(cvode_mem,ny);
    success = boltzmann_check_cvdls_errors(flag,cvode_mem,state,"CVLapackDense");
    break;
    */
  case 2:
    mupper = cvodes_params->mupper;
    mlower = cvodes_params->mlower;
    flag = CVBand(cvode_mem,ny,mupper,mlower);
    success = boltzmann_check_cvdls_errors(flag,cvode_mem,state,"CVBand");
    break;
    /*
  case 3:
    mupper = cvodes_params->mupper;
    mlower = cvodes_params->mlower;
    flag = CVLapackBand(cvode_mem,ny,mupper,mlower);
    success = boltzmann_check_cvdls_errors(flag,cvode_mem,state,"CVLapackBand");
    break;
    */
  case 4:
    flag = CVDiag(cvode_mem);
    success = boltzmann_check_cvdls_errors(flag,cvode_mem,state,"CVDiag");
    break;
    /*
  case 5:
    nnz = cvodes_params->nnz;
    sparsetype = cvodes_params->sparsetype;
    flag = CVKLU(cvode_mem,ny,nnz,sparsetype);
    success = boltzmann_check_cvsls_errors(flag,cvode_mem,state,"CVKLU");
    break;
    */
    /*
  case 6:
    nnz = cvodes_params->nnz;
    num_threads = cvodes_params->num_threads;
    flag = CVSuperLUMT(cvode_mem,num_threads,ny,nnz);
    success = boltzmann_check_cvsls_errors(flag,cvode_mem,state,"CVSuperLUMT");
    break;
    */
  case 7:
    pretype = cvodes_params->pretype;
    maxl    = cvodes_params->maxl;
    flag    = CVSpgmr(cvode_mem,pretype,maxl);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"CVSpgmr");
    break;
  case 8:
    pretype = cvodes_params->pretype;
    maxl    = cvodes_params->maxl;
    flag = CVSpbcg(cvode_mem,pretype,maxl);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"CVSpbcg");
    break;
  case 9:
    pretype = cvodes_params->pretype;
    maxl    = cvodes_params->maxl;
    flag = CVSptfqmr(cvode_mem,pretype,maxl);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"CVSpbcg");
    break;
  }
  return(success);
}
