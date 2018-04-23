#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_check_cvodeset_errors.h"
#include "boltzmann_check_tol_errors.h"
#include "boltzmann_set_cvodes_linear_solver.h"
#include "boltzmann_check_cvspils_errors.h"
#include "boltzmann_check_cvodesens_errors.h"
#include "approximate_ys0.h"
#include "boltzmann_cvodes_psetup.h"
#include "boltzmann_cvodes_psolve.h"
#include "boltzmann_cvodes_jtimes.h"
#include "boltzmann_cvodes_init.h"

int boltzmann_cvodes_init(void *cvode_mem,struct state_struct *state, double *concs) {
  /*
    Initialize the parameters controlling cvodes ode solver.
    Assumes CVodeCreate and CVodeInit have been called.
    Called by: boltzmann_cvodes
    Calls:     N_VNew_Serial,
               N_VMake_Serial,
	       N_VCloneVectorArray_Serial
               CVodeSetErrFile,
               boltzmann_check_cvodeset_errors,
	       CVodeSStolerances,
	       boltzmann_check_tol_errors,
	       CVodeSetUserData,
	       CVodeSetMaxOrder,
	       CVodeSetMaxNumSteps,
	       CVodeSetMaxHnilWarns,
	       CVodeSetStabLimDet,
	       CVodeSetInitStep,
	       CVodeSetMinStep,
	       CVodeSetMaxStep,
	       CVodeSetStopTime,
	       CVodeSetMaxErrTestFails,
	       CVodeSetMaxNonlinIters,
	       CVodeSetMaxConvFails,
	       CVodeSetNonlinConvCoef,
       	       CVSpgmr,
	       CVSpilsSetGStype,
	       CVSpilsSetEpsLin,
	       CVSpilsSetPreconditioner,
	       CVSpilsSetJacTimesVecFn,
	       approximate_ys0,
	       boltzmann_set_cvodes_linear_solver,
	       boltzman_check_cvspils_errors,
	       boltzmann_cvodes_psetup,
	       boltzmann_cvodes_psolve,
	       boltzmann_cvodes_jtime,
	       boltzmann_check_cvodesens_errors
  */
  /*
    Direct cvode error messages to log file, lfp.
  */
  CVSensRhsFn fs;
  struct cvodes_params_struct *cvodes_params;
  N_Vector abs_tol;
  N_Vector y0;
  N_Vector *ys0;
  N_Vector *dys;
  double *abs_tol_data;
  double *p;
  double *pbar;
  double *ke;
  double *ys0v;
  double *ys0vi;
  double reltol;
  double abstol_v;
  double nlscoef;
  double eplifac;
  double hin;
  double hmin;
  double hmax;
  double tstop;
  double dqromax;

  int *plist;

  int success;
  int flag;
  
  int max_ord;
  int mxsteps;

  int maxnef;
  int maxcor;

  int maxncf;
  int pretype;

  int maxl;
  int gstype;

  int ny;
  int i;

  int use_stab_lim_det;
  int mxhnil;

  int ns;
  int ism;

  int dqtype;
  int maxcors;

  int errcons;

  FILE *lfp;
  FILE *efp;

  success = 1;
  lfp     = state->lfp;
  cvodes_params = state->cvodes_params;
  ny            = state->nunique_molecules;
  ke            = state->ke;
  if (success) {
    flag = CVodeSetErrFile(cvode_mem,lfp);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,"ErrFile");
  }
  /*
    Set up tolerances.
  */
  if (success) {
    reltol = cvodes_params->reltol; /* 1e-6 might be a reasonable default. */
    abstol_v = cvodes_params->abstol; /* 1e-12 might be a reasonable default. */
    /*
      Or we could take abs_tol_data from state vector. Hmm?
    */
    abs_tol = N_VNew_Serial(ny);
    /*
    abs_tol_data = N_VGetArrayPointer(abs_tol);
    */
    abs_tol_data = NV_DATA_S(abs_tol);
    for (i=0;i<ny;i++) {
      abs_tol_data[i] = abstol_v;
    }
    flag = CVodeSStolerances(cvode_mem,reltol,abstol_v);
    success = boltzmann_check_tol_errors(flag,cvode_mem,state);
  } /* end if (success) */
  /*
    Set the user data pointer in our case the state struct.
  */
  if (success) {
    flag = CVodeSetUserData(cvode_mem,(void*)state);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,"Userdata");
  } /* end if (success) */
  /*
    Set the order of the lmm to be used.
  */
  if (success) {
    max_ord = cvodes_params->max_ord;
    flag = CVodeSetMaxOrd(cvode_mem,max_ord);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,"MaxOrder");
  }
  /*
    Set the maximum number of internal steps (per cvode call).
    Default 500.
  */
  if (success) {
    mxsteps = cvodes_params->mxsteps;
    flag = CVodeSetMaxNumSteps(cvode_mem,mxsteps);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxNumSteps");
  }
  /*
    Set the maximum number of warnings for zero step size h. Default 10.
  */
  if (success) {
    mxhnil = cvodes_params->mxhnil;
    flag   = CVodeSetMaxHnilWarns(cvode_mem,mxhnil);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxHnilWarns");
  }
  /*
    Set the use stabilit limit detection flag, default is 0. 
  */
  if (success) {
    use_stab_lim_det = cvodes_params->use_stab_lim_det;
    if (use_stab_lim_det) {
      flag = CVodeSetStabLimDet(cvode_mem,TRUE);
    } else {
      flag = CVodeSetStabLimDet(cvode_mem,TRUE);
    }
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "StabLimDet");
  }
  /* 
    Set the initial step size, 0.0 = use default
  */
  if (success) {
    hin = cvodes_params->hin;
    flag = CVodeSetInitStep(cvode_mem,hin);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "InitStep");
  }
  /* 
    Set the minimum step size.
  */
  if (success) {
    hmin = cvodes_params->hmin;
    flag = CVodeSetMinStep(cvode_mem,hmin);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MinStep");
  }
  /* 
    Set the maximum step size.
  */
  if (success) {
    /*
    hmax = cvodes_params->hmax;
    */
    hmax = state->ode_t_final;
    flag = CVodeSetMaxStep(cvode_mem,hmax);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxStep");
  }
  /*
    Set the stop time
  */
  if (success) {
    /*
    tstop = cvodes_params->tstop;
    */
    tstop = state->ode_t_final;
    flag = CVodeSetStopTime(cvode_mem,tstop);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "StopTime");
  }
  /*
    Set the maximum number of error test failures per step. Default 7.
  */
  if (success) {
    maxnef = cvodes_params->maxnef;
    flag = CVodeSetMaxErrTestFails(cvode_mem,maxnef);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxErrTestFails");
  }
  /*
    Set the maximum number of nonlinear sovler iterations per step, Default 3
  */
  if (success) {
    maxcor = cvodes_params->maxcor;
    flag = CVodeSetMaxNonlinIters(cvode_mem,maxcor);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxNonlinIters");
  }
  /*
    Set the maximum number of nonlinear solver convergence failures allowed
    per step. Default 10.
  */
  if (success) {
    maxncf = cvodes_params->maxncf;
    flag   = CVodeSetMaxConvFails(cvode_mem,maxncf);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "MaxConvFails");
  }
  /*
    Set the safety factor in the nonlinear convergences test.
    Default value is .1
  */
  if (success) {
    nlscoef = cvodes_params->nlscoef;
    flag    = CVodeSetNonlinConvCoef(cvode_mem,nlscoef);
    success = boltzmann_check_cvodeset_errors(flag,cvode_mem,state,
					      "NonlinConvCoef");
  }
  /*
    Specify the linear solver.
  */
  if (success) {
    pretype = cvodes_params->pretype;
    maxl    = cvodes_params->maxl;
    flag = CVSpgmr(cvode_mem,pretype,maxl);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"CVSpgmr");
  } /* end if (success) */
  /*
    We need to specify the preconditioner psetup and psolve functions, and
    the jacobian-times vector function.
  */
  /*
    Specify modified Gram-Schmidt
  */
  if (success) {
    gstype = cvodes_params->gstype;
    flag = CVSpilsSetGSType(cvode_mem,gstype);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"SetGSType");
  }
  /*
    Specify factor by which Krylov linear solvers covergence test is reduced
    from Newton iteration test constant. Default 0.05 
  */
  if (success) {
    eplifac = cvodes_params->eplifac;
    flag = CVSpilsSetEpsLin(cvode_mem,eplifac);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,"SetEpsLin");
  }
  /*
    Set the preconditioner routine.
  if (success) {
    flag = CVSpilsSetPreconditioner(cvode_mem,boltzmann_cvodes_psetup,
                                    boltzmann_cvodes_psolve);
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,
                                             "CVSpilsSetPreconditioner");
  }
  */
  /*
    Set the Jacobian vector multpily function.
  */
  if (success) {
    if (state->ode_jacobian_choice == 0) {
      flag = CVSpilsSetJacTimesVecFn(cvode_mem,NULL);
    } else {
      flag = CVSpilsSetJacTimesVecFn(cvode_mem,boltzmann_cvodes_jtimes); 
    }
    success = boltzmann_check_cvspils_errors(flag,cvode_mem,state,
					     "CVSpilsSetJacTimesVecFn");
  }
  if (success) {
    if ((state->compute_sensitivities == 1)  &&
	(state->ode_solver_choice == 1)) {
      /*
	Define the sensitivity problem.
      */
      /*
	Set ns the number of sensitivity parameters,
	= number_reactions.
      */
      ns  = cvodes_params->ns;
      y0  = cvodes_params->y0;
      /*
	Set the sensitivity initial conditions.
      */
      /*
	Allocate space for the ys0 array of ns N_vectors.
      */
      ys0  = N_VCloneVectorArray_Serial(ns,y0);
      if (ys0 == NULL) {
	success =0;
	if (lfp) {
	  fprintf(lfp,"boltzman_cvodes_init: Error allocating space for ys0\n");
	  fflush(lfp);
	}
      }
      if (success) {
	cvodes_params->ys0 = ys0;
	dys = N_VCloneVectorArray_Serial(ns,y0);
	if (dys == NULL) {
	  success =0;
	  if (lfp) {
	    fprintf(lfp,"boltzman_cvodes_init: Error allocating space for ys0\n");
	    fflush(lfp);
	  }
	}
      }
      if (success) {
	cvodes_params->dys = dys;
	ys0v = cvodes_params->ys0v;
	/*
	  Initialize the sensitivy vectors.
	*/
	/*
	  Fill the ys0v array, with dy[i]/dke[j]
	*/
	approximate_ys0(state,concs);
	ys0vi = ys0v;
	for (i=0;i<ns;i++) {
	  ys0[i] = N_VMake_Serial(ny,ys0vi);
	  ys0vi += ny; /* Caution Address artihmetic here. */
	}
	/*
	  Activate the forward sensitivity computation.
	*/
	ism  = CV_STAGGERED;
	/*
	  Use the internal difference quotient sensitivity 
	  right hand side routine.
	*/
	fs   = NULL;
	flag = CVodeSensInit(cvode_mem,ns,ism,fs,ys0);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "Init");
      }
      if (success) {
	/*
	  Set the plist array to be plist[i] = i; i=0;i<number_reactions;
	*/
	plist = cvodes_params->plist;
	for (i=0;i<ns;i++) {
	  plist[i] = i;
	}
	/*
	  Set pbar[i] = abs(ke[i]).
	*/
	pbar = cvodes_params->pbar;
	for (i=0;i<ns;i++) {
	  pbar[i] = abs(ke[i]);
	}
	/*
	  Set the P used by cvodes sensitivity analyses to be state->ke.
	*/
	p    = cvodes_params->p;
	for (i=0;i<ns;i++) {
	  p[i] = ke[i];
	}
	flag = CVodeSetSensParams(cvode_mem,p,pbar,plist);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "SetSensParams");
      }
      if (success) {
	/*
	  Set the differenc quotient strategy, for now using the defaults.
	*/
	dqtype = CV_CENTERED;
	dqromax = 0.0;
	flag = CVodeSetSensDQMethod(cvode_mem,dqtype,dqromax);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "SetSensDQMethod");
      }
      if (success) {
	/*
	  Set the error control strategy to exclude the parameters.
	*/
	errcons = 0;
	flag = CVodeSetSensErrCon(cvode_mem,(booleantype)errcons);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "SetSensErrCon");
      }
      if (success) {
	/*
	  Set the mximum number of nonliear solver iterations for sensitivity 
	  variables per step to its default value.
	*/
	maxcors = 3;
	flag = CVodeSetSensMaxNonlinIters(cvode_mem,maxcors);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "SetSensMaxNonlinIters");
      }
      if (success) {
	/*
	  Set sensitiviy tolerance to be based on toleranfces supplied for
	  states variables and scaling factors in pbar.
	*/
	flag = CVodeSensEEtolerances(cvode_mem);
	success = boltzmann_check_cvodesens_errors(flag,cvode_mem,state,
						   "SensEEtolerances");
	
      }
    } /* end if (ode_solver_choice == 1) */
  } /* end if (success) */
  return(success);
}

