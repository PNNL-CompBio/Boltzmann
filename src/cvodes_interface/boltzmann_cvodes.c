#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_size_jacobian.h"
#include "boltzmann_cvodes_rhs.h"
#include "boltzmann_print_cvodeinit_errors.h"
#include "boltzmann_cvodes_init.h"
#include "boltzmann_print_cvode_error.h"
#include "boltzmann_monitor_ode.h"
#include "boltzmann_print_sensitivities.h"
#include "vec_set_constant.h"
#include  "ode_test_steady_state.h"
#include "boltzmann_cvodes.h"
int boltzmann_cvodes(struct state_struct *state, double *concs) {
  /*
    Invoke the cvodes (from sundials out of llnl) ode solver.
    Called by: ode_solver
    Calls:     N_VMake_Serial, 
               CVodeCreate,
	       CVodeInit,
	       CVode,
	       boltzmann_cvodes_rhs,
	       boltzmann_print_cvodeinit_errors,
	       boltzmann_cvodes_init,
	       boltzmann_print_cvode_error,
	       boltzmann_monitor_ode,
	       fprintf, fflush
  */
  struct cvodes_params_struct *cvodes_params;
  /*
    Allow use of ode_jacobian_choice = 0 to use ode23tb's numerical
    jacobian approximation.
  */
  struct ode23tb_params_struct *ode23tb_params;
  void *cvode_mem;
  /*
    Some kind of nvector declaration here for y0 to be built from
    concs.
  */
  N_Vector y0;
  /*
  N_Vector *ys0;
  */
  double *drfc;
  double *dfdy_a;
  double *dfdy_at;
  double *miter_m;
  double *prec_l;
  double *prec_u;
  double *prec_row;
  double *recip_diag_u;
  double *frow;
  double *srow;
  double *f0;
  double *p;
  double *rp;
  double *pbar;
  double *ys0v;
  double *dfdy;
  double *fac;
  double *thresh;
  double *fdel;
  double *fdiff;
  double *dfdy_tmp;

  int    *dfdy_ia;
  int    *dfdy_ja;

  int    *dfdy_iat;
  int    *dfdy_jat;

  int    *miter_im;
  int    *miter_jm;

  int    *prec_il;
  int    *prec_jl;

  int    *prec_iu;
  int    *prec_ju;

  int    *lindex;
  int    *uindex;

  int    *column_mask;
  int    *sindex;

  int    *plist;

  double t0;
  double tfinal;
  double tout;
  double tret;
  double delt;
  double one_l;
  double zero_l;
  double njthreshold;

  int64_t ode_rxn_view_freq;
  int64_t ode_rxn_view_step;

  int64_t ask_for;

  int success;
  int flag;

  int ny;
  int lmm;

  int iter;
  int itask;

  int print_output;
  int jacobian_choice;

  int num_steps;
  int i;

  int not_done;
  int cvode_mxsteps;

  int num_doubles;
  int num_ints;

  int nnz;
  int nnzm;

  int nnzl;
  int nnzu;

  int drfc_len;
  int ns;

  int done;
  int padi;

  FILE *lfp;
  FILE *efp;
  success = 1;
  one_l   = (int64_t)1;
  zero_l  = (int64_t)0;
  ny                = state->nunique_molecules;
  drfc_len          = state->number_molecules;
  lfp               = state->lfp;
  ns                = state->number_reactions;
  jacobian_choice   = state->ode_jacobian_choice;
  ode_rxn_view_freq = state->ode_rxn_view_freq;
  print_output  = state->print_output;
  cvodes_params = (struct cvodes_params_struct *)state->cvodes_params;
  if (jacobian_choice == 0) {
    ode23tb_params = (struct ode23tb_params_struct *)state->ode23tb_params;
  }
  cvodes_params->ns = ns;
  lmm = cvodes_params->linear_multistep_method;
  iter = cvodes_params->iterative_method;
  cvode_mxsteps = cvodes_params->mxsteps;
  /*
    The cvodes routines will need some space for formation of the jacobian,
    drfc[number_molecules)  jacobian formation workspace.

    doubles:
    drfc    (number_molecules)
    dfdy_a  (nnz <= nunique_molecules*nunique_molecules) 
    dfdy_at (nnz <= nunique_molecules*nunique_molecules) 
    miter_m (nnzm <= nunique_molecules*nunique_molecules) 
    prec_l (nnzl <= nunique_molecules*(nunique_molecules+1)/2)
    prec_u (nnzu <= nunique_molecules*(nunique_molecules+1)/2)
    recip_diag_u (nunique_molecules)
    prec_row (nunique_molecules)
    frow (nunique_molecules)
    srow (nunique_molecules)
    pbar (number_reactions)
    ys0v (nunqiue_molecules * number_reactions)

    ints 
    dfdy_ia(nunqiue_molecules+1)
    dfdy_ja(nnz <= nunique_molecules*nunique_molecules) 
    dfdy_iat(nunqiue_molecules+1)
    dfdy_jat(nnz <= nunique_molecules*nunique_molecules) 
    miter_im (nunique_moleucles+1)
    miter_jm (nnzm <= nunique_molecules*nunique_molecules) 

    prec_il (nunique_molecules+1)
    prec_jl (nnzl <= nunique_molecules*(nunique_molecules+1)/2)
    prec_iu (nunique_molecules+1)
    prec_ju (nnzu <= nunique_molecules*(nunique_molecules+1)/2)

    lindex (unique_molecules)
    uindex (unique_molecules)
    column_mask (unique_molecules)
    sindex (unique_molecules)
    plist  (number_reactions)
    
    Actually though much less space might be needed there 
    might be a formula to compute the needed size of dfdy based on the reaction
    matrix, call it nnz, then the miter_m and miter_jm need at most 
    nnz + nunique_molecules space, and prec_l and prec_jl, prec_u, and prec_ju
    need at most nnz + nunqiue_molecules * (2*fill_level + 1).
    For now space is cheap allocate all the way out.
    Also we will want one scratch vector of ny doubles and 
    2 scratch vectors of ny_ints for forming the factorizations.
    Hmm we also need space for the drfc vector of length number_molecules.
  */
  boltzmann_size_jacobian(state);
  nnz  = cvodes_params->nnz;
  nnzm = cvodes_params->nnzm;
  nnzl = cvodes_params->nnzl;
  nnzu = cvodes_params->nnzu;
  num_doubles = (nnz + nnz + nnzm + nnzl + nnzu + (5*ny) + (3*ns) + (ns*ny));
  num_ints    = (num_doubles + 5*ny + 5 + ns);
  num_doubles = num_doubles + (num_ints + (num_ints & 1))/2;
  drfc_len    = state->number_molecules * 2;
  ask_for = (num_doubles + drfc_len) << 3;
  drfc    = (double *)calloc(one_l,ask_for);
  if (drfc == NULL ) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"boltzmann_cvodes: Error unable to allocate %ld bytes needed for  jacobian workspace\n",ask_for);
      fflush(lfp);
    }
  } else {
    dfdy_a       = &drfc[drfc_len];
    dfdy_at      = &dfdy_a[nnz];
    miter_m  	 = &dfdy_at[nnz];
    prec_l   	 = &miter_m[nnzm];
    prec_u   	 = &prec_l[nnzl];
    prec_row 	 = &prec_u[nnzu];
    recip_diag_u = &prec_row[ny];
    frow         = &recip_diag_u[ny];
    srow         = &frow[ny];
    f0           = &srow[ny];
    p            = &f0[ny];
    rp           = &p[ns];
    pbar         = &rp[ns];
    ys0v         = &pbar[ns];
    dfdy_ia  	 = (int*)&ys0v[ny*ns];
    dfdy_ja  	 = &dfdy_ia[ny+1];
    dfdy_iat     = &dfdy_ja[nnz];
    dfdy_jat     = &dfdy_iat[ny+1];
    miter_im 	 = &dfdy_jat[nnz];
    miter_jm 	 = &miter_im[ny+1];
    prec_il  	 = &miter_jm[nnzm];
    prec_jl  	 = &prec_il[ny+1];
    prec_iu  	 = &prec_jl[nnzl];
    prec_ju      = &prec_iu[ny+1];
    lindex       = &prec_ju[nnzu];
    uindex       = &lindex[ny];
    column_mask  = &uindex[ny];
    sindex       = &column_mask[ny];
    plist        = &sindex[ny];
    cvodes_params->drfc         = drfc;
    cvodes_params->dfdy_a  	= dfdy_a;
    cvodes_params->dfdy_at  	= dfdy_at;
    cvodes_params->miter_m 	= miter_m;
    cvodes_params->prec_l  	= prec_l;
    cvodes_params->prec_u  	= prec_u;
    cvodes_params->prec_row     = prec_row;
    cvodes_params->recip_diag_u = recip_diag_u;
    cvodes_params->frow         = frow;
    cvodes_params->srow         = srow;
    cvodes_params->p            = p;
    cvodes_params->rp           = rp;
    cvodes_params->pbar         = pbar;
    cvodes_params->ys0v         = ys0v;
    cvodes_params->dfdy_ia      = dfdy_ia;
    cvodes_params->dfdy_ja      = dfdy_ja;
    cvodes_params->dfdy_iat     = dfdy_iat;
    cvodes_params->dfdy_jat     = dfdy_jat;
    cvodes_params->miter_im     = miter_im;
    cvodes_params->miter_jm     = miter_jm;
    cvodes_params->prec_il      = prec_il;
    cvodes_params->prec_jl      = prec_jl;
    cvodes_params->prec_iu      = prec_iu;
    cvodes_params->prec_ju      = prec_ju;
    cvodes_params->lindex       = lindex;
    cvodes_params->uindex       = uindex;
    cvodes_params->column_mask  = column_mask;
    cvodes_params->sindex       = sindex;
    cvodes_params->plist        = plist;
  }
  if (success) {
    if (jacobian_choice == 0) {
      num_doubles = ny*ny + 5*ny;
      ask_for = num_doubles << 3;
      dfdy = (double*)calloc(one_l,ask_for);
      if (dfdy != NULL) {
	fac  = (double*)&dfdy[ny*ny];
	thresh = (double*)&fac[ny];
	fdel   = (double*)&thresh[ny];
	fdiff  = (double*)&fdel[ny];
	dfdy_tmp = (double*)&fdiff[ny];
	ode23tb_params->dfdy = dfdy;
	ode23tb_params->fac  = fac;
	ode23tb_params->thresh = thresh;
	ode23tb_params->fdel   = fdel;
	ode23tb_params->fdel   = fdel;
	ode23tb_params->num_jac_first_time = 1;
	ode23tb_params->nf     = 0;
	njthreshold = state->nj_thresh;
	vec_set_constant(ny,thresh,njthreshold);
      } else {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_cvodes: Error unable to allocate %ld bytes needed for numeric approximation to jacobian workspace\n",ask_for);
	  fflush(lfp);
	}
      }
    }
  }
  if (success) {
    if (print_output == 0) {
      ode_rxn_view_freq = 0;
      state->ode_rxn_view_freq = 0;
    }

    /* 
      If tracking net likelihoods and fluxes, set up to
      print the first iteration, 
    */
    ode_rxn_view_step = one_l;
    if (ode_rxn_view_freq > 0) {
      t0 = 0.0;
      boltzmann_monitor_ode(state,t0,concs);
    }
    cvode_mem = CVodeCreate(lmm,iter);
    if (cvode_mem == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_cvodes: Error cvode_mem was null, lmm = %d, iter = %d\n",lmm,iter);
	fflush(lfp);
      }
    }
  }
  /*
    Not quite right here, we need to make concs an "N_vector" type.
    y0 Need to look at sectin 7.1-7.4 for details. Also need to test for
    success from this call, I suppose y0 != NULL
  */
  if (success) {
    y0 = N_VMake_Serial(ny,concs);
    if (y0 == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_cvodes: Error creating nvector y0 of length %d from concs\n",ny);
	fflush(lfp);
      }
    }
  }
  /*
    Initialize the cvode solver specifying the rhs function, t0, and y0.
    This must set the internal problems size, ny from the nvector y0.
  */
  if (success) {
    cvodes_params->y0 = y0;
    t0 = 0.0;
    cvodes_params->cvode_mem = cvode_mem;
    flag = CVodeInit(cvode_mem,boltzmann_cvodes_rhs,t0,y0);
    if (flag < 0) {
      success = 0;
      boltzmann_print_cvodeinit_errors(flag,cvode_mem,state);
    }
  } /* end if (success) */
  /*
    Initialize parameters for the cvodes method from 
    the cvodes_params structure.
  */
  if (success) {
    success = boltzmann_cvodes_init(cvode_mem,state,concs);
  }
  if (success) {
    num_steps = cvodes_params->num_cvode_steps;
    tfinal = state->ode_t_final;
    if (num_steps > 0) {
      itask = CV_NORMAL;
      t0      = 0.0;
      delt    = tfinal/((double)num_steps);
      not_done = 1;
      for (i=0;((i<num_steps) && not_done);i++) {
	tout = delt * (i+1);
	flag = CVode(cvode_mem, tout, y0, &tret, itask);

	if (flag < 0) {
	  success = 0;
	  not_done = 0;
	  boltzmann_print_cvode_error(flag,cvode_mem,state);
	}
	if (success) {
	  if (flag == CV_TSTOP_RETURN) {	
	    not_done = 0;
	    if (lfp) {
	      fprintf(lfp,"Boltzmann_cvodes: CVode returned CV_TSTOP_RETURN\n");
	      fflush(lfp);
	    }
	  } else if (flag == CV_ROOT_RETURN) {
	    not_done = 0;
	    if (lfp) {
	      fprintf(lfp,"Boltzmann_cvodes: CVode returned CV_ROOT_RETURN\n");
	      fflush(lfp);
	    }
	  }
	  /*
	    Do monitoring stuff.
	    See lines 968ff in ode23tb. maybe boltzmann_ode_monitor
	    routine, first need to recover concs from y0.
	    No actually concs should now be up to date.
	  */
	  if (ode_rxn_view_freq > 0) {
	    ode_rxn_view_step -= one_l;
	    if (ode_rxn_view_step == zero_l) {
	      boltzmann_monitor_ode(state,tret,concs);
	      ode_rxn_view_step = ode_rxn_view_freq;
	    }
	  }
	}
	done = ode_test_steady_state(state,ny,concs,f0);
	if (done) {
	  not_done = 0;
	}
      } /* end for (i...) */
      /* end if num_steps > 0) */
    } else {
      /* num_steps == 0) */
      not_done = 1;
      itask = CV_ONE_STEP;
      tout = tfinal;
      while (not_done) {
	flag = CVode(cvode_mem, tout, y0, &tret, itask);
	if (flag < 0) {
	  success = 0;
	  not_done = 0;
	  boltzmann_print_cvode_error(flag,cvode_mem,state);
	}
	if (success) {
	  if (flag == CV_TSTOP_RETURN) {	
	    not_done = 0;
	    if (lfp) {
	      fprintf(lfp,"Boltzmann_cvodes: CVode returned CV_TSTOP_RETURN\n");
	      fflush(lfp);
	    }
	  } else if (flag == CV_ROOT_RETURN) {
	    not_done = 0;
	    if (lfp) {
	      fprintf(lfp,"Boltzmann_cvodes: CVode returned CV_ROOT_RETURN\n");
	      fflush(lfp);
	    }
	  }
	}
	if (success) {
	  /*
	    Do monitoring stuff.
	    See lines 968ff in ode23tb. maybe boltzmann_ode_monitor
	    routine, first need to recover concs from y0.
	    No actually concs should now be up to date.
	  */
	  if (ode_rxn_view_freq > 0) {
	    ode_rxn_view_step -= one_l;
	    if (ode_rxn_view_step == zero_l) {
	      boltzmann_monitor_ode(state,tret,concs);
	      ode_rxn_view_step = ode_rxn_view_freq;
	    }
	  }
	}
	done = ode_test_steady_state(state,ny,concs,f0);
	if (done) {
	  not_done = 0;
	}
      } /* end while not_done */
      
    } /* end else num_cvodes_steps == 0 */
    
  } /* end if (success) */
  if (success) {
    if (state->compute_sensitivities) {
      boltzmann_print_sensitivities(state);
    }
  }
  /*
    Unset the compute_sensitivities field so as to use the original ke vector
    in subsequent boltzmann record-step iterations.
  */
  state->compute_sensitivities = 0;
  /*
    Free cvodes memory.
  */
  return(success);
}
