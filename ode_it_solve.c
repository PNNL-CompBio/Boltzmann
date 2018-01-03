#include "boltzmann_structs.h"

#include "compute_flux_scaling.h"
#include "approximate_fluxes.h"

#include "ode_it_solve.h"
int ode_it_solve(struct state_struct *state,
		 double t,
		 double *y,
		 double *z,
		 double *del,
		 double *znew,
		 double *y_count,
		 double *forward_rxn_likelihoods,
		 double *reverse_rxn_likelihoods,
		 double d,
		 double h,
		 double rtol,
		 double wt,
		 double *rate_p,
		 int *iter_count_p) {
  /*
    Iteratively solve nonlinear  z = h*f(t,v+d*z) where y = v+d*z.
    del, znew, y_count, forward_rxn_likelihoods, 
    and reverse_reaction_likeklihoods are scratch vectors of length ny.    
    computes a result in z, and also modifies y.
    Returns 0 on successful iteration, 1 on a fail.
    Called by: ode23tb
    Calls:     approximate_fluxes, fabs;
  */
  double kappa;
  double errit;
  double eps;
  double flux_scaling;
  double minnrm;
  double cand_nrm;
  double *conc_to_count;
  double recip_wt;
  double newnrm;
  double oldnrm;
  double kappa_rtol;
  double cap;
  double del_scale;
  double cand_rate;
  double p9_old_rate;
  double rate;
  double *dbl_ptr;
  double yi;

  int64_t eps_hex;
  int itfail;
  int max_iter;

  int iter;
  int base_rxn;

  int i;
  int ny;

  int iter_count;
  int padi;

  int nsteps;
  int origin;

  max_iter     = 5;
  kappa        = 0.5;
  itfail       = 0;
  recip_wt     = 1.0/wt;
  base_rxn     = state->base_reaction;
  conc_to_count = state->conc_to_count;
  ny            = state->nunique_molecules;
  kappa_rtol    = kappa * rtol;
  eps_hex       = 0x3CB0000000000000L;
  dbl_ptr       = (double *)&eps_hex;
  eps           = *dbl_ptr;
  minnrm = 0.0;
  rate   = *rate_p;
  iter_count = *iter_count_p;
  for (i=0;i<ny;i++) {
    cand_nrm = 100 * eps * fabs(y[i]) * recip_wt;
    if (cand_nrm > minnrm) {
      minnrm = cand_nrm;
    }
  }
  for (iter = 0; iter<max_iter;iter++) {
    for (i=0;i<ny;i++) {
      y_count[i] = y[i]*conc_to_count[i];
    }
    flux_scaling = compute_flux_scaling(state,y);
    approximate_fluxes(state,y_count,forward_rxn_likelihoods, 
		       reverse_rxn_likelihoods,
		       znew,flux_scaling,base_rxn);
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"ode_it_solve: after call to approximate_fluxes:\n");
      nsteps = 0;
      origin = -1;
      print_concs_fluxes(state,ny,znew,y,y_count,
			 forward_rxn_likelihoods,
			 reverse_rxn_likelihoods,t,h,nsteps,origin);
#endif    

    newnrm = 0.0;
    for (i=0;i<ny;i++) {
      znew[i] = h * znew[i];
      del[i]  = znew[i] - z[i];
      yi      = y[i];
      y[i]    = yi + d * del[i];
      /*
	Here again we need to not let y be negative.
	We do assume that y[i] is nonnegative on input,
	maybe should test that hypothesis.
      */
      if (y[i] < 0.0) {
	y[i] = .5*yi;
	del[i] = -(.5 * yi) / d;
	znew[i] = z[i] + del[i];
      }
      z[i]    = znew[i];
      del_scale = fabs(y[i]);
      if (wt > del_scale) {
	del_scale = wt;
      }
      cand_nrm = fabs(del[i])/del_scale;
      if (cand_nrm > newnrm) {
	newnrm = cand_nrm;
      }
    }
    if (newnrm < minnrm) {
      iter_count = iter + 1;
      break;
    } else {
      if (iter == 0) {
	if (rate > 0) {
	  errit = newnrm * rate/(1.0 - rate);
	  if (errit < .1 * kappa * rtol) {
	    iter_count = iter+1;
	    break;
	  } 
	} else {
	  rate = 0;
	}
      } else {
	if (newnrm > (0.9 * oldnrm)) {
	  itfail = 1;
	  iter_count = iter+1;
	  break;
	} else {
	  p9_old_rate = 0.9 * rate;
	  cand_rate   = newnrm/oldnrm;
	  if (p9_old_rate > cand_rate) {
	    cand_rate = p9_old_rate;
	  }
	  rate = cand_rate;
	  errit = newnrm * rate/(1 - rate);
	  if (errit < kappa * rtol) {
	    iter_count = iter+1;
	    break;
	  } else {
	    if (iter == max_iter) {
	      itfail = 1;
	      iter_count = iter+1;
	      break;
	    } else {
	      cap = 1.0;
	      for (i=iter;i<max_iter-1;i++) {
		cap = cap * rate;
	      }
	      if (kappa_rtol < errit*cap) {
		itfail = 1;
		iter_count = iter+1;
		break;
	      }
	    }
	  }
	}
      }
    }
    oldnrm = newnrm;
  }
  *rate_p = rate;
  *iter_count_p = iter_count;
  return(itfail);
}
	    

