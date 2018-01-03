#include "boltzmann_structs.h"

#include "compute_flux_scaling.h"
#include "approximate_fluxes.h"
#include "ce_approximate_fluxes.h"
#include "blas.h"
#include "lapack.h"
#include "vec_abs.h"
#include "vec_div.h"
#include "vec_mul.h"
#include "vec_max.h"

#include "ode_it_solve.h"
int ode_it_solve(struct state_struct *state,
		 double *miter,
		 int    *ipivot,
		 double t,
		 double *y,
		 double *z,
		 double *del,
		 double *rhs,       
		 double *scratch,   /* used to compute candnrm */
		 double *y_count,
		 double *forward_rxn_likelihoods,
		 double *reverse_rxn_likelihoods,
		 double d,
		 double h,
		 double rtol,
		 double *wt,
		 double *rate_p,
		 int    *iter_count_p) {
  /*
    Iteratively solve nonlinear  Mz = h*f(t,v+d*z) where y = v+d*z.
    del, rhs, scratch,  y_count, forward_rxn_likelihoods, 
    and reverse_reaction_likeklihoods are scratch vectors of length ny.    
    computes a result in z, and also modifies y.
    miter and ipivot contain the LU factorization of I-d*h*dfdy,
    Returns 0 on successful iteration, 1 on a fail.
    Called by: ode23tb
    Calls:     approximate_fluxes, fabs, dgetrs,
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
  double tenth_kappa_rtol;
  double cap;
  double del_scale;
  double cand_rate;
  double p9_old_rate;
  double rate;
  double *dbl_ptr;
  double yi;
  double zi;
  double hundred_eps;
  double local_h;
  double local_d;
  double minus_one;
  double one;
  double rate_ratio;

  int64_t eps_hex;
  int itfail;
  int max_iter;

  int iter;
  int base_rxn;

  int i;
  int ny;

  int nrhs;
  int info;
  
  int nsteps;
  int origin;

  int iter_count;
  int inc1;

  int exit_not;
  int padi;

  char  trans_chars[8];
  char  *trans;
  FILE *lfp;
  FILE *efp;
/*
*/
#define DBG 1
  max_iter     	   = 5;
  kappa        	   = 0.5;
  itfail       	   = 0;
  base_rxn     	   = state->base_reaction;
  conc_to_count    = state->conc_to_count;
  ny               = state->nunique_molecules;
  lfp              = state->lfp;
  kappa_rtol       = kappa * rtol;
  tenth_kappa_rtol = 0.05 * rtol;
  eps_hex          = 0x3CB0000000000000L;
  dbl_ptr          = (double *)&eps_hex;
  eps              = *dbl_ptr;
  inc1             = 1;
  rate             = *rate_p;
  hundred_eps      = 100.0 * eps;
  local_h          = h;
  local_d          = d;
  minus_one        = -1.0;
  one              = 1.0;
  /*
    scratch <- |y|
  */
  vec_abs(&ny,scratch,y);
  /*
    scratch <- 100 * eps * |y|
  */
  dscal_(&ny,&hundred_eps,scratch,&inc1);
  /*
    scratch <- scratch ./ wt
  */
  vec_div(&ny,scratch,scratch,wt);
  /*
    minnrm = 100 * eps * norm(y ./ wt,inf);
  */
  /*
  for (i=0;i<ny;i++) {
    cand_nrm = 100 * eps * fabs(y[i])/wt[i];
    if (cand_nrm > minnrm) {
      minnrm = cand_nrm;
    }
  }
  */
  minnrm = scratch[idamax_(&ny,scratch,&inc1)-1];
  iter_count = 0;
  trans_chars[0] = 'N';
  trans = &trans_chars[0];
  /*
  oldnrm = minnrm;
  */
  exit_not = 1;
  for (iter = 0; ((iter<max_iter) && exit_not);iter++) {
    iter_count = iter + 1;
    /*
      Compute counts to be used in computing concentration derivative.
    for (i=0;i<ny;i++) {
      y_count[i] = y[i]*conc_to_count[i];
    }
      y_count <- y .* conc_to_count
    */
    vec_mul(&ny,y_count,y,conc_to_count);
    /*
      rhs <- y', approximate fluxes.
    */
    
    flux_scaling = compute_flux_scaling(state,y);
    ce_approximate_fluxes(state,y_count,forward_rxn_likelihoods, 
		       reverse_rxn_likelihoods,
		       rhs,flux_scaling,base_rxn);
    /*
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"ode_it_solve: after call to ce_approximate_fluxes:\n");
      nsteps = 0;
      origin = 7;
      print_concs_fluxes(state,ny,rhs,y,y_count,
			 forward_rxn_likelihoods,
			 reverse_rxn_likelihoods,t,h,nsteps,origin);
    }
#endif    
    */
    /*
    approximate_fluxes(state,y_count,forward_rxn_likelihoods, 
		       reverse_rxn_likelihoods,
		       rhs,flux_scaling,base_rxn);
    */
    /*
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"ode_it_solve: after call to approximate_fluxes:, flux_scaling = %le\n",flux_scaling);
      nsteps = 0;
      origin = 7;
      print_concs_fluxes(state,ny,rhs,y,y_count,
			 forward_rxn_likelihoods,
			 reverse_rxn_likelihoods,t,h,nsteps,origin);
    }
#endif
    */
    /*
      for (i=0;i<ny;i++) {
        rhs[i] = h * rhs[i];
        del[i]  = rhs[i] - z[i];
      }
      rhs <- y' * h
    */
    dscal_(&ny,&local_h,rhs,&inc1);
    /*
      rhs <- y' * h - z
      could use a daymx (y<- a*y-x) routine here with y = rhs, a = h, and x = z)
    */
    daxpy_(&ny,&minus_one,z,&inc1,rhs,&inc1);
    /*
      del <- rhs;
    */
    dcopy_(&ny,rhs,&inc1,del,&inc1);
    /*
      Solve Miter * del = rhs for del
    */
    nrhs = 1;
    dgetrs_(trans,&ny,&nrhs,miter,&ny,ipivot,del,&ny,&info);
    if (info != 0) {
      if (lfp) {
	fprintf (lfp,"ode_it_solve: Error nonzero return code from dgetrs was %s\n",info);
	fflush(lfp);
      }
      itfail = 1;
      exit_not = 0;
    } else {
      /*
	dgetrs solve succeeded.
	vec_add could be used here.
	z <- z + del;
      */
      daxpy_(&ny,&one,del,&inc1,z,&inc1);
      /*
	y <- y + d * del
      */
      daxpy_(&ny,&d,del,&inc1,y,&inc1);
      /*
      newnrm = norm(del ./ max(wt,abs(y)),inf);
      scratch <- |y|
      */
      vec_abs(&ny,scratch,y);
      /*
	scratch <- max (|y|,wt))
      */
      vec_max(&ny,scratch,scratch,wt);
      /*
	scratch <- del ./  max(|y|,st)
      */
      vec_div(&ny,scratch,del,scratch);
      newnrm = fabs(scratch[idamax_(&ny,scratch,&inc1)-1]);
      /*
      newnrm = 0.0;
      for (i=0;i<ny;i++) {
        zi      = z[i];
        z[i]    = zi + del[i];
        yi      = y[i];
        y[i]    = yi + d * del[i];

        del_scale = fabs(y[i]);
        if (wt[i] > del_scale) {
      	del_scale = wt[i];
        }
        cand_nrm = fabs(del[i])/del_scale;
        if (cand_nrm > newnrm) {
      	newnrm = cand_nrm;
        }
      } // end for (i...) 
      */
      if (newnrm <= minnrm) {
	exit_not = 0;
      } else {
	if (iter == 0) {
	  /*
	    First iteration.
	  */
	  if (rate > 0.0) {
	    errit = newnrm * rate / ( 1.0 - rate );
	    if (errit <= tenth_kappa_rtol) {
	      exit_not = 0;
	    }
	  } else {
	    rate = 0.0;
	  }
	} else {
	  /*
	    Not first iteration.
	    if newnrm has grown by more than 10 % fail.
	  */
	  if (newnrm > 0.9*oldnrm) {
	    itfail = 1;
	    exit_not = 0;
	  } else {
	    /*
	      newnrm not too large.
	    */
	    rate = 0.9 * rate;
	    errit = newnrm/oldnrm;
	    rate = (rate > errit) ? rate : errit;
	    errit = newnrm * (rate/ (1.0 - rate));
	    if (errit <= kappa_rtol) {
	      /*
		Converged.
	      */
	      exit_not = 0;
	    } else {
	      if (iter == (max_iter - 1)) {
		/* 
		   last iteration and not converged. 
		*/
		itfail   = 1;
		exit_not = 0;
	      } else {
		cap = 1.0;
		for (i=iter;i<max_iter-1;i++) {
		  cap = cap * rate;
		}
		if (kappa_rtol < (errit * cap)) {
		  itfail   = 1;
		  exit_not = 0;
		}
	      } /* end else not last iter */
	    } /* end else not converged */
	  } /* end else newnrm not more than 10% above oldnrm */
	} /* end else not first iteration */
      } /* end else newnrm > minnrm */
    } /* end else dgetrs succeeded */
    oldnrm = newnrm;
  } /* end for (iter ... )*/
  *rate_p = rate;
  *iter_count_p = iter_count;
  return(itfail);
}
	    

