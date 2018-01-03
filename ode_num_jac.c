#include "boltzmann_structs.h"

#include "num_jac_col.h"
#include "blas.h"

#include "ode_num_jac.h"
int ode_num_jac(struct state_struct *state,
		int    first_time, /* 1 for first call, 0 otherwise. */
		double *dfdy, /* ny x ny */
		double t,  /* scalar */
		double *y, /* ny x 1 */
		double *f, /* ny x 1 */
		double *fac, /* ny x 1 */
		double *thresh, /* ny x 1 */
		double *ode_num_jac_scratch, /* 4*ny */
		int64_t *nfcalls /* scalar returned */) {
  /*
    Called by: ode23tb
    Calls:     num_jac_col, sqrt, pow, abs, max, min, dcopy
    thresh[i] = 1.e-6;
    Looks like we need 10 temporary vectors.

    Variable                    TMF
    state                       G*I Uses nunique_molecles, 
                                    base_reaction, conc_to_count fields
				    kf_base_reaction,
				    base_reactants (in num_jac_col)
                                    
    dfdy                        D*O ny by ny jacobian of f w.r.t. y 
    t                           DSI time.
    y                           D*I concentrarions variable, length ny
    f                           D*I vector of fluxes,  length ny
    fac                         D*B history vector of length ny
    thresh                      D*I threshold vector, length ny
    ode_num_jac_scratch         D*W scratch space length 4*ny + 2*nrxns
    *nfcalls                    PSO number of funtion calls to approximate 
                                    fluxes 

  */
  double *y_counts;   /* ny x 1 */
  double *fdiff;        /* ny x 1 */
  double *fdel;   /* ny x 1 */
  double *dfdy_tmp;    /* ny x 1 */
  double *forward_rxn_likelihoods;  /* nrxns x 1 */
  double *reverse_rxn_likelihoods; /* nrxnx x 1 */

  double *conc_to_count; /* from state ny x 1 */
  double *dfdy_colj;    /* pointer into dfdy matrix, not allocated ny x 1 */

  double br;
  double bl;
  double bu;
  double facmin;
  double facmax;
  double eps;
  double sqrt_eps;
  double fourthrt_eps;
  double eighthrt_eps;
  double delj;
  double facj;
  double yj;
  double yscalej;
  double threshj;
  double fj;
  double fscalej; 
  double fscaletmp; 
  double absfdiffmax;
  double absfdelrm;
  double absfvaluerm;
  double infnormdfdy_colj;
  double infnormdfdy_tmp;
  double absfdiffmax_tmp;
  double tmpfac;
  double *dblptr;
  int64_t eps_hex;
  int64_t sqrt_eps_hex;
  int64_t fourthrt_eps_hex;
  int ny;
  int success;

  int j;
  int index1;
  
  int nfc;
  int nrxns;
  
  success = 1;
  index1  = 1;
  ny               = state->nunique_molecules;
  nrxns            = state->number_reactions;
  conc_to_count    = state->conc_to_count;
  eps_hex          = 0x3CB0000000000000L;
  dblptr           = (double*)&eps_hex;
  eps              = *dblptr;
  sqrt_eps_hex     = 0x3E50000000000000L;
  dblptr           = (double*)&sqrt_eps_hex;
  sqrt_eps         = *dblptr;
  fourthrt_eps_hex = 0x3F20000000000000L;
  dblptr           = (double *)&fourthrt_eps_hex;
  fourthrt_eps     = *dblptr;
  eighthrt_eps = sqrt(fourthrt_eps);
  bl           = sqrt_eps * fourthrt_eps;
  br           = eighthrt_eps * bl;
  bu           = fourthrt_eps;
  facmin       = pow(eps,.78);
  facmax       = .1;
  nfc          = 0;

  y_counts      = (double*)ode_num_jac_scratch;
  fdel          = (double*)&y_counts[ny];
  fdiff         = (double*)&fdel[ny];
  dfdy_tmp      = (double*)&fdiff[ny];
  forward_rxn_likelihoods = (double*)&dfdy_tmp[ny];
  reverse_rxn_likelihoods = (double*)&forward_rxn_likelihoods[nrxns];

  
  if (first_time == 1) {
    for (j=0;j<ny;j++) {
      fac[j] = sqrt_eps;
    }
  }
  for (j=0;j<ny;j++) {
    y_counts[j] = y[j] * conc_to_count[j];
  }
  dfdy_colj = dfdy;
  for (j=0;j<ny;j++) {
    yj   = y[j];
    fj   = f[j];
    facj = fac[j];
    threshj = thresh[j];
    /*
    yscalej   = max (abs(yj),threshj);
    */
    yscalej   = abs(yj);
    if (yscalej < threshj) {
      yscalej = threshj;
    }
    delj      = (yj + (facj*yscalej)) - yj;

    while (delj == 0.0) {
      if (facj < facmax) {
	/*
	facj = min(100*facj,facmax);
	*/
	facj = 100*facj;
	if (facj > facmax) {
	  facj = facmax;
	}
	delj = (yj + (facj*yscalej)) - yj;
      } else {
	delj = threshj;
      }
    } /* end while (delj == 0.0 */
    if (fj >= 0.0) {
      delj = abs(delj);
    } else {
      delj = - abs(delj);
    }
    num_jac_col(state,ny,j,y,f,&delj,threshj,
		y_counts,fdel,fdiff,
		forward_rxn_likelihoods,
		reverse_rxn_likelihoods,
		dfdy_colj,
		&absfdiffmax,
		&absfdelrm,
		&absfvaluerm,
		&infnormdfdy_colj);
    nfc += 1;
    /*
    fscalej = max(absfvaluerm,absfdelrm);
    */
    fscalej = absfvaluerm;
    if (absfdelrm > fscalej) {
      fscalej = absfdelrm;
    }
    if ((absfdiffmax == 0.0) || ((absfdelrm > 0.0) && (absfvaluerm > 0.0))) {
      if (absfdiffmax < (br * fscalej)) {
	/*
	tmpfac = min(sqrt(facj),facmax);
	*/
	tmpfac = sqrt(facj);
	if (facmax < tmpfac) {
	  tmpfac = facmax;
	}
	delj   = (yj + (tmpfac * yscalej)) - yj;
	if ((tmpfac != facj) && (delj != 0.0)) {
	  num_jac_col(state,ny,j,y,f,&delj,threshj,
		      y_counts,fdel,fdiff,
		      forward_rxn_likelihoods,
		      reverse_rxn_likelihoods,
		      dfdy_tmp,
		      &absfdiffmax_tmp,
		      &absfdelrm,
		      &absfvaluerm,
		      &infnormdfdy_tmp);
	  nfc += 1;
	  if ((tmpfac * infnormdfdy_tmp) > infnormdfdy_colj) {
	    dcopy(&ny,dfdy_tmp,&index1,dfdy_colj,&index1);
	    /*
	    fscaletmp = max(absfvaluerm,absfdelrm);
	    */
	    fscaletmp = absfvaluerm;
	    if (absfdelrm > fscaletmp) {
	      fscaletmp = absfdelrm;
	    }
	    if (absfdiffmax_tmp <= (bl * fscaletmp)) {
	      /*
	      facj  = min(10*tmpfac,facmax);
	      */
	      facj = 10 * tmpfac;
	      if (facj > facmax) {
		facj = facmax;
	      }
	    } else {
	      if (absfdiffmax_tmp > (bu * fscaletmp)) {
		/*
		facj = max(0.1*tmpfac, facmin);
		*/
		facj = 0.1 * tmpfac;
		if (facmin > facj) {
		  facj = facmin;
		}
	      } else {
		facj = tmpfac;
	      }
	    }
	    fac[j] = facj;
	  } /*end if  ((tmpfac * infnormdfdy_tmp) > infnormdfdy_colj) */
	} /* end if ((tmpfac != facj) && (delj > 0.0)) */
      } /* end if (absfdiffmax < (br * fscalej)) */
      if (absfdiffmax <= bl * fscalej) {
	/*
	fac[j] = min((10*facj),facmax);
	*/
	facj = 10*facj;
	if (facj > facmax) {
	  facj = facmax;
	}
	fac[j] = facj;
      }
      if (absfdiffmax > bu *fscalej) {
	/*
	fac[j] = max(0.1*facj,facmin);
	*/
	facj = 0.1 * facj;
	if (facmin > facj) {
	  facj = facmin;
	}
	fac[j] = facj;
      }
    } /* end if ((absfdiffmax == 0.0) || ((absfdelrm > 0.0) && (absfvaluerm > 0.0))) */
  } /* end for j */
  *nfcalls = nfc;
  return (success);
} 

