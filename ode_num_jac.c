#include "boltzmann_structs.h"

#include "num_jac_col.h"
#include "blas.h"

#include "ode_num_jac.h"
/*
int ode_num_jac(struct state_struct *state,
		int    first_time, // 1 for first call, 0 otherwise. 
		double *dfdy, // ny x ny 
		double t,  // scalar 
		double *y, // ny x 1 
		double *f, // ny x 1 
		double *fac, // ny x 1 
		double *thresh, // ny x 1 
		double *ode_num_jac_scratch, // 4*ny 
		int64_t *nfcalls // scalar returned ) 
*/
int ode_num_jac(struct state_struct *state,
		int    first_time, /* 1 for first call, 0 otherwise. */
		double *dfdy, /* ny x ny */
		double t,  /* scalar */
		double *y, /* ny x 1 */
		double *f, /* ny x 1 */
		double *fac, /* ny x 1 */
		double *thresh, /* ny x 1 */
		double *fdel, /* ny x 1 */
		double *fdiff, /* ny x 1 */
		double *dfdy_tmp, /* ny x 1 */
		int64_t *nfcalls /* scalar returned */) {
  /*
    Called by: ode23tb
    Calls:     num_jac_col, sqrt, pow, fabs, dcopy
    thresh[i] = 1.e-6;
    Looks like we need 10 temporary vectors.

    Variable                    TMF
    state                       G*I Uses nunique_molecles, 
                                    number_reactions
                                    base_reaction, conc_to_count fields
				    kf_base_reaction,
				    base_reactants (in num_jac_col)
                                    
    dfdy                        D*O ny by ny jacobian of f w.r.t. y 
    t                           DSI time.
    y                           D*I concentrarions variable, length ny
    f                           D*I vector of fluxes,  length ny
    fac                         D*B history vector of length ny
    thresh                      D*I threshold vector, length ny
    fdel                        D*W scratch vector length ny
    fdiff                       D*W scratch vector length ny
    dfdy_tmp                    D*W scratch vector length ny
    *nfcalls                    PSO number of funtion calls to approximate 
                                    fluxes 

  */
  /*
    These vectors are now passed in as arcuments.
  double *fdiff;      // ny x 1 
  double *fdel;       // ny x 1 
  double *dfdy_tmp;   // ny x 1 
  */
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
  double absfdelrm_tmp;
  double tmpfac;
  double tmpval;
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

  int rowmax;
  int i;

  FILE *lfp;
  FILE *efp;
  
  success = 1;
  index1  = 1;
  lfp              = state->lfp;
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

  /*
  fdel          = (double*)&y_counts[ny];
  fdiff         = (double*)&fdel[ny];
  dfdy_tmp      = (double*)&fdiff[ny];
  */
  /*
#define DBG 1
  */
#ifdef DBG
  if (lfp) {
    fprintf(lfp,"ode_num_jac_top ------------------------------------\n");
    fflush(lfp);
  }
#endif  
  if (first_time == 1) {
    for (j=0;j<ny;j++) {
      fac[j] = sqrt_eps;
    }
  }
  dfdy_colj = dfdy;
  for (j=0;((j<ny) && success);j++) {
    yj   = y[j];
    fj   = f[j];
    facj = fac[j];
    threshj = thresh[j];
    /*
    yscalej   = max (fabs(yj),threshj);
    */
    yscalej   = fabs(yj);
    if (yscalej < threshj) {
      yscalej = threshj;
    }
    delj      = (yj + (facj*yscalej)) - yj;
#ifdef DBG
    fprintf(lfp,"y[%d]      = %le\n",j,yj);
    fprintf(lfp,"f[%d]      = %le\n",j,fj);
    fprintf(lfp,"thresh[%d] = %le\n",j,threshj);
    fprintf(lfp,"yscale[%d] = %le\n",j,yscalej);
    fprintf(lfp,"fac[%d]    = %le\n",j,facj);
    fprintf(lfp,"del[%d]    = %le\n",j,delj);
    fflush(lfp);
#endif
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
      delj = fabs(delj);
    } else {
      delj = - fabs(delj);
    }
#ifdef DBG
    fprintf(lfp,"fac[%d]    = %le\n",j,facj);
    fprintf(lfp,"del[%d]    = %le\n",j,delj);
    fflush(lfp);
#endif
    success = num_jac_col(state,ny,j,&rowmax,y,f,&delj,threshj,
			  fdel,fdiff,
			  dfdy_colj,
			  &absfvaluerm,
			  &absfdelrm,
			  &absfdiffmax,
			  &infnormdfdy_colj);
    nfc += 1;
#ifdef DBG
    if (lfp) {
      fprintf(lfp,"\n");
      for (i=0;i<ny;i++) {
	fprintf(lfp,"fdel[%d,%d] = %le\n",i,j,fdel[i]);
      }
      fprintf(lfp,"\n");
      for (i=0;i<ny;i++) {
	fprintf(lfp,"fdiff[%d,%d] = %le\n",i,j,fdiff[i]);
      }
      fprintf(lfp,"\n");
      for (i=0;i<ny;i++) {
	fprintf(lfp,"dfdy[%d,%d] = %le\n",i,j,dfdy_colj[i]);
      }
      fprintf(lfp,"\n");
      fprintf(lfp,"rowmax[%d]        = %d\n",j,rowmax);
      fprintf(lfp,"abs(f(rowmax))    = %le\n",absfvaluerm);
      fprintf(lfp,"abs(fdel(rowmax)) = %le\n",absfdelrm);
      fprintf(lfp,"abs(fdiff(rowmax)) = %le\n",absfdiffmax);
      fprintf(lfp,"inf_norm(dfdy[%d]) = %le\n",j,infnormdfdy_colj);
      fflush(lfp);
    }
#endif
    if (success == 0) {
      /*
	Call to num_jac_col failed.
      */
      if (lfp) {
	fprintf(lfp,"ode_num_jac: primary call to num_jac_col failed for j = %d\n",j);
	fflush(lfp);
      }
    } else {
      /* first_call to num_jac_col succeeded.  */
      /*
	fscalej = max(absfvaluerm,absfdelrm);
      */
      fscalej = absfvaluerm;
      if (absfdelrm > fscalej) {
	fscalej = absfdelrm;
      }
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"fscale[%d] = %le\n",j,fscalej);
	fflush(lfp);
      }
#endif    
      /*
	if (j)
      */
      if ((absfdiffmax == 0.0) || ((absfdelrm > 0.0) && (absfvaluerm > 0.0))) {
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," j select true for j = %d\n",j);
	  fflush(lfp);
	}
#endif      
	/*
	  if (k1)
	*/
	if (absfdiffmax <= (br * fscalej)) {
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," k1 select true for j = %d\n",j);
	    fflush(lfp);
	  }
#endif
	  /*
	    tmpfac = min(sqrt(facj),facmax);
	  */
	  tmpfac = sqrt(facj);
	  tmpfac = (tmpfac < facmax) ? tmpfac : facmax ;
	  delj   = (yj + (tmpfac * yscalej)) - yj;
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," tmpfac[%d]  = %le \n",j,tmpfac);
	    fprintf(lfp," new del[%d] = %le \n",j,delj);
	    fflush(lfp);
	  }
#endif	
	  if ((tmpfac != facj) && (delj != 0.0)) {
	    if (fj >= 0.0) {
	      delj = fabs(delj);
	    } else {
	      delj = - fabs(delj);
	    }
#ifdef DBG
	    if (lfp) {
	      fprintf(lfp," Before second num_jac_col\n");
	      fprintf(lfp," del[%d] = %le \n",j,delj);
	      fflush(lfp);
	    }
 #endif	
	    success = num_jac_col(state,ny,j,&rowmax,y,f,&delj,threshj,
				  fdel,fdiff,
				  dfdy_tmp,
				  &absfvaluerm,
				  &absfdelrm_tmp,
				  &absfdiffmax_tmp,
				  &infnormdfdy_tmp);
	    nfc += 1;
	  
#ifdef DBG
	    if (lfp) {
	      fprintf(lfp,"\n");
	      for (i=0;i<ny;i++) {
		fprintf(lfp,"fdel2[%d,%d]    = %le\n",i,j,fdel[i]);
	      }
	      fprintf(lfp,"\n");
	      for (i=0;i<ny;i++) {
		fprintf(lfp,"fdiff2[%d,%d]   = %le\n",i,j,fdiff[i]);
	      }
	      fprintf(lfp,"\n");
	      for (i=0;i<ny;i++) {
		fprintf(lfp,"dfdy_tmp[%d,%d] = %le\n",i,j,dfdy_tmp[i]);
	      }
	      fprintf(lfp,"\n");
	      fprintf(lfp,"rowmax2[%d]       = %d\n",j,rowmax);
	      fprintf(lfp,"abs(f(rowmax))    = %le\n",absfvaluerm);
	      fprintf(lfp,"abs(fdel(rowmax)) = %le\n",absfdelrm_tmp);
	      fprintf(lfp,"abs(fdiff(rowmax)) = %le\n",absfdiffmax_tmp);
	      fprintf(lfp,"inf_norm(dfdy[%d]) = %le\n",j,infnormdfdy_tmp);
	      fflush(lfp);
	    }
#endif
	    if (success == 0) {
	      if (lfp) {
		fprintf(lfp,"ode_num_jac: call to num_jac_col failed for j = %d\n",j);
		fflush(lfp);
	      }
	    } else {
	      /* auxilliary call to num_jac_col succeeded. */
	      if ((tmpfac * infnormdfdy_tmp) > infnormdfdy_colj) {
#ifdef DBG
		if (lfp) {
		  fprintf(lfp," New difference more significant.\n");
		  fflush(lfp);
		}
#endif
		dcopy_(&ny,dfdy_tmp,&index1,dfdy_colj,&index1);
		/*
		  fscaletmp = max(absfvaluerm,absfdelrm_tmp);
		*/
		fscaletmp = (absfvaluerm > absfdelrm_tmp) ? absfvaluerm : absfdelrm_tmp;
		if (absfdiffmax_tmp <= (bl * fscaletmp)) {
		  /*
		    facj  = min(10*tmpfac,facmax);
		  */
		  facj = 10 * tmpfac;
		  facj = (facj < facmax) ? facj : facmax;
		} else {
		  if (absfdiffmax_tmp > (bu * fscaletmp)) {
		    /*
		      facj = max(0.1*tmpfac, facmin);
		    */
		    facj = 0.1 * tmpfac;
		    facj = (facmin > facj) ? facmin : facj ;
		  } else {
		    facj = tmpfac;
		  }
		}
		fac[j] = facj;
#ifdef DBG
		if (lfp) {
		  fprintf(lfp," After new diff more significant fac[%d] = %le\n",
			  j,facj);
		  fflush(lfp);
		}
#endif	 
	      } /*end if  ((tmpfac * infnormdfdy_tmp) > infnormdfdy_colj) */
	    }  /* end else second call to num_jac_col succeeded. */
	  } /* end if ((tmpfac != facj) && (delj != 0.0)) */
	  /* end if (absfdiffmax < (br * fscalej)) (k1) */
	} else {
	  /* 
	    absdiffmax > br * fscale. (j & ~k1)
	  */
	  facj = 10*facj;
	  facj = (facj > facmax) ? facmax : facj;
	  fac[j] = facj;
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," Multiplying fac[%d] by 10 fac[%d] = %le\n",
		    j,j,facj);
	    fflush(lfp);
	  }
#endif	  
	} /* end else absdiffmax > br * fscale. (j & ~k1) */
      } /* end if j */
      if (absfdiffmax > (bu *fscalej)) {
	/*
	  fac[j] = max(0.1*facj,facmin);
	*/
	facj = 0.1 * facj;
	facj = (facmin > facj) ? facmin : facj;
	fac[j] = facj;
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," Dividing fac[%d] by 10 fac[%d] = %le\n",
		  j,j,facj);
	  fflush(lfp);
	}
 #endif
      } /* end if (absfdiffmax > (bu *fscalej)) */
    } /* end else success on first call to num_jac_col */
    dfdy_colj = dfdy_colj + ny; /* Caution Address arithmetic */
  } /* end for j */
  *nfcalls = nfc;
  return (success);
}

