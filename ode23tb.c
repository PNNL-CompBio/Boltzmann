#include "boltzmann_structs.h"

#include "init_base_reactants.h"
#include "approximate_fluxes.h"
#include "compute_flux_scaling.h"
#include "ode_num_jac.h"
#include "ode_it_solve.h"
#include "print_concs_fluxes.h"
#include "blas.h"

#include "ode23tb.h"
int ode23tb (struct state_struct *state, double *counts,
	     double htry, int nonnegative) {

  /*
    NB. This is adapted from the ode23tb matlab code .
    Our mass matrix is the identity matrix.
    species concentrations = y0, species fluxes = f0,
    species concentration = counts .* count_to_conc
    htry is first conc step size, set to zero to set via approximation.
    nonnegative should be set to 1 to specify that concs must be positive.

    Called by: deq_run
    Calls:     init_base_reactants,
               compute_flux_scaling,
               approximate_fluxes, ode_num_jac, ode_it_solve,
               print_concs_fluxes,
	       dnrm2, dgemv, dscal
	       sizeof, calloc, sqrt, pow, max, min, abs

  */
  double *dfdy; /* length nunique_molecules * nunique_molecules */
  double *y0; /* length nunique_molecules */
  double *y;  /* length nunique_molecules */
  double *y1; /* length nunique_molecules */
  double *y2; /* length nunique_molecules */
  double *ynew; /* length nunique_molecules */
  double *yp; /* length nunique_molecules */
  double *y_counts;  /* length nunique_molecules */
  double *y1_counts;  /* length nunique_molecules */
  double *y2_counts; /* length nunique_molecules */
  double *f; /* length nunique_molecules */
  double *f0; /* length nunique_molecules */
  double *f1; /* length nunique_molecules */
  double *z;  /* length nunique_molecules */
  double *z2; /* length nunique_molecules */
  double *z3; /* length nunique_molecules; */
  double *znew; /* length nunique_molecules; */
  double *z_scratch; /* length nunique_molecules */
  double *del_scratch; /* length nunique_molecules */
  double *est1; /* length nunique_molecules */
  double *fac; /* length nunique_molecules */
  double *thresh; /* length nunique_molecules */
  double *delfdelt; /* length nunique_molecules */
  double *ynew_counts; /* length nunique_molecules */
  double *ode_num_jac_scratch; /* length 10*unique_molecules */
  double *forward_rxn_likelihoods; /* length nrxns */
  double *reverse_rxn_likelihoods; /* length nrxns */

  double *count_to_conc;
  double *conc_to_count;
  double *dbl_ptr;
  double t0;
  double t;
  double t2;
  double tdel;
  double tnew;
  double htspan;
  double tfinal;
  double third;
  double alpha;
  double d;
  double gg;
  double c1;
  double c2;
  double c3;
  double p31;
  double p32;
  double p33;
  double root2;
  double eps;
  double sqrt_eps;
  double threshold;
  double h;
  double q;
  double hmin;
  double hmax;
  double tdir;
  double rh;
  double recip_70p;
  double flux_scaling;
  double flux_scaling1;
  double abst;
  double abst_p_h;
  double upper_bound_step;
  double recip_tdel;
  double scalar;
  double recip_cube_root_rtol;
  double wt;
  double normy;
  double normyp;
  double normy2;
  double normynew;
  double err1;
  double err;
  double absh;
  double abslasth;
  double h_ratio;
  double rate;
  double errnn;
  double rtol;
  double norm_delfdelt;
  double err1_est;
  double rtol_o_err;

  int64_t ask_for;
  int64_t one_l;
  int64_t nfevals;
  int64_t nsteps;
  int64_t nfailed;
  int64_t npds;
  int64_t ndecomps;
  int64_t nsolves;
  int64_t nf;
  int64_t eps_hex;
  int64_t sqrt_eps_hex;

  int *base_reactant_indicator; /* length nunique_molecules */
  int *base_reactants;          /* length nunique_molecules */

  int ny;
  int first_time;

  int inc1;
  int not_done;

  int success;
  int base_rxn;

  int need_new_j;
  int jcurrent;

  int one;
  int trans;

  int nrxns;
  int itfail1;

  int itfail2;
  int nofailed;

  int i;
  int unsuccessful_step;

  int tolerance_met;
  int iter_count;

  int nnrejectstep;
  int done;

  int number_base_reaction_reactants;
  int nnreset_znew;


  FILE *lfp;
  FILE *efp;
#define DBG 1
  success = 1;
  t0      = 0.0;
  tfinal  = 10.0;
  tdir    = 1.0;
  one_l   = (int64_t)1;
  rate    = 0.0;
  nsteps   = (int64_t)0;
  nfailed  = nsteps;
  npds     = nsteps;
  ndecomps = nsteps;
  nsolves  = nsteps;
  rtol     = 0.001;
  count_to_conc = state->count_to_conc;
  conc_to_count = state->conc_to_count;
  base_rxn      = state->base_reaction;
  ny            = state->nunique_molecules;
  nrxns         = state->number_reactions;
  lfp           = state->lfp;
  /*
    Allocate double space needed for scratch vectors and matrices.
    actually 26*ny but we'll throw in an extra 7ny for future needs.
  */
  if (success) {
    ask_for = ((ny * 32) + (ny*ny) + (2*nrxns)) * sizeof(double);
    dfdy = (double*)calloc(ask_for,one_l);
    if (dfdy == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"ode23tb: Error could not allocate %ld bytes "
		"for double scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    y0   = &dfdy[ny*ny];
    y    = &y0[ny];
    y1   = &y[ny];
    y2   = &y1[ny];
    ynew = &y2[ny];
    yp   = &ynew[ny];
    y_counts = &yp[ny];
    y1_counts = &y_counts[ny];
    y2_counts = &y1_counts[ny];
    f         = &y2_counts[ny];
    f0        = &f[ny];
    f1        = &f0[ny];
    z         = &f1[ny];
    z2        = &z[ny];
    z3        = &z2[ny];
    znew      = &z3[ny];
    z_scratch = &znew[ny];
    del_scratch = &z_scratch[ny];
    est1        = &del_scratch[ny];
    fac         = &est1[ny];
    thresh      = &fac[ny];
    delfdelt    = &thresh[ny];
    ynew_counts = &delfdelt[ny];
    ode_num_jac_scratch  = &ynew_counts[ny];
    forward_rxn_likelihoods = &ode_num_jac_scratch[4*ny];
    reverse_rxn_likelihoods = &forward_rxn_likelihoods[nrxns];
    /*
      Compute initial concentratiosn from molecule counts.
      Also fill f0 from fluxes input.
    */
    threshold = 1.0e-6;
    for (i=0;i<ny;i++) {
      y0[i] = counts[i] * count_to_conc[i];
      y[i]  = y0[i];
      y_counts[i] = counts[i];
      thresh[i] = threshold;
    }
    /*
      Fill the base_reactant_indicator, and base_reactants vectors,
      and set number_base_reaction_reactants.
    */
    success = init_base_reactants(state);
  }
  if (success) {
    /*
      Compute the scaling flux value.
      Note that flux_scaling is K_f(base_rxn_reaction)*(product of reactant 
      concentrations in base reaction).
    */
    flux_scaling = compute_flux_scaling(state,y0);
    approximate_fluxes(state,y_counts,
		       forward_rxn_likelihoods,reverse_rxn_likelihoods,
		       f0,flux_scaling,base_rxn);
    nfevals = (int64_t)1;
    /*
      Fill the Jacobian. jac is a ny x ny array.
    */
    t0 = 0;
    nf = 0;
    first_time = 1;
    /*
    ode_num_jac(state,dfdy,t0,y0,f0,fac,thresh,&nf,first_time);
    */
    ode_num_jac(state,first_time,
		dfdy,t0,y,f0,fac,thresh,ode_num_jac_scratch, &nf);
    first_time = 0;
    nfevals += nf;
    jcurrent = 1;
    need_new_j = 1-jcurrent;
    t = t0;
    eps_hex       = 0x3CB0000000000000L;
    dbl_ptr       = (double *)&eps_hex;
    eps           = *dbl_ptr;
    sqrt_eps_hex  = 0x3E50000000000000L;
    dbl_ptr       = (double *)&sqrt_eps_hex;
    sqrt_eps      = *dbl_ptr;
    third     = 1.0/3.0;
    root2     = sqrt(2.0);
    alpha     = 2.0 - root2;
    d         = alpha * 0.5;
    gg        = root2 * 0.25;
    c1        = (alpha - 1.0) * third;
    c2        = third;
    c3        = -alpha * third;
    p31       = 1.5 + root2;
    p32       = 2.5 + (root2 + root2);
    p33       = -(6.0 + (4.5*root2));
    inc1      = 1;
    recip_70p = 1.0/0.7;
    /*
    recip_cube_root_rtol = 1.0 / (pow(rtol,third));
    */
    recip_cube_root_rtol = 10.0; 
    for (i=0;i<ny;i++) {
      yp[i] = f0[i];
    }
    hmin      = 0.0;
    hmax      = 0.1;
    normy     = dnrm2(&ny,y,&inc1);
    normyp    = dnrm2(&ny,yp,&inc1);
    need_new_j = 0;
    htspan  = tfinal - t0;
    tdir = 1.0;

    if (htry == 0.0) {
      /*
	Compute an initial step size h using yp = y'(t)
      */
      /*
      wt = max(normy,threshold);
      */
      wt = normy;
      if (threshold > wt) {
	wt = threshold;
      }
      rh = (recip_70p * (normyp/wt)) * recip_cube_root_rtol;
      /*
      absh = min(hmax, htspan);
      */
      absh = htspan;
      if (absh > hmax) {
	absh = hmax;
      }
      if ((absh * rh) > 1.0) {
	absh = 1.0/rh;
      }
      /*
      absh = max(absh,hmin);
      */
      if (hmin > absh) {
	absh = hmin;
      }

      /*
	Estimate error of first order Taylor seris, 0.5*h^2*y''(t)
        and use rule of thumb to select step sizze for second order method.
      */
      h = tdir * absh;
      abst = abs(t);
      abst_p_h = abs(t+h);
      /*
      upper_bound_step = sqrt_eps * max(abst,abst_p_h);
      */
      upper_bound_step = abst_p_h;
      if (abst > upper_bound_step) {
	upper_bound_step = abst;
      }
      upper_bound_step = upper_bound_step * sqrt_eps;

      /*
      tdel = (t + (tdir * min(upper_bound_step,absh))) - t;
      */
      tdel = upper_bound_step;
      if (absh < tdel) {
	tdel = absh;
      }
      tdel = tdir * tdel;
      tdel = (t + tdel) - t;
      /*
	f1   = flux at t+tdel, y + tdel*yp
	So we need to first recompute concentrations y = y + tdel*yp;
	conc = count/volume
      */
    /*
      Now the flux concentration needs likelihoods which needs counts,
      so we need to convert the y1 to counts.
    */
      for (i=0;i<ny;i++) {
	y1[i] = y[i] + tdel * yp[i];
	y1_counts[i] = y1[i] * conc_to_count[i];
      }
      flux_scaling1 = compute_flux_scaling(state,y1);
      approximate_fluxes(state,y1_counts,
			 forward_rxn_likelihoods,reverse_rxn_likelihoods,
			 f1,flux_scaling1,base_rxn);
      nfevals += 1;
      recip_tdel = 1.0/tdel;
      for (i=0;i<ny;i++) {
	/*
	dfdt[i] = (f1[i] - f0[i]) * recip_tdel;
	delfdelt[i] = dfdt[i];
	*/
	delfdelt[i] = (f1[i] - f0[i]) * recip_tdel;
      }
      /*
	delfdelt = dfdt + dfdy * yp
      */
      trans = 0;
      scalar = 1.0;
      dgemv(&trans,&ny,&ny,&scalar,dfdy,&ny,yp,&inc1,&scalar,delfdelt,&inc1);
  
      norm_delfdelt = dnrm2(&ny,delfdelt,&inc1);
      rh = (recip_70p * (norm_delfdelt/wt)) * recip_cube_root_rtol;
      /*
      absh = min(hmax,htspan);
      */
      absh = htspan;
      if (hmax < absh) {
	absh = hmax;
      }
      if ((absh *rh) > 1.0) {
	absh = 1.0/rh;
      }
      /*
      absh = max(absh,hmin);
      */
      if (hmin > absh) {
	absh = hmin;
      }
    } else {
      /*
      htry is not 0.0 but specified.
      absh = min(hmax,max(hmin,htry));
      */
      absh = htry;
      if (hmin > absh) {
	absh = hmin;
      }
      if (hmax < absh) {
	absh = hmax;
      }
    }
    h = tdir * absh;
    for (i=0;i<ny;i++) {
      z[i] = h*yp[i];
    }
#ifdef DBG
    if (lfp) {
      fprintf(lfp," After first call to approximate_fluxes\n");
      print_concs_fluxes(state,ny,f0,y,forward_rxn_likelihoods,reverse_rxn_likelihoods,t0,h);
    }
#endif
    /*
      Main loop.
    */
    not_done = 1;
    while (not_done) {
      hmin = 16.0*eps*abs(t);
      abslasth = absh;
      /*
      absh = min(hmax, max(hmin,absh));
      */
      if (absh < hmin) {
	absh = hmin;
      }
      if (hmax < absh) {
	absh = hmax;
      }
      h = tdir * absh;
      if ((1.1*absh) >= abs(tfinal - t)) {
	h = tfinal - t;
	absh = abs(h);
	not_done = 0;
      }
      if (absh != abslasth) {
	h_ratio = absh/abslasth;
	dscal(&ny,&h_ratio,z,&inc1);
      }
      nofailed = 1;
      unsuccessful_step = 1;
      tolerance_met     = 1;
      /*
	Loop for advancing one step.
      */
      while (unsuccessful_step && tolerance_met) {
	/*
	wt = max(normy,threshold);
	*/
	wt = normy;
	if (threshold > wt) {
	  wt = threshold;
	}
	if (need_new_j) {
	  /*
	    Recompute the jacobian at current y, and f (yp)
	    compute f0. at current y.
	  */
	  flux_scaling = compute_flux_scaling(state,y);
	  approximate_fluxes(state,y_counts,
			     forward_rxn_likelihoods,reverse_rxn_likelihoods,
			     f0,flux_scaling,base_rxn);
	  /*
	    Comput new Jacobian, dfdy.
	  */
	  ode_num_jac(state,first_time,dfdy,t,y,f0,fac,thresh,
		      ode_num_jac_scratch, &nf);
	  nfevals += (nf + 1);
	  npds += 1;
	  jcurrent = 1;
	  need_new_j = 0;
	}
	/*
	  Stage 1.
	*/
	t2 = t + alpha*h;
	for (i=0;i<ny;i++) {
	  y2[i] = y[i] + (alpha * z[i]);
	  y2_counts[i] = y2[i] * conc_to_count[i];
	  z2[i] = z[i];
	}
	iter_count = 0;
	itfail1 = ode_it_solve(state,t2,y2,z2,del_scratch,z_scratch,
			       y2_counts,forward_rxn_likelihoods,
			       reverse_rxn_likelihoods,d,h,rtol,wt,
			       &rate, &iter_count);
	nfevals += iter_count;
	itfail2 = 0;
	if (itfail1 == 0) {
	  /* Stage 2. */
	  normy2 = dnrm2(&ny,y2,&inc1);
	  /*
	  wt     = max(wt,normy2)
	  */
	  if (normy2 > wt) {
	    wt = normy2;
	  }
	  tnew = t + h;
	  if (done) {
	    tnew = tfinal;
	  }
	  for (i=0;i<ny;i++) {
	    znew[i] = (p31 * z[i]) + (p32*z2[i]) + (p33 * (y2[i]-y[i]));
	    ynew[i] = y[i] + (gg * (z[i] + z2[i])) + (d*znew[i]);
	    ynew_counts[i] = ynew[i] * conc_to_count[i];
	  }
	  iter_count = 0;
	  itfail2 = ode_it_solve(state,tnew,ynew,znew,del_scratch,z_scratch,
				 ynew_counts,forward_rxn_likelihoods,
				 reverse_rxn_likelihoods,d,h,rtol,wt,
				 &rate, &iter_count);
	  nfevals += iter_count;
}
	if (itfail1 || itfail2) {
	  nofailed = 0;
	  nfailed += 1;
	  if (jcurrent) {
	    if (absh <= hmin) {
	      if (lfp) {
		fprintf(lfp,"ode23tb: Error integration tolerance not met, t= %le, hmin = %le\n",t,hmin);
		fflush(lfp);
	      }
	      not_done          = 0;
	      unsuccessful_step = 0;
	      tolerance_met = 0;
	    } else {
	      abslasth = absh;
	      /*
	      absh = max(0.3*absh,hmin);
	      */
	      absh = 0.3 * absh;
	      if (hmin > absh) {
		absh = hmin;
	      }
	      h = tdir * absh;
	      h_ratio = absh/abslasth;
	      dscal(&ny,&h_ratio,z,&inc1);
	      done = 0;
	      not_done = 1;
	    }
	  } else {
	    need_new_j = 1 - jcurrent;
	  }
	} else {
	  /* stage 1 and 2 succeeded, estimate local truncation error.*/
	  normynew = dnrm2(&ny,ynew,&inc1);
	  /*
	  wt = max(wt,normynew);
	  */
	  if (normynew > wt) {
	    wt = normynew;
	  }
	  err1 = 0.0;
	  for (i=0;i<ny;i++) {
	    err1_est = (c1 * z[i]) + (c2 * z2[i]) + (c3 * znew[i]); 
	    est1[i] = err1_est;
            /*
	    err1 = max(abs(est1[i]/wt),err1);
	    */
	    err1_est = abs(err1_est/wt);
	    if (err1_est > err1) {
	      err1 = err1_est;
	    }
	  }
	  err = err1/16.0;
	  nnrejectstep = 0;
	  errnn = 0.0;
	  if (nonnegative) {
	    if (err <= rtol) {
	      for (i=0;i<ny;i++) {
		if (ynew[i] < 0.0) {
		  errnn += (ynew[i] * ynew[i]);
		}
	      } 
	      errnn = sqrt(errnn)/wt;
	      if (errnn > rtol) {
		err = errnn;
		nnrejectstep = 1;
	      }
	    }
	  } /* end if nonnegative */
	  if (err > rtol) {
	    /*
	      Step failed.
	    */
	    nfailed = nfailed + 1;
	    if (absh <= hmin) {
	      if (lfp) {
		fprintf(lfp,"ode23tb: Error Integration tol not met, t= %le, hmin = %le\n",t,hmin);
		fflush(lfp);
	      }
	      not_done       =   0;
	      tolerance_met  =   0;
	      unsuccessful_step = 0;
	    } else {
	      /* step succeeded. */
	      nofailed = 0;
	      abslasth = absh;
	      if (nnrejectstep) {
		/*
		absh = max(hmin,0.5*absh);
		*/
		absh = absh * .5;
		if (hmin > absh) {
		  absh = hmin;
		}
	      } else {
		/*
		absh = max(hmin, abslasth * max(0.1,0.7*pow(rtol/err,third)));
		*/
		rtol_o_err = rtol/err;
		absh = 0.7 * pow(rtol_o_err,third);
		if (.1 > absh) {
		  absh = .1;
		}
		absh = absh * abslasth;
		if (hmin > absh) {
		  absh = hmin;
		}
	      }
	      h = tdir * absh;
	      h_ratio = absh/abslasth;
	      for (i=0;i<ny;i++) {
		z[i] *= h_ratio;
	      }
	      not_done = 1;
	      unsuccessful_step = 1;
	    }
	  } else {
	    /* Successful step get out of inner loop. */
	    unsuccessful_step = 0;
	  }
	} /* end else iterative solves succeeded. */
      } /* end while (unsuccessful_step && tolerance_met) */
      nsteps = nsteps + 1;
      nnreset_znew = 0;
      if (nonnegative) {
	for (i=0;i<ny;i++) {
	  z3[i] = znew[i];
	  if (ynew[i] < 0.0) {
	    ynew[i] = 0.0;
	    z3[i] = 0.0;
	    nnreset_znew = 1;
	  }
	}
	if (nnreset_znew) {
	  normynew = dnrm2(&ny,ynew,&inc1);
	}
      }
#ifdef DBG
      if (lfp) {
	print_concs_fluxes(state,ny,f0,ynew,forward_rxn_likelihoods,
			   reverse_rxn_likelihoods,tnew,h);
      }
#endif
      if (nnreset_znew) {
	for (i=0;i<ny;i++) {
	  znew[i] = z3[i];
	}
      }
      t = tnew;
      for (i=0;i<ny;i++) {
	y[i] = ynew[i];
	z[i] = znew[i];
	y_counts[i] = y[i] * conc_to_count[i];
      }
      if (not_done) {
	normy = normynew;
	jcurrent = 0;
	need_new_j = 1;
	if (nofailed) {
	  q = pow((err/rtol),third);
	  h_ratio = hmax/absh;
	  if (0.7 < (q * h_ratio)) {
	    h_ratio = 0.7/q;
	  }
	  /*
	    h_ratio = min(5,max(0.2,h_ratio));
	  */
	  if (.2 > h_ratio) {
	    h_ratio = .2;
	  }
	  if (h_ratio > 5.0) {
	    h_ratio = 5.0;
	  }
	  if (abs(h_ratio - 1) > 0.2) {
	    absh = h_ratio *absh;
	    dscal(&ny,&h_ratio,z,&inc1);
	  }
	}
      }
    } /* end while (not_done) */
  } /* end if success - allocation succeeded */
  for(i=0;i<ny;i++) {
    counts[i] = y_counts[i];
  }
  return (success);
}

