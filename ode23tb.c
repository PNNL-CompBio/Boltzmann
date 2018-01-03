#include "boltzmann_structs.h"

#include "init_base_reactants.h"
#include "approximate_delta_concs.h"
/*#include "ce_approximate_delta_concs.h" */
/*#include "lr_approximate_delta_concs.h" */
#include "compute_flux_scaling.h"
#include "ode_num_jac.h"
#include "ode_it_solve.h"
#include "ode_print_concs.h"
#include "blas.h"
#include "lapack.h"
/*
#define DBG 1
*/
#ifdef DBG 
#include "print_concs_fluxes.h"
#endif

#include "ode23tb.h"
int ode23tb (struct state_struct *state, double *counts,
	     double htry, int nonnegative, int normcontrol,
	     int print_concs, int choice) {

  /*
    NB. This is adapted from the ode23tb matlab code .
    Our mass matrix is the identity matrix.
    species concentrations = y0, species concentration changes wrt time = f0,
    species concentration = counts .* count_to_conc
    htry is first conc step size, set to zero to set via approximation.
    nonnegative should be set to 1 to specify that concs must be positive.

    choice is not used by ode23tb but is a method selector for its
    parent, ode_solver, and might be used in printing debug information.
    

    Called by: ode_solver
    Calls:     init_base_reactants,
               compute_flux_scaling,
               approximate_fluxes, ode_num_jac, ode_it_solve,
               print_concs_fluxes,
	       dnrm2_, dgemv_, dscal_, idamax_
	       sizeof, calloc, sqrt, pow, fabs, dgetrf_, dgetrs_

  */
  double *dfdy; /* length nunique_molecules * nunique_molecules */
  double *miter; /* length nunique_molecules * nunique_molecules */
  double *y0; /* length nunique_molecules */
  double *y;  /* length nunique_molecules */
  double *y1; /* length nunique_molecules */
  double *y2; /* length nunique_molecules */
  double *ynew; /* length nunique_molecules */
  double *yp; /* length nunique_molecules */
  double *absy; /* length nunique_molecules */
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
  double *it_solve_rhs; /* length nunique_molecules */
  double *it_solve_del; /* length nunique_molecules */
  double *it_solve_scr; /* length nunique_molecules */
  double *est1; /* length nunique_molecules */
  double *est2; /* length nunique_molecules */
  double *fac; /* length nunique_molecules */
  double *thresh; /* length nunique_molecules */
  double *delfdelt; /* length nunique_molecules */
  double *ynew_counts; /* length nunique_molecules */
  double *wti; /* length nunique_molecules */
  double *pivot; /* length nunique_molecules */
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
  double njthreshold;
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
  double wt2;
  double normy;
  double normyp;
  double normy2;
  double normynew;
  double err1;
  double err2;
  double err;
  double absh;
  double abslasth;
  double h_ratio;
  double rate;
  double errnn;
  double rtol;
  double norm_delfdelt;
  double err1_est;
  double err2_est;
  double rtol_o_err;
  double min_conc;
  double yp_o_wt;
  double max_yp_o_wt;
  double norm_delfdelt_o_wt;
  double norm_delfdelt_o_wti;
  double dzero;
  double mdh;

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
  int *ipivot;                  /* overlaid on pivot space */

  int ny;
  int first_time;

  int inc1;
  int not_done;

  int success;
  int base_rxn;

  int need_new_j;
  int need_new_m;

  int need_new_lu;
  int i;

  int jcurrent;
  int mcurrent;


  int one;
  int nrhs;

  int nrxns;
  int nofailed;

  int itfail1;
  int itfail2;

  int unsuccessful_step;
  int tolerance_met;

  int iter_count;
  int nnrejectstep;

  int number_base_reaction_reactants;
  int nnreset_znew;

  int origin;
  int done;

  int nysq;
  int info;

  int ode_solver_choice;
  int delta_concs_choice;

  char  trans_chars[8];
  char  *trans;

  FILE *lfp;
  FILE *efp;
  success = 1;
  t0      = 0.0;
  /*
  tfinal  = 10.0;
  */
  tfinal  = state->ode_t_final;
  tdir    = 1.0;
  one_l   = (int64_t)1;
  rate    = 0.0;
  nsteps   = (int64_t)0;
  nfailed  = nsteps;
  npds     = nsteps;
  ndecomps = 0;
  nsolves  = 0;
  rtol     = 0.001;
  dzero    = 0.0;
  count_to_conc = state->count_to_conc;
  conc_to_count = state->conc_to_count;
  base_rxn      = state->base_reaction;
  ny            = state->nunique_molecules;
  delta_concs_choice = (int)state->delta_concs_choice;
  nrxns         = state->number_reactions;
  min_conc      = state->min_conc;
  lfp           = state->lfp;
  nysq          = ny * ny;
  trans_chars[0] = 'N';
  trans_chars[1] = 'T';
  trans_chars[3] = 'C';
  trans = &trans_chars[0];
  /*
    Allocate double space needed for scratch vectors and matrices.
    actually 31*ny but we'll throw in an extra 9ny for future needs.
  */
  if (success) {
    ask_for = ((ny * 40) + (2*ny*ny) + (2*nrxns)) * sizeof(double);
    dfdy = (double*)calloc(ask_for,one_l);
    if (dfdy == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"ode23tb: Error could not allocate %lld bytes "
		"for double scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    miter     = &dfdy[nysq];
    y0        = &miter[nysq];
    y         = &y0[ny];
    y1        = &y[ny];
    y2        = &y1[ny];
    ynew      = &y2[ny];
    yp        = &ynew[ny];
    absy      = &yp[ny];
    y_counts  = &absy[ny];
    y1_counts = &y_counts[ny];
    y2_counts = &y1_counts[ny];
    f         = &y2_counts[ny];
    f0        = &f[ny];
    f1        = &f0[ny];
    z         = &f1[ny];
    z2        = &z[ny];
    z3        = &z2[ny];
    znew      = &z3[ny];
    it_solve_rhs = &znew[ny];
    it_solve_del = &it_solve_rhs[ny];
    it_solve_scr = &it_solve_del[ny];
    est1        = &it_solve_scr[ny];
    est2        = &est1[ny];
    fac         = &est2[ny];
    thresh      = &fac[ny];
    delfdelt    = &thresh[ny];
    ynew_counts = &delfdelt[ny];
    wti         = &ynew_counts[ny];
    pivot       = &wti[ny];
    ipivot      = (int*)pivot;
    ode_num_jac_scratch  = &pivot[ny];
    forward_rxn_likelihoods = &ode_num_jac_scratch[5*ny];
    reverse_rxn_likelihoods = &forward_rxn_likelihoods[nrxns];
    /*
      We need to move from stochastic counts to contiuous counts
      so as not to have zero concentrations.
    */
    /*
      Compute initial concentratiosn from molecule counts.
      Also fill f0 from fluxes input.
    */
    threshold = 1.0e-3;
    njthreshold = 1.0e-6;
    
    
    for (i=0;i<ny;i++) {
      y0[i] = counts[i] * count_to_conc[i];
      y[i]  = y0[i];
      y_counts[i] = counts[i];
      thresh[i] = njthreshold;
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
    approximate_delta_concs(state,y_counts,
			    forward_rxn_likelihoods,reverse_rxn_likelihoods,
			    f0,flux_scaling,base_rxn,
			    delta_concs_choice);
#ifdef DBG
    if (lfp) {
      fprintf(lfp," ode23tb: After first call to approximate_delta_concs, flux_scaling = %le\n",flux_scaling);
      origin = 0; 
      t = 0;
      h = 0;
      print_concs_fluxes(state,ny,f0,y,y_counts,
			 forward_rxn_likelihoods,
			 reverse_rxn_likelihoods,t,h,nsteps,origin);
    }
#endif
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
    mcurrent = 1;
    need_new_m = 0;
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
    /*
    recip_70p = 1.0/0.7;
    to match matlab we use
    */
    recip_70p = 1.43;
    /*
    recip_cube_root_rtol = 1.0 / (pow(rtol,third));
    */
    recip_cube_root_rtol = 10.0; 
    dcopy_(&ny,f0,&inc1,yp,&inc1);
    /*
    for (i=0;i<ny;i++) {
      yp[i] = f0[i];
    }
    */
    hmin      = 0;
    hmax      = 1;
    if (normcontrol) {
      normy     = dnrm2_(&ny,y,&inc1);
      normyp    = dnrm2_(&ny,yp,&inc1);
    } else {
      normy     = fabs(y[idamax_(&ny,y,&inc1)-1]);
      normyp    = fabs(yp[idamax_(&ny,yp,&inc1)-1]);
    }
    jcurrent   = 1;
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
      if (normcontrol) {
	wt2 = normy;
	if (threshold > wt2) {
	  wt2 = threshold;
	}
	for (i=0;i<ny;i++) {
	  wti[i] = wt2;
	}
	rh = (recip_70p * (normyp/wt2)) * recip_cube_root_rtol;
      } else {
	max_yp_o_wt = 0.0;
	for (i=0;i<ny;i++) {
	  wt = fabs(y[i]);
	  if (threshold > wt) {
	    wt = threshold;
	  }
	  wti[i] = wt;
	  yp_o_wt = fabs(yp[i])/wt;
	  if (yp_o_wt > max_yp_o_wt) {
	    max_yp_o_wt = yp_o_wt;
	  }
	}
	rh = recip_70p * max_yp_o_wt * recip_cube_root_rtol;
      }
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
      abst = fabs(t);
      abst_p_h = fabs(t+h);
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
	f1   = flux at t+tdel, y
	conc = count/volume
      */
      /*
	Now the flux concentration needs likelihoods which needs counts,
	so we need to convert the y1 to counts.
      */
      flux_scaling1 = compute_flux_scaling(state,y1);
      approximate_delta_concs(state,y_counts,
			 forward_rxn_likelihoods,reverse_rxn_likelihoods,
			      f1,flux_scaling1,base_rxn,delta_concs_choice);
#ifdef DBG
      if (lfp) {
	fprintf(lfp,"ode23tb: After call to approximate_delta_concs in f1, flux_scaling = %le\n",flux_scaling);
	origin = 1; 
	print_concs_fluxes(state,ny,f1,y1,y1_counts,
			   forward_rxn_likelihoods,
			   reverse_rxn_likelihoods,t,h,nsteps,origin);
      }
#endif
      nfevals += 1;
      recip_tdel = 1.0/tdel;
      for (i=0;i<ny;i++) {
	/*
	dfdt[i] = (f1[i] - f0[i]) * recip_tdel;
	delfdelt[i] = dfdt[i];
	delfdelt[i] = (f1[i] - f0[i]) * recip_tdel;
	*/
	delfdelt[i] = 0.0;
      }
      /*
	delfdelt = dfdt + dfdy * yp
      */
      scalar = 1.0;
      dgemv_(trans,&ny,&ny,&scalar,dfdy,&ny,yp,&inc1,&scalar,delfdelt,&inc1);
  
      if (normcontrol) {
	norm_delfdelt = dnrm2_(&ny,delfdelt,&inc1);
	norm_delfdelt_o_wt = norm_delfdelt/wt2;
      } else {
	norm_delfdelt_o_wt = 0.0;
	for (i=0;i<ny;i++) {
	  norm_delfdelt_o_wti = fabs(delfdelt[i])/wti[i];
	  if (norm_delfdelt_o_wti > norm_delfdelt_o_wt) {
	    norm_delfdelt_o_wt = norm_delfdelt_o_wti;
	  }
	}
      }
      rh = (recip_70p * sqrt(.5*norm_delfdelt_o_wt)) * recip_cube_root_rtol;
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
    /*
      Main loop.
    */
    need_new_lu = 1;
    not_done = 1;
    done     = 0;
    info = 0;
    while (not_done) {
      hmin = 16.0*eps*fabs(t);
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
      if ((1.1*absh) >= fabs(tfinal - t)) {
	h = tfinal - t;
	absh = fabs(h);
	not_done = 0;
	done     = 1;
      }
      if (absh != abslasth) {
	h_ratio = absh/abslasth;
	dscal_(&ny,&h_ratio,z,&inc1);
	need_new_lu = 1;
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
	if (normcontrol) {
	  wt2 = normy;
	  if (threshold > wt2) {
	    wt2 = threshold;
	  }
	  for (i=0;i<ny;i++) {
	    wti[i] = wt2;
	  }
	} else {
	  for (i=0;i<ny;i++) {
	    wt = fabs(y[i]);
	    if (threshold > wt) {
	      wt = threshold;
	    }
	    wti[i] = wt;
	  }
	}
	if (need_new_j) {
	  /*
	    Recompute the jacobian at current y, and f (yp)
	    compute f0. at current y.
	  */
	  flux_scaling = compute_flux_scaling(state,y);
	  approximate_delta_concs(state,y_counts,
			     forward_rxn_likelihoods,reverse_rxn_likelihoods,
				  f0,flux_scaling,base_rxn,delta_concs_choice);
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp,"ode23tb: after new Jacobian call to approximate_delta_concs, flux_scaling = %le\n",flux_scaling);
	    origin = 2; 
	    print_concs_fluxes(state,ny,f0,y,y_counts,
			       forward_rxn_likelihoods,
			       reverse_rxn_likelihoods,t,h,nsteps,origin);
	  }
#endif

	  /*
	    Compute new Jacobian, dfdy.
	  */
	  ode_num_jac(state,first_time,dfdy,t,y,f0,fac,thresh,
		      ode_num_jac_scratch, &nf);
	  nfevals += (nf + 1);
	  npds += 1;
	  jcurrent = 1;
	  need_new_j = 0;
	  need_new_lu = 1;
	}
	if (need_new_m) {
	  mcurrent = 1;
	  need_new_m = 0;
	  need_new_lu = 1;
	}
	if (need_new_lu) {
	  /*
	    This could be its own subroutine say ode_build_miter:
	  */
	  /*
	    Set miter to be all 0.
	  */
	  dscal_(&nysq,&dzero,miter,&inc1);
	  /*
	    Set miter to be the identity matrix.
	  */
	  for (i=0;i<nysq;i += (ny+1)) {
	    miter[i] = 1.0;
	  }
	  /*
	    set miter to be I - (d*h)*dfdy
	  */
	  mdh = 0.0 - (d*h);
	  daxpy_(&nysq,&mdh,dfdy,&inc1,miter,&inc1);
	  /*
	    Now we need to factor miter with say dgetrf (lapack routines)
	    overwrites miter with its factorization, generating
	    a pivot vector, ipivot,
	    as arguments to ode_it_solve, 
	  */
	  info = 0;
	  dgetrf_(&ny,&ny,miter,&ny,ipivot,&info);
	  if (info != 0) {
	    if (lfp) {
	      fprintf(lfp,"ode23tb: dgetrf_ call failed with info = %d\n",
		      info);
	      fflush(lfp);
	    }
	    break;
	  }
	  ndecomps = ndecomps + 1;
	  rate = 0.0;
	  need_new_lu = 0;
	}
	/*
	  Stage 1.
	*/
	t2 = t + alpha*h;
	for (i=0;i<ny;i++) {
	  y2[i] = y[i] + (alpha * z[i]);
	  if (y2[i] < min_conc) {
	    if (lfp) {
	      fprintf(lfp,"ode23tb: Warning y2[%d] was < 0, setting to %le\n",
		      i,min_conc);
	    } else {
	      fprintf(lfp,"ode23tb: Warning y2[%d] was < %le, setting to %le\n",
		      i,min_conc,min_conc);
	    }
	    y2[i] = min_conc;
	  }
	  y2_counts[i] = y2[i] * conc_to_count[i];
	  z2[i] = z[i];
	}
	/*
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," Stage1 y2, y2_counts, and z2 before ode_it_solve\n");
	  origin = 3; 
	  print_concs_fluxes(state,ny,z2,y2,y2_counts,
			     forward_rxn_likelihoods,
			     reverse_rxn_likelihoods,t,h,nsteps,origin);
	}
#endif
	*/
	iter_count = 0;
	itfail1 = ode_it_solve(state,miter,ipivot,t2,y2,z2,
			       it_solve_del,it_solve_rhs,it_solve_scr,
			       y2_counts,forward_rxn_likelihoods,
			       reverse_rxn_likelihoods,d,h,rtol,wti,
			       &rate, &iter_count);
	/*
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," Stage1 y2, y2_counts, and z2 after ode_it_solve\n");
	  origin = 4; 
	  print_concs_fluxes(state,ny,z2,y2,y2_counts,
			     forward_rxn_likelihoods,
			     reverse_rxn_likelihoods,t,h,nsteps,origin);
	}
#endif
	*/
	nfevals += iter_count;
	nsolves += iter_count;
	itfail2 = 0;
	if (itfail1 == 0) {
	  /* Stage 2. */
	  if (normcontrol) {
	    normy2 = dnrm2_(&ny,y2,&inc1);
	    /*
	      wt     = max(wt,normy2)
	    */
	    if (normy2 > wt2) {
	      wt2 = normy2;
	      for (i=0;i<ny;i++) {
		wti[i] = wt2;
	      }
	    }
	  } else {
	    for (i=0;i<ny;i++) {
	      wt = fabs(y2[i]);
	      if (wt > wti[i]) {
		wti[i] = wt;
	      }
	    }
	  }
	  tnew = t + h;
	  if (done) {
	    tnew = tfinal;
	  }
	  for (i=0;i<ny;i++) {
	    znew[i] = (p31 * z[i]) + (p32*z2[i]) + (p33 * (y2[i]-y[i]));
	    ynew[i] = y[i] + (gg * (z[i] + z2[i])) + (d*znew[i]);
	    if (ynew[i] < min_conc) {
	      if (lfp) {
		if (ynew[i] < 0) {
		  fprintf(lfp,"ode23tb: Warning ynew[%d] was < 0, setting to %le\n",
			  i,min_conc);
		} else {
		  fprintf(lfp,"ode23tb: Warning ynew[%d] was < %le, setting to %le\n",
			  i,min_conc,min_conc);
		}
	      }
	      ynew[i] = min_conc;
	    }
	    ynew_counts[i] = ynew[i] * conc_to_count[i];
	  }
	  iter_count = 0;
	  /*
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," Stage2 ynew, ynew_counts, and znew "
		    "before ode_it_solve\n");
	    origin = 5; 
	    print_concs_fluxes(state,ny,znew,ynew,ynew_counts,
			       forward_rxn_likelihoods,
			       reverse_rxn_likelihoods,t,h,nsteps,origin);
	  }
#endif
	  */
	  itfail2 = ode_it_solve(state,miter,ipivot,tnew,ynew,znew,
				 it_solve_del,it_solve_rhs,it_solve_scr,
				 ynew_counts,forward_rxn_likelihoods,
				 reverse_rxn_likelihoods,d,h,rtol,wti,
				 &rate, &iter_count);
	  nfevals += iter_count;
	  nsolves += iter_count;
	  /*
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," Stage2 ynew, ynew_counts, and znew "
		    "after ode_it_solve\n");
	    origin = 6; 
	    print_concs_fluxes(state,ny,znew,ynew,ynew_counts,
			       forward_rxn_likelihoods,
			       reverse_rxn_likelihoods,t,h,nsteps,origin);
	  }
#endif
	  */
	} /* end if (itfail1 == 0) */
	if (itfail1 || itfail2) {
	  nofailed = 0;
	  nfailed += 1;
	  if (jcurrent & mcurrent) {
	    if (absh <= hmin) {
	      if (lfp) {
		fprintf(lfp,"ode23tb: Error integration tolerance not met, t= %le, hmin = %le\n",t,hmin);
		fflush(lfp);
	      }
	      not_done          = 0;
	      done              = 1;
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
	      dscal_(&ny,&h_ratio,z,&inc1);
	      need_new_lu = 1;
	      not_done = 1;
	      done = 0;
	    }
	  } else {
	    need_new_j = 1 - jcurrent;
	    need_new_m = 1 - mcurrent;
	  }
	} else {
	  /* stage 1 and 2 succeeded, estimate local truncation error.*/
	  if (normcontrol) {
	    normynew = dnrm2_(&ny,ynew,&inc1);
	  /*
	  wt2 = max(wt,normynew);
	  */
	    if (normynew > wt2) {
	      wt2 = normynew;
	    }
	    for (i=0;i<ny;i++) {
	      wti[i] = wt2;
	    }
	  } else {
	    normynew = 0.0;
	    for (i=0;i<ny;i++) {
	      wt = fabs(ynew[i]);
	      if (wt > normynew) {
		normynew= wt;
	      }
	      if (wt > wti[i]) {
		wti[i] = wt;
	      }
	    }
	  }
	  err1 = 0.0;
	  for (i=0;i<ny;i++) {
	    err1_est = (c1 * z[i]) + (c2 * z2[i]) + (c3 * znew[i]); 
	    est1[i] = err1_est;
            /*
	    err1 = max(fabs(est1[i]/wt),err1);
	    */
	    err1_est = fabs(err1_est/wti[i]);
	    if (err1_est > err1) {
	      err1 = err1_est;
	    }
	  }
	  /*
	    Here we need to solve miter * est2 = est1
	    With calls to something like dgetrs_
	  */
	  nrhs = 1;
	  dcopy_(&ny,est1,&inc1,est2,&inc1);
	  dgetrs_(trans,&ny,&nrhs,miter,&ny,ipivot,est2,&ny,&info);
	  if (info != 0) {
	    if (lfp) {
	      fprintf(lfp,"ode23tb: dgetrs_ call failed with info = %d\n",
		      info);
	      fflush(lfp);
	    }
	    break;
	  }
	  nsolves += 1;
	  err2 = 0.0;
	  for (i=0;i<ny;i++) {
	    err2_est = fabs(est2[i]/wti[i]);
	    if (err2_est > err2) {
	      err2 = err2_est;
	    }
	  }
	  err = err1/16.0;
	  if (err2 > err) {
	    err = err2;
	  }
	  nnrejectstep = 0;
	  errnn = 0.0;
	  if (nonnegative) {
	    if (err <= rtol) {
	      for (i=0;i<ny;i++) {
		if (ynew[i] < 0.0) {
		  errnn += (ynew[i] * ynew[i]);
		}
	      }
	      errnn = sqrt(errnn);
	      if (normcontrol) {
		errnn = errnn/wt2;
	      } else {
		errnn  = errnn/threshold;
	      }
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
	      done           =   1;
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
	      need_new_lu = 1;
	      not_done = 1;
	      done     = 0;
	      unsuccessful_step = 1;
	    }
	  } else {
	    /* Successful step get out of inner loop. */
	    unsuccessful_step = 0;
	  }
	} /* end else iterative solves succeeded. */
      } /* end while (unsuccessful_step && tolerance_met) */
      if (info != 0) {
	break;
      }

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
	  if (normcontrol) {
	    normynew = dnrm2_(&ny,ynew,&inc1);
	  } else {
	    normynew = fabs(ynew[idamax_(&ny,ynew,&inc1)-1]);
	  }
	}
      }
      /* 
	Advance the integration one step. 
      */
      if (nnreset_znew) {
	for (i=0;i<ny;i++) {
	  znew[i] = z3[i];
	}
      }
      t = tnew;
      for (i=0;i<ny;i++) {
	y[i] = ynew[i];
	z[i] = znew[i];
	if (y[i] < min_conc) {
	  if (lfp) {
	    if (y[i] < 0) {
	      fprintf(lfp,"ode23tb: Warning y[%d] was < 0, setting to %le\n",
		      i,min_conc);
	    } else {
	      fprintf(lfp,"ode23tb: Warning y[%d] was < %le, setting to %le\n",
		      i,min_conc,min_conc);
	    }
	  }
	  y[i] = min_conc;
	}
	y_counts[i] = y[i] * conc_to_count[i];
      }
      if (print_concs) {
	ode_print_concs(state,t,y);
      }
      if (not_done) {
	if (normcontrol) {
	  normy = normynew;
	} 
	jcurrent = 0;
	mcurrent = 1;
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
	  if (fabs(h_ratio - 1) > 0.2) {
	    absh = h_ratio *absh;
	    dscal_(&ny,&h_ratio,z,&inc1);
	    need_new_lu = 1;
	  }
	}
      }
    } /* end while (not_done) MAIN loop */
  } /* end if success - allocation succeeded */
  for(i=0;i<ny;i++) {
    counts[i] = y_counts[i];
  }
  return (success);
}

