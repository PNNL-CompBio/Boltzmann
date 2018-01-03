#include "boltzmann_structs.h"

#include "approximate_delta_concs.h"
#include "ode_num_jac.h"
#include "ode_it_solve.h"
#include "ode_print_concs.h"
#include "ode_print_dconcs.h"
#include "blas.h"
#include "lapack.h"
#include "ode23tb_normyp_o_wt.h"
#include "ode23tb_limit_h.h"
#include "ode23tb_init_wt.h"
#include "ode23tb_update_wt.h"
#include "vec_set_constant.h"
#include "ode23tb_build_factor_miter.h"
#include "ode23tb_max_abs_ratio.h"
#include "ode23tb_nonneg_err.h"
#include "ode23tb_enforce_nonneg.h"
#include "boltzmann_monitor_ode.h"
#include "boltzmann_size_jacobian.h"
/*
#include "update_rxn_likelihoods.h"
#include "get_counts.h"
#include "ode_print_lklhds.h"
#include "ode_print_bflux.h"
#include "compute_net_likelihoods.h"
#include "compute_net_lklhd_bndry_flux.h"
#include "print_net_likelihood_header.h"
#include "print_net_lklhd_bndry_flux_header.h"
#include "print_net_likelihoods.h"
#include "print_net_lklhd_bndry_flux.h"
*/
/*
#define DBG 1
*/

#ifdef DBG 
#include "print_concs_dconcs.h"
#endif

#include "ode23tb.h"
int ode23tb (struct state_struct *state, double *concs) {

  /*
    NB. This is adapted from the ode23tb matlab code .
    Our mass matrix is the identity matrix.
    species concentrations = y, species concentration changes wrt time = f0,
    species concentration = counts .* count_to_conc
    htry is first conc step size, set to zero to set via approximation.
    nonnegative should be set to 1 to specify that concs must be positive.

    choice is not used by ode23tb but is a method selector for its
    parent, ode_solver, and might be used in printing debug information.
    

    Called by: ode_solver
    Calls:     approximate_delta_concs, ode_num_jac, ode_it_solve,
               ode_print_concs, ode_norm_yp_o_wt, ode23tb_limit_h,
	       ode23tb_init_wt, od23tb_update_wt, vec_set_constant,
	       ode23tb_build_factor_miter, ode23tb_max_abs_ration,
	       ode23tb_nonneg_err, ode23tb_enforce_nonneg,
	       boltzmann_monitor_ode, boltzmann_size_jacobian,
	       dcopy_, dnrm2_, dgemv_, dscal_, idamax_
	       sizeof, calloc, sqrt, pow, fabs, dgetrf_, dgetrs_

  */
  struct ode23tb_params_struct *ode23tb_params;
  struct cvodes_params_struct *cvodes_params;
  double *dfdy; /* length nunique_molecules * nunique_molecules */
  double *miter; /* length nunique_molecules * nunique_molecules */
  double *y;  /* length nunique_molecules */
  double *y2; /* length nunique_molecules */
  double *ynew; /* length nunique_molecules */
  double *yp; /* length nunique_molecules */
  double *f0; /* length nunique_molecules */
  /*double *f1;*/ /* length nunique_molecules */
  double *z;  /* length nunique_molecules */
  double *z2; /* length nunique_molecules */
  double *znew; /* length nunique_molecules; */
  double *it_solve_rhs; /* length nunique_molecules */
  double *it_solve_del; /* length nunique_molecules */
  double *it_solve_scr; /* length nunique_molecules */
  double *est; /* length nunique_molecules */
  double *fac; /* length nunique_molecules */
  double *thresh; /* length nunique_molecules */
  double *delfdelt; /* length nunique_molecules */
  double *wt; /* length nunique_molecules */
  double *pivot; /* length nunique_molecules */
  double *net_lklhd_bndry_flux;    /* length unique_molecules.*/
  double *fdel; /* length unique_molecules */                
  double *fdiff; /* length unique_molecules */
  double *dfdy_tmp; /* length unique_molecules */

  double *drfc; /* length 2*number_molecules. */
  double *dfdy_row; /* length ny */
  double *dfdy_a; /* length nnz */
  double *dfdy_at; /* length nnz */

  double *forward_rxn_likelihoods; /* length nrxns */
  double *reverse_rxn_likelihoods; /* length nrxns */
  /*double *net_likelihoods;        */ /* length nrxns */


  double *count_to_conc; /* length unique_molecules */
  double *conc_to_count; /* length unique_molecules */
  double *counts;
  double *dbl_ptr;

  double t0;
  double t;
  double t2;
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
  /*
  double t1;
  double tdel;
  double abst;
  double abst_p_h;
  double upper_bound_step;
  double recip_tdel;
  */
  double scalar;
  double recip_cube_root_rtol;
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
  double rtol;
  double norm_delfdelt;
  double errnn_scale;
  double rtol_o_err;
  double min_conc;
  double norm_delfdelt_o_wt;
  double dzero;
  double htry;
  double znew_coeff[4];
  double ynew_coeff[5];
  double est_coeff[4];

  int64_t ask_for;
  int64_t one_l;
  int64_t zero_l;
  int64_t nfevals;
  int64_t nsteps;
  int64_t nfailed;
  int64_t npds;
  int64_t ndecomps;
  int64_t nsolves;
  int64_t nf;
  int64_t eps_hex;
  int64_t sqrt_eps_hex;
  int64_t ode_rxn_view_freq;
  int64_t ode_rxn_view_step;

  int *dfdy_ia;
  int *dfdy_iat;
  int *dfdy_ja;
  int *dfdy_jat;
  int *column_mask;

  int *ipivot;                  /* overlaid on pivot space */

  int    ifour;
  int    ifive;

  int ny;
  int first_time;

  int inc1;
  int not_done;

  int success;
  int two_ny;  
  /*
  int nl_success;
  */

  int need_new_j;
  int need_new_m;

  int need_new_lu;
  int i;

  int jcurrent;
  int mcurrent;


  int done;
  int nrhs;

  int nrxns;
  int nofailed;

  int itfail1;
  int itfail2;

  int unsuccessful_step;
  int tolerance_met;

  int iter_count;
  int nnrejectstep;

#ifdef DBG
  int origin;
  int padj;
#endif

  int nysq;
  int info;

  int delta_concs_choice;
  int nnreset_znew;

  int print_output;
  int ierr;
  
  int normcontrol;
  int nonnegative;

  int number_molecules;
  int nnz;

  int num_ints;
  int ode_jacobian_choice;

  int num_doubles;
  int print_concs;

  int drfc_len;
  int ia_len;

  char  trans_chars[8];
  char  *trans;



  FILE *lfp;
  FILE *ode_dconcs_fp;

  success = 1;
  ode23tb_params = state->ode23tb_params;
  cvodes_params  = state->cvodes_params;
  ode_jacobian_choice = state->ode_jacobian_choice;
  t0      = 0.0;
  tnew    = t0;
  /*
  tfinal  = 10.0;
  */
  tfinal      = state->ode_t_final;
  print_concs = state->print_ode_concs;
  tdir    = 1.0;
  one_l   = (int64_t)1;
  zero_l  = (int64_t)0;
  rate    = 0.0;
  nsteps   = (int64_t)0;
  nfailed  = nsteps;
  npds     = nsteps;
  nonnegative = 1;
  normcontrol = 0;
  htry        = 0.0;
  ode23tb_params->nonnegative = nonnegative;
  ode23tb_params->normcontrol = normcontrol;
  ode23tb_params->htry        = htry;
  /*
    If we want to control these parameters externally we could
    set them in read_params:
    nonnegative = ode23tb_params->nonnegative;
    normcontrol = ode23tb_params->normcontrol;
    htry        = ode23tb_params->htry;
  */
  ndecomps = 0;
  nsolves  = 0;
  rtol     = 0.001;
  dzero    = 0.0;
  nnreset_znew  = 0;
  count_to_conc = state->count_to_conc;
  conc_to_count = state->conc_to_count;
  ny            = state->nunique_molecules;
  print_output  = state->print_output;
  two_ny        = ny + ny;
  nysq          = ny * ny;
  delta_concs_choice = (int)state->delta_concs_choice;
  nrxns         = state->number_reactions;
  /*
  min_conc      = state->min_conc;
  */
  min_conc      = 0.0;
  lfp           = state->lfp;
  ode_dconcs_fp = state->ode_dconcs_fp;
  ode_rxn_view_freq = state->ode_rxn_view_freq;
  if (print_output == 0) {
    ode_rxn_view_freq = 0;
    state->ode_rxn_view_freq = 0;
  }
  trans_chars[0] = 'N';
  trans_chars[1] = 'T';
  trans_chars[3] = 'C';
  trans = &trans_chars[0];
  ifour         = 4;
  ifive         = 5;
  inc1          = 1;

  dfdy         		  = NULL;
  miter        		  = NULL;
  y            		  = NULL;
  z            		  = NULL;
  y2           		  = NULL;
  z2           		  = NULL;
  znew         		  = NULL;
  ynew         		  = NULL;
  yp           		  = NULL;
  f0           		  = NULL;
  it_solve_rhs 		  = NULL;
  it_solve_del 		  = NULL;
  it_solve_scr 		  = NULL;
  est          		  = NULL;
  fac          		  = NULL;
  thresh       		  = NULL;
  delfdelt     		  = NULL;
  wt           		  = NULL;
  pivot        		  = NULL;
  ipivot       		  = NULL;    
  net_lklhd_bndry_flux 	  = NULL;
  fdel                 	  = NULL;
  fdiff                	  = NULL;
  dfdy_tmp             	  = NULL;
  forward_rxn_likelihoods = NULL;
  reverse_rxn_likelihoods = NULL;
  counts                  = NULL;
  drfc                    = NULL;
  dfdy_row                = NULL;
  dfdy_a                  = NULL;
  dfdy_at                 = NULL;
  dfdy_ia                 = NULL;
  dfdy_iat                = NULL;
  dfdy_ja                 = NULL;
  dfdy_jat                = NULL;
  column_mask             = NULL;
  /*
    Allocate double space needed for scratch vectors and matrices.
    actually 31*ny but we'll throw in an extra 9 ny for future needs.
  */
  if (success) {
    ask_for = ((ny * 40) + (2*ny*ny) + (3*nrxns)) * sizeof(double);
    dfdy = (double*)calloc(ask_for,one_l);
    ode23tb_params->dfdy = dfdy;
    if (dfdy == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"ode23tb: Error could not allocate %ld bytes "
		"for scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    miter     = &dfdy[nysq];
    y         = &miter[nysq];
    z         = &y[ny];
    y2        = &z[ny];
    z2        = &y2[ny];
    znew      = &z2[ny];
    ynew      = &znew[ny];
    yp        = &ynew[ny];
    f0        = &yp[ny];
    it_solve_rhs = &f0[ny];
    it_solve_del = &it_solve_rhs[ny];
    it_solve_scr = &it_solve_del[ny];
    est         = &it_solve_scr[ny];
    fac         = &est[ny];
    ode23tb_params->fac = fac;
    thresh      = &fac[ny];
    ode23tb_params->thresh = thresh;
    delfdelt    = &thresh[ny];
    wt          = &delfdelt[ny];
    pivot       = &wt[ny];
    ipivot      = (int*)pivot;
    net_lklhd_bndry_flux = &pivot[ny];
    fdel                  = &net_lklhd_bndry_flux[ny];
    ode23tb_params->fdel  = fdel;
    fdiff                 = &fdel[ny];
    ode23tb_params->fdiff = fdiff;
    dfdy_tmp              = &fdiff[ny];
    ode23tb_params->dfdy_tmp = dfdy_tmp;
    forward_rxn_likelihoods = state->ode_forward_lklhds;
    reverse_rxn_likelihoods = state->ode_reverse_lklhds;
    counts                  = state->ode_counts;
    /*
    net_likelihoods         = &dfdy_tmp[ny];
    */
    /*
      We need to move from stochastic counts to contiuous counts
      so as not to have zero concentrations.
    */
    /*
      Compute initial concentratiosn from molecule counts.
      Also fill f0 from fluxes input.
    */
    threshold = 1.0e-3;
    njthreshold = state->nj_thresh; /* 1.0e-6; */
    ode23tb_params->threshold = threshold;
    ode23tb_params->njthreshold = njthreshold;
    dcopy_(&ny,concs,&inc1,y,&inc1);
    vec_set_constant(ny,thresh,njthreshold);
    t= t0;
    if (print_concs) {
      ode_print_concs(state,t,y);
    }
  }
  if (success) {
    if (ode_jacobian_choice != 0) {
      /* 
	 We need to allocate vectors pointed to by cvodes_params,
	 for use in approximate_delta_concs.
      */
      drfc_len = state->number_molecules*2;
      /*
	Set the cvodes_params nnz field.
      */
      boltzmann_size_jacobian(state);
      nnz = cvodes_params->nnz;
      num_doubles = drfc_len + nnz + nnz + ny;
      /*
	Num ints is only 2*nnz + 3*ny + 2, but we have 
	five entities and we want each of them to have an
	even number of elements, hence the + 7.
      */
      num_ints    = 2*nnz + 3*ny + 7;
      num_ints    = num_ints + (num_ints & 1);
      ask_for     = (num_doubles << 3) + (num_ints << 2);

      drfc = (double *)calloc(one_l,ask_for);
      if (drfc == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"ode23tb: Error could not allocate %ld bytes "
		  "for double scratch space.\n",ask_for);
	  fflush(lfp);
	}
      } else {
	dfdy_row = (double *)&drfc[drfc_len];
	dfdy_a   = (double *)&dfdy_row[ny];
	dfdy_at  = (double *)&dfdy_a[nnz];
	dfdy_ia  = (int *)&dfdy_at[nnz];
	if (ny & 1) {
	  ia_len = ny + 1;
	} else {
	  ia_len = ny + 2;
	}
	dfdy_iat = (int*)&dfdy_ia[ia_len];
	nnz      = nnz + (nnz &1);
	dfdy_ja  = (int*)&dfdy_iat[ia_len];
	dfdy_jat = (int*)&dfdy_ja[nnz];
	column_mask = (int*)&dfdy_jat[nnz];
	cvodes_params->drfc     = drfc;
	cvodes_params->prec_row = dfdy_row;
	cvodes_params->dfdy_a   = dfdy_a;
	cvodes_params->dfdy_at  = dfdy_at;
	cvodes_params->dfdy_ia  = dfdy_ia;
	cvodes_params->dfdy_iat = dfdy_iat;
	cvodes_params->dfdy_ja  = dfdy_ja;
	cvodes_params->dfdy_jat = dfdy_jat;
	cvodes_params->column_mask = column_mask;
      }
    }
  } 
  if (success) {
    approximate_delta_concs(state,y,f0,delta_concs_choice);
    if (ode_rxn_view_freq > 0) {
      ode_print_dconcs(state,t,f0);
    }
    nfevals = (int64_t)1;
    /*
      Fill the Jacobian. jac is a ny x ny array.
    */
    t0 = 0;
    nf = 0;
    ode23tb_params->num_jac_first_time = 1;
    ode23tb_params->nf                 = nf;
    first_time = 1;
    ode_num_jac(state,first_time,
		dfdy,t0,y,f0,fac,thresh,
		fdel,
		fdiff,
		dfdy_tmp,
		&nf);
    first_time = 0;
    ode23tb_params->nf                 = nf;
    ode23tb_params->num_jac_first_time = 0;
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
    znew_coeff[0] = -p33;
    znew_coeff[1] =  p31;
    znew_coeff[2] =  p33;
    znew_coeff[3] =  p32;
    ynew_coeff[0] =  1.0;
    ynew_coeff[1] =  gg;
    ynew_coeff[2] =  0.0;
    ynew_coeff[3] =  gg;
    ynew_coeff[4] =  d;
    est_coeff[0] =  c1;
    est_coeff[1] =  0.0;
    est_coeff[2] =  c2;
    est_coeff[3] =  c3;

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
      ode23tb_init_wt(normcontrol, ny,normy,threshold, 
		      y, wt);
      /*
        these would be used to set h for tdel, but dfdy does not
	depend on t, hence tdel is not needed.

      ode23tb_init_h(normcontrol, ny, normy, normyp, threshold, recip_70p,
		     recip_cube_root_rtol, htspan, hmin, hmax, tdir, 
		     yp, wt,  &h, &absh)
       		     
      */
      /*
	Estimate error of first order Taylor seris, 0.5*h^2*y''(t)
        and use rule of thumb to select step size for second order method.
      */
      /*
      abst = fabs(t);
      abst_p_h = fabs(t+h);
      upper_bound_step = sqrt_eps * max(abst,abst_p_h);
      upper_bound_step = abst_p_h;
      if (abst > upper_bound_step) {
	upper_bound_step = abst;
      }
      upper_bound_step = upper_bound_step * sqrt_eps;
      */

      /*
      tdel = (t + (tdir * min(upper_bound_step,absh))) - t;
      tdel = upper_bound_step;
      if (absh < tdel) {
	tdel = absh;
      }
      tdel = tdir * tdel;
      tdel = (t + tdel) - t;
      */
      /*
	f1   = flux at t+tdel, y
	conc = count/volume
      approximate_delta_concs(state,y,f1,delta_concs_choice);
      if (ode_rxn_view_freq > 0) {
	t1 = t + tdel;
	ode_print_dconcs(state,t1,f1);
      }
      nfevals += 1;
      recip_tdel = 1.0/tdel;
      */
	/*
      for (i=0;i<ny;i++) {
	dfdt[i] = (f1[i] - f0[i]) * recip_tdel;
	delfdelt[i] = dfdt[i];
	delfdelt[i] = (f1[i] - f0[i]) * recip_tdel;
      }
	*/
      /*
	Because f1 = f0 as approxiomate_delta_concs has nodependency on t,
	and y has not changed value dfdt is 0.
      */
      vec_set_constant(ny,delfdelt,dzero);
      /*
	delfdelt = dfdt + dfdy * yp
      */
      scalar = 1.0;
      dgemv_(trans,&ny,&ny,&scalar,dfdy,&ny,yp,&inc1,&dzero,delfdelt,&inc1);

      norm_delfdelt = dnrm2_(&ny,delfdelt,&inc1);
  
      norm_delfdelt_o_wt = ode23tb_normyp_o_wt(normcontrol,ny,norm_delfdelt,
					       delfdelt,wt);
      rh = (recip_70p * sqrt(.5*norm_delfdelt_o_wt)) * recip_cube_root_rtol;
      
      /*
      absh = min(hmax,htspan);
      */
      absh = htspan;
      ode23tb_limit_h(hmax,hmin,rh,tdir,&absh,&h);
    } else {
      /*
      htry is not 0.0 but specified.
      absh = min(hmax,max(hmin,htry));
      */
      rh = 1.0/hmax;
      absh = htry;
      ode23tb_limit_h(hmax,hmin,rh,tdir,&absh,&h);
    }
    dcopy_(&ny,yp,&inc1,z,&inc1);
    dscal_(&ny,&h,z,&inc1);
    /*
    for (i=0;i<ny;i++) {
      z[i] = h*yp[i];
    }
    */
    /*
      Set up to print net likelihoods and net-likelihood-boundary fluxes
    */
    if (ode_rxn_view_freq > 0) {
      /* 
	If tracking net likelihoods and fluxes, set up to
	print the first iteration, 
      */
      ode_rxn_view_step = one_l;
    }
    /*
      Main loop.
    */
    need_new_lu = 1;
    not_done = 1;
    done     = 0;
    info = 0;
    rh = 1.0/hmax;
    while (not_done) {
      hmin = 16.0*eps*fabs(t);
      abslasth = absh;
      /*
      absh = min(hmax, max(hmin,absh));
      */
      ode23tb_limit_h(hmax,hmin,rh,tdir,&absh,&h);
      /*
	Check for last step.
      */
      if ((1.1*absh) >= fabs(tfinal - t)) {
	h = tfinal - t;
	absh = fabs(h);
	not_done = 0;
	done     = 1;
      }
      /*
	Rescale z, if timestep changed.
      */
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
#ifdef DBG
	fprintf(ode_dconcs_fp,"top of inner loop\n");
	fflush(ode_dconcs_fp);
#endif
	ode23tb_init_wt(normcontrol, ny,normy,threshold, 
			y, wt);
	if (need_new_j) {
	  /*
	    Recompute the jacobian at current y, and f (yp)
	    compute f0 and dfdy. at current y. f0 is set here.
	    This might becom ode23tb_new_dfdy
	  */
	  approximate_delta_concs(state,y,f0,delta_concs_choice);
	  if (ode_rxn_view_freq > 0) {
#ifdef DBG
	    fprintf(ode_dconcs_fp,"Inner_loop,New j needed, f0 recomputed\n");
	    fflush(ode_dconcs_fp);
#endif
	    ode_print_dconcs(state,t,f0);
	    
	  }
	  /*
	    Compute new Jacobian, dfdy.
	  ode_num_jac(state,first_time,dfdy,t,y,f0,fac,thresh,
		      ode_num_jac_scratch, &nf); dfdy is set here.
	  */
	  ode_num_jac(state,first_time,dfdy,t,y,f0,fac,thresh,
		      fdel,
		      fdiff,
		      dfdy_tmp,
		      &nf);
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
	    This could be its own subroutine say ode23tb_build_factor_miter:
	  */
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Inner loop needs new_lu\n");
#endif
	  success = ode23tb_build_factor_miter(ny,nysq,d,h,dfdy,
					       miter,ipivot,&info,lfp);
	  if (success == 0) {
	    break;
	  }
	  ndecomps = ndecomps + 1;
	  rate = 0.0;
	  need_new_lu = 0;
	} /* end if need_new_lu */
	/*
	  Stage 1.
	*/
#ifdef DBG
	fprintf(ode_dconcs_fp,"Inner loop Stage 1\n");
#endif

	
	t2 = t + alpha*h;

	/*
	  The following copy copies y in to y2 and z into z2,
	  takes advantage of the fact that the vectors 
	  y,z,y2,z2 are stored consecutively in that order (see
	  above partitioning of workspace.
	*/
	dcopy_(&two_ny,y,&inc1,y2,&inc1);
	daxpy_(&ny,&alpha,z,&inc1,y2,&inc1);
	/*
	for (i=0;i<ny;i++) {
	  y2[i] = y[i] + (alpha * z[i]);
	  z2[i] = z[i];
	}
	*/
	/*
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," Stage1 y2, y2_counts, and z2 before ode_it_solve\n");
	  origin = 3; 
	  print_concs_dconcs(state,ny,z2,y2,t,h,nsteps,origin);
	}
#endif
	*/
#ifdef DBG
	fprintf(ode_dconcs_fp,"Calling ode_it_solve in stage 1 to modify y2,z2\n");
	fflush(ode_dconcs_fp);
#endif

	iter_count = 0;
	itfail1 = ode_it_solve(state,miter,ipivot,t2,y2,z2,
			       it_solve_del,it_solve_rhs,it_solve_scr,
			       d,h,rtol,wt,&rate, &iter_count);
	/*
#ifdef DBG
	if (lfp) {
	  fprintf(lfp," Stage1 y2, y2_counts, and z2 after ode_it_solve\n");
	  origin = 4; 
	  print_concs_dconcs(state,ny,z2,y2,t,h,nsteps,origin);
	}
#endif
	*/
	nfevals += iter_count;
	nsolves += iter_count;
	itfail2 = 0;
	if (itfail1 == 0) {
	  /* Stage 2. */
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Inner loop Stage 2\n");
#endif
	  normy2 = dnrm2_(&ny,y2,&inc1);
	  ode23tb_update_wt(normcontrol, ny,normy2,y2, wt);

	  tnew = t + h;
	  if (done) {
	    tnew = tfinal;
	  }
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Stage 2, Computing new ynew,znew\n");
	  fflush(ode_dconcs_fp);
#endif
	  /*
	    The following loop may be done via gemv operations
	    if the vectors y,z,y2,z2,ynew,znew are stored consecutively.
	    ynew and znew cannotbe computed simultaneously as
	    ynew depends on znew.

	    znew = -p33 * y + p31 * z + p33 * y2 + p32 * z2 = 

                              -p33
                               p31
	       (y,z,y2,z2) * ( p32)
	                       p33

                                        
                                       1.0 
                                        gg 
            ynew = (y,z,y2,z2,znew) * ( 0   )
                                        gg     
				        d 	      
	  */                                  
	  scalar = 1.0;
	  //vec_set_constant(ny,znew,dzero);
	  dgemv_(trans,&ny,&ifour,&scalar,y,&ny,znew_coeff,&inc1,
		 &dzero,znew,&inc1);
	  //vec_set_constant(ny,ynew,dzero);
	  dgemv_(trans,&ny,&ifive,&scalar,y,&ny,ynew_coeff,&inc1,
		 &dzero,ynew,&inc1);
	  
	  iter_count = 0;
	  /*
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," Stage2 ynew, ynew_counts, and znew "
		    "before ode_it_solve\n");
	    origin = 5; 
	    print_concs_dconcs(state,ny,znew,ynew,t,h,nsteps,origin);
	  }
#endif
	  */
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Stage 2, calling ode_it_solve to modify ynew,znew\n");
	  fflush(ode_dconcs_fp);
#endif
	  itfail2 = ode_it_solve(state,miter,ipivot,tnew,ynew,znew,
				 it_solve_del,it_solve_rhs,it_solve_scr,
				 d,h,rtol,wt,&rate,&iter_count);
	  nfevals += iter_count;
	  nsolves += iter_count;
	  /*
#ifdef DBG
	  if (lfp) {
	    fprintf(lfp," Stage2 ynew, ynew_counts, and znew "
		    "after ode_it_solve\n");
	    origin = 6; 
	    print_concs_dconcs(state,ny,znew,ynew,t,h,nsteps,origin);
	  }
#endif
	  */
	} /* end if (itfail1 == 0) */
	if (itfail1 || itfail2) {
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Ode_it_solve failure, itfail1 == %d, itfail2 == %d\n",itfail1, itfail2);
	  fflush(ode_dconcs_fp);
#endif
	  nofailed = 0;
	  nfailed += 1;
	  if (jcurrent & mcurrent) {
	    if (absh <= hmin) {
	      if (lfp) {
		fprintf(lfp,"ode23tb: Error integration tolerance not met, t= %le, hmin = %le\n",t,hmin);
		fflush(lfp);
	      }
	      tolerance_met = 0;
	      /*
		Set flags to break out of inner and main loop.
	      */
	      not_done          = 0;
	      done              = 1;
	      unsuccessful_step = 0;
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
#ifdef DBG
	  fprintf(ode_dconcs_fp,"Inner loop stage1 and stage2 converged\n");
	  fflush(ode_dconcs_fp);
#endif
	  if (normcontrol) {
	    normynew = dnrm2_(&ny,ynew,&inc1);
	  }
	  ode23tb_update_wt(normcontrol,ny,normynew,ynew,wt);
	  /*
	    est = c1 *z + c2 * z2 + c3*znew

                                        c1
                                         0
	          (z,y2,z2,znew) *     ( c2 )
		                        c3
	  */
	  //vec_set_constant(ny,est,dzero);
	  dgemv_(trans,&ny,&ifour,&scalar,z,&ny,est_coeff,&inc1,
		 &dzero,est,&inc1);

	  err1 = ode23tb_max_abs_ratio(ny,est,wt);
	  /*
	    Here we need to solve miter * est2 = est1
	    With calls to something like dgetrs_
	  */
	  nrhs = 1;
	  /*
	  dcopy_(&ny,est1,&inc1,est2,&inc1);
	  */
	  dgetrs_(trans,&ny,&nrhs,miter,&ny,ipivot,est,&ny,&info);
	  if (info != 0) {
	    if (lfp) {
	      fprintf(lfp,"ode23tb: dgetrs_ call to compute est, failed with info = %d\n",
		      info);
	      fflush(lfp);
	    }
	    break;
	  }
	  nsolves += 1;
	  err2 = ode23tb_max_abs_ratio(ny,est,wt);
	  err = err1/16.0;
	  if (err2 > err) {
	    err = err2;
	  }
	  /*
	    Evaluate the error from nonnegativity violations
	  */
	  if (nonnegative) {
	    if (normcontrol) {
	      errnn_scale = 1.0/wt[0];
	    } else {
	      errnn_scale  = 1.0/threshold;
	    }
	    ode23tb_nonneg_err(ny,ynew,rtol,errnn_scale,
			       &err,&nnrejectstep);
	  }
	  if (err > rtol) {
	    /*
	      Step failed.
	    */
#ifdef DBG
	    fprintf(ode_dconcs_fp,"Inner loop step failed err = %ld\n",err);
	    fflush(ode_dconcs_fp);
#endif
	    nfailed = nfailed + 1;
	    if (absh <= hmin) {
	      if (lfp) {
		fprintf(lfp,"ode23tb: Error Integration tol not met, t= %le, hmin = %le\n",t,hmin);
		fflush(lfp);
	      }
	      tolerance_met  =   0;
	      /*
		set flags to break out of loop.
	      */
	      not_done       =   0;
	      done           =   1;
	      unsuccessful_step = 0;
	    } else {
	      /* reduce step size */
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
	      dscal_(&ny,&h_ratio,z,&inc1);
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
		
#ifdef DBG
      fprintf(ode_dconcs_fp," After inner_loop, info = %d\n",info);
      fflush(ode_dconcs_fp);
#endif
      if (not_done) {
	nsteps = nsteps + 1;
	nnreset_znew = 0;
	if (nonnegative) {
	  ode23tb_enforce_nonneg(ny,normcontrol,ynew,znew,&normynew);
	}
	/* 
	  Advance the integration one step. 
	*/
#ifdef DBG
	fprintf(ode_dconcs_fp,"Advancing integration copy ynew,znew to y,z\n");
	fflush(ode_dconcs_fp);
#endif
	t = tnew;
	dcopy_(&ny,ynew,&inc1,y,&inc1);
	dcopy_(&ny,znew,&inc1,z,&inc1);
	if (ode_rxn_view_freq > 0) {
	  ode_rxn_view_step -= one_l;
	  if (ode_rxn_view_step == zero_l) {
	    boltzmann_monitor_ode(state,t,y);
	    /*
	    ode_print_concs(state,t,y);
	    get_counts(ny,y,conc_to_count,counts);

	    ierr = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
					  reverse_rxn_likelihoods);

	    ode_print_lklhds(state,t,forward_rxn_likelihoods,
			     reverse_rxn_likelihoods);
	    approximate_delta_concs(state,y,f0,delta_concs_choice);
	    ode_print_dconcs(state,t,f0);
	    */
	    ode_rxn_view_step = ode_rxn_view_freq;
	  }
	}
	if (normcontrol) {
	  normy = normynew;
	} 
	jcurrent = 0;
	mcurrent = 1;
	if (nofailed) {
	  /*
	    Modify step size if inner steps succeeded.
	  */
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
	} /* end if nofailed */
      } /* end if (not_done) */
#ifdef DBG
      fprintf(ode_dconcs_fp,"Bottom of main_loop\n");
      fflush(ode_dconcs_fp);
#endif      
    } /* end while (not_done) MAIN loop */
  } /* end if success - allocation succeeded */

  dcopy_(&ny,y,&inc1,concs,&inc1);
  if (dfdy != NULL) {
    free(dfdy);
  }
  return (success);
}

