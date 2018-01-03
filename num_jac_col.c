#include "boltzmann_structs.h"

#include "approximate_delta_concs.h"
/*
#define DBG 1 
*/
#ifdef DBG
#include "print_concs_fluxes.h"
#endif

#include "num_jac_col.h"

int num_jac_col(struct state_struct *state,
		int ny, int j,
		int    *rowmax_p,
		double *y,
		double *f,
		double *delj_p,
		double threshj,
		double *fdel,
		double *fdiff,
		double *dfdy_colj,
		double *absfvaluerm_p,
		double *absfdelrm_p,
		double *absfdiffmax_p,
		double *infnormdfdy_colj_p) {
  /*
    Compute column j of the jacobian matrix dfdy,
    given a delta for elemant j of y, delj, a 
    theshhold for elemntj, the conc_to_count ratio for element j,
    the count vector, the flux vector  f and concetration vector y,
    two scractch

    Called by: ode_num_jac
    Calls:     compute_flux_scaling, approximate_delta_concs, fabs

    Arguments          TMF    Descriptin
    state              G*I    Boltzmann state for passing to approxmate_fluxes(f)
    ny                 ISI    length of f, y, fdel, fdiff, dfdy_colj
                              vectors,

    j                  ISI    integer in 0:ny-1 specifying which y element
                              f is to be differentiated against

    y                  D*I    concentrations vector.

    f                  D*I    flux vector,

    delj_p             D*I    value of del[j] from calling routine,
                              if yj + del[j] < 0 we set del[j] to be
			      -.5*yj and ydelj to be .5*yj

    threshj            DIS    value of thresh[j] from calling routine.


    fdel               D*W    scratch vector of length ny used to compute
                              flux at ydel.

    fdiff              D*W    scratch vector of length ny used to compute
                              flux at ydel minus flux at y.

    dfdy_colj          D*O    output vector of length ny corresponding to
                              the j'th column of dfdy.
   
    absfdiffmax_p      D*O    scalar that is the inf norm of fdiff,
 
    absfvaluerm_p      D*O    scalar that is the fabs(flux) at location (rowmax)
                              of absfdiffmax

    absfdelrm_p        D*O    scalar that is the fabs(fdel) at location (rowmax)
                              of absfdiffmax

    infnormdfdy_colj_p D*O    scalar that is the inf norm of dfdy_colj
  */
  struct molecule_struct *molecules;
  struct molecule_struct *moleculej;
  double recip_delj;
  double yj;
  double delj;
  double y_counts_savej;
  double ydelj;
  double conc_scalej;
  double absdfdy_coljk;
  double absfdiffmax;
  double absfvaluerm;
  double absfdelrm;
  double infnormdfdy_colj;
  double *conc_to_count;
  double zero;
#ifdef DBG
  double t0;
  double h;
  int    origin;
  int    nsteps;
#endif

  int    base_rxn;
  int    k;

  int    success;
  int    rowmax;

  int    variable;
  int    choice;

  FILE   *lfp;
  FILE   *efp;

  success = 1;
  base_rxn      = state->base_reaction;
  conc_to_count = state->conc_to_count;
  lfp           = state->lfp;
  molecules     = state->sorted_molecules;
  choice        = (int)state->delta_concs_choice;
  delj          = *delj_p;
  moleculej     = (struct molecule_struct *)&molecules[j];
  variable      = moleculej->variable;
  zero          = 0.0;
  if ((j < 0) || (j >= ny)) {
    success = 0;
  } else {
    if ((delj == 0.0) || (variable == 0)) {
      /*
	Fixed concentration .
      */
      for (k=0;k<ny;k++) {
	fdel[k]  = zero;
	fdiff[k] = zero;
	dfdy_colj[k] = zero;
      }
      *absfdiffmax_p = zero;
      *absfdelrm_p   = zero;
      *absfvaluerm_p = zero;
      *infnormdfdy_colj_p = zero;
    } else {
      yj = y[j];
      conc_scalej    = conc_to_count[j];
      ydelj          = yj + delj;
      /*
	Now here we want to be careful, to keep y non-negativ.
	we assume that yj is non-negative on input (maybe ought to
	check that and print message if not.
      */
      if (ydelj < 0.0) {
	ydelj = .5 *yj;
	delj  = -ydelj;
	*delj_p = delj;
      }
      if (ydelj < 0.0) {
	ydelj = 0.0;
      }
      /*
	evalueate flux at ydelj
      */
      y[j] = ydelj;
      approximate_delta_concs(state,y,fdel,choice);
#ifdef DBG
      if (lfp) {
	origin = 1000+j; 
	t0 = 0.0;
	h  = 0.0;
	nsteps = 0;
	print_concs_fluxes(state,ny,fdel,y,t0,h,nsteps,origin);
      }
#endif
      /*
	Restore y vector.
      */
      y[j]           = yj;
      if (delj != 0) {
	recip_delj = 1.0/delj;
	for (k=0;k<ny;k++) {
	  fdiff[k] = fdel[k] - f[k];
	  dfdy_colj[k] = fdiff[k] * recip_delj;
	}
      } else {
	for (k=0;k<ny;k++) {
	  dfdy_colj[k] = 0.0;
	}
      }
      infnormdfdy_colj = fabs(dfdy_colj[0]);
      rowmax           = 0;
      for (k=1;k<ny;k++) {
	absdfdy_coljk = fabs(dfdy_colj[k]);
	if (absdfdy_coljk > infnormdfdy_colj) {
	  infnormdfdy_colj = absdfdy_coljk;
	  rowmax = k;
	}
      }
      *rowmax_p = rowmax;
      absfdiffmax    = fabs(fdiff[rowmax]);
      absfdelrm      = fabs(fdel[rowmax]);
      absfvaluerm    = fabs(f[rowmax]);
      *absfdiffmax_p = absfdiffmax;
      *absfdelrm_p   = absfdelrm;
      *absfvaluerm_p = absfvaluerm;
      *infnormdfdy_colj_p = infnormdfdy_colj;
    } /* end else variable concentration. */
  } /* end else valid j */
  return(success);
}
