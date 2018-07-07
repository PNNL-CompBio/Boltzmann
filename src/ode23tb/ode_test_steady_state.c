#include "boltzmann_structs.h"
#include "blas.h"
#include "gradient.h"
#include "ode_test_steady_state.h"
int ode_test_steady_state(struct state_struct *state,
			  int ny,
			  double *y,
			  double *f) {
  /*
    Test whethter or not steady state has been reached by testing the
    sized of the derivative vector f agains the concentration vector, y,
    size or against a set threshold. Compute f from y with 
    gradient routine.
    Returns 1 if steady state has been reached according to the
    criteria, 0 if not.

    if ode_stop_style = 0, always returns 0 - don't stop at steady state,
      integrate till final time.
     
    if ode_stop_style = 1, then 

      if (ode_stop_rel is 1) {
         if (||f||_a < threshold * ||y||_a) {
	    done = 1;
	 }
      } else {
         if (||f||_a < threshold) 
	    done = 1;
         }
      }

      where ||.||_a specifies the infinity, 1, or 2 norm based on 
      ode_stop_norm = 0,1, or 2.
   
    Called by ode23tb, boltzmann_cvodes
    Calls: gradient, idamax_, dnrm2_, fabs
  */
  double ode_stop_thresh;
  double fnorm;
  double ynorm;
  int gradient_choice;
  int ode_stop_rel;

  int ode_stop_norm;
  int ode_stop_style;

  int done;
  int incx;

  int i;
  int padi;

  ny                 = state->nunique_molecules;
  gradient_choice    = state->gradient_choice;
  ode_stop_thresh    = state->ode_stop_thresh;
  ode_stop_rel       = state->ode_stop_rel;
  ode_stop_norm      = state->ode_stop_norm;
  ode_stop_style     = state->ode_stop_style;
  incx               = 1;

  done = 0;
  if (ode_stop_style) {
    gradient(state,y,f,gradient_choice);
    if (ode_stop_norm == 1) {
      /*
	Use the 1 norm.
      */
      fnorm = 0.0;
      for (i=0;i<ny;i++) {
	fnorm += fabs(f[i]);
      }
    } else if (ode_stop_norm == 2) {
      /*
	Use the 2 norm.
      */
      fnorm = dnrm2_(&ny,f,&incx);
    } else {
      /*
	Use the infinity norm - note idamax_ is a fortran function and
	returns a fortran index from which we must subtract 1.
      */
      fnorm = fabs(f[idamax_(&ny,f,&incx)-1]);
    }
    if (ode_stop_rel) {
      if (ode_stop_norm == 1) {
        /*
      	Use the 1 norm.
        */
        ynorm = 0.0;
        for (i=0;i<ny;i++) {
	  ynorm += fabs(y[i]);
        }
      } else if (ode_stop_norm == 2) {
        /*
      	Use the 2 norm.
        */
        ynorm = dnrm2_(&ny,y,&incx);
      } else {
        /*
      	Use the infinity norm
        */
        ynorm = fabs(y[idamax_(&ny,y,&incx)-1]);
      } 
      ode_stop_thresh = ode_stop_thresh * ynorm;
    } /* end if ode_stop_rel */
    if (fnorm < ode_stop_thresh) {
      done = 1;
    }
  } /* end if ode_stop_style */
  return(done); 
}
