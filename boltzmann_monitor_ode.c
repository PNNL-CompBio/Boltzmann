#include "boltzmann_structs.h"
#include "ode_print_concs.h"
#include "get_counts.h"
#include "update_rxn_likelihoods.h"
#include "ode_print_lklhds.h"
#include "approximate_delta_concs.h"
#include "ode_print_dconcs.h"
#include "boltzmann_monitor_ode.h"
void boltzmann_monitor_ode(struct state_struct *state,
			   double time,
			   double *concs) {
  /*
    Print out the concentrations, likelihoods, and concentrations derivative
    for and ode step.
    Called by boltzmann_cvodes, ode23tb
  */
  double *counts;
  double *forward_rxn_likelihoods;
  double *reverse_rxn_likelihoods;
  double *conc_to_count;
  double *f;
  int ny;
  int ierr;
  int delta_concs_choice;
  int padi;
  counts                  = state->ode_counts;
  forward_rxn_likelihoods = state->ode_forward_lklhds;
  reverse_rxn_likelihoods = state->ode_reverse_lklhds;
  f                       = state->ode_f;
  ny                      = state->nunique_molecules;
  conc_to_count           = state->conc_to_count;
  delta_concs_choice      = state->delta_concs_choice;
  ode_print_concs(state,time,concs);
  get_counts(ny,concs,conc_to_count,counts);
  ierr = update_rxn_likelihoods(state,counts,forward_rxn_likelihoods,
				reverse_rxn_likelihoods);

  ode_print_lklhds(state,time,forward_rxn_likelihoods,
		   reverse_rxn_likelihoods);
  approximate_delta_concs(state,concs,f,delta_concs_choice);
  ode_print_dconcs(state,time,f);
}
  
