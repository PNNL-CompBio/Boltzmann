/* deq_run.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:


Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/
#include "boltzmann_structs.h"
/*
#define DBG_DEQ_RUN 1
*/
#include "alloc7.h"
#include "init_base_reactants.h"
#include "init_relative_rates.h"
#include "ode_print_concs_header.h"
#include "ode_print_grad_header.h"
#include "ode_print_bflux_header.h"
#include "ode_print_lklhd_header.h"
#include "print_net_likelihood_header.h"
#include "print_net_lklhd_bndry_flux_header.h"
#include "get_counts.h"
/*
#include "fill_flux_pieces.h"
#include "ode23tb.h"
*/
#include "ode_solver.h"

#include "deq_run.h"
int deq_run(struct state_struct *state) {
  /*
    Run the boltzmann simulations, to be called after
    boltzmann_init has been called.
    Alters the current_counts field of state with an
    estimate of the steady state concentrations.
  

    Called by: deq, boltzmann_run
    Calls:     alloc7,
	       init_base_reactants,
	       init_relative_rates,
	       ode_solver
  */
  struct molecule_struct *molecules;
  struct molecule_struct *molecule;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *current_counts;
  double *counts;
  double *bndry_flux_counts;
  double *activities;
  double *no_op_likelihood;
  double *reactant_term;
  double *product_term;
  double *concs;
  double *conc_to_count;
  double *count_to_conc;
  double *dg0s;
  double *free_energy;
  double htry;
  double min_conc;
  int64_t *rxn_fire;
  int64_t i;
  int64_t n_warmup_steps;
  int64_t n_record_steps;
  int64_t one_l;
  int64_t zero_l;


  int64_t rxn_view_step;
  int64_t rxn_view_pos;

  int64_t rxn_view_freq;
  int64_t rxn_view_hist_length;

  int64_t lklhd_view_step;
  int64_t lklhd_view_freq;

  int64_t choice_view_step;
  int64_t choice_view_freq;

  int64_t count_view_step;
  int64_t count_view_freq;

  int success;
  int number_reactions;

  int number_reactions_t2;
  int number_reactions_t2_p1;

  int unique_molecules;
  int padi;

  int print_output;
  int noop_rxn;

  int nonnegative;
  int cindex;

  int normcontrol;
  int print_ode_concs;

  int solver_choice;
  int ode_rxn_view_freq;


  FILE *lfp;
  FILE *ode_kq_fp;
  success = 1;
  one_l   = (int64_t)1;
  zero_l  = (int64_t)0;

  n_warmup_steps    	 = state->warmup_steps;
  n_record_steps    	 = state->record_steps;
  number_reactions       = (int)state->number_reactions;
  unique_molecules     	 = (int)state->nunique_molecules;
  current_counts         = state->current_counts;
  bndry_flux_counts  	 = state->bndry_flux_counts;
  activities        	 = state->activities;
  rxn_fire          	 = state->rxn_fire;
  no_op_likelihood  	 = state->no_op_likelihood;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  print_output           = (int)state->print_output;
  conc_to_count          = state->conc_to_count;
  count_to_conc          = state->count_to_conc;
  min_conc               = state->min_conc;
  number_reactions_t2    = number_reactions << 1;
  number_reactions_t2_p1 = number_reactions_t2 + 1;
  rxn_view_freq        	 = state->rxn_view_freq;
  rxn_view_hist_length 	 = state->rxn_view_hist_length;
  lklhd_view_freq        = state->lklhd_view_freq;
  count_view_freq        = state->count_view_freq;
  molecules              = state->sorted_molecules;
  compartments           = state->sorted_compartments;
  print_ode_concs        = state->print_ode_concs;
  ode_rxn_view_freq      = state->ode_rxn_view_freq;
  solver_choice          = (int)state->ode_solver_choice;
  rxn_view_pos         	 = zero_l;
  choice_view_freq       = lklhd_view_freq;
  rxn_view_step        	 = one_l;
  count_view_step        = one_l;
  lklhd_view_step 	 = one_l;
  choice_view_step       = one_l;
  lfp                    = state->lfp;
  noop_rxn               = number_reactions + number_reactions;
  normcontrol            = 0;
  print_ode_concs        = print_ode_concs && print_output;
  state->print_ode_concs = print_ode_concs;
  /*
    We need to allocate space needed by the derivative approximation routine
    ode_counts (length = num_species)
    reactant_term vector (length = num_rxns), 
    product_term_vector (length = num_rxns),  
    ode_forward_lklhds (length = num_rxns)
    ode_reverse_lklhds (length = num_rxns)
    kf_rel (length = num_rxns),
    kr_rel (length = num_rxns)
    This is now called in run_init.
  if (success) {
    success = alloc7(state);
  }
  */
  /*
    At some point we might want to pull these initializations into
    run_init as well.
    Initialize the free_energy to be the delta_g0.
  */
  if (success) {
    reactant_term = state->reactant_term;
    product_term  = state->product_term;
    concs          = state->ode_concs;
    dg0s = state->dg0s;
    free_energy  = state->free_energy;
    for (i=0;i<state->number_reactions;i++) {
      free_energy[i] = dg0s[i];
    }
  }
  /*
    Fill the base_reactant_indicator, and base_reactants vectors,
    and set number_base_reaction_reactants.
  */
  if (success) {
    success = init_base_reactants(state);
  }
  /*
    Fill the relative forward and reverse reaction rate vectors,
    kf_rel, and kr_rel. For use in the lr6_gradient routine.
  */
  if (success) {
    if (state->gradient_choice == 6) {
      success = init_relative_rates(state);
    }
  }
  /*
    Compute concentrations from counts;
  */
  counts = current_counts;
  if (success) {
    molecule    = molecules;
    for(i=0;i<unique_molecules;i++) {
      molecule  = (struct molecule_struct *)&molecules[i];
      cindex = molecule->c_index;
      compartment = (struct compartment_struct *)&compartments[cindex];
      concs[i] = counts[i] * count_to_conc[i];
    }
  }
  /*
    Convert counts to continuous setting if 0 to avoid 0 concentrations.
  */
  /*
  for (i=0;i<unique_molecules;i++) {
    if (counts[i] <= 0.0) {
      counts[i] = min_conc * conc_to_count[i];
    }
  }
  */
  htry = 0.0;
  nonnegative = 1;
  if (success) {
    if (print_output) {
      if (ode_rxn_view_freq > 0) {
	ode_print_concs_header(state);
	ode_print_lklhd_header(state);
	ode_print_grad_header(state);
	ode_print_bflux_header(state);
	ode_kq_fp = fopen(state->ode_kq_file,"w");
	state->ode_kq_fp = ode_kq_fp;
	print_net_likelihood_header(state);
	print_net_lklhd_bndry_flux_header(state);
      }
    }
    success = ode_solver(state,concs,solver_choice);
  }
  /*
  j = 1;
  print_counts(state,j);
  */
  /*
    Convert concs to counts,
  */
  get_counts(unique_molecules,concs,conc_to_count,counts);
  if (state->no_round_from_deq == (int64_t)0) {
    /*
      Convert counts back to integers,
   */
    for (i=0;i<unique_molecules;i++) {
      counts[i] = (double)((int64_t)(counts[i] + 0.5));
    }
  }
  return(success);
}
