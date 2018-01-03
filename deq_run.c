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

#include "flatten_state.h"
#include "alloc4.h"
#include "form_molecules_matrix.h"
#include "alloc7.h"
#include "update_rxn_log_likelihoods.h"
#include "fill_flux_pieces.h"
#include "ode23tb.h"

#include "deq_run.h"
int deq_run(struct state_struct *state) {
  /*
    Run the boltzmann simulations, to be called after
    boltzmann_init has been called.

    Called by: deq
    Calls:     flatten_state,
               alloc4,
	       form_molelcules_matrix,
	       alloc7,
	       update_rxn_log_likelihoods
	       fill_flux_pieces,
	       ode23tb
  */ 
  struct state_struct *nstate;
  struct molecules_matrix_struct *molecules_matrix;
  struct molecule_struct *molecules;
  struct molecule_struct *molecule;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  int64_t *transpose_work;
  int64_t choice;
  double dchoice;
  double uni_multiplier;
  double vall;
  double rvall;
  double scaling;
  double dg_forward;
  double sum_likelihood;
  double r_sum_likelihood;
  double entropy;
  double scaled_likelihood;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *current_counts;
  double *counts;
  double *future_counts;
  double *bndry_flux_counts;
  double *activities;
  double *no_op_likelihood;
  double *flux_vector;
  double *flux_jacobian;
  double *reactant_term;
  double *product_term;
  double *p_over_r;
  double *r_over_p;
  double *concs;
  int    *rxn_fire;
  char   *cmpt_string;
  double *dg0s;
  double *free_energy;
  double htry;
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

  int64_t fe_view_step;
  int64_t fe_view_freq;

  int success;
  int number_reactions;

  int number_reactions_t2;
  int number_reactions_t2_p1;

  int rxn_choice;
  int unique_molecules;

  int print_output;
  int noop_rxn;

  int rxn_no;
  int j;

  int base_reaction;
  int fill_jacobi;

  int nonnegative;
  int cindex;


  FILE *lfp;
  success = 1;
  one_l   = (int64_t)1;
  zero_l  = (int64_t)0;
  nstate = state;
  state->workspace_base = NULL;
  success = flatten_state(state,&nstate);
  n_warmup_steps    	 = state->warmup_steps;
  n_record_steps    	 = state->record_steps;
  number_reactions       = (int)state->number_reactions;
  unique_molecules     	 = (int)state->nunique_molecules;
  current_counts         = state->current_counts;
  future_counts          = state->future_counts;
  bndry_flux_counts  	 = state->bndry_flux_counts;
  activities        	 = state->activities;
  rxn_fire          	 = state->rxn_fire;
  no_op_likelihood  	 = state->no_op_likelihood;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  print_output           = (int)state->print_output;
  number_reactions_t2    = number_reactions << 1;
  number_reactions_t2_p1 = number_reactions_t2 + 1;
  rxn_view_freq        	 = state->rxn_view_freq;
  rxn_view_hist_length 	 = state->rxn_view_hist_length;
  lklhd_view_freq        = state->lklhd_view_freq;
  count_view_freq        = state->count_view_freq;
  fe_view_freq           = state->fe_view_freq;
  molecules              = state->sorted_molecules;
  compartments           = state->sorted_cmpts;
  rxn_view_pos         	 = zero_l;
  choice_view_freq       = lklhd_view_freq;
  rxn_view_step        	 = one_l;
  count_view_step        = one_l;
  lklhd_view_step 	 = one_l;
  choice_view_step       = one_l;
  lfp                    = state->lfp;
  noop_rxn               = number_reactions + number_reactions;
  /*
    Allocate space for forming the molecules matrix.
  */
  if (success) {
    success = alloc4(state,&molecules_matrix,&transpose_work);
  }
  /*
    Compute the molecules matrix.
  */
  if (success) {
    success = form_molecules_matrix(state,molecules_matrix,transpose_work);
  }
  /*
    We need to allocate space for a flux vector ( lenth = num_species)
    The Jacobian of the flux vector (length = num_species * num_species);
    the reactant_term vector (length = num_rxns), 
    product_term_vector (length = num_rxns),  p_over_r, r_over_p both of 
    length num_rxns.
  */
  if (success) {
    success = alloc7(state);
  }
  /*
    Initialize the free_energy to be the delta_g0.
  */
  if (success) {
    flux_vector = state->flux_vector;
    flux_jacobian = state->flux_jacobian;
    reactant_term = state->reactant_term;
    product_term  = state->product_term;
    p_over_r       = state->p_over_r;
    r_over_p       = state->r_over_p;
    concs          = state->concs;
    dg0s = state->dg0s;
    free_energy  = state->free_energy;
    for (i=0;i<state->number_reactions;i++) {
      free_energy[i] = dg0s[i];
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
      concs[i] = compartment->recip_volume * counts[i];
    }
  }
  /*
    Compute the reaction likelihoods: forward_rxn_likelihood, 
    and reverse_rxn_likelihood fields of state..
  success = update_rxn_log_likelihoods(state);
  */
  /*
  base_reaction = 0;
  fill_jacobi = 1;
  if (success) {
    success = fill_flux_pieces(state,molecules_matrix,base_reaction,
			       fill_jacobi);
  }
  */

  htry = 0.1;
  nonnegative = 1.0;
  if (success) {
    success = ode23tb(state,counts,htry,nonnegative);
  }
  j = 1;
  print_counts(state,j);
  return(success);
}
