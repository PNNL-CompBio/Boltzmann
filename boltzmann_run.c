/* boltzmann_run.c
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
#define DBG_BOLTZMANN_RUN 1
*/

#include "flatten_state.h"
#include "update_rxn_log_likelihoods.h"
#include "choose_rxn.h"
#include "compute_delta_g_forward_entropy_free_energy.h"
#include "print_concentrations.h"
#include "print_likelihoods.h"
#include "save_likelihoods.h"
#include "print_free_energy.h"
#include "print_boundary_flux.h"
#include "print_restart_file.h"
#include "print_reactions_view.h"

#include "boltzmann_run.h"
int boltzmann_run(struct state_struct *state) {
  /*
    Run the boltzmann simulations, to be called after
    boltzmann_init has been called.

    Called by: boltzmann/client
    Calls:     update_rxn_log_likelihoods,
	       choose_rxn,
	       print_restart_file
	       print_reactions_view
  */
  struct state_struct *nstate;
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
  double *current_concentrations;
  double *future_concentrations;
  double *bndry_flux_concs;
  double *activities;
  double *no_op_likelihood;
  int    *rxn_fire;
  char   *cmpt_string;
  double *dg0s;
  double *free_energy;

  int success;
  int number_reactions;

  int number_reactions_t2;
  int number_reactions_t2_p1;

  int i;
  int j;

  int n_warmup_steps;
  int n_record_steps;

  int rxn_choice;
  int unique_molecules;

  int rxn_view_step;
  int rxn_view_pos;

  int rxn_view_freq;
  int rxn_view_hist_length;

  int lklhd_view_step;
  int lklhd_view_freq;

  int conc_view_step;
  int conc_view_freq;

  int print_output;
  int noop_rxn;

  FILE *lfp;
  success = 1;
  nstate = state;
  state->workspace_base = NULL;
  success = flatten_state(state,&nstate);
  n_warmup_steps    	 = (int)state->warmup_steps;
  n_record_steps    	 = (int)state->record_steps;
  number_reactions       = (int)state->number_reactions;
  unique_molecules     	 = (int)state->nunique_molecules;
  current_concentrations = state->current_concentrations;
  future_concentrations  = state->future_concentrations;
  bndry_flux_concs  	 = state->bndry_flux_concs;
  activities        	 = state->activities;
  rxn_fire          	 = state->rxn_fire;
  no_op_likelihood  	 = state->no_op_likelihood;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  print_output           = (int)state->print_output;
  number_reactions_t2    = number_reactions << 1;
  number_reactions_t2_p1 = number_reactions_t2 + 1;
  rxn_view_freq        	 = (int)state->rxn_view_freq;
  rxn_view_hist_length 	 = (int)state->rxn_view_hist_length;
  lklhd_view_freq        = (int)state->lklhd_view_freq;
  conc_view_freq         = (int)state->conc_view_freq;
  rxn_view_pos         	 = 0;
  rxn_view_step        	 = 1;
  conc_view_step         = 1;
  lklhd_view_step 	 = 1;
  lfp                    = state->lfp;
  noop_rxn               = number_reactions + number_reactions;
  /*
    Initialize the free_energy to be the delta_g0.
  */
  if (success) {
    dg0s = state->dg0s;
    free_energy  = state->free_energy;
    for (i=0;i<state->number_reactions;i++) {
      free_energy[i] = dg0s[i];
    }
  }
  if (print_output && lfp) {
    fprintf(lfp,
	    "\nWarmup_step rxn_choice forward_likelihood "
	    "reverse_likelihood\n");
  }
  for (i=0;i<n_warmup_steps;i++) {
    /*
      Compute the reaction likelihoods: forward_rxn_likelihood, 
      and reverse_rxn_likelihood fields of state..
    */
    success = update_rxn_log_likelihoods(state);
    /*
      Choose a reaction by computing the partial sums of the reaction 
      likelihoods and then using a uniform random number generator to pick one
      with probability proportional to the relative size of the reaction
      likelihood ratio. A second step called the metropolis method is 
      employed to allow reactions that use the last reactant molecules or
      produce the first product molecules to fire.
      This call call also updates the future_concentrations vector.
    */
    rxn_choice = choose_rxn(state,&r_sum_likelihood);
    if (rxn_choice < 0) break;
    if (print_output && lfp) {
      if (rxn_choice == noop_rxn) {
	fprintf(lfp,"%d\tnone\n",i);
      } else {
	fprintf(lfp,"%d\t%d\t%le\t%le\n",i,rxn_choice,forward_rxn_likelihood[rxn_choice],
		reverse_rxn_likelihood[rxn_choice]);
	fflush(lfp);
      }
    }
    /*
      Copy the future concentrations, resulting from the reaction firing
      to the current concentrations.
    */
    for (j=0;j<unique_molecules;j++) {
      current_concentrations[j] = future_concentrations[j];
    }
    /*
      Doug thinks we can remove these calls.
    success = update_rxn_log_likelihoods(state);
    success = compute_delta_g_forward_entropy_free_energy(state,
							  &dg_forward,
							  &entropy);
    */
  } /* end for(i...) */
  if (success) {
    success = update_rxn_log_likelihoods(state);
  }
  /* 
    Data collection phase (recording).
  */
  if (success) {
    /*
      Set the rxn_fire counts to 0.
    */
    if (print_output) {
      if (lfp) {
	fprintf(lfp,
	    "\nRecord_step rxn_choice forward_likelihood "
	    "reverse_likelihood\n");
      }
      for (i=0;i<number_reactions_t2;i++) {
	rxn_fire[i] = 0;
      }
    }
    /*
      Initialize the boundary fluxes to 0.
    */
    for (i=0;i<unique_molecules;i++) {
      bndry_flux_concs[i] = 0.0;
    }
    for (i=0;i<n_record_steps;i++) {
      /*
	Choose a reaction setting the future_concentrations field of
	the state structure and counting the number of times this
	reaction in this direction has fired.
      */
      rxn_choice = choose_rxn(state,&r_sum_likelihood);
      if (rxn_choice < 0) break;
      if (print_output) {
	if (lfp) {
	  if (rxn_choice == noop_rxn) {
	    fprintf(lfp,"%d\tnone\n",i);
	  } else {
	    fprintf(lfp,"%d\t%d\t%le\t%le\n",i,rxn_choice,forward_rxn_likelihood[rxn_choice],
		    reverse_rxn_likelihood[rxn_choice]);
	    fflush(lfp);
	  }
	}
	if (rxn_choice <= number_reactions_t2) {
	  rxn_fire[rxn_choice] += 1;
	}
      }
      /*
	Copy the future concentrations, the result of the reaction firing
	to the current concentrations.
      */
      for (j=0;j<unique_molecules;j++) {
	current_concentrations[j] = future_concentrations[j];
      }
      /*
	Compute the reaction likelihoods and their logarithms
	in the forward_rxn_likelihood, reverse_rxn_likelihood 
	forward_rxn_log_likelihood_ratio and 
	reverse_rxn_log_likelihood_ratio fields
	based on the current_concentrations field of state.
      */
      success = update_rxn_log_likelihoods(state);
      /*
	Compute the dg_forward and entropy values and free_energy field of
	the state structure at the current concentations, and 
	rxn_likelihoods.
      */
      success = compute_delta_g_forward_entropy_free_energy(state,
							    &dg_forward,
							    &entropy,
							    i);
      if (print_output) {
	/* 
	  print the concentrations. 
	*/
	conc_view_step = conc_view_step - 1;
	if ((conc_view_step <= 0) || (i == (n_record_steps-1))) {
	  print_concentrations(state,i);
	  conc_view_step = conc_view_freq;
	}
	/* 
	  print the entropy, dg_forward and the reaction likelihoods, 
	*/
	lklhd_view_step = lklhd_view_step - 1;
	if ((lklhd_view_step <= 0) || (i == (n_record_steps-1))) {
	  print_likelihoods(state,entropy,dg_forward,i) ;
	  lklhd_view_step = lklhd_view_freq;
	}
	if (rxn_view_freq > 0) {
	  rxn_view_step = rxn_view_step - 1;
	  /*
	    Save the likelihoods on a per reaction basis for later 
	    printing to the rxns.view file.
	  */
	  if ((rxn_view_step <= 0) || (i == (n_record_steps-1))) {
	    no_op_likelihood[rxn_view_pos] = r_sum_likelihood;
	    save_likelihoods(state,rxn_view_pos);
	    rxn_view_step = rxn_view_freq;
	    rxn_view_pos  += 1;
	  }
	}
	/*
	  If user has requested them, print out free energies as well.
	*/
	if (state->free_energy_format > (int64_t)0) {
	  print_free_energy(state,i);
	}
      } /* end if (print_output) */
    } /* end for(i...) */
    if (print_output) {
      if (state->num_fixed_concs > (int64_t)0) {
	print_boundary_flux(state);
      } /* end if (state->num_fixed_concs ...) */
      if (success) {
	success = print_restart_file(state);
      }
      if (success) {
	if (rxn_view_freq > 0) {
	  success = print_reactions_view(state);
	}
      }
    }
    state->entropy = entropy;
    state->dg_forward = dg_forward;
  }
  return(success);
}
