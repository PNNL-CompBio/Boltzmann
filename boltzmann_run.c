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
#include "update_rxn_log_likelihoods.h"
#include "choose_rxn.h"
#include "deq_run.h"
#include "compute_delta_g_forward_entropy_free_energy.h"
#include "print_counts.h"
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
    boltzmann_init has been called to allocate and initialize state structure.

    Called by: boltzmann/client
    Calls:     update_rxn_log_likelihoods,
	       choose_rxn,
               deq_run,
	       compute_delta_g_forward_entropy_free_energy
	       print_counts,
	       print_concs,
	       print_likelihoods,
	       print free_energy,
	       print boundary_flux,
	       print_restart_file
	       print_reactions_view
  */
  double dg_forward;
  double r_sum_likelihood;
  double entropy;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *current_counts;
  double *future_counts;
  double *bndry_flux_counts;
  double *activities;
  double *no_op_likelihood;
  double *dg0s;
  double *free_energy;
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

  int64_t fe_view_step;
  int64_t fe_view_freq;

  int64_t use_deq;

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

  int ierr;
  int padi;

  FILE *lfp;
  success = 1;
  one_l   = (int64_t)1;
  zero_l  = (int64_t)0;
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
  use_deq                = state->use_deq;
  rxn_view_pos         	 = zero_l;
  choice_view_freq       = lklhd_view_freq;
  rxn_view_step        	 = one_l;
  count_view_step        = one_l;
  lklhd_view_step 	 = one_l;
  choice_view_step       = one_l;
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
  if ((print_output > 1) && lfp) {
    fprintf(lfp,
	    "\nWarmup_step rxn_choice forward_likelihood "
	    "reverse_likelihood\n");
  }
  if (use_deq == zero_l) {
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
	This call also updates the future_counts vector.
      */
      rxn_choice = choose_rxn(state,&r_sum_likelihood);
      if (rxn_choice < 0) break;
      if ((print_output > 1) && lfp) {
	if (choice_view_freq > zero_l) {
	  choice_view_step = choice_view_step - one_l;
	  if ((choice_view_step <= zero_l) || (i == (n_warmup_steps-one_l))) {
	    if (rxn_choice == noop_rxn) {
	      fprintf(lfp,"%lld\tnone\n",i);
	    } else {
	      if (rxn_choice < number_reactions) {
		fprintf(lfp,"%lld\t%d\t%le\t%le\n",i,rxn_choice,forward_rxn_likelihood[rxn_choice],
			reverse_rxn_likelihood[rxn_choice]);
	      } else {
		rxn_no = rxn_choice - number_reactions;
		fprintf(lfp,"%lld\t%d\t%le\t%le\n",i,rxn_choice,reverse_rxn_likelihood[rxn_no],
			forward_rxn_likelihood[rxn_no]);
	      
	      }
	    }
	    fflush(lfp);
	    choice_view_step = choice_view_freq;
	  }
	}
      }
      /*
	Copy the future counts, resulting from the reaction firing
	to the current counts.
      */
      for (j=0;j<unique_molecules;j++) {
	current_counts[j] = future_counts[j];
      }
      /*
	Doug thinks we can remove these calls.
	success = update_rxn_log_likelihoods(state);
	success = compute_delta_g_forward_entropy_free_energy(state,
	                                                     &dg_forward,
							     &entropy);
      */
    } /* end for(i...) */
    i = -2;
  } else {
    /*
      Use ode solver to move from initial concentrations to
      steady state.
    */
    state->print_ode_concs = 0;
    success = deq_run(state);
    i = -1;
  }
  print_counts(state,i);
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
    if (print_output > 1) {
      if (lfp) {
	fprintf(lfp,
	    "\nRecord_step rxn_choice forward_likelihood "
	    "reverse_likelihood\n");
      }
      for (i=0;i<number_reactions_t2;i++) {
	rxn_fire[i] = (int64_t)0;
      }
    }
    /*
      Initialize the boundary fluxes to 0.
    */
    for (i=0;i<unique_molecules;i++) {
      bndry_flux_counts[i] = 0.0;
    }
    rxn_view_step    = one_l;
    count_view_step  = one_l;
    lklhd_view_step  = one_l;
    fe_view_step     = one_l;
    choice_view_step = one_l;
    for (i=0;i<n_record_steps;i++) {
      /*
	Choose a reaction setting the future_counts field of
	the state structure and counting the number of times this
	reaction in this direction has fired.
      */
      rxn_choice = choose_rxn(state,&r_sum_likelihood);
      if (rxn_choice < 0) break;
      if (print_output > 1) {
	if (lfp) {
	  if (choice_view_freq > zero_l) {
	    choice_view_step = choice_view_step - one_l;
	    if ((choice_view_step <= zero_l) || (i == (n_record_steps-one_l))) {
	      if (rxn_choice == noop_rxn) {
		fprintf(lfp,"%lld\tnone\n",i);
	      } else {
		if (rxn_choice < number_reactions) {
		  fprintf(lfp,"%lld\t%d\t%le\t%le\n",i,rxn_choice,forward_rxn_likelihood[rxn_choice],
			  reverse_rxn_likelihood[rxn_choice]);
		} else {
		  rxn_no = rxn_choice - number_reactions;
		  fprintf(lfp,"%lld\t%d\t%le\t%le\n",i,rxn_choice,reverse_rxn_likelihood[rxn_no],
			  forward_rxn_likelihood[rxn_no]);
		  
		}
	      }
	      fflush(lfp);
	      choice_view_step = choice_view_freq;
	    }
	  }
	}
      }
      if (print_output) {
	if (rxn_choice <= number_reactions_t2) {
	  rxn_fire[rxn_choice] += (int64_t)1;
	}
      }
      /*
	Copy the future counts, the result of the reaction firing
	to the current counts.
      */
      for (j=0;j<unique_molecules;j++) {
	current_counts[j] = future_counts[j];
      }
      /*
	Compute the reaction likelihoods and their logarithms
	in the forward_rxn_likelihood, reverse_rxn_likelihood 
	forward_rxn_log_likelihood_ratio and 
	reverse_rxn_log_likelihood_ratio fields
	based on the current_counts field of state.
      */
      success = update_rxn_log_likelihoods(state);
      /*
	Compute the dg_forward and entropy values and free_energy field of
	the state structure at the current counts, and 
	rxn_likelihoods.
      */
      success = compute_delta_g_forward_entropy_free_energy(state,
							    &dg_forward,
							    &entropy,
							    i);
      if (print_output) {
	/* 
	  print the counts. 
	*/
	count_view_step = count_view_step - one_l;
	if ((count_view_step <= zero_l) || (i == (n_record_steps-1))) {
	  print_counts(state,i);
	  count_view_step = count_view_freq;
	}
	/* 
	  print the entropy, dg_forward and the reaction likelihoods, 
	*/
	if (lklhd_view_freq > zero_l) {
	  lklhd_view_step = lklhd_view_step - one_l;
	  if ((lklhd_view_step <= zero_l) || (i == (n_record_steps-one_l))) {
	    print_likelihoods(state,entropy,dg_forward,i) ;
	    lklhd_view_step = lklhd_view_freq;
	  }
	}
	if (rxn_view_freq > zero_l) {
	  rxn_view_step = rxn_view_step - one_l;
	  /*
	    Save the likelihoods on a per reaction basis for later 
	    printing to the rxns.view file.
	  */
	  if ((rxn_view_step <= zero_l) || (i == (n_record_steps-one_l))) {
	    no_op_likelihood[rxn_view_pos] = r_sum_likelihood;
	    save_likelihoods(state,rxn_view_pos);
	    rxn_view_step = rxn_view_freq;
	    rxn_view_pos  += one_l;
	  }
	}
	/*
	  If user has requested them, print out free energies as well.
	*/
	if (state->free_energy_format > zero_l) {
	  if (fe_view_freq > zero_l) {
	    fe_view_step = fe_view_step - one_l;
	    if ((fe_view_step <= zero_l) || (i == (n_record_steps-one_l))) {
	      print_free_energy(state,i);
	      fe_view_step = fe_view_freq;
	    }
	  }
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
	if (rxn_view_freq > zero_l) {
	  success = print_reactions_view(state);
	}
      }
    }
    state->entropy = entropy;
    state->dg_forward = dg_forward;
  }
  return(success);
}
