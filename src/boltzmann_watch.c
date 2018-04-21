#include "boltzmann_structs.h"
#include "print_rxn_choice.h"
#include "print_counts.h"
#include "print_likelihoods.h"
#include "save_likelihoods.h"
#include "print_free_energy.h"
#include "boltzmann_watch.h"
void boltzmann_watch(struct state_struct *state,
		     int64_t *choice_view_step_p,
		     int64_t *count_view_step_p,
		     int64_t *lklhd_view_step_p,
		     int64_t *rxn_view_step_p,
		     int64_t *fe_view_step_p,
		     int64_t *rxn_view_pos_p,
		     double  r_sum_likelihood,
		     double  dg_forward,
		     double  entropy,
		     int64_t iter,
		     int64_t rxn_choice) {
  /*
    Print per iteration information called only if print_output = 1.
    Called by: boltzmann_run
    Calls:     print_rxn_choice, print_counts, print_likelihoods, 
               save_likelihoods, print_free_energy, fprintf, fflush
  */
  double *no_op_likelihood;
  int64_t choice_view_step;
  int64_t choice_view_freq;
  int64_t count_view_step;
  int64_t count_view_freq;
  int64_t lklhd_view_step;
  int64_t lklhd_view_freq;
  int64_t rxn_view_step;
  int64_t rxn_view_freq;
  int64_t rxn_view_pos;
  int64_t fe_view_step;
  int64_t fe_view_freq;
  int64_t zero_l;
  int64_t one_l;

  zero_l = (int64_t)0;
  one_l  = (int64_t)1;

  count_view_freq  = state->count_view_freq;
  lklhd_view_freq  = state->lklhd_view_freq;
  rxn_view_freq    = state->rxn_view_freq;
  fe_view_freq     = state->fe_view_freq;  
  choice_view_freq = lklhd_view_freq;
  no_op_likelihood = state->no_op_likelihood;
  
  choice_view_step = *choice_view_step_p;
  count_view_step  = *count_view_step_p;
  lklhd_view_step  = *lklhd_view_step_p;
  rxn_view_step    = *rxn_view_step_p;
  fe_view_step     = *fe_view_step_p;
  rxn_view_pos     = *rxn_view_pos_p;
  /*
    For first iteration no reaction has been chosen yet.
  */
  if (iter > zero_l) {
    
    if (choice_view_freq > zero_l) {
      choice_view_step = choice_view_step - one_l;
      if (choice_view_step <= zero_l) {
	/* 
	   time to print reaction choice info .
	*/
	print_rxn_choice(state,iter,rxn_choice);
	choice_view_step = choice_view_freq;
      }
    } /* end if (choice_view_freq > 0) */
  } /* end if iter > 0 */
  /* 
    print the counts. 
  */
  if (count_view_freq > zero_l) {
    count_view_step = count_view_step - one_l;
    if (count_view_step <= zero_l)  {
      print_counts(state,iter);
      count_view_step = count_view_freq;
    }
  }
  /* 
     print the entropy, dg_forward and the reaction likelihoods, 
  */
  if (lklhd_view_freq > zero_l) {
    lklhd_view_step = lklhd_view_step - one_l;
    if (lklhd_view_step <= zero_l) {
      print_likelihoods(state,entropy,dg_forward,iter) ;
      lklhd_view_step = lklhd_view_freq;
    }
  }
  /*
    Save transpose of likeliihoods if rxn_view_freq > 0
  */
  if (rxn_view_freq > zero_l) {
    rxn_view_step = rxn_view_step - one_l;
    /*
      Save the likelihoods on a per reaction basis for later 
      printing to the rxns.view file.
    */
    if (rxn_view_step <= zero_l) {
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
      if (fe_view_step <= zero_l) {
	print_free_energy(state,iter);
	fe_view_step = fe_view_freq;
      }
    }
  }
  *choice_view_step_p = choice_view_step;
  *count_view_step_p  =	count_view_step;
  *lklhd_view_step_p  =	lklhd_view_step; 
  *rxn_view_step_p    =	rxn_view_step; 
  *fe_view_step_p     =	fe_view_step;   
  *rxn_view_pos_p     = rxn_view_pos;
}
