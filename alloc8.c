#include "boltzmann_structs.h"
#include "alloc8.h"
int alloc8(struct state_struct *state) {
  /*
    Allocate workspace vectors. This used to happen in the
    flatten_state routine.

    Allocates space for and sets
      future_counts,
      free_energy,
      forward_rxn_likelihood,
      reverse_rxn_likelihood,
      forward_rxn_log_likelihood_ratio,
      reverse_rxn_log_likelihood_ratio,
      rxn_likelihood_ps
      
    pointers in state

    Called by: run_init
    Calls:     calloc, fprintf, fflush
  */
  double  *future_counts;
  double  *free_energy;
  double  *forward_rxn_likelihood;
  double  *reverse_rxn_likelihood;
  double  *forward_rxn_log_likelihood_ratio;
  double  *reverse_rxn_log_likelihood_ratio;
  double  *rxn_likelihood_ps;
  
  int64_t usage;
  int64_t number_reactions;
  int64_t number_species;
  int64_t align_len;
  int64_t align_mask;
  int64_t data_pad;
  int64_t ask_for;
  int64_t one_l;

  int     success;
  int     padi;

  FILE    *lfp;
  FILE    *efp;
  
  align_len            = state->align_len;
  align_mask           = state->align_mask;
  number_reactions     = state->number_reactions;
  number_species       = state->nunique_molecules;
  lfp                  = state->lfp;
  success              = 1;
  usage = (int64_t)0;
  one_l = (int64_t)1;
  
  ask_for = number_species * sizeof(double);
  data_pad = (align_len - (ask_for & align_mask)) & align_mask;
  ask_for += data_pad;
  usage   += ask_for;
  future_counts = (double *)calloc(one_l,ask_for);
  if (future_counts == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"alloc8: Error unable to allocate %lld bytes for future_counts\n", ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    state->future_counts = future_counts;
    ask_for = number_reactions * sizeof(double);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    free_energy = (double *)calloc(one_l,ask_for);
    if (free_energy == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
	"alloc8: Error unable to allocate %lld bytes for free_energy\n", 
		ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->free_energy = free_energy;
    usage+= ask_for;
    forward_rxn_likelihood = (double *)calloc(one_l,ask_for);
    if (forward_rxn_likelihood == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
    "alloc8: Error unable to allocate %lld bytes for forward_rxn_likelihood\n",
		ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->forward_rxn_likelihood = forward_rxn_likelihood;
    usage += ask_for;
    reverse_rxn_likelihood = (double *)calloc(one_l,ask_for);
    if (reverse_rxn_likelihood == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
    "alloc8: Error unable to allocate %lld bytes for reverse_rxn_likelihood",
		ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->reverse_rxn_likelihood = reverse_rxn_likelihood;
    usage += ask_for;
    forward_rxn_log_likelihood_ratio = (double*)calloc(one_l,ask_for);
    if (forward_rxn_log_likelihood_ratio == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc8: Error unable to allocate %lld bytes for forward_rxn_log_likelihood_ratio\n", ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->forward_rxn_log_likelihood_ratio = forward_rxn_log_likelihood_ratio;
    usage += ask_for;
    reverse_rxn_log_likelihood_ratio = (double*)calloc(one_l,ask_for);
    if (reverse_rxn_log_likelihood_ratio == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc8: Error unable to allocate %lld bytes for reverse_rxn_log_likelihood_ratio\n", ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->reverse_rxn_log_likelihood_ratio = reverse_rxn_log_likelihood_ratio;
    ask_for = (number_reactions + number_reactions + 2) * sizeof(double);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    rxn_likelihood_ps = (double *)calloc(one_l,ask_for);
    if (rxn_likelihood_ps == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc8: Error unable to allocate %lld bytes for rxn_likelihood_ps\n", ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->rxn_likelihood_ps = rxn_likelihood_ps;
  }
  state->usage += usage;
  return(success);
}
