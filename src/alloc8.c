#include "boltzmann_structs.h"
#include "alloc8.h"
/*
int alloc8(struct state_struct *state, int direction, void *workspace, void **next_addr_p) {
*/
int alloc8(struct state_struct *state) {
  /*
    Allocate workspace vectors. This used to happen in the
    flatten_state routine.
    If direction is 0, allocate space and set pointers to it in the state
    structure. 
    If diretion is 1, set pointers into the address starting at workspace,
    and set *next_addr_p to be the first byte after the allocated space.

    Allocates space for and sets the following pointers in state.
      future_counts,
      free_energy,
      forward_rxn_likelihood,
      reverse_rxn_likelihood,
      forward_rxn_log_likelihood_ratio,
      reverse_rxn_log_likelihood_ratio,
      rxn_likelihood_ps
      
    Called by: run_init, boltzmann_flatten_state
    Calls:     calloc, fprintf, fflush
  */
  /*
  void    *next_addr;
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
  int64_t number_molecules;
  int64_t align_len;
  int64_t align_mask;
  int64_t data_pad;
  int64_t ask_for;
  int64_t one_l;
  int64_t run_workspace_bytes;

  int     success;
  int     padi;

  FILE    *lfp;
  FILE    *efp;
  
  align_len            = state->align_len;
  align_mask           = state->align_mask;
  number_reactions     = state->number_reactions;
  number_molecules     = state->nunique_molecules;
  lfp                  = state->lfp;
  run_workspace_bytes  = state->run_workspace_bytes;
  /*
  next_addr            = workspace;
  */
  success              = 1;
  usage = (int64_t)0;
  one_l = (int64_t)1;
  
  ask_for = number_molecules * sizeof(double);
  data_pad = (align_len - (ask_for & align_mask)) & align_mask;
  ask_for += data_pad;
  usage   += ask_for;
  run_workspace_bytes += ask_for;

  /*
  if (direction == 0) {
  */
    future_counts = (double *)calloc(one_l,ask_for);
    if (future_counts == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc8: Error unable to allocate %ld bytes for future_counts\n", ask_for);
	fflush(lfp);
      }
    } else {
      state->future_counts = future_counts;
    }
  /*
  } else {
    state->future_counts = (double*)next_addr;
    next_addr += ask_for;
  }
  */
  if (success) {
    ask_for = number_reactions * sizeof(double);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    run_workspace_bytes += ask_for;

    /*
    if (direction == 0) {
    */
      free_energy = (double *)calloc(one_l,ask_for);
      if (free_energy == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,
		  "alloc8: Error unable to allocate %ld bytes for free_energy\n", 
		  ask_for);
	  fflush(lfp);
	}
      } else {
	state->free_energy = free_energy;
      }
    /*
    } else {
      state->free_energy = (double*)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      forward_rxn_likelihood = (double *)calloc(one_l,ask_for);
      if (forward_rxn_likelihood == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,
		  "alloc8: Error unable to allocate %ld bytes for forward_rxn_likelihood\n",
		  ask_for);
	  fflush(lfp);
	}
      } else {
	state->forward_rxn_likelihood = forward_rxn_likelihood;
      }
    /*
    } else {
      state->forward_rxn_likelihood = (double*)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      reverse_rxn_likelihood = (double *)calloc(one_l,ask_for);
      if (reverse_rxn_likelihood == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,
    "alloc8: Error unable to allocate %ld bytes for reverse_rxn_likelihood",
		  ask_for);
	  fflush(lfp);
	}
      } else {
	state->reverse_rxn_likelihood = reverse_rxn_likelihood;
      }
    /*
    } else {
      state->reverse_rxn_likelihood = (double *)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      forward_rxn_log_likelihood_ratio = (double*)calloc(one_l,ask_for);
      if (forward_rxn_log_likelihood_ratio == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc8: Error unable to allocate %ld bytes for forward_rxn_log_likelihood_ratio\n", ask_for);
	  fflush(lfp);
	}
      } else {
	state->forward_rxn_log_likelihood_ratio = forward_rxn_log_likelihood_ratio;
      }
    /*
    } else {
      state->forward_rxn_log_likelihood_ratio = (double*)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      reverse_rxn_log_likelihood_ratio = (double*)calloc(one_l,ask_for);
      if (reverse_rxn_log_likelihood_ratio == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc8: Error unable to allocate %ld bytes for reverse_rxn_log_likelihood_ratio\n", ask_for);
	  fflush(lfp);
	}
      } else {
	state->reverse_rxn_log_likelihood_ratio = reverse_rxn_log_likelihood_ratio;
      }
    /*
    } else {
      state->reverse_rxn_log_likelihood_ratio = (double*)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    ask_for = (number_reactions + number_reactions + 2) * sizeof(double);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      rxn_likelihood_ps = (double *)calloc(one_l,ask_for);
      if (rxn_likelihood_ps == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc8: Error unable to allocate %ld bytes for rxn_likelihood_ps\n", ask_for);
	  fflush(lfp);
	}
      } else {
	state->rxn_likelihood_ps = rxn_likelihood_ps;
      }
    /*
    } else {
      state->rxn_likelihood_ps = (double*)next_addr;
      next_addr += ask_for;
    }
    */
  }
  state->usage += usage;
  state->run_workspace_bytes = run_workspace_bytes;
  /*
  *next_addr_p = next_addr;
  */
  return(success);
}
