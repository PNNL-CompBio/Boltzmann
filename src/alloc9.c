#include "boltzmann_structs.h"
#include "alloc9.h"
/*
int alloc9(struct state_struct *state,
	   int direction,
	   void *workspace,
	   void **next_addr_p) {
*/
int alloc9(struct state_struct *state) {
  /*
    Allocate workspace needed for printing reaction view pieces in
    the record loop of boltzmann_run.
    Allocates space for and sets
         no_op_likelihood, 
         rxn_view_likelihoods,
	 rev_rxn_view_likelihoods,
	 rxn_fire,
	 rxn_mat_row 
    pointers in state

    Called by: run_init, boltzmann_flatten_alloc1
    Calls:     calloc, fprintf, fflush
  */
  /*
  void    *next_addr;
  */
  double  *no_op_likelihood;
  double  *rxn_view_likelihoods;
  double  *rev_rxn_view_likelihoods;
  double  *rxn_mat_row;
  int64_t *rxn_fire;
  int64_t usage;
  int64_t rxn_view_hist_length;
  int64_t number_reactions;
  int64_t number_species;
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
  rxn_view_hist_length = state->rxn_view_hist_length;
  number_reactions     = state->number_reactions;
  number_species       = state->nunique_molecules;
  lfp                  = state->lfp;
  run_workspace_bytes  = state->run_workspace_bytes;
  /*
  next_addr            = workspace;
  */
  success              = 1;
  usage = (int64_t)0;
  one_l = (int64_t)1;
  
  ask_for = rxn_view_hist_length * sizeof(double);
  data_pad = (align_len - (ask_for & align_mask)) & align_mask;
  ask_for += data_pad;
  usage   += ask_for;
  run_workspace_bytes += ask_for;
  /*
  if (direction == 0) {
  */
    no_op_likelihood = (double *)calloc(one_l,ask_for);
    if (no_op_likelihood == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc9: Error unable to allocate %ld bytes for no_op_likelihood, turning off printing\n", ask_for);
	fflush(lfp);
      }
      state->print_output = 0;
    } else {
      state->no_op_likelihood = no_op_likelihood;
    } 
  /*
  } else {
    state->no_op_likelihood = (double*)next_addr;
    next_addr += ask_for;
  }
  */
  if (success) {
    ask_for = rxn_view_hist_length * sizeof(double) * number_reactions;
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      rxn_view_likelihoods = (double *)calloc(one_l,ask_for);
      if (rxn_view_likelihoods == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc9: Error unable to allocate %ld bytes for rxn_view_likelihoods, turning off printing\n", ask_for);
	  fflush(lfp);
	}
	state->print_output = 0;
      } else {
	state->rxn_view_likelihoods = rxn_view_likelihoods;
      }
    /*
    } else {
      state->rxn_view_likelihoods = (double *)next_addr;
      next_addr += ask_for;
    }      
    */
  }
  if (success) {

    usage+= ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      rev_rxn_view_likelihoods = (double *)calloc(one_l,ask_for);
      if (rev_rxn_view_likelihoods == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc9: Error unable to allocate %ld bytes for rev_rxn_view_likelihoods, turning off printing\n", ask_for);
	  fflush(lfp);
	}
	state->print_output = 0;
      } else {
	state->rev_rxn_view_likelihoods = rev_rxn_view_likelihoods;
      }
    /*
    } else {
      state->rev_rxn_view_likelihoods = (double *)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    ask_for = (number_reactions + 1) * 2 * sizeof(int64_t);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      rxn_fire = (int64_t*)calloc(one_l,ask_for);
      if (rxn_fire == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc9: Error unable to allocate %ld bytes for rxn_fire, turning off printing\n", ask_for);
	  fflush(lfp);
	}
	state->print_output = 0;
      } else {
	state->rxn_fire = rxn_fire;
      }
    /*
    } else {
      state->rxn_fire = (int64_t *)next_addr;
      next_addr += ask_for;
    }
    */
  }
  if (success) {
    ask_for = number_species * sizeof(double);
    data_pad = (align_len - (ask_for & align_mask)) & align_mask;
    ask_for += data_pad;
    usage += ask_for;
    run_workspace_bytes += ask_for;
    /*
    if (direction == 0) {
    */
      rxn_mat_row = (double*)calloc(one_l,ask_for);
      if (rxn_mat_row == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"alloc9: Error unable to allocate %ld bytes for rxn_mat_row, turning off printing\n", ask_for);
	  fflush(lfp);
	}
	state->print_output = 0;
      } else {
	state->rxn_mat_row = rxn_mat_row;
      }
    /*
    } else {
      state->rxn_mat_row = (int*)next_addr;
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
