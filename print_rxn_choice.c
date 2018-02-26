#include "boltzmann_structs.h"
#include "print_rxn_choice.h"
void print_rxn_choice(struct state_struct *state,
		      int64_t iter,
		      int rxn_choice) {
  /*
    Print the reaction choice and likelihood information to the logfile.
    Called by: bwarmp_run, boltzmann_watch
    Calls      fprintf,fflush
  */
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  int number_reactions;
  int noop_rxn;
  int rxn_no;
  int padi;
  FILE *lfp;
  FILE *efp;
  lfp              = state->lfp;
  if (lfp) {
    number_reactions = (int)state->number_reactions;
    noop_rxn = number_reactions + number_reactions;
    forward_rxn_likelihood = state->forward_rxn_likelihood;
    reverse_rxn_likelihood = state->reverse_rxn_likelihood;
    if (rxn_choice == noop_rxn) {
      fprintf(lfp,"reaction_choice: %ld\tnone\n",iter);
    } else {
      if (rxn_choice < number_reactions) {
	fprintf(lfp,"%ld\t%d\t%le\t%le\n",iter,rxn_choice,forward_rxn_likelihood[rxn_choice],
		reverse_rxn_likelihood[rxn_choice]);
      } else {
	rxn_no = rxn_choice - number_reactions;
	fprintf(lfp,"%ld\t%d\t%le\t%le\n",iter,rxn_choice,reverse_rxn_likelihood[rxn_no],
		forward_rxn_likelihood[rxn_no]);
      }
    }
    fflush(lfp);
  }
}
