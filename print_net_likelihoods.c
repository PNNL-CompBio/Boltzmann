#include "boltzmann_structs.h"

#include "print_net_likelihoods.h"
void print_net_likelihoods(struct state_struct *state, double time) {
  /*
    open the net likelihood file and print a header
    with the reaction titles
    Caleld by: ode23tb
    Calls: fopen, fprint, fflush
  */
  struct rxn_struct *reaction;
  int64_t number_reactions;
  double  *net_likelihood;

  int  i;
  int  padi;
  FILE *net_lklhd_fp;
  FILE *lfp;

  net_lklhd_fp      = state->net_lklhd_fp;
  number_reactions  = state->number_reactions;
  net_likelihood    = state->net_likelihood;
  if (net_lklhd_fp != NULL) {
    fprintf(net_lklhd_fp,"%le",time);
    for (i=0;i<number_reactions;i++) {
      fprintf(net_lklhd_fp,"\t%le",net_likelihood[i]);
    }
    fprintf(net_lklhd_fp,"\n");
    fflush(net_lklhd_fp);
  } /* end if (net_lklhd_fp != NULL) */
}
