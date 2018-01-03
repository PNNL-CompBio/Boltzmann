#include "boltzmann_structs.h"

#include "print_net_likelihood_header.h"
void print_net_likelihood_header(struct state_struct *state) {
  /*
    open the net likelihood file and print a header
    with the reaction titles
    Caleld by: ode23tb
    Calls: fopen, fprint, fflush
  */
  struct rxn_struct *reaction;
  int64_t number_reactions;
  char *net_lklhd_file;
  char *rxn_title_text;
  char *title;

  int  i;
  int  padi;
  FILE *net_lklhd_fp;
  FILE *lfp;

  number_reactions  = state->number_reactions;
  rxn_title_text   = state->rxn_title_text;
  reaction         = state->reactions;
  net_lklhd_file = state->net_lklhd_file;
  net_lklhd_fp  = fopen(net_lklhd_file,"w");
  if (net_lklhd_fp != NULL) {
    state->net_lklhd_fp = net_lklhd_fp;
    fprintf (net_lklhd_fp," Reaction Net Likelihoods\n");
    fprintf (net_lklhd_fp,"time\\reaction");
    for (i=0;i<number_reactions;i++) {
      title = (char*)&rxn_title_text[reaction->title];
      fprintf(net_lklhd_fp,"\t%s",title);
      reaction += 1; /* Caution address arithmetic */
    }
    fprintf(net_lklhd_fp,"\n");
    fflush(net_lklhd_fp);
  } /* end if (net_lklhd_fp != NULL) */
}

