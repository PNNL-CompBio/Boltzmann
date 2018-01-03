#include "boltzmann_structs.h"

#include "print_net_lklhd_bndry_flux_header.h"
void print_net_lklhd_bndry_flux_header(struct state_struct *state) {
  /*
    open the net likelihood boundary flux file and print a header
    with the fixed concentration species (and their compartments.
    Called by: ode23tb
    Calls: fopen, fprint, fflush
  */
  struct molecule_struct *molecule;
  struct compartment_struct *sorted_cmpts;
  struct compartment_struct *cur_cmpt;
  int64_t nunique_molecules;
  char *nl_bndry_flx_file;
  char *molecules_text;
  char *compartment_text;
  char *cmpt_string;
  char *molecule_str;

  int i;
  int ci;

  FILE *nl_bndry_flx_fp;
  FILE *lfp;

  nunique_molecules = state->nunique_molecules;
  molecule          = state->sorted_molecules;
  sorted_cmpts      = state->sorted_cmpts;
  molecules_text    = state->molecules_text;
  compartment_text  = state->compartment_text;
  nl_bndry_flx_file = state->nl_bndry_flx_file;
  nl_bndry_flx_fp  = fopen(nl_bndry_flx_file,"w");
  if (nl_bndry_flx_fp != NULL) {
    state->nl_bndry_flx_fp = nl_bndry_flx_fp;
    fprintf (nl_bndry_flx_fp," Net Likelihood Boundary Fluxes\n");
    fprintf (nl_bndry_flx_fp,"time\\species");
    for (i=0;i<nunique_molecules;i++) {
      if (molecule->variable == 0) {
	molecule_str = (char *)&molecules_text[molecule->string];
	ci = molecule->c_index;
	if (ci > 0) {
	  cur_cmpt = (struct compartment_struct *)&(sorted_cmpts[ci]);
	  cmpt_string = (char *)&compartment_text[cur_cmpt->string];
	  fprintf(nl_bndry_flx_fp,"\t%s:%s",cmpt_string,molecule_str);
	} else {
	  fprintf(nl_bndry_flx_fp,"\t%s",molecule_str);
	}
      }
      molecule += 1; /* Caution address arithmetic; */
    } /* end for (i...) */
    fprintf(nl_bndry_flx_fp,"\n");
    fflush(nl_bndry_flx_fp);
  } /* end if (nl_bndry_flx_fp != NULL) */
}

