#include "boltzmann_structs.h"

#include "ode_print_bflux_header.h"
void ode_print_bflux_header(struct state_struct *state) {
  /*
    open the ode boundary flux file and print a header
    with the fixed concentration species (and their compartments.
    Called by: deq_run
    Calls: fopen, fprint, fflush
  */
  struct molecule_struct *molecule;
  struct compartment_struct *sorted_cmpts;
  struct compartment_struct *cur_cmpt;
  int64_t nunique_molecules;
  char *ode_bflux_file;
  char *molecules_text;
  char *compartment_text;
  char *cmpt_string;
  char *molecule_str;

  int i;
  int ci;

  FILE *ode_bflux_fp;
  FILE *lfp;

  nunique_molecules = state->nunique_molecules;
  molecule          = state->sorted_molecules;
  sorted_cmpts      = state->sorted_cmpts;
  molecules_text    = state->molecules_text;
  compartment_text  = state->compartment_text;
  ode_bflux_file    = state->ode_bflux_file;
  ode_bflux_fp      = fopen(ode_bflux_file,"w");
  if (ode_bflux_fp != NULL) {
    state->ode_bflux_fp = ode_bflux_fp;
    fprintf (ode_bflux_fp," Ode Boundary Fluxes\n");
    fprintf (ode_bflux_fp,"time\\species");
    for (i=0;i<nunique_molecules;i++) {
      if (molecule->variable == 0) {
	molecule_str = (char *)&molecules_text[molecule->string];
	ci = molecule->c_index;
	if (ci > 0) {
	  cur_cmpt = (struct compartment_struct *)&(sorted_cmpts[ci]);
	  cmpt_string = (char *)&compartment_text[cur_cmpt->string];
	  fprintf(ode_bflux_fp,"\t%s:%s",cmpt_string,molecule_str);
	} else {
	  fprintf(ode_bflux_fp,"\t%s",molecule_str);
	}
      }
      molecule += 1; /* Caution address arithmetic; */
    } /* end for (i...) */
    fprintf(ode_bflux_fp,"\n");
    fflush(ode_bflux_fp);
  } /* end if (ode_bflux_fp != NULL) */
}

