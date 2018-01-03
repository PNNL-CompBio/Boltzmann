#include "boltzmann_structs.h"

#include "print_net_lklhd_bndry_flux.h"
void print_net_lklhd_bndry_flux(struct state_struct *state, 
				double *net_lklhd_bndry_flux,
				double time) {
  /*
    Print the net likelihood boundary fluxes.
    Called by: ode23tb
    Calls: fopen, fprint, fflush
  */
  struct molecule_struct *molecule;
  int64_t nunique_molecules;

  int i;
  int padi;

  FILE *nl_bndry_flx_fp;
  FILE *lfp;

  nunique_molecules = state->nunique_molecules;
  molecule  = state->sorted_molecules;
  nl_bndry_flx_fp  = state->nl_bndry_flx_fp;
  if (nl_bndry_flx_fp != NULL) {
    fprintf(nl_bndry_flx_fp,"%le",time);
    for (i=0;i<nunique_molecules;i++) {
      if (molecule->variable == 0) {
	fprintf(nl_bndry_flx_fp,"\t%le",net_lklhd_bndry_flux[i]);
      }
      molecule += 1; /* Caution address arithmetic; */
    } /* end for (i...) */
    fprintf(nl_bndry_flx_fp,"\n");
    fflush(nl_bndry_flx_fp);
  } /* end if (nl_bndry_flux_fp != NULL) */
}

