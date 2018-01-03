#include "boltzmann_structs.h"

#include "print_concs_fluxes.h"
int print_concs_fluxes(struct state_struct *state,int ny,
		       double *fluxes, double *concs, 
		       double *forward_rxn_likelihoods,
		       double *reverse_rxn_likelihoods,
		       double time, double h) {
/* 
  Debugging routine to print concentraions and fluxes for each molecule.
  Called by: ode23tb 
  Calls  fprintf
*/
  struct molecule_struct *molecule;  
  char   *molecules_text;
  FILE *lfp;
  FILE *efp;
  int  i;
  int  success;
  int  nrxns;
  int  padi;
  success = 1;
  lfp =  state->lfp;
  molecule = state->sorted_molecules;
  molecules_text = state->molecules_text;
  nrxns          = state->number_reactions;
  if (lfp) {
    fprintf(lfp," time = %le, h = %le\n"
	    "Molecule\tindex\tconcentration\tflux\n",time,h);
    for (i=0;i<ny;i++) {
      fprintf(lfp,"%s\t%d\t%le\t%le\n",
	      (char*)&molecules_text[molecule->string],
	      i,concs[i],fluxes[i]);
      molecule   += 1; /* Caution address arithmetic here. */
    }
    fprintf(lfp,"rxn_n\tfrwrd_lklhd\trvrs_lklhd\n");
    for (i=0;i<nrxns;i++) {
      fprintf(lfp,"%d\t%le\t%le\n",i,forward_rxn_likelihoods[i],
	      reverse_rxn_likelihoods[i]);
    }
    fflush(lfp);
  }
  
  return(success);
}
