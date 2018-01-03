#include "boltzmann_structs.h"

#include "molecules_lookup.h"

#include "translate_regulation_metabolites.h"
int translate_regulation_metabolites(struct state_struct *state) {
  /*
    Translate the reg_species entries from pointers into the regulation_text
    array to species position, overwriting the reg_species entry with the
    metabolite location in the concentration array.
    Called by: species_init.
    Calls:     molecules_lookup, fprintf, fflush
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;
  struct reactions_matrix_struct *rxns_matrix;
  int64_t *reg_species;
  int64_t *compartment_indices;
  char    *regulation_text;
  char    *metabolite;
  
  int success;
  int max_regs_per_rxn;

  int num_rxns;
  int reg_pos;

  int reg_base;
  int rmet_num;

  int i;
  int j;

  int lc;
  int rc;

  int met_pos;
  int met_num;

  FILE *lfp;
  FILE *efp;
  success = 1;
  if (state->use_regulation) {
    lfp                 = state->lfp;
    max_regs_per_rxn    = (int)state->max_regs_per_rxn;
    num_rxns            = (int)state->number_reactions;
    reactions           = state->reactions;
    regulation_text     = state->regulation_text;
    reg_species         = state->reg_species;
    rxns_matrix         = state->reactions_matrix;
    compartment_indices = rxns_matrix->compartment_indices;

    reg_pos          = 0;
    reg_base         = 0;
    reaction         = reactions;
    for (i=0;i<num_rxns;i++) {
      lc = reaction->left_compartment;
      rc = reaction->right_compartment;
      lc = compartment_indices[lc];
      rc = compartment_indices[rc];
      for (j=0; ((j<max_regs_per_rxn) && success); j++) {
	reg_pos = reg_base + j;
	met_pos = reg_species[reg_pos];
	if (met_pos >=0) {
	  /*
	    This regulation exists.
	  */
	  metabolite = (char*)&regulation_text[met_pos];
	  met_num = molecules_lookup(metabolite,lc,state);
	  if (lc != rc) {
	    /*
	      Left and Right compartments are different.
	    */
	    rmet_num = molecules_lookup(metabolite,rc,state);
	    if (met_num < 0) {
	      /*
		Left compartment did not have metabolite.
	      */
	      met_num = rmet_num;
	    } else {
	      /*
		Left compartment had metabolite. 
	      */
	      if (rmet_num >= 0) {
		if (lfp) {
		  fprintf(lfp,"translate_regulation_metabolites: WARNING: "
			  "metabolite %s is in two different compartments "
			  "for reaction %d using the left compartment value\n",
			  metabolite,i);
		  fflush(lfp);
		} /* end if (lfp) */
	      } /* end if both compartments had metabolite. */
	    } /* end else left compartment had metabolite. */
	  } /* end if left and right compartments differ. */
	  if (met_num >= 0) {
	    /*
	      Metabolite found.
	    */
	    reg_species[reg_pos] = (int64_t)met_num;
	  } else {
	    /* 
	      Metabolite not found.
	    */
	    success = 0;
	    if (lfp) {
	      fprintf(lfp,"translate_regulation_metabolites: Error, metabolite %s, not found in compartment %d\n",metabolite,lc);
	      fflush(lfp);
	    }
	  } /* end else metabolite not found. */
	} /* end if this regulation exists. */
      } /* end for (j...) */
      /*
        Caution address arithmetic on the following line.
      */
      reg_base += max_regs_per_rxn;
      reaction += 1;
    } /* end for (i...) */ 
  } /* end if using regulation */
  return(success);
}
