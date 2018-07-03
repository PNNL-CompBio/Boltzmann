#include "boltzmann_structs.h"

#include "molecules_lookup.h"
#include "compartment_lookup.h"
#include "find_colon.h"


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
  char    *compartment_name;

  int num_rxns;
  int max_regs_per_rxn;

  int reg_pos;
  int reg_base;

  int i;
  int j;

  int met_pos;
  int met_num;

  int colon_loc;
  int ci;

  int success;
  int padi;

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
      for (j=0; ((j<max_regs_per_rxn) && success); j++) {
	reg_pos = reg_base + j;
	met_pos = reg_species[reg_pos];
	if (met_pos >=0) {
	  /*
	    This regulation exists.
	  */
	  metabolite = (char*)&regulation_text[met_pos];
	  /*
	    In storing the metabolite if it had a compartment that
	    compartment is still attached via a :compartment_name.
	    Extract that compartment name, look it up, and use that
	    for the compartment instead of inheriting it from
	    the reaction characteristics (left or right compartment).
	  */
	  colon_loc = find_colon(metabolite);
	  ci = 0;
	  if (colon_loc >=0) {
	    compartment_name = (char*)&metabolite[colon_loc+1];
	    ci = compartment_lookup(compartment_name,state);
	    metabolite[colon_loc] = '\0';
	    if (ci < 0) {
	      ci = 0;
	      if (lfp) {
		fprintf(lfp,"translate_regulation_metabolites: Warning: compartment %s not found, using default global compartment.\n",compartment_name);
		fflush(lfp);
	      }
	    }
	  }
	  met_num = molecules_lookup(metabolite,ci,state);
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
	      fprintf(lfp,"translate_regulation_metabolites: Error, metabolite %s, not found in compartment %d\n",metabolite,ci);
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
