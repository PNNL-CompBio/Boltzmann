#include "boltzmann_structs.h"

#include "init_base_reactants.h"
int init_base_reactants(struct state_struct *state) {
  /*
    Initialize the base_reactant_indicator vector (length nunique_molecules),
    and the base_reactants vector length number_base_reaction_reactants,
    and the number_base_reaction_reactants value.
    Called by: ode23tb
  */
  struct reaction_struct *reactions; 
  struct reaction_struct *reaction; 
  struct reactions_matrix_struct *reactions_matrix; 
  int64_t *molecules_indices;
  int64_t *coefficients;
  int64_t *rxn_ptrs;
  int *base_reactant_indicator;
  int *base_reactants;
  int ny;
  int i;

  int nrxns;
  int number_base_reaction_reactants;

  int success;
  int r_count;

  int moleculei;
  int base_rxn;

  FILE *lfp;
  FILE *efp;
  success = 1;
  lfp       = state->lfp;
  reactions = state->reactions;
  nrxns     = state->number_reactions;
  ny        = state->nunique_molecules;
  /*
    Allocate integer space 
    We need a base_reactant_indicator vector of length ny (number of species)
    that is 1 if the species is a reactant (on the left side) of the
    base reaction and 0 otherwise.
    We needed a base_reactants vector length ny (though we won't use all 
    of them) that stores the species indices of species that are reactants
    in the base reaction (list of the indices of base_reactant_indicator
    vector that are 1), 
  */
  base_rxn = (int)state->base_reaction;
  /* 
    The following two vectors are allocated in alloc7.
  */
  base_reactant_indicator = state->base_reactant_indicator;
  base_reactants          = state->base_reactants;
  if ((base_rxn < 0) || (base_rxn >= nrxns) ) {
    success = 0;
    if (lfp) {
      fprintf(lfp,"init_base_reactants: Error base_rxn = %d, must be in [0, %d)\n",base_rxn,nrxns);
      fflush(lfp);
    }
  }
  if (success) {
    reaction  = (struct reaction_struct *)&reactions[base_rxn];
    number_base_reaction_reactants = reaction->num_reactants;
    reactions_matrix  = state->reactions_matrix;
    molecules_indices = reactions_matrix->molecules_indices;
    coefficients      = reactions_matrix->coefficients;
    rxn_ptrs          = reactions_matrix->rxn_ptrs;
    state->number_base_reaction_reactants = number_base_reaction_reactants;
    for (i=0;i<ny;i++) {
      base_reactant_indicator[i] = 0;
    }
    r_count = 0;
    for (i=rxn_ptrs[base_rxn];i<rxn_ptrs[base_rxn+1];i++) {
      if (coefficients[i] < 0) {
	moleculei = molecules_indices[i];
	base_reactant_indicator[moleculei] = 1;
	base_reactants[r_count] = moleculei;
	r_count += 1;
      }
    }
  }
  return(success);
}
