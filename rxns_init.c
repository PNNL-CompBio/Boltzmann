#include "boltzmann_structs.h"
#include "parse_reactions_file.h"
#include "sort_compartments.h"
#include "unique_compartments.h"
#include "translate_compartments.h"
#include "sort_molecules.h"
#include "unique_molecules.h"

#include "rxns_init.h"
int rxns_init(struct state_struct *state) {
/*
    Read in the reactions file, building a species dictionary, 
    of sorted species and compartments.
    
    Called by: boltzmann_init_core
    Calls:     parse_reactions_file,
	       sort_compartments,
	       unique_compartments,
	       translate_compartments,
	       sort_molecules,
	       unique_molecules,
  */
  int success;
  int padi;
  /*
    Read reactions file to count molecules and reactions
    Then allocate space then read for real.
  */
  success = parse_reactions_file(state,state->reaction_file);
  /*
    First we need to sort the compartments.
  */
  if (success) {
    success = sort_compartments(state->unsorted_cmpts,
				state->sorted_cmpts,
				state->compartment_text,
				state->number_compartments);
  }
  /*
    Then we extract the unique compartments.
    and fill the compartment_indices vector of the rxns_matrix structure.
  */
  if (success) {
    success = unique_compartments(state);
  }
  /*
    Now we need to assign the proper compartment numbers to the 
    unsorted molecules, using the compartment_indices field
    of the rxns_matrix structure.
  */
  if (success) {
    success = translate_compartments(state);
  }
  /*
    Now we need to sort the molecules, by compartment and name.
  */
  if (success) {
    success = sort_molecules(state->unsorted_molecules,
			     state->sorted_molecules,
			     state->molecules_text,
			     state->number_molecules);
  }
  /*
    Then we extract the unique molecules and set the
    molecules_indices field of the rxn_matrix struct.
    Also set the solvent_pos field of state.
  */
  if (success) {
    success = unique_molecules(state);
  }
  return(success);
}
