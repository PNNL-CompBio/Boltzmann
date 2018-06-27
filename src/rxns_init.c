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
#define DBG 1
*/
#ifdef DBG
  struct molecule_struct *molecules;
  struct molecule_struct *molecule;
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  char   *molecules_text;
  char   *compartment_text;
  char   *compartment_name;
  char   *molecule_name;

  int number_molecules;
  int number_compartments;

  int i;
  int ci;

  int mi;
  int padj;

  FILE *lfp;
  FILE *efp;
#endif
  /*
    Read reactions file to count molecules and reactions
    Then allocate space then read for real.
  */
  success = parse_reactions_file(state,state->reaction_file);
  /*
    Lets have a look at the unsorted molecules and unsorted compartments.
    number_molecules and number_compartments fields of state are
    set in size_rxns_file call.
  */
#ifdef DBG
  lfp = state->lfp;
  molecules = state->unsorted_molecules;
  compartments = state->unsorted_cmpts;
  molecules_text = state->molecules_text;
  compartment_text = state->compartment_text;
  number_molecules = state->number_molecules;
  number_compartments = state->number_compartments;
  if (lfp) {
    fprintf(lfp,"Raw molecule:compartment data, number_molecules = %ld, number_compartments = %ld, number_reactions = %ld\n",state->number_molecules,state->number_compartments,state->number_reactions);
    fprintf(lfp,"i\tm_index\tc_index\tmolecule:compartment\n");
    molecule = molecules;
    for (i=0;i<number_molecules;i++) {
      mi = molecule->m_index;
      ci = molecule->c_index;
      compartment = (struct compartment_struct *)&compartments[ci];
      molecule_name = (char*)&molecules_text[molecule->string];
      compartment_name=(char*)&compartment_text[compartment->string];
      if (ci > 0) {
	fprintf(lfp,"%d\t%d\t%d\t%s:%s\n",
		i,mi,ci,molecule_name,compartment_name);
      } else {
	fprintf(lfp,"%d\t%d\t%d\t%s\n",
		i,mi,ci,molecule_name);
      }
      molecule += 1;/* Caution address arithmetic here */
    }
    compartment = compartments;
    fprintf(lfp,"\ni\tc_index\tcompartment\n");
    for (i=0;i<number_compartments;i++) {
      ci = compartment->c_index;
      if (ci > 0) {
	compartment_name=(char*)&compartment_text[compartment->string];
	fprintf(lfp,"%d\t%d\t%s\n",i,ci,compartment_name);
      }
      compartment += 1;/* Caution address arithmetic here */
    }
    fflush(lfp);
  }
#endif 
  /*
    First we need to sort the compartments.
  */
  if (success) {
    success = sort_compartments(state->unsorted_cmpts,
				state->sorted_compartments,
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
