#include "boltzmann_structs.h"
#include "create_output_filenames.h"
#include "open_output_files.h"
#include "sbml_to_boltzmann.h"
#include "size_rxns_file.h"

#include "io_size_init.h"
int io_size_init(struct state_struct *state) {
/*
  Set the output filenames, open the outputfiles, 
  and size the reactions file.
    
    Called by: boltzmann_init_core
    Calls:     create_output_filenames,
	       open_output_files,
	       sbml_to_boltzman,
	       size_rxns_file

  */
  int success;
  int print_output;
  
  FILE *lfp;
  FILE *efp;

  success = create_output_filenames(state);

  if (success) {
    print_output = (int)state->print_output;
    if (print_output) {
      success = open_output_files(state);
    }
  }
  /*
    Here is where the sbml file processing will need to go.
  */
  if (success) {
    if (state->sbml_file[0] != '\0') {
      success = sbml_to_boltzmann(state);
    }
  }
  if (success) {
    success = size_rxns_file(state,state->reaction_file);
  }
#ifdef DBG
  lfp = state->lfp;
  if (lfp) {
    fprintf(lfp,"io_size_init: after size_rxns_file, success = %d\n",success);
    fprintf(lfp,"io_size_init: rxns = %ld, cmpts = %ld, molecules = %ld, "
	    "reaction_file_length = %ld\n",
	    state->number_reactions,state->number_compartments,
	    state->number_molecules,
	    state->reaction_file_length);
    fprintf(lfp,"io_size_init: molecules_len = %ld, reaction_title_len = %ld, "
	    "pathway_len = %ld, compartment_len = %ld\n",
	    state->molecules_len,
	    state->rxn_title_len,
	    state->pathway_len,
	    state->compartment_len);
    fflush(lfp);
  }
#endif
  return(success);
}
