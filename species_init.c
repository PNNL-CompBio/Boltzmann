#include "boltzmann_structs.h"

#include "set_compartment_ptrs.h"
#include "translate_regulation_metabolites.h"
#include "read_compartment_sizes.h"
#include "read_initial_concentrations.h"
#include "check_initial_concentrations.h"
#include "set_count_trans.h"

#include "species_init.h"
int species_init(struct state_struct *state) {
/*
    Initialize the conncentrations and compartment
    sizes for a reactions file.
    
    Called by: boltzmann_init_core
    Calls:     set_compartment_ptrs,
               translate_regulation_metabolites,
	       read_compartment_sizes,
	       read_intial_concentrations,
	       check_intial_concentrations
  */
  int success;
  int padi;
  success = 1;
  /*
    Now in order to enable molecule/compartment lookup for
    read_initial_concentrations we need to set the compartment 
    pointers, these are pointers into the list of sorted molecules
    All the molecules within a compartment are adjacent in the
    sorted molecules list. This call just sets pointers to the
    first molecule in each compartment, and one past the last
    molecule so that molecules in compartment i in the sorted molecules list
    are in positions compartment_ptrs[i]:compartment_ptrs[i+1]-1 
    inclusive.
  */
  if (success) {
    success = set_compartment_ptrs(state);
  }
  if (success) {
    success = translate_regulation_metabolites(state);
  }
  /*
    If we are going to read in a list of compartment sizes that should
    happen here after the compartments have been set, and before
    the initinal concentrations are read in.
  */
  if (success) {
    if (state->compartment_file[0] != '\0') {
      success = read_compartment_sizes(state);
    } 
  }
  /*
    Read initial concentrations, convert them to counts,
    and print them to the counts output file.
  */
  if (success) {
    success = read_initial_concentrations(state);
  }
  if (success) {
    success = check_initial_concentrations(state);
  }
  if (success) {
    /*
      Set the count_to_conc and conc_to_count vectors, one entry per unique
      molecule based on its compartment size.
    */
    success = set_count_trans(state);
  }
  return(success);
}
