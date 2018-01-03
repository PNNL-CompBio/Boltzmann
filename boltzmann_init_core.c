#include "boltzmann_structs.h"
#include "io_size_init.h"
#include "alloc2.h"
#include "rxns_init.h"
#include "alloc3.h"
#include "species_init.h"
#include "alloc4.h"
#include "form_molecules_matrix.h"
#include "energy_init.h"
#include "echo_inputs.h"
#include "run_init.h"

#include "boltzmann_init_core.h"
int boltzmann_init_core(struct state_struct *state) {
/*
    Initialize the reactions and data structures for boltzmann, after
    a state variable has been allocated and the parameters read.
    
    Called by: boltzmann_init, boltzmann_boot
    Calls:     io_size_init,
	       alloc2,
	       rxns_init,
	       alloc3,
	       species_init
	       alloc4
	       form_molecules_matrix
	       energy_init
	       echo_inputs
	       run_init
  */
  int success;
  int print_output;

  print_output = state->print_output;
  success = io_size_init(state);
  /*
    At this point in time we can compute how much space the aligned
    reaction titles, pathway descriptions, compartments and molecules
    verbage will take.
    We want an uppercase version of the molecules, so we need two 
    molecules copies. Also we have an upperbound of the number of molecules,
    state->number_molecules, and we can allocate space for the molecules 
    sorting.
  */
  if (success) {
    success = alloc2(state);
  }
  /*
    Read reactions file to count molecules and reactions,
    Sort the compartments, remove duplicates and build a 
    translation table, sort the molecules (species) and 
    remove duplicates.
  */
  if (success) {
    success = rxns_init(state);
  }
  /*
    Now we need to allocate space for the counts, concentrations,
    and read in the intial concentrations converting them to counts.
  */
  if (success) {
    success = alloc3(state);
  }
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
    success = species_init(state);
  }
  /*
    Compute the molecules matrix.
  */
  if (success) {
    success = alloc4(state,&state->molecules_matrix,&state->transpose_workspace);
  }
  if (success) {
    success = form_molecules_matrix(state,state->molecules_matrix,state->transpose_workspace);
  }
  /*
    Compute the reaction energies of formation if called for.
  */
  if (success) {
    success = energy_init(state);
  }
  if (success) {
    if (print_output) {
      success = echo_inputs(state);
    }
  }
  if (success) {
    success = run_init(state);
  }
  return(success);
}
