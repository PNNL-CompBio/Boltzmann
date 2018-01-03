#include "boltzmann_structs.h"
#include "echo_params.h"
#include "echo_reactions_file.h"
#include "print_molecules_dictionary.h"
#include "print_dg0_ke.h"
/*
#include "print_reactions_matrix.h"
*/

#include "echo_inputs.h"
int echo_inputs(struct state_struct *state) {
/*
    Echo the paramters, reactions, molecules dictionary,
    delta g0's, and the ke's to output files.
    
    Called by: boltzmann_init_core
    Calls:     echo_params,
	       echo_reactions_file,
	       print_molecules_dictionary,
	       print_dg0_ke,
	       print_reactions_matrix,

  */
  struct formation_energy_struct *formation_energies;
  struct vgrng_state_struct *vgrng_state;
  struct vgrng_state_struct *vgrng2_state;
  double *dg0s;
  double *free_energy;
  double *activities;
  int64_t vgrng_start;
  int64_t i;

  int success;
  int vgrng_start_steps;

  int print_output;
  int padi;
  
  FILE *bndry_flux_fp;
  FILE *lfp;

  lfp = state->lfp;
  /* 
    Echo params to the log file.
  */
  success = echo_params(lfp,state);
  if (success) {
    /*
      create the rxns.echo file.
    */
    success = echo_reactions_file(state);
  }
  if (success) {
    /*
      create the rxns.dict file.
    */
    success = print_molecules_dictionary(state);
  }
  if (success) {
    /*
      Create the rxns.dg0ke file.
    */
    success = print_dg0_ke(state);
  }
  if (success) {
    /*
      Create the rxns.mat file. print_reactions_matrix needs
      rxn_mat_row which is only set by flatten_state.
    success = print_reactions_matrix(state);
    */
  }
  return(success);
}
