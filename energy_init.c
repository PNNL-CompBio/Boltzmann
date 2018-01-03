#include "boltzmann_structs.h"

#include "compute_standard_energies.h"
#include "compute_ke.h"
#include "zero_solvent_coefficients.h"
#include "compute_kss.h"

#include "energy_init.h"
int energy_init(struct state_struct *state) {
/*
    Initialize the delta g0's, kinetic energies, and kinetic
    steady state coefficients, and zero out the solvent coefficients
    in the reactiosn matrix.
    
    Called by: boltzmann_init_core
    Calls:     compute_standard_energies,
	       compute_ke,
	       zero_solvent_coefficients,
	       compute_kss
  */
  int success;
  int padi;
  /*
    Compute the reaction energies of formation if called for.
  */
  success = 1;
  if (state->use_pseudoisomers) {
    success = compute_standard_energies(state);
  }
  /*
    Compute the reaction ke's.
  */
  if (success) {
    success = compute_ke(state);
  }
  /*
    Zero out the solvent molecule coefficients in the reactions matrix.
  */
  if (success) {
    success = zero_solvent_coefficients(state);
  }
  /*
    Compute the reaction kss's.
  */
  if (success) {
    success = compute_kss(state);
  }
  return(success);
}
