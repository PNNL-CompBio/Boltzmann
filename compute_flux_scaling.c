#include "boltzmann_structs.h"

#include "compute_flux_scaling.h"
double compute_flux_scaling(struct state_struct *state, double *concs) {
  /*
    if flux_scaling field of state is 0, 
      Compute the kf_base_reaction * 
      (product of base_reaciont reactant concentrations).
    else 
      use the value in state.

    Called by: num_jac_col, ode23tb
    Uses kf_base_reaction, number_base_reaction_reactants
         and base_reactants fields of state.
    Concs is an input vector of length nunique_species with the 
    concentrations in them.
  */
  double kf_base;
  double flux_scaling;
  int *base_reactants;
  int nbr;
  int i;
  int species;
  int padi;
  if (state->flux_scaling == 0) {
    nbr = state->number_base_reaction_reactants;
    flux_scaling = state->kf_base_reaction;
    base_reactants = state->base_reactants;
    for (i=0;i<nbr;i++) {
      species = base_reactants[i];
      flux_scaling = flux_scaling * concs[species];
    }
  } else {
    flux_scaling = state->flux_scaling;
  }
  return (flux_scaling);
}
