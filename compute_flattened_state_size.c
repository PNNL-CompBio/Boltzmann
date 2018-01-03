#include "boltzmann_structs.h"
#include "compute_flattened_state_size.h"
int64_t compute_flattened_state_size(struct state_struct *state) {
  /*
    Compute the size of the flattened state.
    Called by: boltzmann_flatten_alloc, boltzmann_boot
  */
  int64_t fsize;
  int nr;
  int nm;
  int nc;
  int nz;
  nr = state->number_reactions;
  nm = state->nunique_molecules;
  nc = state->nunique_compartments;
  nz = state->number_molecules;
  fsize = (int64_t)276 + (34*nr) + (15*nm) + (8*nc) + (5*nz);
  return(fsize);
}  
