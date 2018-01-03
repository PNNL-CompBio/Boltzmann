/* boltzmann_global_to_local_fluxes.c
*******************************************************************************
boltzmann

Pacific Northwest National Laboratory, Richland, WA 99352.

Copyright (c) 2010 Battelle Memorial Institute.

Publications based on work performed using the software should include 
the following citation as a reference:


Licensed under the Educational Community License, Version 2.0 (the "License"); 
you may not use this file except in compliance with the License. 
The terms and conditions of the License may be found in 
ECL-2.0_LICENSE_TERMS.TXT in the directory containing this file.
        
Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR 
CONDITIONS OF ANY KIND, either express or implied. See the License for the 
specific language governing permissions and limitations under the License.
******************************************************************************/

#include "boltzmann_structs.h"
#include "flatten_super_state.h"
#include "flatten_state.h"

#include "boltzmann_global_to_local_fluxes.h"
int boltzmann_global_to_local_fluxes(struct super_state_struct *super_state,
				    double *global_fluxes,
				    struct state_struct *local_state) {
  /*
    This routine takes a pointer to the superstate, a pointer to a vector
    of global fluxes, and a pointer to a local_state_struct
    and sets the flux values in the local state_struct from the
    global fluxes vector, using the molecules map in the super_state 
    struct.
    Called by: User (part of the API)
  */
  struct super_state_pointers_struct ssps;
  struct super_state_pointers_struct *super_state_pointers;
  struct state_struct *local_statep;
  double *local_fluxes;
  int64_t *molecule_map_starts;
  int64_t *molecule_map;
  int64_t *molecule_maps;
  int64_t state_index;
  int64_t nunique_molecules;
  int64_t i;
  int success;
  int padi;
  super_state_pointers = &ssps;
  success = flatten_super_state(super_state,super_state_pointers);
  if (success) {
    local_statep = local_state;
    success = flatten_state(local_state,&local_statep);
  }
  if (success) {
    state_index         = local_state->agent_type;
    local_fluxes        = local_state->bndry_flux_counts;
    nunique_molecules   = local_state->nunique_molecules;
    molecule_map_starts = super_state_pointers->molecule_map_starts;
    molecule_maps       = super_state_pointers->molecule_map;
    molecule_map        = &molecule_maps[molecule_map_starts[state_index]];
    for (i=0;i<nunique_molecules;i++) {
      local_fluxes[i]   = global_fluxes[molecule_map[i]];
    }
  }
  return (success);
}
