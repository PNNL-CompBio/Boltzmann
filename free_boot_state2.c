/* free_boot_state2.c
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
/*
#define DBG_FREE_BOOT_STATE2 1
*/
#include "free_boot_state2.h"
int free_boot_state2(struct state_struct *state) {
  /*
    Free space allocated to fields of the state_struct by alloc2, alloc3,
    alloc4, and alloc5 calls.
    Called by: boltzmann_init_core
    Calls:     free.
  */
  struct rxn_matrix_struct *reactions_matrix; 
  struct molecules_matrix_struct *molecules_matrix;
  struct formation_energy_struct *formation_energies;
  int success;
  int padi;
  success = 1;
  if (state) {
    /*
      Free allocations made in alloc2
    */
    if (state->rxn_title_text) {
      free(state->rxn_title_text);
    }
    if (state->reactions) {
      free(state->reactions);
    }
    reactions_matrix = state->reactions_matrix;
    if (reactions_matrix) {
      if (reactions_matrix->rxn_ptrs) {
	free(reactions_matrix->rxn_ptrs);
      }
      if (reactions_matrix->molecules_indices) {
	free(reactions_matrix->molecules_indices);
      }
      if (reactions_matrix->compartment_indices) {
	free(reactions_matrix->compartment_indices);
      }
      if (reactions_matrix->coefficients) {
	free(reactions_matrix->coefficients);
      }
      if (reactions_matrix->solvent_coefficients) {
	free(reactions_matrix->solvent_coefficients);
      }
      if (reactions_matrix->text) {
	free(reactions_matrix->text);
      }
      free(state->reactions_matrix);    
    }
    /*
      Need to be carefule with unsorted_molecules and unsorted_cmpts as
      they may have exchanged values with their sorted analogs.
    */
    if (state->unsorted_molecules) {
      if (state->sorted_molecules < state->unsorted_molecules) {
	free(state->sorted_molecules);
      } else {
	free(state->unsorted_molecules);
      }
    }
    if (state->unsorted_cmpts) {
      if (state->sorted_cmpts < state->unsorted_cmpts) {
        free(state->sorted_cmpts);
      } else {
	free(state->unsorted_cmpts);
      }
    }
    if (state->activities) {
      free(state->activities);
    }
    if (state->reg_constant) {
      free(state->reg_constant);
    }
    if (state->reg_exponent) {
      free(state->reg_exponent);
    }
    if (state->reg_species) {
      free(state->reg_species);
    }
    if (state->reg_drctn) {
      free(state->reg_drctn);
    }
    if (state->vgrng_state) {
      free(state->vgrng_state);
    }
    if (state->vgrng2_state) {
      free(state->vgrng2_state);
    }
    /*
      Free fields allocated in alloc3
    */
    if (state->compartment_ptrs) {
      free(state->compartment_ptrs);
    }
    if (state->current_counts) {
      free(state->current_counts);
    }
    if (state->bndry_flux_counts) {
      free(state->bndry_flux_counts);
    }
    if (state->net_lklhd_bndry_flux) {
      free(state->net_lklhd_bndry_flux);
    }
    if (state->net_likelihood) {
      free(state->net_likelihood);
    }
    if (state->count_to_conc) {
      free(state->count_to_conc);
    }
    if (state->conc_to_count) {
      free(state->conc_to_count);
    }
    if (state->dg0s) {
      free(state->dg0s);
    }
    if (state->ke) {
      free(state->ke);
    }
    if (state->kss) {
      free(state->kss);
    }
    if (state->kss_e_val) {
      free(state->kss_e_val);
    }
    if (state->kss_u_val) {
      free(state->kss_u_val);
    }
    if (state->molecule_dg0tfs) {
      free(state->molecule_dg0tfs);
    }
    if (state->molecule_probabilities) {
      free(state->molecule_probabilities);
    }
    if (state->molecule_chemical_potentials) {
      free(state->molecule_chemical_potentials);
    }
    /*
      Free fields allocated in alloc4
    */
    molecules_matrix = state->molecules_matrix;
    if (molecules_matrix) {
      if (molecules_matrix->molecules_ptrs) {
	free(molecules_matrix->molecules_ptrs);
      }
      if (molecules_matrix->rxn_indices) {
	free(molecules_matrix->rxn_indices);
      }
      /*
      if (molecules_matrix->transpose_work) {
	free(molecules_matrix->transpose_work);
      }
      */
      free(state->molecules_matrix);
    }
    /*
      Free variables allocated in alloc5
    */
    /*
    formation_energies = state->formation_energies;
    if (formation_energies) {
      if (formation_energies->pseudoisomers) {
	free(formation_energies->pseudoisomers);
      }
      if (formation_energies->pseudoisomer_strings) {
	free(formation_energies->pseudoisomer_strings);
      }
      free(state->formation_energies);
    }
    */
    /*
      Free workspace if it existed
    */
    if (state->workspace_base) {
      free(state->workspace_base);
    }
  }
  return(success);
}

