/* compute_standard_energies.c
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

#include "size_pseudoisomer_file.h"
#include "alloc5.h"
#include "parse_pseudoisomer_dg0f_file.h"
#include "alloc6.h"
#include "compute_molecule_dg0tfs.h"
#include "compute_molecular_partition_probability.h"
#include "compute_chemical_potential.h"
#include "compute_reaction_dg0.h"
#include "unalloc6.h"

#include "compute_standard_energies.h"
int compute_standard_energies(struct state_struct *state,
		       struct formation_energy_struct **formation_energies_p) {
  struct formation_energy_struct *formation_energies;
  /*
    Called by: boltzmann_init, boltzmann_boot
    Calls size_pseudoisomer_file,
          alloc5,
	  parse_pseudoisomer_dg0f_file
	  alloc6
	  compute_molecule_dg0tfs
	  compute_molecular_partiton_probability
	  compute_chemical_potential
	  compute_reaction_dg0
	  unalloc6
  */
  int64_t num_pseudoisomers;
  int64_t length_pseudoisomer_strings;
  int success;
  success = 1;
  if (*formation_energies_p == NULL) {
    /*
      Formation_energies_p is null means the pseudoisomer data base has
      not yet been initialized.
      First determine the size of the pseudoisomers file.
    */
    success = size_pseudoisomer_file(state,
				     &num_pseudoisomers,
				     &length_pseudoisomer_strings);
    if (success) {
      /*
	Now allocate the formation_energies struct and its fields 
	needed to parse and sort the pseudoisomer file.
      */
      success = alloc5(num_pseudoisomers,
		       length_pseudoisomer_strings,
		       state->align_len,
		       &state->usage,
		       &formation_energies);
    }
    if (success) {
      /*
	Read in formation energies from pseudoisomer file.
      */
      *formation_energies_p = formation_energies;
      formation_energies->num_pseudoisomers = num_pseudoisomers;
      formation_energies->length_pseudoisomer_strings = 
	length_pseudoisomer_strings;
      success = parse_pseudoisomer_dg0f_file(formation_energies->pseudoisomers,
					     formation_energies->pseudoisomer_strings,
					     state->pseudoisomer_file,
					     num_pseudoisomers,
					     state->align_len);
    } 
    if (success) {
      /*
	Sort the pseudoisomer data by json_cpd_name, priority and dg0f value.
      success = sort_pseudoisomers(formation_energies);
      */
    }
    /*
    if (success  && formation_energies->pseudoisomer_scratch) {
      free(formation_energies->pseudoisomer_scratch);
    }
    */
  } else {
    formation_energies = *formation_energies_p;
  }
  if (success) {
    /*
      Now we need to allocate space for computing the reaction dg0fs.
      First we need to fill the relevant fields of the formation_energies
      struct from the state struct.
    */
    formation_energies->reactions           = state->reactions;
    formation_energies->reactions_matrix    = state->reactions_matrix;
    formation_energies->sorted_molecules    = state->sorted_molecules;
    formation_energies->current_counts      = state->current_counts;
    formation_energies->number_reactions    = state->number_reactions;
    formation_energies->nunique_molecules   = state->nunique_molecules;
    formation_energies->align_len           = state->align_len;
    formation_energies->align_mask          = state->align_mask;
    formation_energies->max_filename_len    = state->max_filename_len;
    formation_energies->use_pseudoisomers   = state->use_pseudoisomers;
    formation_energies->print_output        = state->print_output;


    formation_energies->ideal_gas_r         = state->ideal_gas_r;
    formation_energies->temp_kelvin         = state->temp_kelvin;
    formation_energies->rt                  = state->rt;
    formation_energies->m_r_rt              = state->m_r_rt;
    formation_energies->m_rt                = state->m_rt;
    formation_energies->cals_per_joule      = state->cals_per_joule;
    formation_energies->joules_per_cal      = state->joules_per_cal;
    formation_energies->ph                  = state->ph;
    formation_energies->ionic_strength      = state->ionic_strength;
    formation_energies->molecules_text      = state->molecules_text;
    formation_energies->molecule_dg0tfs     = state->molecule_dg0tfs;
    formation_energies->molecule_probabilities = state->molecule_probabilities;
    formation_energies->molecule_chemical_potentials = 
      state->molecule_chemical_potentials;

    formation_energies->log_fp              = state->lfp;

    success = alloc6(formation_energies,&state->usage);
  }
  if (success) {
    success = compute_molecule_dg0tfs(formation_energies);
  }
  if (success) {
    success = compute_molecular_partition_probability(formation_energies);
  }
  if (success) {
    success = compute_chemical_potential(formation_energies);
  } 
  if (success) {
    success = compute_reaction_dg0(formation_energies);
  }
  if (success) {
    success = unalloc6(formation_energies);
  }
  return(success);
}
