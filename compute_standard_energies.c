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
/*
#include "alloc6.h"
*/
#include "compute_molecule_dg0tfs.h"
/*
#include "compute_molecular_partition_probability.h"
#include "compute_chemical_potential.h"
*/
#include "compute_reaction_dg0.h"
#include "unalloc6.h"


#include "compute_standard_energies.h"
int compute_standard_energies(struct state_struct *state) {
  /*
    Called by: enery_init
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
  struct pseudoisomer_struct *pseudoisomers;
  char *pseudoisomer_strings;
  int64_t num_pseudoisomers;
  int64_t length_pseudoisomer_strings;
  void   *pointers[4];
  int success;
  int padi;
  /*
  int nu_molecules;
  FILE *lfp;
  FILE *efp;
  lfp = state->lfp;
  nu_molecules = state->nunique_molecules;
  */
  /*
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
		     pointers);
  }
  if (success) {
    /*
      Read in formation energies from pseudoisomer file.
    */
    pseudoisomers = (struct pseudoisomer_struct *)pointers[0];
    pseudoisomer_strings = (char *)pointers[1];
    success = parse_pseudoisomer_dg0f_file(pseudoisomers,
					   pseudoisomer_strings,
					   state->pseudoisomer_file,
					   num_pseudoisomers,
					   state->align_len);
  } 
  /*
  if (success) {
    success = alloc6(nu_molecules,
                     (void**)&pointers[2],&state->usage,lfp);
  }
  if (success) {
    //Sort the pseudoisomer data by json_cpd_name, priority and dg0f value.
    success = sort_pseudoisomers(num_pseudoisomers,pointers);
  }
  */
  if (success) {
    success = compute_molecule_dg0tfs(state,
				      pseudoisomers,
				      pseudoisomer_strings,
				      num_pseudoisomers);
  }
  if (success) {
    success = unalloc6(2,pointers);
    /*
    success = unalloc6(2,(void**)&pointers[2]);
    */
  }
  /*
    // These two calls are coupled as compute_chemical_potential
    // needs the molecule_probabilities.
    // However molecule_probabilities and molecule_chemical_potentials
    // are not currently used anywhwere else so we comment these out.
  if (success) {
    success = compute_molecular_partition_probability(state);
  }
  if (success) {
    success = compute_chemical_potential(state);
  } 
  */
  if (success) {
    success = compute_reaction_dg0(state);
  }
  return(success);
}
