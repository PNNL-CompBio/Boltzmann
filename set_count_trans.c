/* set_count_trans.c
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

#include "set_count_trans.h"
int set_count_trans(struct state_struct *state) {
  /*
    Fill the count_to_conc and conc_to_count vectors with
    appropriate multipliers for use in rxn_likelihood and deq stuff.
    Called by: species_init
    Calls:
  */
  struct compartment_struct *compartments;
  struct compartment_struct *compartment;
  struct molecule_struct *cur_molecule;
  double *count_to_conc;
  double *conc_to_count;

  int nzm;
  int cni;
  int success;
  int j;

  success = 1;
  cur_molecule  = state->sorted_molecules;
  compartments  = state->sorted_cmpts;
  nzm           = state->nunique_molecules;
  count_to_conc = state->count_to_conc;
  conc_to_count = state->conc_to_count;

  for (j = 0; j < nzm; j++) {
    cni              = cur_molecule->c_index;
    compartment      = &compartments[cni];
    count_to_conc[j] = compartment->count_to_conc;
    conc_to_count[j] = compartment->conc_to_count;
    cur_molecule     += 1; /* Caution address arithmetic */
  }
  return(success);
}
