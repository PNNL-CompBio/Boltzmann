/* set_compartment_ptrs.c
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

#include "set_compartment_ptrs.h"
int set_compartment_ptrs(struct state_struct *state) {
  /*
    Set pointers to the first molecule in each compartment
    Called by: species_init, rxn_map_init
    Calls:     
  */
  int64_t *cmpt_ptrs;
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *cur_molecule;
  int nzc;
  int nzm;

  int success;
  int j;

  int fni;
  int cni;

  int i;
  int padi;

  success = 1;
  nzc            = state->nunique_compartments;
  nzm            = state->nunique_molecules;
  cmpt_ptrs      = state->compartment_ptrs;
  sorted_molecules = state->sorted_molecules;
  cur_molecule     = sorted_molecules;
  cmpt_ptrs[0] = 0;
  cni = 0;
  if (nzc < 2) {
    cmpt_ptrs[1] = nzm;
  } else {
    i = 1;
    for (j = 0; j < nzm; j++) {
      fni  = cur_molecule->c_index;
      if (fni != cni) {
	cmpt_ptrs[i] = j;
	i += 1;
	cni = fni;
      }
      cur_molecule += 1; /* Caution address arithmetic */
    }
    cmpt_ptrs[i] = nzm;
  }
  return(success);
}
