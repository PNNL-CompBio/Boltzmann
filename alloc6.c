/* alloc6.c
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

#include "alloc6.h"
int alloc6(struct formation_energy_struct *formation_energies,
	   int64_t *usage) {
  /*
    Allocate space for the formation_energies struct fields needed
    after parsing the pseudoisomers file is complete. This space 
    depends on fields from state structure.
    Called by: formation_energy_rxn_dg0fs
    Calls:     calloc
  */
  int64_t ask_for;
  int64_t nu_molecules;
  int64_t enu_molecules;
  int64_t one_l;
  int success;
  success = 1;
  one_l   = (int64_t)1;
  nu_molecules = formation_energies->nunique_molecules;
  /*
    allocate space for ordering molecule names.
    not used for now.
  */
  /*
  if (success) {
    enu_molecules = nu_molecules + (nu_molecules & one_l);
    ask_for = enu_molecules * ((int64_t)sizeof(int));
    *usage += ask_for;
    formation_energies->sorted_molecule_order = (int*)calloc(one_l,ask_for);
    if (formation_energies->sorted_molecule_order == NULL) {
      fprintf(stderr,
	      "alloc6: Unable to allocate %ld bytes for "
	      "sorted_molecule_order\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }    
  if (success) {
    ask_for += ask_for;
    *usage += ask_for;
    formation_energies->order_scratch = (int*)calloc(one_l,ask_for);
    if (formation_energies->order_scratch == NULL) {
      fprintf(stderr,
	      "alloc6: Unable to allocate %ld bytes for "
	      "order_scratch\n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  */
  return(success);
}
