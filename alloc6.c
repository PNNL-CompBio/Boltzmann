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
int alloc6(int nu_molecules,
	   void **pointers,
	   int64_t *usage,
	   FILE *lfp) {
  /*
    Allocate space for the formation_energies struct fields needed
    after parsing the pseudoisomers file is complete. This space 
    depends on fields from state structure.
    Called by: compute_standard_energies
    Calls:     calloc
  */
  /*
  */
  /*
  int64_t enu_molecules;
  int64_t ask_for;
  int64_t one_l;
  int     *sorted_molecule_order;
  int     *order_scratch;
  */
  int success;
  int padi;
  success = 1;
  /*
    allocate space for ordering molecule names.
    not used for now.
  */
  /*
  one_l   = (int64_t)1;
  enu_molecules = nu_molecules + (nu_molecules & one_l);
  ask_for = enu_molecules * ((int64_t)sizeof(int));
  *usage += ask_for;
  sorted_molecule_order = (int*)calloc(one_l,ask_for);
  if (sorted_molecule_order == NULL) {
    success = 0;
    if (lfp) {
      fprintf(lfp,
              "alloc6: Unable to allocate %ld bytes for "
              "sorted_molecule_order\n",
	      ask_for);
      fflush(lfp);
    }
  }    
  if (success) {
    pointers[0] = (void*)sorted_molecule_order;
    ask_for += ask_for;
    *usage += ask_for;
    order_scratch = (int*)calloc(one_l,ask_for);
    if (order_scratch == NULL) {
      success = 0;
      if (lfp) {
         fprintf(lfp,
  	         "alloc6: Unable to allocate %ld bytes for "
	         "order_scratch\n",
	         ask_for);
         fflush(lfp);
	 }
      }
    }
  }
  if (success) {
     pointers[1] = (void*)order_scratch;
  }
  */
  return(success);
}
