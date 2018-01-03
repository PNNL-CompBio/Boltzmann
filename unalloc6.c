/* unalloc6.c
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
int unalloc6(struct formation_energy_struct *formation_energies) {
  /*
    free space allocated by alloc6 routine.
    Called by: compute_standard_energies
    Calls:     free
  */
  int64_t ask_for;
  int64_t nu_molecules;
  int64_t enu_molecules;
  int64_t one_l;
  int success;
  success = 1;
  one_l   = (int64_t)1;
  /*
  free(formation_energies->molecule_dg0tfs);
  free(formation_energies->molecule_probabilities);
  free(formation_energies->molecule_chemical_potentials);
  free(formation_energies->sorted_molecule_order);
  free(formation_energies->order_scratch);
  */
  return (success);
}
