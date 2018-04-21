/* sbml_merge_species_trans.c
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

#include "sbml_merge_species_trans.h"
void sbml_merge_species_trans(int64_t l1, char **species1, char **trans1,
			      int64_t l2, char **species2, char **trans2,
			      char **species3, char **trans3) {
  /*
    Merge two sorted list species and coresponding translations.
    Called by: sbml_sort_species_trans
    Calls:     strcmp
  */
  int64_t j;
  int64_t j1;
  int64_t j2;
  int64_t j3;
  int64_t n;
  j  = (int64_t)0;
  j1 = j;
  j2 = j;
  j3 = j;
  n  = l1 + l2;
  for (j3=(int64_t)0;j3 < n;j3++) {
    if (strcmp(species1[j1],species2[j2]) <= 0) {
      species3[j3] = species1[j1];
      trans3[j3] = trans1[j1];
      j1++;
      if (j1 == l1) {
	for (j = j2;j<l2;j++) {
	  j3++;
	  species3[j3] = species2[j];
	  trans3[j3] = trans2[j];
	}
	break;
      }
    } else {
      species3[j3] = species2[j2];
      trans3[j3] = trans2[j2];
      j2 ++;
      if (j2 == l2) {
	for (j = j1;j<l1;j++) {
	  j3++;
	  species3[j3] = species1[j];
	  trans3[j3] = trans1[j];
	}
	break;
      }
    }
  }
}
