/* merge_sorted_strings.c
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

#include "merge_sorted_strings.h"
void merge_sorted_strings(int64_t l1, char **list1,
			  int64_t l2, char **list2,
			  char **list3) {
  /*
    Merge two sorted lists of strings.
    Called by: sort_json_ids
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
    if (strcmp(list1[j1],list2[j2]) <= 0) {
      list3[j3] = list1[j1];
      j1++;
      if (j1 == l1) {
	for (j = j2;j<l2;j++) {
	  j3++;
	  list3[j3] = list2[j];
	}
	break;
      }
    } else {
      list3[j3] = list2[j2];
      j2 ++;
      if (j2 == l2) {
	for (j = j1;j<l1;j++) {
	  j3++;
	  list3[j3] = list1[j];
	}
	break;
      }
    }
  }
}
