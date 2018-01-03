/* sort_global_compartments.c
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
#include "merge_compartments.h"

#include "sort_global_compartments.h"
int sort_global_compartments(struct istring_elem_struct **unsorted_compartments,
			     struct istring_elem_struct **sorted_compartments,
			     int64_t *compartment_map_indices,
			     char *compartment_text,
			     int n) {
  /*
    Sort presorted compartment lists from different reactions files.
    Called by: boltzman_boot
    Calls    : merge_compartments

    n is the number reaction files or sorted compartment lists to be merged.
    compartment_map_indices[i] is the index of the first compartment for
    the i'th reaction file. compartment_map_indices[n] is the global number
    of compartments.
    This is essentially sort_compartments starting with presorted lists of 
    different lengths than the doubling step size.
    for now we use a simple merge sort, with strcmp as the
    comparator function. In the future might want to bin sort
    on the first character, then sort by length and then
    alphabetically for same length entities so as to only
    do strcmp's for strings starting with the same character, and
    of the same length.
    
  
  */
  struct istring_elem_struct *u_compartments;
  struct istring_elem_struct *s_compartments;
  struct istring_elem_struct *temp;
  int global_compartments;
  int success;
  int step;

  int l1;
  int l2;

  int ln;
  int j;
  
  int list1_first;
  int list2_first;

  success = 1;
  u_compartments = *unsorted_compartments;
  s_compartments = *sorted_compartments;
  global_compartments = compartment_map_indices[n];
  if (n == 1) {
    /*
      Only 1 Global compartment nothing to be done.
    */
    *sorted_compartments = u_compartments;
  } else {
    step = 1;
    for (step = 1; step < n; step += step) {
      for(j=0;j<(n-step);j+=(step + step)) {
	list1_first = compartment_map_indices[j];
	list2_first = compartment_map_indices[j+step];
	l1 = list2_first - list1_first;
	l2 = global_compartments - list2_first;
	if ((j + step + step) < n) {
	  l2 = compartment_map_indices[j+step+step] - list2_first;
	}
	merge_compartments((struct istring_elem_struct *)&u_compartments[list1_first],
		       (struct istring_elem_struct *)&u_compartments[list2_first],
			   (struct istring_elem_struct *)&s_compartments[list1_first],
			   compartment_text,l1,l2);
      }
      /* Now if the last group is <= step they just need to be copied
	 into the sorted list.
      */ 
      ln = n & (step + step - 1);
      if (ln <= step) {
	list1_first = compartment_map_indices[n-ln];
	l1 = global_compartments - list1_first;
	for (j = list1_first;j<global_compartments;j++) {
	  s_compartments[j].string   = u_compartments[j].string;
	  s_compartments[j].m_index  = u_compartments[j].m_index;
	  s_compartments[j].c_index  = u_compartments[j].c_index;
	  s_compartments[j].g_index  = u_compartments[j].g_index;
	  s_compartments[j].variable = u_compartments[j].variable;
	}
      }
      temp               = s_compartments;
      s_compartments     = u_compartments;
      u_compartments     = temp;
    }
  }
  *sorted_compartments   = u_compartments;
  *unsorted_compartments = s_compartments;
  return(success);
}
