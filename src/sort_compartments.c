/* sort_compartments.c
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

#include "sort_compartments.h"
int sort_compartments(struct compartment_struct *unsorted_compartments,
		      struct compartment_struct *sorted_compartments,
		      char *compartment_text,
		      int n) {
  /*
    Sort the uppercase compartment names.
    Called by: rxns_init, rxn_map_init, sbml_process_list_of_compartments
    Calls    : merge_compartments, memmove

    for now we use a simple merge sort, with strcmp as the
    comparator function. In the future might want to bin sort
    on the first character, then sort by length and then
    alphabetically for same length entities so as to only
    do strcmp's for strings starting with the same character, and
    of the same length.
  
  */
  struct compartment_struct *u_compartments;
  struct compartment_struct *s_compartments;
  struct compartment_struct *temp;
  struct compartment_struct ies;
  int64_t move_size;
  int64_t e_size;

  int success;
  int step;

  int l1;
  int l2;

  int ln;
  int j;
  success = 1;
  u_compartments = unsorted_compartments;
  s_compartments = sorted_compartments;
  e_size = (int64_t) sizeof(ies);
  if (n > 2) {
    for (step = 1; step < n; step += step) {
      for(j=0;j<(n-step);j = j + step + step) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) l2 = step;
	merge_compartments((struct compartment_struct *)&u_compartments[j],
		       (struct compartment_struct *)&u_compartments[j+step],
			   (struct compartment_struct *)&s_compartments[j],
			   compartment_text,l1,l2);
      }
      /* Now if the last group is <= step they just need to be copied
	 into the sorted list.
      */ 
      ln = n & (step + step - 1);
      if (ln <= step) {
	if (ln > 0) {
	  move_size = ((int64_t)ln) * e_size;
	  j = n-ln;
	  memcpy((void*)&s_compartments[j],(void*)&u_compartments[j],move_size);
	}
	/*
	for (j = n-ln;j<n;j++) {
	  s_compartments[j].volume   = u_compartments[j].volume;
	  s_compartments[j].recip_volume   = u_compartments[j].recip_volume;
	  s_compartments[j].ntotal_exp     = u_compartments[j].ntotal_exp;
	  s_compartments[j].ntotal_opt     = u_compartments[j].ntotal_opt;
	  s_compartments[j].conc_to_count  = u_compartments[j].conc_to_count;
	  s_compartments[j].count_to_conc  = u_compartments[j].count_to_conc;
	  s_compartments[j].string   = u_compartments[j].string;
	  s_compartments[j].c_index  = u_compartments[j].c_index;
	  s_compartments[j].g_index  = u_compartments[j].g_index;
	}
	*/
      }
      temp               = s_compartments;
      s_compartments     = u_compartments;
      u_compartments     = temp;
    }
  }
  if (u_compartments != sorted_compartments) {
    move_size = ((int64_t)n) * e_size;
    memcpy((void*)s_compartments,(void*)u_compartments,move_size);
  }
  /*
  *sorted_compartments   = u_compartments;
  *unsorted_compartments = s_compartments;
  */
  return(success);
}
