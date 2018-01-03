/* sort_global_molecules.c
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
#include "merge_molecules.h"

#include "sort_global_molecules.h"
int sort_global_molecules(struct istring_elem_struct **unsorted_molecules,
			     struct istring_elem_struct **sorted_molecules,
			     int64_t *molecule_map_indices,
			     char *molecule_text,
			     int n) {
  /*
    Sort presorted molecule lists from different reactions files.
    Called by: boltzman_boot
    Calls    : merge_molecules

    n is the number reaction files or sorted molecule lists to be merged.
    molecule_map_indices[i] is the index of the first molecule for
    the i'th reaction file. molecule_map_indices[n] is the global number
    of molecules.
    This is essentially sort_molecules starting with presorted lists of 
    different lengths than the doubling step size.
    for now we use a simple merge sort, with strcmp as the
    comparator function. In the future might want to bin sort
    on the first character, then sort by length and then
    alphabetically for same length entities so as to only
    do strcmp's for strings starting with the same character, and
    of the same length.
    
  
  */
  struct istring_elem_struct *u_molecules;
  struct istring_elem_struct *s_molecules;
  struct istring_elem_struct *temp;
  int global_molecules;
  int success;
  int step;

  int l1;
  int l2;

  int ln;
  int j;

  int list1_first;
  int list2_first;
  success = 1;
  u_molecules = *unsorted_molecules;
  s_molecules = *sorted_molecules;
  global_molecules = molecule_map_indices[n];
  if (n == 1) {
    /*
      Only 1 Global molecule nothing to be done.
    */
    *sorted_molecules = u_molecules;
  } else {
    step = 1;
    for (step = 1; step < n; step += step) {
      for(j=0;j<(n-step);j+=(step + step)) {
	list1_first = molecule_map_indices[j];
	list2_first = molecule_map_indices[j+step];
	l1 = list2_first - list1_first;
	l2 = global_molecules - list2_first;
	if ((j + step + step) < n) {
	  l2 = molecule_map_indices[j+step+step] - list2_first;
	}
	merge_molecules((struct istring_elem_struct *)&u_molecules[list1_first],
		       (struct istring_elem_struct *)&u_molecules[list2_first],
			   (struct istring_elem_struct *)&s_molecules[list1_first],
			   molecule_text,l1,l2);
      }
      /* Now if the last group is <= step they just need to be copied
	 into the sorted list.
      */ 
      ln = n & (step + step - 1);
      if (ln <= step) {
	list1_first = molecule_map_indices[n-ln];
	l1 = global_molecules - list1_first;
	for (j = list1_first;j<global_molecules;j++) {
	  s_molecules[j].string   = u_molecules[j].string;
	  s_molecules[j].m_index  = u_molecules[j].m_index;
	  s_molecules[j].c_index  = u_molecules[j].c_index;
	  s_molecules[j].g_index  = u_molecules[j].g_index;
	  s_molecules[j].variable = u_molecules[j].variable;
	}
      }
      temp               = s_molecules;
      s_molecules     = u_molecules;
      u_molecules     = temp;
    }
  }
  *sorted_molecules   = u_molecules;
  *unsorted_molecules = s_molecules;
  return(success);
}
