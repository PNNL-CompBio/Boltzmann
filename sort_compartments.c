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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"
#include "merge_compartments.h"

#include "sort_compartments.h"
int sort_compartments(struct istring_elem_struct **unsorted_compartments,
		      struct istring_elem_struct **sorted_compartments,
		      int n) {
  /*
    Sort the uppercase istring names.
    Called by: boltzman_init
    Calls    : merge_compartments

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

  int success;
  int step;

  int l1;
  int l2;

  int ln;
  int j;
  success = 1;
  u_compartments = *unsorted_compartments;
  s_compartments = *sorted_compartments;
  
  if (n > 2) {
    for (step = 1; step < n; step += step) {
      for(j=0;j<(n-step);j = j + step + step) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) l2 = step;
	merge_compartments((struct istring_elem_struct *)&u_compartments[j],
		       (struct istring_elem_struct *)&u_compartments[j+step],
	               (struct istring_elem_struct *)&s_compartments[j],l1,l2);
      }
      /* Now if the last group is <= step they just need to be copied
	 in to the sorted list.
      */ 
      ln = n & (step + step - 1);
      if (ln <= step) {
	for (j = n-ln;j<n;j++) {
	  s_compartments[j].string = u_compartments[j].string;
	  s_compartments[j].m_index  = u_compartments[j].m_index;
	  s_compartments[j].c_index  = u_compartments[j].c_index;
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
