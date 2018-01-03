/* sort_molecules.c
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
#include "merge_molecules.h"

#include "sort_molecules.h"
int sort_molecules(struct istring_elem_struct **unsorted_istrings,
		   struct istring_elem_struct **sorted_istrings,
		   int n) {
  /*
    Sort the uppercase istring names.
    Called by: boltzman_init
    Calls    : merge_istrings

    for now we use a simple merge sort, with strcmp as the
    comparator function. In the future might want to bin sort
    on the first character, then sort by length and then
    alphabetically for same length entities so as to only
    do strcmp's for strings starting with the same character, and
    of the same length.
  
  */
  struct istring_elem_struct *u_istrings;
  struct istring_elem_struct *s_istrings;
  struct istring_elem_struct *temp;

  int success;
  int step;

  int l1;
  int l2;

  int ln;
  int j;
  success = 1;
  u_istrings = *unsorted_istrings;
  s_istrings = *sorted_istrings;
  
  if (n > 2) {
    for (step = 1; step < n; step += step) {
      for(j=0;j<(n-step);j = j + step + step) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) l2 = step;
	merge_molecules((struct istring_elem_struct *)&u_istrings[j],
			(struct istring_elem_struct *)&u_istrings[j+step],
			(struct istring_elem_struct *)&s_istrings[j],l1,l2);
      }
      /* Now if the last group is <= step they just need to be copied
	 in to the sorted list.
      */ 
      ln = n & (step + step - 1);
      if (ln <= step) {
	for (j = n-ln;j<n;j++) {
	  s_istrings[j].string = u_istrings[j].string;
	  s_istrings[j].m_index  = u_istrings[j].m_index;
	  s_istrings[j].c_index  = u_istrings[j].c_index;
	}
      }
      temp                  = s_istrings;
      s_istrings            = u_istrings;
      u_istrings            = temp;
    }
  }
  *sorted_istrings   = u_istrings;
  *unsorted_istrings = s_istrings;
  return(success);
}
