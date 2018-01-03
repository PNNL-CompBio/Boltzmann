/* sbml_find_string.c
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

#include "sbml_find_string.h"
int sbml_find_string (char* string, char **strings, int num_strings) {
  /*
    Return the index of string in the strings array of length num_strings
    strings array is assumed to be sorted in strcmp order.
    If the string is not found, -1 is returned.
    Called by: sbml_generate_init_conc_line, 
               sbml_process_species_reference_tag
  */
  int index;
  int i;
  int left;
  int right;
  int mid;
  int crslt;
  char *last;
  char *first;
  first = strings[0];
  last = strings[num_strings-1];
  index = -1;
  crslt = strcmp(string,first);
  if (crslt == 0) {
    index = 0;
  } else {
    if (crslt > 0) {
      left = 0;
      crslt = strcmp(string,last);
      if (crslt == 0) {
	index = num_strings-1;
      } else {
	if (crslt < 0) {
	  right = num_strings -1;
	  mid = (left + right) >> 1;
	  while (mid != left) {
	    crslt = strcmp(string,strings[mid]);
	    if (crslt == 0) {
	      index = mid;
	      break;
	    } else {
	      if (crslt < 0) {
		right = mid;
	      } else {
		left = mid;
	      }
	    }
	    mid = (left + right) >> 1;
	  }
	}
      }
    }
  }
  return(index);
}
