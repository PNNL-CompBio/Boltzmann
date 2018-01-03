/* sbml_lookup_reaction_attribute.c
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

#include "sbml_lookup_reaction_attribute.h"

int sbml_lookup_reaction_attribute(char *key) {
/*
  Return position  of  key in following list.
    0 compartment
    1 fast,
    2 id,
    3 reversible
  Returns -1 if key is not in list.

  Called by: sbml_process_list_of_reactions
  Calls:     strcmp
*/
  int cmp0;
  int cmp1;

  int cmp2;
  int cmp3;

  int cmp4;
  int cmp5;

  int tag;
  int padi;

  cmp0 = strcmp(key,"compartment");
  tag = -1;
  if (cmp0 == 0) {
    tag = 0;
  } else {
    cmp1 = strcmp(key,"delta_g0");
    if (cmp1 == 0) {
      tag = 1;
    } else {
      cmp2 = strcmp(key,"delta_g0_units");
      if (cmp2 == 0) {
	tag = 2;
      } else {
	cmp3 = strcmp(key,"fast");
	if (cmp3 == 0) {
	    tag = 3;
	} else {
	  cmp4 = strcmp(key,"id");
	  if (cmp4 == 0) {
	    tag = 4;
	  } else {
	    cmp5 = strcmp(key,"reversible");
	    if (cmp5 == 0) {
	      tag = 5;
	    }
	  }
	}
      }
    }
  }
  return (tag);
}
