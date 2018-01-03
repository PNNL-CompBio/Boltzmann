/* sbml_lookup_species_attribute.c
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

#include "sbml_lookup_species_attribute.h"

int sbml_lookup_species_attribute(char *key) {
/*
  Return position  of  key in following list.
    0 bondaryCondition
    1 compartment
    2 constant
    3 hasOnlySubstanceUnits
    4 id
    5 initialAmount
    6 substanceUnits
  Returns -1 if key is not in list.

  Called by: sbml_process_list_of_species
  Calls:     strcmp
*/
  int cmp0;
  int cmp1;
  int cmp2;
  int cmp3;
  int cmp4;
  int cmp5;
  int cmp6;
  int tag;
  cmp0 = strcmp(key,"bondaryCondition");
  cmp6 = strcmp(key,"substanceUnits");
  tag = -1;
  if (cmp0 == 0) {
    tag = 0;
  } else {
    if (cmp6 == 0) {
      tag = 6;
    } else {
      cmp3 = strcmp(key,"hasOnlySubstanceUnits");
      if (cmp3 == 0) {
	tag = 3;
      } else {
	if (cmp3 < 0) {
	  cmp1 = strcmp(key,"compartment");
	  if (cmp1 == 0) {
	    tag = 1;
	  } else {
	    cmp2 = strcmp(key,"constant");
	    if (cmp2 == 0) {
	      tag = 2;
	    }
	  }
	} else {
	  cmp4 = strcmp(key,"id");
	  if (cmp4 == 0) {
	    tag = 4;
	  } else {
	    cmp5 = strcmp(key,"initialAmount");
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
