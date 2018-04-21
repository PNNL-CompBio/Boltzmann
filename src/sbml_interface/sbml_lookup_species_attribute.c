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
    6 initialConcentration
    7 name
    8 substanceUnits
  Returns -1 if key is not in list.

  Called by: sbml_parse_species_key_value
  Calls:     strcmp
*/
  int tag;
  int fc;
  int bfc;
  int cfc;
  int hfc;
  int ifc;
  int nfc;
  int sfc;
  tag = -1;
  fc =  (int)key[0];
  bfc = (int)'b';
  cfc = (int)'c';
  hfc = (int)'h';
  ifc = (int)'i';
  nfc = (int)'n';
  sfc = (int)'s';
  if (fc < ifc) {
    if (fc <= bfc) {
      if (strcmp(key,"bondaryCondition") == 0) {
	tag = 0;
      }
    } else {
      if (fc == cfc) {
	if (strcmp(key,"compartment") == 0) {
	  tag = 1;
	} else {
	  if (strcmp(key,"constant") == 0) {
	    tag = 2;
	  }
	}
      } else {
	if (fc == hfc) {
	  if (strcmp(key,"hasOnlySubstanceUnits") == 0) {
	    tag = 3;
	  }
	}
      }
    }
  } else {
    if (fc >= sfc) {
      if (strcmp(key,"substanceUnits") == 0) {
	tag = 8;
      } 
    } else {
      if (fc == ifc) {
	if (strcmp(key,"id") == 0) {
	  tag = 4;
	} else {
	  if (strcmp(key,"initialAmount") == 0) {
	    tag = 5;
	  } else {
	    if (strcmp(key,"initialConcentration") == 0) {
		tag = 6;
	    }
	  }
	}
      } else {
	if (fc == nfc) {
	  if (strcmp(key,"name") == 0) {
	    tag = 7;
	  }
	}
      }
    }
  }
  return(tag);
}

