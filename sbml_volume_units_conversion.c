/* sbml_process_list_of_compartments.c
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

#include "sbml_volume_units_conversion.h"
double sbml_volume_units_conversion(char *units, FILE *error_fp) {
  /*
    Compute the conversion factor for units into liters.
    Called by: sbml_process_list_of_compartments
    Calls:     strcmp, fprintf,fflushf
  */
  double multiplier;
  if (strcmp(units,"litre") == 0) {
    multiplier = 1.0;
  } else {
    if (strcmp(units,"liter") == 0) {
      multiplier = 1.0;
    } else {
      if (strcmp(units,"millilitre") == 0) {
	multiplier = 1.0e-3;
      }
      else {
	if (strcmp(units,"milliliter") == 0) {
	  multiplier = 1.0e-3;
	} else {
	  if (strcmp(units,"microlitre") == 0) {
	    multiplier = 1.0e-6;
	  } else {
	    if (strcmp(units,"microliter") == 0) {
	      multiplier = 1.0e-6;
	    } else {
	      fprintf(error_fp,"sbml_volume_units_conversion: Error "
		      "invalid units field was %s\n",units);
	      multiplier = 1.0;
	    }
	  }
	}
      }
    }
  }
  return (multiplier);
}
