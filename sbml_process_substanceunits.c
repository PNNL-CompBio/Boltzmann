/* sbml_process_substanceunits.c
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

#include "sbml_process_substanceunits.h"
void sbml_process_substanceunits(char *units, double *multiplier_p, 
				 FILE *error_fp) {
  /*
    Determine the intiial concentration multiplier based on the
    units.
    Called by: sbml_process_list_of_species
    Calls:     fprintf, fflush, strcmp
  */
  double mole;
  double millimole;
  double micromole;
  double nanomole;
  double picomole;
  double femtomole;
  double multiplier;
  mole = 1.0;
  millimole = 1.e-3;
  micromole = 1.e-6;
  nanomole  = 1.e-9;
  picomole  = 1.e-12;
  femtomole = 1.e-15;
  if (strcmp(units,"mole") == 0)  {
    multiplier = mole;
  } else {
    if (strcmp(units,"millimole") == 0) {
      multiplier = millimole;
    } else {
      if (strcmp(units,"micromole") == 0) {
	multiplier = micromole;
      } else {
	if (strcmp(units,"nanomole") == 0) {
	  multiplier = nanomole;
	} else {
	  if (strcmp(units,"picomole") == 0) {
	    multiplier = picomole;
	  } else {
	    if (strcmp(units,"femtomole") == 0) {
	      multiplier = femtomole;
	    } else {
	      fprintf(error_fp,"sbml_process_substanceunits: "
		      "Error, unrecognized unit for ammount "
		      "was %s\n",units);
	      fflush(error_fp);
	      multiplier = mole;
	    }
	  }
	}
      }
    }
  }
  *multiplier_p = multiplier;
}
