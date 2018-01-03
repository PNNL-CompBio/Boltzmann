/* check_initial_concentrations.c
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

#include "check_initial_concentrations.h"
int check_initial_concentrations(struct state_struct *state) {
  /*
    Check that all molecules in the rxns.dat file had entries in the
    intial concentrations file. Code is mostly lifteed from
    print_dictionary.c
    Called by: species_init
    Calls:     fprintf,fflush (intrinsic)
  */
  struct molecule_struct *molecules;
  struct compartment_struct *compartments;
  struct compartment_struct *cur_cmpt;
  char *molecules_text;
  char *compartment_text;
  char *molecule_string;
  char *cmpt_string;

  int i;
  int oi;

  int ci;
  int success;

  int nunique_molecules;
  int padi;
  FILE *lfp;
  
  success = 1;
  lfp = state->lfp;
  nunique_molecules = state->nunique_molecules;
  molecules         = state->sorted_molecules;
  compartments      = state->sorted_cmpts;
  molecules_text    = state->molecules_text;
  compartment_text  = state->compartment_text;
  oi          = -1;
  for (i=0;i<nunique_molecules;i++) {
    ci = molecules->c_index;
    molecule_string    = (char *)&molecules_text[molecules->string];
    if (ci != oi) {
      oi = ci;
      if (ci > 0) {
	cur_cmpt = (struct compartment_struct *)&(compartments[ci]);
	cmpt_string = (char *)&compartment_text[cur_cmpt->string];
      }
    }
    if (molecules->solvent == 0) {
      if (molecules->variable == -1) {
	/*
	  We have a molecule that did not get an intial concentration
	  in read_intial_concentrations and is not a solvent molecule.
	*/
	if (lfp) {
	  if (ci > 0) {
	    fprintf(lfp,"Error molecule %d %s:%s did not have an intial "
		    "concentration.\n",i,molecule_string,cmpt_string);
	    fflush(lfp);
	  } else {
	    fprintf(lfp,"Error molecule %d %s did not have an intial "
		    "concentration.\n",i,molecule_string);
	    fflush(lfp);
	  }
	}
	success = 0;
      }
    }
    molecules += 1; /* Caution address arithmetic. */
  }
  return(success);
}
