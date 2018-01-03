/* read_initial_concentrations.c
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

#include "boltzmann_structs.h"

#include "molecules_lookup.h"

#include "read_initial_concentrations.h"
int read_initial_concentrations(struct state_struct *state) {
  /*
    Read the concs.in file for initial concentrations, and
    set the concentrations array.
    Called by: boltzmann_init
    Calls:     molecules_lookup
               fopen, fgets, fclose, fprintf, fflush (intrinsic)
  */
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *molecule;
  double  conc;
  double *concs;
  int64_t molecules_buff_len;
  int64_t one_l;
  int success;
  int nzr;

  int nu_molecules;
  int i;

  int nscan;
  int variable;

  int si;
  int ci;

  
  char *molecules_buffer;
  char *molecule_name;
  char *compartment_name;
  char *variable_c;
  char vc[2];
  char *fgp;
  
  FILE *conc_fp;
  nu_molecules = state->unique_molecules;
  molecules_buff_len = state->rxn_buff_len;
  molecules_buffer   = state->rxn_buffer;
  molecule_name      = state->molecule_name;
  compartment_name   = state->compartment_name;
  sorted_molecules  = state->sorted_molecules;
  concs            = state->current_concentrations;
  success = 1;
  one_l = (int64_t)1;
  variable_c = (char *)&vc[0];
  for (i=0;i<nu_molecules;i++) {
    concs[i] = -1.0;
  }
  conc_fp = fopen(state->init_conc_file,"r");
  if (conc_fp) {
    while (!feof(conc_fp)) {
      fgp = fgets(molecules_buffer,molecules_buff_len,conc_fp);
      if (fgp) {
	nscan = sscanf(molecules_buffer,"%s:%s %le %1s",
		       molecule_name, compartment_name, &conc,variable_c);
	if (nscan >= 3) {
	  /*
	    A compartment was specified.
	  */
	  ci = compartment_lookup(compartment_name,state);
	  variable = 1;
	  if (nscan == 4) {
	    vc[0] = vc[0] & 95;
	    if (strncmp(variable_c,"C",one_l) == 0) {
	      variable = 0;
	    }
	  }
	} else {
	  ci = -1;
	  nscan = sscanf(molecules_buffer,"%s %le %1s",molecule_name,
			 &conc,variable_c);
	  variable = 1;
	  if (nscan == 3) {
	    vc[0] = vc[0] & 95;
	    if (strncmp(variable_c,"C",one_l) == 0) {
	      variable = 0;
	    }
	  } else {
	    if (nscan < 2) {
	      fprintf(stderr,"read_initial_concentrations: Error "
		      "poorly formated line was\n%s\n",molecules_buffer);
	      fflush(stderr);
	      success = 0;
	    }
	  }
	}
	if (success) {
	  if (strcmp(molecule_name,"*") == 0) {
	    state->default_initial_conc = conc;
	  } else {
	    si = molecules_lookup(molecule_name,ci,state);
	    if ((si >=0) && si < nu_molecules) {
	      concs[si] = conc;
	      molecule = (struct istring_elem_struct *)&sorted_molecules[si];
	      molecule->variable = variable;
	    } else {
	      fprintf(stderr,"read_initial_concentrations: Error "
		      "unrecognized molecules in conc.in was %s\n",
		      molecule_name);
	      fflush(stderr);
	    }
	  }
	}
      }
    }
    fclose(conc_fp);
  } else {
    fprintf(stderr,
	    "read_initial_concentrations: Warning unable to open concs.in ."
	    " Will use %le for all molecules\n",state->default_initial_conc);
    fflush(stderr);
  }
  if (success) {
    conc = state->default_initial_conc;
    for (i=0;i<nu_molecules;i++) {
      if (concs[i] < 0.0) {
	concs[i] = conc;
      }
    }
  }
  return(success);
}
