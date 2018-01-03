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

#include "species_lookup.h"

#include "read_initial_concentrations.h"
int read_initial_concentrations(struct state_struct *state) {
  /*
    Read the concs.in file for initial concentrations, and
    set the concentrations array.
    Called by: boltzmann
    Calls:     species_lookup
  */
  double  conc;
  double *concs;
  int64_t species_buff_len;
  int success;
  int nzr;

  int nu_species;
  int i;

  int nscan;
  int ws;

  int nws;
  int ll;

  int si;
  int padi;
  char *species_buffer;
  char *species_name;
  char *fgp;
  
  FILE *conc_fp;
  nu_species = state->unique_species;
  species_buff_len = state->rxn_buff_len;
  species_buffer   = state->rxn_buffer;
  species_name     = state->species_name;
  concs            = state->concentrations;
  success = 1;
  for (i=0;i<nu_species;i++) {
    concs[i] = -1.0;
  }
  conc_fp = fopen("concs.in","r");
  if (conc_fp) {
    while (!feof(conc_fp)) {
      fgp = fgets(species_buffer,species_buff_len,conc_fp);
      if (fgp) {
	ll = strlen(species_buffer);
	nscan = sscanf(species_buffer,"%s %le",species_name,&conc);
	if (nscan == 2) {
	  if (strcmp(species_name,"*") == 0) {
	    state->default_initial_conc = conc;
	  } else {
	    si = species_lookup(species_name,state);
	    if ((si >=0) && si < nu_species) {
	      concs[si] = conc;
	    } else {
	      fprintf(stderr,"read_initial_concentrations: Error "
		      "unrecognized species in conc.in was %s\n",
		      species_name);
	      fflush(stderr);
	    }
	  }
	}
      }
    }
  } else {
    fprintf(stderr,
	    "read_initial_concentrations: Warning unable to open concs.in ."
	    " Will use %le for all species\n",state->default_initial_conc);
    fflush(stderr);
  }
  if (success) {
    conc = state->default_initial_conc;
    for (i=0;i<nu_species;i++) {
      if (concs[i] < 0.0) {
	concs[i] = conc;
      }
    }
  }
  return(success);
}
