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
#include "compartment_lookup.h"
#include "upcase.h"

#include "read_initial_concentrations.h"
int read_initial_concentrations(struct state_struct *state) {
  /*
    Read the concs.in file for initial concentrations, and
    set the concentrations array.
    Called by: boltzmann_init
    Calls:     molecules_lookup
               compartment_lookup,
	       upcase,
               fopen, fgets, fclose, fprintf, fflush (intrinsic)
  */
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *molecule;
  double  conc;
  double *concs;
  double *bndry_flux_concs;
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

  int mol_len;
  int cmpt_len;

  int num_fixed_concs;
  int padi;
  
  char *molecules_buffer;
  char *molecule_name;
  char *compartment_name;
  char *variable_c;
  char vc[2];
  char *fgp;
  
  FILE *conc_fp;
  nu_molecules       = state->unique_molecules;
  molecules_buff_len = state->rxn_buff_len;
  molecules_buffer   = state->rxn_buffer;
  molecule_name      = state->molecule_name;
  compartment_name   = state->compartment_name;
  sorted_molecules   = state->sorted_molecules;
  concs              = state->current_concentrations;
  bndry_flux_concs   = (double *)state->bndry_flux_concs;
  success = 1;
  one_l = (int64_t)1;
  variable_c = (char *)&vc[0];
  for (i=0;i<nu_molecules;i++) {
    concs[i] = -1.0;
    bndry_flux_concs[i] = -1;
  }
  num_fixed_concs = 0;
  conc_fp = fopen(state->init_conc_file,"r");
  if (conc_fp) {
    while (!feof(conc_fp)) {
      fgp = fgets(molecules_buffer,molecules_buff_len,conc_fp);
      if (fgp) {
	nscan = sscanf(molecules_buffer,"%s %le %1s",
		       molecule_name, &conc,variable_c);
	variable = 1;
	if (nscan == 3) {
	  /*
	    A variable or constant specifier was givven.
	  */
	  vc[0] = vc[0] & 95;
	  if (strncmp(variable_c,"C",one_l) == 0) {
	    variable = 0;
	    num_fixed_concs += 1;
	  }
	}
	mol_len = strlen(molecule_name);
	compartment_name = molecule_name;
	for (i=0;i<mol_len-1;i++) {
	  if (molecule_name[i] == ':') {
	    compartment_name = (char *)&molecule_name[i+1];
	    molecule_name[i] = '\0';
	    cmpt_len = mol_len - i - 1;
	    mol_len = i;
	    break;
	  }
	}
	if (compartment_name == molecule_name) {
	  ci = -1;
	} else {
	  upcase(cmpt_len,compartment_name,compartment_name);
	  ci = compartment_lookup(compartment_name,state);
	}
	if (nscan >= 2) {
	  upcase(mol_len,molecule_name,molecule_name);
	  si = molecules_lookup(molecule_name,ci,state);
	  if ((si >=0) && si < nu_molecules) {
	    concs[si] = conc;
	    molecule = (struct istring_elem_struct *)&sorted_molecules[si];
	    molecule->variable = variable;
	  } else {
	    fprintf(stderr,"read_initial_concentrations: Error "
		    "unrecognized molecule in conc.in was %s\n",
		    molecule_name);
	    fflush(stderr);
	    success = 0;
	    break;
	  }
	} else {
	  fprintf(stderr,"read_initial_concentrations: Error "
		  "poorly formated line was\n%s\n",molecules_buffer);
	  fflush(stderr);
	  success = 0;
	}
      } /* end if (fgp) */
    } /* end while (!feof(conc_fp)) */
    fclose(conc_fp);
  } else {
    fprintf(stderr,
	    "read_initial_concentrations: Warning unable to open %s\n",
	    state->init_conc_file);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    conc = state->default_initial_conc;
    for (i=0;i<nu_molecules;i++) {
      if (concs[i] < 0.0) {
	concs[i] = conc;
      }
    }
    for (i=0;i<nu_molecules;i++) {
      bndry_flux_concs[i] = concs[i];
    }
    state->num_fixed_concs = num_fixed_concs;
    /*
      Print the initial concentrations to the concentrations output file.
    */
    if (state->print_output) {
      if (state->concs_out_fp) {
	fprintf(state->concs_out_fp,"init");
	for (i=0;i<nu_molecules;i++) {
	  fprintf(state->concs_out_fp,"\t%le",concs[i]);
	}
	fprintf(state->concs_out_fp,"\n");
      } else {
	fprintf(stderr,
		"read_initial_concentrations: Error concs_out_fp not open\n");
	fflush(stderr);
	success = 0;
      }
    }
  }
  return(success);
}
