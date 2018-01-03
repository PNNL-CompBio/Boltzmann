/* print_boundary_flux.c
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

#include "print_boundary_flux.h"
void print_boundary_flux(struct state_struct *state) {
  /*
    Print out the boundary fluxes, one molecule:compartment flux per line.
    
    Called by: boltzmann_run
    Calls:     fprintf

    Arguments:
    Name           TMF       Descripton
    state          G*I       state structure. No fields are modified.
                             Used fields are sorted_molecules, 
			     bndry_flux_fp, sorted_cmpts,
			     unqiue_molecules and bndry_flux_concs.
  */
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *molecule;
  struct molecule_struct *sorted_cmpts;
  struct molecule_struct *cur_cmpt;
  double *bndry_flux_concs;
  char *molecules_text;
  char *compartment_text;
  char *cmpt_string;
  char *molecule_str;
  FILE *bndry_flux_fp;

  int ci;
  int oi;
  
  int j;
  int unique_molecules;

  sorted_molecules = state->sorted_molecules;
  bndry_flux_fp    = state->bndry_flux_fp;
  sorted_cmpts     = state->sorted_cmpts;
  unique_molecules = state->nunique_molecules;
  bndry_flux_concs = state->bndry_flux_concs;
  molecules_text   = state->molecules_text;
  compartment_text = state->compartment_text;
  if (bndry_flux_fp) {
    fprintf(bndry_flux_fp,"final flux\n");
    molecule    = sorted_molecules;
    cmpt_string = NULL;
    oi = -1;
    for (j=0;j<unique_molecules;j++) {
      ci = molecule->c_index;
      if (ci != oi) {
	oi = ci;
	if (ci > 0) {
	  cur_cmpt = (struct molecule_struct *)&(sorted_cmpts[ci]);
	  cmpt_string = (char *)&compartment_text[cur_cmpt->string];
	} else {
	  cmpt_string = NULL;
	}
      }
      if (molecule->variable == 0) {
	molecule_str = (char *)&molecules_text[molecule->string];
	if (ci > 0) {
	  fprintf(bndry_flux_fp,"%s:%s\t%le\n",molecule_str,
		  cmpt_string,bndry_flux_concs[j]);
	} else {
	  fprintf(bndry_flux_fp,"%s\t%le\n",molecule_str,
		  bndry_flux_concs[j]);
	}
      }
      molecule += 1; /* Caution address arithmetic.*/
    } /* end for (j...) */
  }
  return;
}
