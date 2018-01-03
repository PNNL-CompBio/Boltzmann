/* print_counts.c
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

#include "print_counts.h"
void print_counts(struct state_struct *state, int64_t step) {
  /* 
    print the molecule counts and or concentrations.
    if concs_or_counts is 1 or 3, the current counts are printed.
    if it is 2 or 3 the concentrations are printed - each to their own file.

    Called by boltzmann_run, deq

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are unique_molecules,
					    current_concentrations,
					    counts_out_fp,
					    concs_out_fp,
					    concs_or_counts;
					    sorted_molecules;
                           no fields of state are modified.

    step          JSI      eight byte integer step number, -1 for initial step.
                  
    
  */
  struct molecule_struct *cur_molecule;
  double *current_counts;
  double *count_to_conc;
  double conc;
  
  int unique_molecules;
  int j;
  int concs_or_counts;
  int padi;

  FILE *counts_out_fp;
  FILE *concs_out_fp;
  counts_out_fp          = state->counts_out_fp;
  concs_out_fp           = state->concs_out_fp;
  unique_molecules       = state->nunique_molecules;
  current_counts         = state->current_counts;
  cur_molecule           = state->sorted_molecules;
  concs_or_counts        = (int)state->concs_or_counts;
  count_to_conc          = state->count_to_conc;
  if (concs_or_counts & 1) {
    if (counts_out_fp) {
      switch(step) {
      case -3:
	fprintf(counts_out_fp,"init");
	break;
      case -2:
	fprintf(counts_out_fp,"awm");
	break;
      case -1:
	fprintf(counts_out_fp,"adeq");
	break;
      default:
	fprintf(counts_out_fp,"%ld",step);
      }
      for (j=0;j<unique_molecules;j++) {
	if (j != solvent_pos) {
	  fprintf(state->counts_out_fp,"\t%le",current_counts[j]);
	}
      }
      fprintf(state->counts_out_fp,"\n");
    }
  }
  if (concs_or_counts & 2) {
    if (concs_out_fp) {
      switch(step) {
      case -3:
	fprintf(concs_out_fp,"init");
	break;
      case -2:
	fprintf(concs_out_fp,"awm");
	break;
      case -1:
	fprintf(concs_out_fp,"adeq");
	break;
      default:
	fprintf(concs_out_fp,"%ld",step);
      }
      for (j=0;j<unique_molecules;j++) {
	if ((cur_molecule->solvent == 0) || (cur_molecule->variable == 1)) {
	  conc = current_counts[j] * count_to_conc[j];
	  fprintf(state->concs_out_fp,"\t%le",conc);
	}
	cur_molecule += 1; /* caution address arithmetic.*/
      }
      fprintf(state->concs_out_fp,"\n");
    }
  }
  return;
}
