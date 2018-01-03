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
    print the molecule counts.
    Prints out the current counts field of the state structure
    in a tab delimited row terminated by a newline.

    Called by boltzmann_run, deq

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are unique_molecules,
					    current_concentrations,
					    counts_out_fp
                           no fields of state are modified.

    step          JSI      eight byte integer step number, -1 for initial step.
                  
    
  */
  double *current_counts;
  int unique_molecules;
  int j;
  int solvent_pos;
  int padi;

  FILE *counts_out_fp;
  counts_out_fp          = state->counts_out_fp;
  unique_molecules       = state->nunique_molecules;
  current_counts         = state->current_counts;
  solvent_pos            = (int)state->solvent_pos;
  if (counts_out_fp) {
    if (step < (int64_t)0) {
      fprintf(counts_out_fp,"init");
    } else {
      fprintf(counts_out_fp,"%lld",step);
    }
    for (j=0;j<unique_molecules;j++) {
      if (j != solvent_pos) {
	fprintf(state->counts_out_fp,"\t%le",current_counts[j]);
      }
    }
    fprintf(state->counts_out_fp,"\n");
  }
  return;
}
