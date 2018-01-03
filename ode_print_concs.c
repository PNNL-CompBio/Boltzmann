/* ode_print_concs.c
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

#include "ode_print_concs.h"
void ode_print_concs(struct state_struct *state, double time, double *concs) {
  /* 
    print the molecule concentraions
    Prints out the current concentrations field of the state structure
    in a tab delimited row terminated by a newline.

    Called by: ode23tb

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are unique_molecules,
					    sorted_molecules,
					    ode_concs_fp
                           no fields of state are modified.

    step          JSI      eight byte integer step number, -1 for initial step.
                  
    
  */
  struct molecule_struct *cur_molecule;
  int unique_molecules;
  int j;

  int print_concs_or_counts;
  int do_concs;

  int do_counts;
  int padi;

  FILE *ode_concs_fp;
  FILE *ode_counts_fp;

  ode_concs_fp           = state->ode_concs_fp;
  ode_counts_fp          = state->ode_counts_fp;
  print_concs_or_counts  = state->print_concs_or_counts;
  unique_molecules       = state->nunique_molecules;
  cur_molecule           = state->sorted_molecules;

  do_concs = 0;
  do_counts = 0;
  if (print_concs_or_counts & 2) {
    if (ode_concs_fp) {
      do_concs = 1;
    }
  }
  if (print_concs_or_counts & 1) {
    if (ode_counts_fp) {
      do_counts = 1;
    }
  }
  if (do_concs) {
    fprintf(ode_concs_fp,"%le",time);
  }
  if (do_counts) {
    fprintf(ode_counts_fp,"%le",time);
  }
  if (do_concs || do_counts) {
    for (j=0;j<unique_molecules;j++) {
      if ((cur_molecule->solvent == 0) || (cur_molecule->variable == 1)) {
	if (do_concs) {
	  fprintf(ode_concs_fp,"\t%le",concs[j]);
	}
	if (do_counts) {
	  fprintf(ode_counts_fp,"\t%le",concs[j]);
	}
      }
      cur_molecule += 1; /* caution address arithmetic.*/
    }
    if (do_concs) {
      fprintf(ode_concs_fp,"\n");
      fflush(ode_concs_fp);
    }
    if (do_counts) {
      fprintf(ode_counts_fp,"\n");
      fflush(ode_counts_fp);
    }
  } /* end if (do_concs || do_counts) */
}
