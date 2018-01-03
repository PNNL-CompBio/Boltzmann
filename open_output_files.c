/* open_output_files.c
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

#include "open_output_files.h"
int open_output_files(struct state_struct *state) {
  /*
    Open the log, counts, rxn_likelihoods, free_energy, and 
    bndry_flux output files.
    Called by: boltzmann_init
    Calls:     fopen, fprintf, fflush.
  */
  int success;
  int padi;
  success = 1;
  if (success) {
    if (state->log_file) {
      /*
	Open the log file.
      */
      state->lfp = fopen(state->log_file,"w");
      if (state->lfp == NULL) {
	fprintf(stderr,
		"open_output_files: unable to open log_file, %s, quitting.\n",
		state->log_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->counts_out_file) {
      /*
	Open the counts output file.
      */
      state->counts_out_fp = fopen(state->counts_out_file,"w");
      if (state->counts_out_fp == NULL) {
	fprintf(stderr,
		"output_output_files: unable to open counts_out_file, %s, quitting.\n",
		state->counts_out_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->rxn_lklhd_file) {
      /*
	Open the likelihoods output file.
      */
      state->rxn_lklhd_fp = fopen(state->rxn_lklhd_file,"w");
      if (state->rxn_lklhd_fp == NULL) {
	fprintf(stderr,
		"open_output_files unable to open rxn_lklhd_file, %s, quitting.\n",
		state->rxn_lklhd_file);
	fflush(stderr);
	success = 0;
      }
    }
  }
  if (success) {
    if (state->free_energy_format > (int64_t)0) {
      /*
	Open the free energy output file.
      */
      if (state->free_energy_file) {
	state->free_energy_fp = fopen(state->free_energy_file,"w");
	if (state->free_energy_fp == NULL) {
	  fprintf(stderr,
		  "open_output_files: unable to open free_energy_file, %s, quitting.\n",
		  state->free_energy_file);
	  fflush(stderr);
	  success = 0;
	}
      }
    }
  }
  if (success) {
    if (state->bndry_flux_file) {
      /*
	Open the boundary flux output file.
      */
      state->bndry_flux_fp = fopen(state->bndry_flux_file,"w");
      if (state->bndry_flux_fp == NULL) {
	fprintf(stderr,
		"open_output_files: unable to open bndry_flux_file, %s, quitting.\n",
		state->bndry_flux_file);
	fflush(stderr);
	success = 0;
      } 
    }
  }
  return (success);
}
