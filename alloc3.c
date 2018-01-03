/* alloc3.c
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

#include "alloc3.h"

int alloc3(struct state_struct *state) {
  /*
    Allocate space for the molecule_name field for reading
    in intial concentrations, the concentrations vector,
    and the molecules matrix and its fields. The molecules 
    matrix is a transpose of the reactions matris.
    Called by: boltzmann_init
    Calls:     calloc, fprintf, fflush (intrinsic)
  */
  /*
  struct molecules_matrix_struct sms;
  struct molecules_matrix_struct *molecules_matrix;
  */
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t align_len;
  int64_t align_mask;
  int nu_molecules;
  int max_molecule_len;

  int success;
  int nrxns;
  int rxn_view_freq;
  int rxn_view_hist_length;

  int nzr;
  int max_compartment_len;
  success = 1;
  one_l      		  = (int64_t)1;
  usage      		  = state->usage;
  align_mask 		  = state->align_mask;
  align_len  		  = state->align_len;
  nu_molecules            = (int)state->unique_molecules;
  nzr                     = (int)state->number_molecules;
  nrxns                   = (int)state->number_reactions;
  max_molecule_len        = (int)state->max_molecule_len + 1;
  max_compartment_len     = (int)state->max_compartment_len + 1;
  rxn_view_freq           = (int)state->rxn_view_freq;
  /*
    Allocate space for compartment pointers in the sorted molecules list -
    length is unique_compartments + 1;
  */
  if (success) {
    ask_for = ((int64_t)state->unique_compartments + one_l) * 
      ((int64_t)sizeof(int64_t));
    state->compartment_ptrs = (int64_t*)calloc(one_l,ask_for);
    if (state->compartment_ptrs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "compartment_ptrs field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }

  if (success) {
    /*
      Allocate space for the current concentrations buffer.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->current_concentrations = (double *)calloc(one_l,ask_for);
    if (state->current_concentrations == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "current_concentrations field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    /*
      Allocate space for the boudary flux concentrations buffer.
    */
    ask_for = ((int64_t)nu_molecules) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->bndry_flux_concs = (double *)calloc(one_l,ask_for);
    if (state->bndry_flux_concs == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "bndry_flux_concs field\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->dg0s = (double *)calloc(one_l,ask_for);
    if (state->dg0s == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->dg0s field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
    
  }
  if (success) {
    ask_for = ((int64_t)nrxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->ke = (double *)calloc(one_l,ask_for);
    if (state->ke == NULL) {
      fprintf(stderr,"alloc3: Error unable to allocate %ld bytes for "
	      "state->ke field.\n",ask_for);
      fflush(stderr);
      success = 0;
    } 
  }
  state->usage = usage;
  return(success);
}
