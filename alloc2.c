/* alloc2.c
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

#include "alloc2.h"
int alloc2(struct state_struct *state) {
  /*
    Allocate space for the text strings to store all of the
    information within the reactions file. Also allocate space
    for the reaction structure and meta data.
    Called by: boltzmann
    Calls:     calloc, fprintf (intrinsic)
  */
  struct rxn_struct rs;
  struct rxn_struct *reactions;
  struct rxn_matrix_struct rms;
  struct rxn_matrix_struct *reactions_matrix;
  struct istring_elem_struct ises;
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *unsorted_molecules;
  struct istring_elem_struct *sorted_cmpts;
  struct istring_elem_struct *unsorted_cmpts;
  int64_t usage;
  int64_t rxn_title_space;
  int64_t pathway_space;
  int64_t compartment_space;
  int64_t molecules_space;
  int64_t align_len;
  int64_t align_mask;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t molecules_len;
  int64_t ask_for;
  int64_t one_l;
  int64_t nze;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *molecules_text;
  char *raw_molecules_text;
  int num_rxns;
  int num_molecules;
  int num_cmpts;
  int success;
  success    = 1;
  usage      = state->usage;
  one_l      = (int64_t)1;
  align_mask = state->align_mask;
  align_len  = state->align_len;
  num_rxns   = state->number_reactions;
  num_molecules = state->number_molecules;
  num_cmpts   = state->number_compartments;
  rxn_title_space =  state->rxn_title_len + ((int64_t)num_rxns) * align_len;
  rxn_title_space += align_len - (rxn_title_space & align_mask);
  pathway_space   =  state->pathway_len + ((int64_t)num_rxns) * align_len;
  pathway_space   += align_len - (pathway_space & align_mask);
  if (state->number_compartments > 0) {
    compartment_space = state->compartment_len +
                        ((int64_t)num_rxns) * align_len;
    compartment_space += align_len - (compartment_space & align_mask);
  } else {
    compartment_space = 0;
  }
  molecules_space  = state->molecules_len + ((int64_t)num_molecules) * align_len;
  molecules_space  += align_len - (molecules_space & align_mask);
  ask_for = rxn_title_space + pathway_space + compartment_space +
    molecules_space + molecules_space;
  usage += ask_for;
  rxn_title_text = (char *) calloc(one_l,ask_for);
  if (rxn_title_text) {
    state->rxn_title_space   = rxn_title_space;
    state->pathway_space     = pathway_space;
    state->compartment_space = compartment_space;
    state->molecules_space     = molecules_space;
    /*
      Caution address arithmetic follows.
    */
    pathway_text = rxn_title_text + rxn_title_space;
    compartment_text = pathway_text + pathway_space;
    molecules_text     = compartment_text + compartment_space;
    raw_molecules_text = molecules_text + molecules_space;
    state->rxn_title_text   = rxn_title_text;
    state->pathway_text     = pathway_text;
    state->compartment_text = compartment_text;
    state->molecules_text     = molecules_text;
    state->raw_molecules_text = raw_molecules_text;
  } else {
    fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	    "for text strings in core.\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    /*
      Allocate space for the reaction meta data.
    */
    ask_for = ((int64_t)num_rxns+one_l) * ((int64_t)sizeof(rs));
    usage   += ask_for;
    state->reactions = (struct rxn_struct *)calloc(one_l,ask_for);
    if (state->reactions == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reaction meta data \n",
	      ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    /*
      Allocate space for the reaction matrix.
    */
    ask_for = (int64_t)sizeof(rms);
    usage   += ask_for;
    reactions_matrix  = (struct rxn_matrix_struct *)calloc(one_l,ask_for);
    if (reactions_matrix) {
      state->reactions_matrix = reactions_matrix;
    } else {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reaction meta data \n",
	      ask_for);
      success = 0;
    }
  }
  /*
    Allocate space for the reaction matrix arrays, each having nze elements
    where nze = num_molecules + num_rxns + 1;
  */
  if (success) {
    nze = (int64_t)num_molecules + (int64_t)num_rxns + one_l;
    ask_for = nze * ((int64_t)sizeof(int64_t));
    usage += ask_for;
    reactions_matrix->rxn_ptrs = (int64_t*)calloc(one_l,ask_for);
    if (reactions_matrix->rxn_ptrs == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reactions_matrix->rxn_ptrs\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    usage += ask_for;
    reactions_matrix->molecules_indices = (int64_t*)calloc(one_l,ask_for);
    if (reactions_matrix->molecules_indices == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reactions_matrix->molecules_indices\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    usage += ask_for;
    reactions_matrix->compartment_indices = (int64_t*)calloc(one_l,ask_for);
    if (reactions_matrix->compartment_indices == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reactions_matrix->compartment_indices\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  
  if (success) {
    usage += ask_for;
    reactions_matrix->coefficients = (int64_t*)calloc(one_l,ask_for);
    if (reactions_matrix->coefficients == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reactions_matrix->coefficients\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = nze * ((int64_t)sizeof(char *));
    usage += ask_for;
    reactions_matrix->text = (char **)calloc(one_l,ask_for);
    if (reactions_matrix->text == NULL) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for reactions_matrix->text\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  /*
    Allocate space to store the sorted molecules pointers,
    and scratch space for sorting the molecules pointers.
  */
  if (success) {
    ask_for = ((int64_t)num_molecules) * ((int64_t)sizeof(ises));
    ask_for = ask_for << 1;
    usage += ask_for;
    state->unsorted_molecules = (struct istring_elem_struct *)calloc(one_l,ask_for);
    if (state->unsorted_molecules) {
      /*
	Caution address arithmetic follows.
	state->sorted_molecules = &state->unsorted_molecules[num_molecules];
      */
      state->sorted_molecules = state->unsorted_molecules + num_molecules;
    } else {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for sorted and unsorted molecules pointers\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = ((int64_t)num_cmpts) * ((int64_t)sizeof(ises));
    ask_for = ask_for << 1;
    usage += ask_for;
    state->unsorted_cmpts = (struct istring_elem_struct *)calloc(one_l,ask_for);
    if (state->unsorted_molecules) {
      /*
	Caution address arithmetic follows.
	state->sorted_cmpts = &state->unsorted_cmpts[num_cmpts];
      */
      state->sorted_cmpts = state->unsorted_cmpts + num_cmpts;
    } else {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for sorted and unsorted molecules pointers\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  } /* end if (success) */
  if (success) {
    ask_for = ((int64_t)num_cmpts);
    ask_for += (ask_for & one_l);
    ask_for *= ((int64_t)sizeof(int));
    usage += ask_for;
    state->cmpts_map = (int*)calloc(one_l,ask_for);
    if (state->cmpts_map == 0) {
      fprintf(stderr,"alloc2: Error, unable to allocate %ld bytes of space "
	      "for cmpts_map.\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = ((int64_t)num_rxns) * ((int64_t)sizeof(double));
    usage += ask_for;
    state->activities= (double*)calloc(one_l,ask_for);
    if (state->activities == NULL) {
      fprintf(stderr,"alloc2: Error unable to allocate %ld bytes for "
	      "state->activities field.\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }


  state->usage = usage;
  return(success);
}
