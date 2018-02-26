/* sbml_alloc1.c
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

#include "sbml_alloc1.h"

int sbml_alloc1(struct sbml2bo_struct *state) {
  /*
    Allocate space for sorting the sbml compartments and
    storing the compartment text.
    Called by: sbml_to_boltzmann, sbml2bo
    Calls:     calloc
  */
  struct compartment_struct compartment;
  struct compartment_struct *compartments;
  int64_t ask_for;
  int64_t one_l;
  int64_t trans_pos;
  char *specid_2_json_strings;
  char **species_trans;
  char **sort_species_trans_scratch;

  int     max_compartment_len;
  int     success;

  int     num_cmpts;
  int     num_species;

  int     i;
  int     j;

  int     max_specid_len;
  int     max_species_len;

  FILE    *lfp;
  FILE    *error_fp;
  success = 1;
  one_l   = (int64_t)1;
  max_compartment_len = state->max_compartment_len;
  num_cmpts           = state->num_cmpts;
  num_species         = state->num_species;
  max_specid_len      = state->max_specid_len;
  max_species_len     = state->max_species_len;
  lfp         = state->log_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  ask_for = (max_compartment_len + state->alignment) * num_cmpts;
  state->compartment_text = (char*)calloc(one_l,ask_for);

  if (state->compartment_text == NULL) {
    success = 0;
    fprintf(error_fp,"sbml_alloc1: Error unable to allocate %ld bytes for "
	    "compartment_text\n",ask_for);
    fflush(error_fp);
  }
  if (success) {
    ask_for = ((int64_t)(num_cmpts + num_cmpts)) * sizeof(compartment);
    compartments = (struct compartment_struct *)calloc(one_l,ask_for);
    if (compartments == NULL) {
      success = 0;
      fprintf(error_fp,"sbml_alloc1: Error unable to allocate %ld bytes for "
	      "sorted_compartment work space\n",ask_for);
      fflush(error_fp);
    }
  }
  if (success) {
    state->sorted_compartments = compartments;
    state->unsorted_compartments = (struct compartment_struct*)&(compartments[num_cmpts]);
    ask_for = (num_species*(max_specid_len+max_species_len));
    specid_2_json_strings = (char *)calloc(one_l,ask_for);
    if (specid_2_json_strings == NULL) {
      success = 0;
      fprintf(error_fp,"sbml_alloc1: Error unable to allocate %ld bytes for "
	      "specid_2_json_strings\n",ask_for);
      fflush(error_fp);
    }
  }
  if (success) {
    state->specid_2_json_strings = specid_2_json_strings;
    ask_for = ((int64_t)(num_species << 1)) * sizeof(char *);
    species_trans = (char **)calloc(one_l,ask_for);
    if (species_trans == NULL) {
      success = 0;
      fprintf(error_fp,"sbml_alloc1: Error unable to allocate %ld bytes for "
	      "species_trans\n",ask_for);
      fflush(error_fp);
    }
  }
  if (success) {
    state->spec_ids = species_trans;
    state->translations = (char **)&species_trans[num_species];
    if (state->num_kegg_ids < (num_species << 2)) {
      ask_for = ask_for + ask_for;
      sort_species_trans_scratch = (char **)calloc(one_l,ask_for);
      if (sort_species_trans_scratch == NULL) {
	success = 0;
	fprintf(error_fp,"sbml_alloc1: Error unable to allocate %ld bytes for "
		"sort_species_trans_scratch\n",ask_for);
	fflush(error_fp);
      } else {
	state->sort_species_trans_scratch = sort_species_trans_scratch;
      }
    } else {
      state->sort_species_trans_scratch = state->sort_json_ids_scratch;
    }
  }
  if (success) {
    trans_pos = (int64_t)0;
    j = 0;
    for (i=0;i<num_species;i++) {
      species_trans[j] = (char*)&specid_2_json_strings[trans_pos];
      trans_pos += max_specid_len;
      species_trans[j+1] = (char*)&specid_2_json_strings[trans_pos];
      trans_pos += max_species_len;
      j += 2;
    }
  }
  return(success);
}
