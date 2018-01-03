/* sbml_alloc0.c
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

#include "sbml_alloc0.h"

int sbml_alloc0(struct sbml2bo_struct **state_p) {
  /*
    Allocate the sbml2bo_struct and its file name fields.
    Called by: sbml_to_boltzmann, sbml2bo
    Calls:     calloc
  */
  struct sbml2bo_struct state_instance;
  struct sbml2bo_struct *state;
  int64_t ask_for;
  int64_t one_l;
  char *sbml_file;
  char *concs_in_file;
  char *rxns_dat_file;
  char *cmpts_dat_file;
  char *ms2js_file;
  char *kg2js_file;
  char *id_name_file;
  char *log_file;
  int     max_file_name_len;
  int     success;
  int     num_files;
  int     padi;
  success = 1;
  /*
    Currently we have only seven files, but will allocate
    space for an extra one.
  */
  num_files = 8;
  one_l   = (int64_t)1;
  max_file_name_len = 128;
  ask_for = (int64_t)sizeof(state_instance);
  state = (struct sbml2bo_struct *)calloc(one_l,ask_for);
  if (state == NULL) {
    fprintf(stderr,
	    "sbml_alloc0: Error, unable to allocate %lld bytes for "
	    "sbml2bo_state\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    *state_p = state;
    state->num_files = num_files;
    state->file_name_len = max_file_name_len - 1;
    ask_for = num_files * max_file_name_len;
    sbml_file = (char *)calloc(one_l,ask_for);
    if (sbml_file == NULL) {
      fprintf(stderr,
	    "sbml_alloc0: Error, unable to allocate %lld bytes for "
	    "file names\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    concs_in_file  = (char*)&sbml_file[max_file_name_len];
    rxns_dat_file  = (char*)&concs_in_file[max_file_name_len];
    cmpts_dat_file = (char*)&rxns_dat_file[max_file_name_len];
    ms2js_file     = (char*)&cmpts_dat_file[max_file_name_len];
    kg2js_file     = (char*)&ms2js_file[max_file_name_len];
    id_name_file   = (char*)&kg2js_file[max_file_name_len];
    log_file       = (char*)&id_name_file[max_file_name_len];
    state->sbml_file      = sbml_file;
    state->concs_in_file  = concs_in_file;
    state->rxns_dat_file  = rxns_dat_file;
    state->cmpts_dat_file = cmpts_dat_file;
    state->ms2js_file     = ms2js_file;
    state->kg2js_file     = kg2js_file;
    state->id_name_file   = id_name_file;
    state->log_file       = log_file;
    state->num_reactions  = 0;
    state->num_species    = 0;
    state->num_cmpts      = 0;
    state->alignment      = 16;
    state->align_mask     = state->alignment-1;
    state->max_compartment_len = 1024;
    state->max_specid_len      = 1024;
    state->max_species_len     = 1024;
  }
  return(success);
}
