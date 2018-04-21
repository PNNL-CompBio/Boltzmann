/* sbml_alloc2.c
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

#include "sbml_alloc2.h"

int sbml_alloc2(struct sbml2bo_struct *state,
		int64_t num_modelseed_ids,
		int64_t length_ms2js_strings,
		int64_t num_kegg_ids,
		int64_t length_kg2js_strings) {
  /*
    Allocate the sbml2bo_struct and its file name fields.
    Called by: sbml_to_boltzmann, sbml2bo
    Calls:     calloc
  */
  struct t2js_struct t2js_instance;
  struct t2js_struct *ms2js_data;
  struct t2js_struct *kg2js_data;
  int64_t ask_for;
  int64_t one_l;
  char **ms_ids;
  char **kegg_ids;
  char **json_ids;
  char *ms2js_strings;
  char *kg2js_strings;
  int     success;
  int     padi;
  FILE *lfp;
  success = 1;
  one_l   = (int64_t)1;
  lfp     = state->log_fp;
  ask_for = (int64_t)sizeof(t2js_instance);
  ms2js_data = (struct t2js_struct*)calloc(one_l,ask_for);
  if (lfp == NULL) {
    lfp = stderr;
  }
  if (ms2js_data == NULL) {
    success = 0;
  }
  if (success) {
    state->ms2js_data = ms2js_data;
    ask_for = (int64_t)sizeof(t2js_instance);
    kg2js_data = (struct t2js_struct*)calloc(one_l,ask_for);
    if (kg2js_data == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for kg2js_data\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    state->kg2js_data = kg2js_data;
    ask_for = length_ms2js_strings;
    ms2js_strings = (char*)calloc(one_l,ask_for);
    if (ms2js_strings == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for ms2js_strings\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    ms2js_data->strings = ms2js_strings;
    ask_for = num_modelseed_ids * (sizeof(char *) << 1);
    ms_ids  = (char **)calloc(one_l,ask_for);
    if (ms_ids == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for ms_ids\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    ms2js_data->dictionary_ids = ms_ids;
    ms2js_data->json_ids = (char **)&ms_ids[num_modelseed_ids]; 
    ms2js_data->num_ids = num_modelseed_ids;
    ms2js_data->length_strings = length_ms2js_strings;
  }
  if (success) {
    ask_for = length_kg2js_strings;
    kg2js_strings = (char*)calloc(one_l,ask_for);
    if (kg2js_strings == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for kg2js_strings\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    kg2js_data->strings = kg2js_strings;
    ask_for = num_kegg_ids * (sizeof(char *) << 1);
    kegg_ids = (char **)calloc(one_l,ask_for);
    if (kegg_ids == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for kegg_ids\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    kg2js_data->dictionary_ids = kegg_ids;
    kg2js_data->json_ids = (char **)&kegg_ids[num_kegg_ids]; 
    kg2js_data->num_ids = num_kegg_ids;
    kg2js_data->length_strings = length_kg2js_strings;
    /*
      We need a list of sorted json id's for the json id's matching a kegg id.
    */
    ask_for = (int64_t)((num_kegg_ids << 1) + num_kegg_ids) * sizeof(char *);
    json_ids = (char **)calloc(one_l,ask_for);
    if (json_ids == NULL) {
      success = 0;
      fprintf(lfp,"sbml_alloc2: Error unable to allocate %ld bytes for json_ids and sort_json_scratch\n",ask_for);
      fflush(lfp);
    }
  }
  if (success) {
    state->json_ids = json_ids;
    state->sort_json_ids_scratch = (char**)&json_ids[num_kegg_ids];
  }
  return(success);
}
