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
		int64_t length_ms2js_strings) {
  /*
    Allocate the sbml2bo_struct and its file name fields.
    Called by: sbml_to_boltzmann, sbml2bo
    Calls:     calloc
  */
  struct ms2js_struct ms2js_instance;
  struct ms2js_struct *ms2js_data;
  int64_t ask_for;
  int64_t one_l;
  char **ms_ids;
  char **js_ids;
  char *sbml_file;
  char *concs_in_file;
  char *rxns_dat_file;
  char *cmpts_dat_file;
  char *ms2js_file;
  char *log_file;
  char *ms2js_strings;
  int     max_file_name_len;
  int     success;
  int     num_files;
  int     padi;
  success = 1;
  one_l   = (int64_t)1;
  ask_for = (int64_t)sizeof(ms2js_instance);
  ms2js_data = (struct ms2js_struct*)calloc(one_l,ask_for);
  if (ms2js_data == NULL) {
    success = 0;
  }
  if (success) {
    state->ms2js_data = ms2js_data;
    ask_for = length_ms2js_strings;
    ms2js_strings = (char*)calloc(one_l,ask_for);
    if (ms2js_strings == NULL) {
      success = 0;
    }
  }
  if (success) {
    ms2js_data->ms2js_strings = ms2js_strings;
    ask_for = num_modelseed_ids * (sizeof(char *) << 1);
    ms_ids  = (char **)calloc(one_l,ask_for);
    if (ms_ids == NULL) {
      success = 0;
    }
  }
  if (success) {
    ms2js_data->ms_ids = ms_ids;
    ms2js_data->js_ids = (char **)&ms_ids[num_modelseed_ids]; /* Caution address arithmetic */
    ms2js_data->num_modelseed_ids = num_modelseed_ids;
    ms2js_data->length_ms2js_strings = length_ms2js_strings;
  }
  return(success);
}
