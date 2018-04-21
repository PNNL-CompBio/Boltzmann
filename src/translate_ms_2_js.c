/* translate_ms_2_js.c
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

#include "translate_ms_2_js.h"
char *translate_ms_2_js (struct ms2js_struct *ms2js_data, char *sid) {
  /*
    Convert a modelseed id, sid to a json id.
    Called by: sbml_process_list_of_species
    Calls:     strcmp
  */
  char **ms_ids;
  char **js_ids;
  char *js_id;
  char *ms_id;
  int i;
  int low;

  int high;
  int mid;

  int num_modelseed_ids;
  int cmp;

  int last;
  int padi;

  last = ms2js_data->num_modelseed_ids - 1;
  ms_ids            = ms2js_data->ms_ids;
  js_ids            = ms2js_data->js_ids;
  ms_id = *ms_ids;
  js_id = NULL;
  cmp = strcmp(ms_id,sid);
  if (cmp >= 0) {
    if (cmp == 0) {
      js_id = *js_ids;
    } else {
      ms_id = ms_ids[last];
      cmp = strcmp(ms_id,sid);
      if (cmp <= 0) {
	if (cmp == 0) {
	  js_id = js_ids[last];
	} else {
	  low = 0;
	  high = last;
	  while ((high - low) > 1) {
	    mid = (high + low) >> 1;
	    ms_id = ms_ids[mid];
	    cmp = strcmp(ms_id,sid);
	    if (cmp == 0) {
	      js_id = js_ids[mid];
	      break;
	    } else {
	      if (cmp < 0) {
		high = mid;
	      } else {
		low = mid;
	      }
	    }
	  }
	}
      }
    }
  }
  return(js_id);
}
