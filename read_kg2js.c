/* read_kg2js.c
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

#include "count_ntb.h"

#include "read_kg2js.h"

int read_kg2js(struct sbml2bo_struct *sbml_state) {
  /*
    read the moduleseed to json id file translations, assumed to be
    sorted on moduleseed id's.
    Called by: sbml_to_boltzmann.
    Calls:     count_ntb, strcpy, fopen, fgets, fclose
  */
  struct t2js_struct *kg2js_data;
  FILE *kg2js_fp;
  char *kg2js_file;
  char *kg2js_string;
  char **kegg_ids;
  char **js_ids;
  int64_t num_kegg_ids;
  int64_t i;
  int64_t alignment;
  int64_t align_mask;
  int64_t pads;
  int64_t extra_len;
  int success;
  int line_len;
  int nc;
  int id_len;
  char *id_string;
  char line_buffer[1024];
  char *line;
  success = 1;
  line_len = 1024;
  kg2js_data = sbml_state->kg2js_data;
  kg2js_file = sbml_state->kg2js_file;
  kg2js_fp = fopen(kg2js_file,"r");
  alignment = sbml_state->alignment;
  align_mask = sbml_state->align_mask;
  num_kegg_ids = kg2js_data->num_ids;
  if (kg2js_fp == NULL) {
    success = 0;
  }
  if (success) {
    kegg_ids = kg2js_data->dictionary_ids;
    js_ids = kg2js_data->json_ids;
    line = (char*)&line_buffer[0];
    kg2js_string = kg2js_data->strings;
    for (i=0;i<num_kegg_ids;i++) {
      fgets(line,line_len,kg2js_fp);
      id_string = line;
      nc = count_ntb(id_string);
      if (nc > 0) {
	id_string[nc] = 0;
	*kegg_ids = kg2js_string;
	kegg_ids += 1; /* Caution address arithmetic */
	strcpy(kg2js_string,id_string);
	id_len = nc + 1;
	extra_len = (int64_t)id_len & align_mask;
	pads      = id_len + ((alignment - extra_len) & align_mask);
	id_string += (nc+1);
	kg2js_string += pads; /* Caution address arithmetic */
	nc = count_ntb(id_string);
	if (id_string[nc-1] == '\n') {
	  nc = nc-1;
	}
	id_string[nc] = 0;
	*js_ids = kg2js_string;
	js_ids += 1;  /* Caution address arithmetic */
	strcpy(kg2js_string,id_string);
	id_len = nc + 1;
	extra_len = (int64_t)id_len & align_mask;
	pads      = id_len + ((alignment - extra_len) & align_mask);
	kg2js_string += pads; /* Caution address arithmetic */
      }
    }
    fclose (kg2js_fp);
  }
  return (success);
}
