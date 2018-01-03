/* sort_json_ids.c
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

#include "merge_sorted_strings.h"

#include "sort_json_ids.h"
int sort_json_ids(struct sbml2bo_struct *state) {
  /*
    Sort the json id's from a kg2js data struct (kg2js_data)
    List of sorted string pointers is returned in state->json_ids.
    Called by: sbml2bo, sbml_to_boltzmann
    Calls:     merge_sorted_strings
  */
  struct t2js_struct *kg2js_data;
  int success;
  int padi;
  int64_t num_kegg_ids;
  int64_t n;
  int64_t i;
  int64_t j;
  int64_t step;
  int64_t step2;
  int64_t l1;
  int64_t l2;
  int64_t ln;
  int64_t offset;
  int64_t noffset;
  char **cur_list;
  char **next_list;
  char **tmp_list;
  char **kjs_ids;
  char **json_ids;
  char **sort_scratch;
  success = 1;
  json_ids = state->json_ids;
  sort_scratch = state->sort_json_ids_scratch;
  num_kegg_ids = state->num_kegg_ids;
  kg2js_data = state->kg2js_data;
  kjs_ids    = kg2js_data->dictionary_ids;
  for (i=0;i<num_kegg_ids;i++) {
    sort_scratch[i] = kjs_ids[i];
  }
  n = num_kegg_ids;
  if (n > 1) {
    cur_list = (char**)&sort_scratch[0];
    next_list = (char**)&sort_scratch[n];
    for (step=1;step < n;step += step) {
      step2 = step + step;
      if (step2 >= n) {
	next_list = json_ids;
      }
      for (j=0;j<(n - step); j = j + step2) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) l2 = step;
	merge_sorted_strings(l1,&cur_list[j],
			     l2,&cur_list[j+step],
			     &next_list[j]);
      }
      ln = n & (step2 - 1);
      if (ln <= step) {
	for (j=n-ln;j<n;j++) {
	  next_list[j] = cur_list[j];
	}
      }
      tmp_list = cur_list;
      cur_list = next_list;
      next_list = tmp_list;
    }
  } else {
    for (i=0;i<n;i++) {
      json_ids[i] = kjs_ids[i];
    }
  }
  return(success);
}

