/* sbml_sort_species_trans.c
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

#include "sbml_merge_species_trans.h"

#include "sbml_sort_species_trans.h"
int sbml_sort_species_trans(int num_species,
			    char **spec_ids, char **translations,
			    char **sort_species_trans_scratch) {
  /*
    Sort specid,translation(to json_id) string pointer pairs on
    the specid field.
    Called by: parse_sbml
    Calls:     sbml_merge_species_trans
  */


  int success;
  int padi;
  int64_t n;
  int64_t i;
  int64_t j;
  int64_t step;
  int64_t step2;
  int64_t l1;
  int64_t l2;
  int64_t ln;
  char **cur_spec_ids;
  char **next_spec_ids;
  char **tmp_spec_ids;
  char **cur_trans_ids;
  char **next_trans_ids;
  char **tmp_trans_ids;
  char **species_scratch;
  char **trans_scratch;

  success = 1;
  n       = num_species;
  species_scratch = sort_species_trans_scratch;
  trans_scratch   = &sort_species_trans_scratch[(n << 1)];
  for (i=0;i<n;i++) {
    species_scratch[i] = spec_ids[i];
    trans_scratch[i] = translations[i];
  }
  if (n > 1) {
    cur_spec_ids = (char**)&species_scratch[0];
    next_spec_ids = (char**)&species_scratch[n];
    cur_trans_ids = (char**)&trans_scratch[0];
    next_trans_ids = (char**)&trans_scratch[n];
    for (step=1;step < n;step += step) {
      step2 = step + step;
      if (step2 >= n) {
	next_spec_ids = spec_ids;
	next_trans_ids = translations;
      }
      for (j=0;j<(n - step); j = j + step2) {
	l1 = step;
	l2 = n - j - step;
	if (l2 > step) l2 = step;
	sbml_merge_species_trans(l1,&cur_spec_ids[j],&cur_trans_ids[j],
				 l2,&cur_spec_ids[j+step],
				 &cur_trans_ids[j+step],
				 &next_spec_ids[j],
				 &next_trans_ids[j]);
      }
      ln = n & (step2 - 1);
      if (ln <= step) {
	for (j=n-ln;j<n;j++) {
	  next_spec_ids[j]  = cur_spec_ids[j];
	  next_trans_ids[j] = cur_trans_ids[j];
	}
      }
      tmp_spec_ids = cur_spec_ids;
      cur_spec_ids = next_spec_ids;
      next_spec_ids = tmp_spec_ids;
      tmp_trans_ids = cur_trans_ids;
      cur_trans_ids = next_trans_ids;
      next_trans_ids = tmp_trans_ids;
    }
  } 
  return(success);
}

