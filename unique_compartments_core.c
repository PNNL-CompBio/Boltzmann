/* unique_compartments_core.c
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

#include "unique_compartments_core.h"
int unique_compartments_core(int nzr,
			     struct istring_elem_struct *sorted_cmpts,
			     char   *compartment_text,
			     int64_t *compartment_indices,
			     int64_t *nunique_compartments,
			     int64_t *sum_compartment_len,
			     int64_t align_len, 
			     int64_t align_mask) {

  /*
    Remove duplicates from the sorted_compartments list
    and set the compartment_indices fields in the
    reactions_matrix appropriately.
    Called by: unique_compartments, boltzmann_boot
    Calls:     strcmp (intrinsic)
  */
  struct istring_elem_struct *srtd_cmpts;
  struct istring_elem_struct *cur_cmpt;
  struct istring_elem_struct *ucmpts_next;
  char *sstring;
  char *cstring;
  int64_t compartment_len;
  int64_t cmpt_size;
  int64_t pad_size;

  int i;
  int padi;

  int success;
  int nu_cmpts;

  success = 1;
  compartment_len = (int64_t)0;
  srtd_cmpts = sorted_cmpts;
  /* loop over sorted compartments. */
  nu_cmpts = 0;
  
  if (nzr > 1) {
    if (srtd_cmpts->string >= 0) {
      compartment_indices[srtd_cmpts->c_index] = nu_cmpts;
      cstring = (char*)&compartment_text[srtd_cmpts->string];
      cmpt_size = strlen(cstring);
      pad_size = (align_len  - (cmpt_size & align_mask)) & align_mask;
      if (cmpt_size == 0) {
	pad_size = 16;
	compartment_indices[srtd_cmpts->c_index] = 0;
      }
      compartment_len += (pad_size + cmpt_size);
    } else {
      compartment_indices[srtd_cmpts->c_index] = 0;
    }
    cur_cmpt = srtd_cmpts;
    srtd_cmpts += 1; /* Caution address arithmetic. */
    ucmpts_next  = srtd_cmpts;
  }
  for (i=1;i<nzr;i++) {
    if (srtd_cmpts->string>=0) {
      sstring = (char*)&compartment_text[srtd_cmpts->string];
      if ((cur_cmpt->string < 0) ||
	  (strcmp(sstring,cstring) != 0)) {
	nu_cmpts += 1;
	cur_cmpt = srtd_cmpts;
	cstring  = sstring;
	cmpt_size = strlen(cstring) + 1;
	pad_size = (align_len  - (cmpt_size & align_mask)) & align_mask;
	compartment_len += (int64_t)(pad_size + cmpt_size);
	ucmpts_next->string = cur_cmpt->string;
	ucmpts_next += 1; /* Caution address arithmetic. */
      }
      compartment_indices[srtd_cmpts->c_index] = nu_cmpts;
    } else {
      compartment_indices[srtd_cmpts->c_index] = 0;
    }
    srtd_cmpts += 1; /* Caution address arithmetic. */
  }
  *nunique_compartments = nu_cmpts + 1;
  *sum_compartment_len = compartment_len;
  return(success);
}
