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
int unique_compartments_core(int num_compartments,
			     struct compartment_struct *sorted_cmpts,
			     char   *compartment_text,
			     int    *cmpt_tracking,
			     int64_t *nunique_compartments,
			     int64_t *sum_compartment_len,
			     int64_t align_len) {

  /*
    Remove duplicates from the sorted_compartments list
    and set the compartment_indices fields in the
    reactions_matrix appropriately.
    Called by: unique_compartments, boltzmann_boot
    Calls:     strcmp (intrinsic)

    Arguments            TMF
    num_compartments     I0I Number of compartments as counted by 
                             parse_reactions_file = state->number_compartments)

 		     
    sorted_cmpts         G*I List of compartment structures sorted 
                             alphabetically by sort_compartments
		     
    compartment_text     C*I Character array containing compartment names.
    

    cmpt_tracking        I*W Integer of vector of length nzr
                             Used to map unsorted compartment numbers
			     to sorted compartment numbers.

    nunique_compartments L*O Number of unique compartments including the
                             default global compartment.

    sum_compartment_len  L*O Sum of the compartment lengths. Only
                             used by boltzmann_boot.

    align_len            L0I String alignment length  = state->align_len

    align_mask           

  */
  struct compartment_struct *compartment;
  struct compartment_struct *next_ucmpt;
  char *sstring;
  char *cstring;
  int64_t compartment_len;
  int64_t cmpt_size;
  int64_t pad_size;
  int64_t align_mask;
  int64_t string_offset;
  

  int i;
  int ci;

  int success;
  int nu_cmpts;

  success = 1;
  align_mask = align_len - (int64_t)1;

  compartment_len = (int64_t)0;

  compartment = sorted_cmpts;
  /*
    srtd_cmpts = sorted_cmpts;
  */
  /* loop over sorted compartments. */
  nu_cmpts = 0;
  cstring = "";
  /*
    We always have an default no-name global compartment.
    (see parse_reactions_file.c)
  */
  cmpt_tracking[0] = 0;
  compartment_len = align_len;
  /*
    cur_cmpt = srtd_cmpts;
    srtd_cmpts   += 1; // Caution address arithmetic. 
  */
  compartment  += 1; /* Caution address arithmetic. */
  next_ucmpt   =  compartment;
  for (i=1;i<num_compartments;i++) {
    string_offset = compartment->string;
    ci = compartment->c_index;
    if (string_offset > 0) {
      sstring = (char*)&compartment_text[string_offset];
      if ((strcmp(sstring,cstring) != 0)) {
	nu_cmpts += 1;
	cstring  = sstring;
	cmpt_size = strlen(cstring) + 1;
	pad_size = (align_len  - (cmpt_size & align_mask)) & align_mask;
	compartment_len += (int64_t)(pad_size + cmpt_size);
	next_ucmpt->string = string_offset;
	next_ucmpt->c_index = nu_cmpts;
	next_ucmpt+= 1; /* Caution address arithmetic. */
      }
    }
    cmpt_tracking[ci] = nu_cmpts;
    compartment  += 1; /* Caution address arithmetic. */
  } /* end for (i...) */
  *nunique_compartments = nu_cmpts + 1;
  *sum_compartment_len = compartment_len;
  return(success);
}
