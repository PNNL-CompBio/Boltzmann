/* unique_molecules_core.c
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

#include "unique_molecules_core.h"
int unique_molecules_core(int nzr,
			  struct molecule_struct *sorted_molecules,
			  char *molecules_text,
			  int64_t *molecules_map,
			  int64_t *nunique_molecules,
			  int64_t *sum_molecule_len,
			  int64_t align_len,
			  int64_t align_mask) {
  /*
    Remove duplicates from the sorted_molecules list
    and set the column_indices fields in the
    reactions_matrix appropriately.
    Called by: unique_molecules, boltzmann_boot
    Calls:     strcmp (intrinsic)
  */
  struct molecule_struct *cur_molecule;
  struct molecule_struct *umolecules_next;
  char *cstring;
  char *sstring;
  int64_t molecule_len;
  int64_t m_size;
  int64_t pad_size;
  
  int i;
  int success;

  int nu_molecules;
  int padi;

  int ni;
  int cni;

  success = 1;
  molecule_len = (int64_t)0;
  /* loop over sorted molecules. */
  nu_molecules = 0;
  molecules_map[sorted_molecules->m_index] = nu_molecules;
  cur_molecule = sorted_molecules;
  cstring = NULL;
  if (cur_molecule->string >= 0) {
    cstring = (char *)&molecules_text[cur_molecule->string];
    m_size = strlen(cstring);
    pad_size = (align_len - (m_size & align_mask)) & align_mask;
    molecule_len += (m_size + pad_size);
  }
  /*
    This translation has already been done in translate_compartments.
  cur_molecule->c_index = compartment_indices[cur_molecule->c_index];
  */
  cni = sorted_molecules->c_index;
  sorted_molecules += 1; /* Caution address arithmetic. */
  umolecules_next  = sorted_molecules;
  for (i=1;i<nzr;i++) {
    ni = sorted_molecules->c_index;
    sstring = NULL;
    if (sorted_molecules->string >= 0) {
      sstring = (char *)&molecules_text[sorted_molecules->string];
    }
    if ((ni != cni)  ||
	(strcmp(sstring,cstring) != 0)) {
      cstring = sstring;
      m_size = strlen(cstring) + 1;
      pad_size = (align_len - (m_size & align_mask)) & align_mask;
      molecule_len += (int64_t)(m_size + pad_size);
      nu_molecules += 1;
      cur_molecule = sorted_molecules;
      molecules_map[sorted_molecules->m_index] = nu_molecules;
      umolecules_next->string = cur_molecule->string;
      umolecules_next->m_index = nu_molecules;
      umolecules_next->c_index = ni;
      cni = ni;
      umolecules_next += 1; /* Caution address arithmetic. */
    } else {
      molecules_map[sorted_molecules->m_index] = nu_molecules;
    }
    sorted_molecules += 1; /* Caution address arithmetic. */
  }
  *nunique_molecules = nu_molecules + 1;
  *sum_molecule_len = molecule_len;

  return(success);
}
