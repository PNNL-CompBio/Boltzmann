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
int unique_molecules_core(int number_molecules,
			  struct molecule_struct *sorted_molecules,
			  char *molecules_text,
			  char *solvent_string,
			  int64_t *molecules_indices,
			  int64_t *compartment_indices,
			  int64_t *nunique_molecules,
			  int64_t *sum_molecule_len,
			  int     *solvent_pos,
			  int64_t align_len,
			  int64_t align_mask) {
  /*
    Remove duplicates from the sorted_molecules list
    and set the column_indices fields in the
    reactions_matrix appropriately.

    The m_index field of the sorted molecules structure contains its
    original position in the reactions_file of each molecule.
    This is used to put the new molecule number and compartment number
    into the molecules_indices and compartment_indices vectors of
    the reactions matrix.
    The compartment number, c_index field, in the sorted molecules 
    has already been translated to the its reduces range in the
    sorted unique compartments by the previous call to translate_compartments.

    Called by: unique_molecules, global_merge_molecules
    Calls:     strcmp (intrinsic)
  */
  struct molecule_struct *molecule;
  struct molecule_struct *umolecules_next;
  char *prev_molecule_name;
  char *curr_molecule_name;
  int64_t molecule_len;
  int64_t m_size;
  int64_t pad_size;
  
  int i;
  int success;

  int nu_molecules;
  int molecule_index;

  int curr_cmpt_index;
  int prev_cmpt_index;

  success = 1;
  molecule_len = (int64_t)0;
  /* loop over sorted molecules. */
  nu_molecules = 0;
  *solvent_pos  = -1;

  molecule = sorted_molecules;
  
  molecule_index = molecule->m_index;
  prev_cmpt_index = molecule->c_index;
  molecules_indices[molecule_index] = nu_molecules;
  compartment_indices[molecule_index] = prev_cmpt_index;
  prev_molecule_name = NULL;
  if (molecule->string >= 0) {
    prev_molecule_name = (char *)&molecules_text[molecule->string];
    m_size = strlen(prev_molecule_name);
    pad_size = (align_len - (m_size & align_mask)) & align_mask;
    molecule_len += (m_size + pad_size);
    if (strcmp(prev_molecule_name,solvent_string) == 0) {
      molecule->solvent = 1;
      *solvent_pos          = 0;
    }
  }
  /*
    This translation has already been done in translate_compartments.
  molecule->c_index = cmpt_tracking[molecule->c_index];
  */
  /*
    Advance the molecule pointer down the sorted_molecule array.
  */
  molecule += 1; /* Caution address arithmetic. */
  umolecules_next  = molecule;
  for (i=1;i<number_molecules;i++) {
    curr_cmpt_index = molecule->c_index;
    molecule_index =  molecule->m_index;
    curr_molecule_name = NULL;
    if (molecule->string >= 0) {
      curr_molecule_name = (char *)&molecules_text[molecule->string];
    }
    if ((curr_cmpt_index != prev_cmpt_index)  ||
	(strcmp(curr_molecule_name,prev_molecule_name) != 0)) {
      prev_molecule_name = curr_molecule_name;
      prev_cmpt_index    = curr_cmpt_index;
      m_size = strlen(curr_molecule_name) + 1;
      pad_size = (align_len - (m_size & align_mask)) & align_mask;
      molecule_len += (int64_t)(m_size + pad_size);
      nu_molecules += 1;

      umolecules_next->string = molecule->string;
      umolecules_next->m_index = nu_molecules;
      umolecules_next->c_index = curr_cmpt_index;
      if (strcmp(solvent_string,curr_molecule_name) == 0) {
	umolecules_next->solvent = 1;
	if (curr_cmpt_index == 0) {
	  *solvent_pos = nu_molecules;
	}
      } else {
	umolecules_next->solvent = 0;
      }
      umolecules_next += 1; /* Caution address arithmetic. */
    } 
    /*
      Set the position dependent molecules_indices and 
      compartment_indices fields of the reaction matrix.
    */
    molecules_indices[molecule_index] = nu_molecules;
    compartment_indices[molecule_index] = curr_cmpt_index;
    /*
      Advance the molecule pointer down the sorted_molecule array.
    */
    molecule += 1; /* Caution address arithmetic. */
  }
  *nunique_molecules = nu_molecules + 1;
  *sum_molecule_len = molecule_len;
  return(success);
}
