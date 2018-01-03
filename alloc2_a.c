/* alloc2_a.c
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

#include "alloc2_a.h"
int alloc2_a(struct state_struct *state, int setup) {
  /*
    Allocate space for the text strings to store all of the
    information within the reactions file.

    if Setup = 2, raw_molecules_text is not allocated.

      reaction_titles_length,
      pathway_text_length,
      compartment_text_length,
      rxn_title_text,
      molecule_text_length,
      regulation_text_length,
      rxn_title_text,
      pathway_text,
      compartment_text,
      regulation_text,
      molecules_text,
      raw_molecules_text,

    Called by: alloc2,
    Calls:     calloc, fprintf, fflush (intrinsic)
  */
  int64_t usage;
  int64_t reaction_titles_length;
  int64_t pathway_text_length;
  int64_t compartment_text_length;
  int64_t molecule_text_length;
  int64_t regulation_text_length;
  int64_t align_len;
  int64_t align_mask;
  int64_t ask_for;
  int64_t one_l;
  int64_t nze;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *regulation_text;
  char *molecules_text;
  char *raw_molecules_text;
  int64_t num_rxns;
  int64_t num_molecules;
  int64_t num_cmpts;
  int64_t num_regulations;
  int64_t max_regs_per_rxn;
  int64_t double_size;
  int64_t int64_t_size;
  int success;
  int padi;
  success = 1;
  one_l              = (int64_t)1;
  align_mask 	     = state->align_mask;
  align_len  	     = state->align_len;
  num_rxns   	     = state->number_reactions;
  num_molecules      = state->number_molecules;
  num_cmpts          = state->number_compartments;
  max_regs_per_rxn   = state->max_regs_per_rxn;
  usage              = state->usage;
  num_regulations    = max_regs_per_rxn * num_rxns;
  if (setup == 1) {
    reaction_titles_length    =  state->reaction_titles_length + num_rxns * align_len;
    reaction_titles_length    += ((align_len - (reaction_titles_length & align_mask)) & align_mask);
    pathway_text_length      =  state->pathway_text_length + (num_rxns * align_len);
    pathway_text_length      += (align_len - (pathway_text_length & align_mask));
    if (state->number_compartments > 0) {
      compartment_text_length = state->compartment_text_length +
	(num_rxns * align_len);
      compartment_text_length += ((align_len - (compartment_text_length & align_mask)) & align_mask);
    } else {
      compartment_text_length = 0;
    }
    molecule_text_length  = state->molecule_text_length + (num_molecules * align_len);
    molecule_text_length  += ((align_len - (molecule_text_length & align_mask)) & align_mask);
    regulation_text_length  = state->regulation_text_length + (num_regulations * align_len);
    regulation_text_length  += (align_len - (regulation_text_length & align_mask));
    ask_for = reaction_titles_length + pathway_text_length + 
      compartment_text_length + molecule_text_length + 
      molecule_text_length + regulation_text_length;
    /* end if setup == 1 */
  } else {
    reaction_titles_length   =  state->reaction_titles_length;
    pathway_text_length      =  state->pathway_text_length;
    compartment_text_length  =  state->compartment_text_length;
    molecule_text_length     =  state->molecule_text_length;
    regulation_text_length   = state->regulation_text_length;
    ask_for = reaction_titles_length + pathway_text_length + 
      compartment_text_length + molecule_text_length + 
      regulation_text_length;
  }
  usage += ask_for;
  rxn_title_text = (char *) calloc(one_l,ask_for);
  if (rxn_title_text) {
    if (setup == 1) {
      state->reaction_titles_length   = reaction_titles_length;
      state->pathway_text_length      = pathway_text_length;
      state->compartment_text_length  = compartment_text_length;
      state->molecule_text_length     = molecule_text_length;
      state->regulation_text_length   = regulation_text_length;
    }
    /*
      Caution address arithmetic follows.
    */
    pathway_text              = rxn_title_text + reaction_titles_length;
    compartment_text          = pathway_text + pathway_text_length;
    regulation_text           = compartment_text + compartment_text_length;
    molecules_text     	      = regulation_text + regulation_text_length;
    state->rxn_title_text     = rxn_title_text;
    state->pathway_text       = pathway_text;
    state->compartment_text   = compartment_text;
    state->regulation_text    = regulation_text;
    state->molecules_text     = molecules_text;
    if (setup == 1) {
      raw_molecules_text      = molecules_text + molecule_text_length;
      state->raw_molecules_text = raw_molecules_text;
    }
  } else {
    fprintf(stderr,"alloc2_a: Error, unable to allocate %lld bytes of space "
	    "for text strings in core.\n",ask_for);
    fflush(stderr);
    success = 0;
  }
  state->usage = usage;
  return(success);
}
