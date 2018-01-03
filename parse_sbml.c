/* parse_sbml.c
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

#include "sbml_find_section.h"
#include "sbml_process_list_of_compartments.h"
#include "sbml_process_list_of_species.h"
#include "sbml_sort_species_trans.h"
#include "sbml_process_list_of_reactions.h"

#include "parse_sbml.h"

int parse_sbml(struct sbml2bo_struct *state) {
  /*
    Parse an sbml file generating concs.in, rxnx.dat and cmpts.dat files.
    Called by: sbml2bo, sbml_to_boltzmann
    Calls: sbml_find_section,
           sbml_proces_list_of_compartments,
           sbml_proces_list_of_species,
           sbml_proces_list_of_reactions,
	   calloc, fopen, fclose, fflush
  */
  char sbml_buffer_c[2048];
  char *sbml_buffer;
  char **spec_ids;
  char **translations;
  char **sort_species_trans_scratch;

  int64_t ask_for;
  int64_t one_l;

  int in_list_of_compartments;
  int in_list_of_species;

  int in_list_of_reactions;
  int success;

  int in_species_tag;
  int sbml_buffer_len;

  int num_species;
  int padi;

  FILE *lfp;
  FILE *sbml_fp;
  FILE *concs_fp;
  FILE *cmpts_fp;
  FILE *rxns_fp;
  FILE *id_name_fp;
  FILE *error_fp;

  /*
    Allocation and local pointers.
  */
  sbml_buffer = (char *)&(sbml_buffer_c[0]);
  sbml_buffer_len = 2048;
  num_species = state->num_species;
  lfp         = state->log_fp;
  if (lfp == NULL) {
    error_fp = stderr;
  } else {
    error_fp = lfp;
  }
  success     = 1;
  one_l   = (int64_t)1;
  if (success) {
    sbml_fp = fopen(state->sbml_file,"r");
    if (sbml_fp == NULL) {
      fprintf(error_fp,"parse_sbml: unable to open sbml file, %s\n",
	      state->sbml_file);
      fflush(error_fp);
      success = 0;
    }
  }
  if (success) {
    state->sbml_fp = sbml_fp;
    concs_fp = fopen(state->concs_in_file,"w");
    if (concs_fp == NULL) {
      fprintf(error_fp,"parse_sbml: unable to open concs.in file, %s\n",
	      state->concs_in_file);
      fflush(error_fp);
      success = 0;
    }
  }
  if (success) {
    state->concs_fp = concs_fp;
    cmpts_fp = fopen(state->cmpts_dat_file,"w");
    if (cmpts_fp == NULL) {
      fprintf(error_fp,"parse_sbml: unable to open cmpts.dat file, %s\n",
	      state->cmpts_dat_file);
      fflush(error_fp);
      success = 0;
    }
  }
  if (success) {
    state->cmpts_fp = cmpts_fp;
    rxns_fp = fopen(state->rxns_dat_file,"w");
    if (rxns_fp == NULL) {
      fprintf(error_fp,"parse_sbml: unable to open rxns.dat file, %s\n",
	      state->rxns_dat_file);
      fflush(error_fp);
      success = 0;
    }
  }
  if (success) {
    state->rxns_fp = rxns_fp;
    id_name_fp = fopen(state->id_name_file,"w");
    if (id_name_fp == NULL) {
      fprintf(error_fp,"parse_sbml: unable to open _names.lis file, %s\n",
	      state->id_name_file);
      fflush(error_fp);
      success = 0;
    }
  }
  if (success) {
    state->id_name_fp = id_name_fp;
    success = sbml_find_section(sbml_fp,sbml_buffer,
				sbml_buffer_len,
				"<listOfCompartments>");
  }
  if (success) {
    success = sbml_process_list_of_compartments(sbml_fp,sbml_buffer,
					   sbml_buffer_len,state);
  } 
  if (success) {
    fclose(cmpts_fp);
    success = sbml_find_section(sbml_fp,sbml_buffer,
				sbml_buffer_len,"<listOfSpecies>");
  }
  if (success) {
    success = sbml_process_list_of_species(sbml_fp,sbml_buffer,
					  sbml_buffer_len,state);
  }
  if (success) {
    /*
      Now we need to sort the species id, carrying along their
      json_ids with them.
      sbml_process_list_of_species generated a list of pairs of 
      string pointers in the species_trans array. The i'th species
      id is pointed to by species_trans[i+i] and the i'th species
      transaltion is pointed to by species_trans[i+i+1].
      The strings themselves are stored in specid_2_json_strings array,
      also set in sbml_process_list_of_species.
    */
    spec_ids = state->spec_ids;
    translations = state->translations;
    sort_species_trans_scratch = state->sort_species_trans_scratch;
    success = sbml_sort_species_trans(num_species,
				      spec_ids,translations,
				      sort_species_trans_scratch);

  }
  if (success) {
    fclose(id_name_fp);
    fclose(concs_fp);
    success = sbml_find_section(sbml_fp,sbml_buffer,
				sbml_buffer_len,"<listOfReactions>");
  }
  if (success) {
    success = sbml_process_list_of_reactions(sbml_fp,sbml_buffer,
					  sbml_buffer_len,state);
  }
  if (success) {
    fclose(rxns_fp);
  }
  return(success);
}
