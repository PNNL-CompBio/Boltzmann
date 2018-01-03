/* boltzmann_flatten_aux_data.c
*******************************************************************************
Boltzmann

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
#include "boltzmann_set_filename_ptrs.h"
#include "boltzmann_flatten_aux_data.h"
int boltzmann_flatten_aux_data(struct state_struct *state, 
			      void **aux_data_p,
			      int direction) {
/*
  This routine flattens the filenames and text fields of the
  boltzmann state struct and their lengths into one contigues 
  string. if direcion = 0, and writes it to the
  file pointed to by state->aux_data_file.
  If direction = 1, it extracts the lengths, allocates a vector for
  the fields and sets the pointers in the state struct and reads
  the string data from the state->aux_data_file.

  This data is not needed by boltzmann run (unless
  printing is turned on: state->print_output = 1)
  but may be helpful in mapping molecules and compartment names.
 
  format of the string in the file:
  
  if direction is 0, a flattened string is returned with all of the
  string data and its lengths: the first 128 characters are
  8 16 character fields storing the asci representation of the
  following lengths:
           total_length_in_bytes,
	   number_of_file_names,
	   filename_length,
	   rxn_titles_length,
	   pathway_text_length,
	   compartment_text_length,
	   molecule_text_length,
	   regulation_text_length

	   The lengths are in bytes.
	   The filename offset is always 128,
	   rxn_titles_offset is 128 + (filename_length * num_files)
	   pathway_text_offset = rxn_titles_offset + rxn_titles_length
	   compartment_text_offset = pathway_text_offset + pathway_text_length,
	   molecule_text_offset    = compartment_text_offset + 
	                               compartment_text_length,
	   regulation_text_offset  = molecule_text_offset + molecule_text_length
	   total_length = regulation_text_offset + regulation_text_length.

  if direction is 1 the string pointers in state are
  set to point to the appropriate places in aux_data.

  Called by: boltzmann_save_aux_data, boltzmann_load_aux_data
  Calls    : boltzmann_set_filename_ptrs,
             memset, memcpy, calloc, fprintf, fflush, sprintf, sscanf

*/
  int64_t filenames_length;
  int64_t reaction_titles_length;
  int64_t pathway_text_length;
  int64_t compartment_text_length;
  int64_t molecule_text_length;
  int64_t regulation_text_length;
  int64_t aux_length;
  int64_t string_length;
  int64_t one_l;

  char *aux_data;
  char *strpos;
  char *filenames;
  char *ffilenames;
  char *fstrings;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *molecules_text;
  char *regulation_text;
  int num_files;
  int filename_length;

  int success;
  int blank;

  int nr;
  int padi;

  FILE *lfp;
  FILE *efp;
  success = 1;
  one_l = (int64_t)1;
  lfp = state->lfp;
  if (direction == 0) {
    num_files               = state->num_files;
    filename_length         = state->max_filename_len;
    reaction_titles_length  = state->reaction_titles_length;
    pathway_text_length     = state->pathway_text_length;
    compartment_text_length = state->compartment_text_length;
    molecule_text_length    = state->molecule_text_length;;
    regulation_text_length  = state->regulation_text_length;  
    aux_length = 128 + (num_files * filename_length) + 
      reaction_titles_length + pathway_text_length + compartment_text_length +
      molecule_text_length + regulation_text_length;
    aux_data = (char *)calloc(one_l,aux_length);
    if (aux_data == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,
		"boltzmann_flatten_aux_data unable to allocate %ld bytes "
		"for auxilliary data strings\n",aux_length);
	fflush(lfp);
      }
    }
    if (success) {
      blank = 32;
      memset(aux_data,blank,aux_length);
      strpos = (char*)&aux_data[0];
      sprintf(strpos,"%ld",aux_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%d",num_files);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%d",filename_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%ld",reaction_titles_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%ld",pathway_text_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%ld",compartment_text_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%ld",molecule_text_length);
      strpos += 16; // Caution address arithmetic.
      sprintf(strpos,"%ld",regulation_text_length);


      filenames      = (char *)&aux_data[128];
      filenames_length = num_files * filename_length;
      memcpy(filenames,state->params_file,filenames_length);
      // Caution address arithmetic.
      rxn_title_text = filenames + filenames_length;
      memcpy(rxn_title_text,state->rxn_title_text,reaction_titles_length);
      // Caution address arithmetic.
      pathway_text = rxn_title_text + reaction_titles_length; 
      memcpy(pathway_text,state->pathway_text,pathway_text_length);
      // Caution address arithmetic.
      compartment_text = pathway_text + pathway_text_length;
      memcpy(compartment_text,state->compartment_text,compartment_text_length);
      // Caution address arithmetic.
      molecules_text = compartment_text + compartment_text_length;
      memcpy(molecules_text,state->molecules_text,molecule_text_length);
      // Caution address arithmetic.
      regulation_text = molecules_text + molecule_text_length;
      memcpy(regulation_text,state->regulation_text,regulation_text_length);
      *aux_data_p = (void *)aux_data;
    }
  } else {
    /* 
       direction is 1.
    */
    aux_data = (char*)*aux_data_p;
    strpos = (char*)&aux_data[0];
    nr = sscanf(strpos,"%ld",&aux_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%d",&num_files);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%d",&filename_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%ld",&reaction_titles_length);    
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%ld",&pathway_text_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%ld",&compartment_text_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%ld",&molecule_text_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    strpos += 16; // Caution address arithmetic.
    nr = sscanf(strpos,"%ld",&regulation_text_length);
    if (nr != 1) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_data error scanning text lengths\n");
	fflush(lfp);
      }
    }
    if (success) {
      state->num_files = num_files;
      state->max_filename_len = filename_length;
      state->reaction_titles_length = reaction_titles_length;
      state->pathway_text_length = pathway_text_length;
      state->compartment_text_length = compartment_text_length;
      state->molecule_text_length = molecule_text_length;
      state->regulation_text_length = regulation_text_length;
      filenames_length = num_files * filename_length;

      /*
	Here we need to allocate separate space for 
	filenames, and the other text struct, so that
	free works.
      */
      ffilenames      = (char *)&aux_data[128];
      filenames      = (char*)calloc(one_l,filenames_length);
      if (filenames == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_flatten_aux_data: Error on a load, unable to allocate %ld bytes for filenames\n",filenames_length);
	  fflush(lfp);
	}
      }
      if(success) { 
	memcpy(filenames,ffilenames,filenames_length);
	state->params_file      = filenames;
	boltzmann_set_filename_ptrs(state);
	
	string_length  = reaction_titles_length + pathway_text_length +
	compartment_text_length + molecule_text_length +
	regulation_text_length;

	rxn_title_text = (char *)calloc(one_l,string_length);
	if (rxn_title_text == NULL) {
	  success = 0;
	  if (lfp) {
	    fprintf(lfp,"boltzmann_flatten_aux_data: Error on a load, unable to allocate %ld bytes for string data\n",string_length);
	    fflush(lfp);
	  }
	}
	if (success) {
	  //* Caution address arithmetic next5 statments.
  fstrings = ffilenames + filenames_length;
	  memcpy(rxn_title_text,fstrings,string_length);


	  pathway_text = rxn_title_text + reaction_titles_length; 
	  compartment_text = pathway_text + pathway_text_length;
	  molecules_text = compartment_text + compartment_text_length;
	  regulation_text = molecules_text + molecule_text_length;
	  state->rxn_title_text   = rxn_title_text;
	  state->pathway_text     = pathway_text;
	  state->compartment_text = compartment_text;
	  state->molecules_text   = molecules_text;
	  state->regulation_text  = regulation_text;
	}
      }
    }
  }
  return(success);
}
    
