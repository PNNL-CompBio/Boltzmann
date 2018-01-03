/* boltzmann_boot.c
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
#include "alloc0.h"
#include "read_params.h"
#include "open_output_files.h"
#include "vgrng_init.h"
#include "echo_params.h"
#include "size_rxns_list.h"
#include "size_rxns_file.h"
#include "alloc2.h"
#include "parse_reactions_file.h"
#include "echo_reactions_file.h"
#include "sort_compartments.h"
#include "unique_compartments.h"
#include "unique_compartments_core.h"
#include "translate_compartments.h"
#include "sort_molecules.h"
#include "unique_molecules.h"
#include "unique_molecules_core.h"
#include "print_molecules_dictionary.h"
#include "alloc3.h"
#include "set_compartment_ptrs.h"
#include "read_initial_concentrations.h"
/*
#include "form_molecules_matrix.h"
*/
#include "compute_standard_energies.h"
#include "compute_ke.h"
#include "compute_kss.h"
#include "zero_solvent_coefficients.h"
#include "print_rxn_likelihoods_header.h"
#include "print_free_energy_header.h"
#include "flatten_state.h"
#include "free_boot_state2.h"
#include "sort_global_compartments.h"
#include "sort_global_molecules.h"
#include "count_ws.h"
#include "count_nws.h"
#include "boltzmann_mmap_superstate.h"
#include "flatten_super_state.h"
/*
#define DBG_BOLTZMANN_BOOT  
*/
#include "boltzmann_boot.h"
int boltzmann_boot(char *param_file_name, 
		   struct super_state_struct **super_statep) {
  /*
    Initialize the reactions and data structures for boltzmann.
    This is a variation on boltzmann_init in that it processes
    multiple input reaction files, and builds a single output file of the
    form 
    meta_data[0] 8-byte number of reaction files.
    meta_data[1] 8-byte total  length in bytes.
    meta_data[2] 8-byte page_size in bytes, default val is 4K = 4096
    meta_data[3] 8-byte number of pages = (total_length >> log2(page_size)) + 
                 (page_size - (total_length & (page_size-1))) & (page_size-1))
    meta_data[4] 8-byte align_len for strings, default is 16.
    meta_data[5] 8-byte align_mask             default is 15.
    meta_data[6] 8-byte global number of molecules.
    meta_data[7] 8-byte unique_global_molecules
    meta_data[8] 8-byte molecule_text_length (in bytes)
    meta_data[9] 8-byte global number of compartments. 
    meta_data[10] 8-byte unique_global_compartments 
    meta_data[11] 8-byte compartment_text_length (in bytes)
    meta_data[12] 8 byte state_offsets_sizes_offset (in 8byte words)
    meta_data[13] 8 byte molecule_map_starts_offset; (in 8byte words)
    meta_data[14] 8-byte molecule_map_offset (in 8byte words)
    meta_data[15] 8-byte molecule_names_offset (in 8byte words)
    meta_data[16] 8-byte compartment_map_offset (in 8byte words)
    meta_data[17] 8-byte compartment_names_offset (in 8byte words)
    meta_data[18] 8-byte molecules_text_offset (in bytes)
    meta_data[19] 8-byte compartments_text_offset (in bytes)

    state_offsets_sizes[number_of_reaction_files]
       offset in bytes, length pairs of the reaction states.

    molecule_map_starts [number_of_reaction_files+1]
      These are indexes to the first element per reaction file 
      for the molecles_map vector.

    molecule_map           [global_number_of_molecules]
    compartment_map        [global_number_of_molecules]

    molecule_names     [unqiue_global_number_of_molecules]
       These will be character offsets relative to the
       beginning of the molecules_text address.
    
    compartment_names  [unique_global_number_of_compartments]
       These will be character offsets relative to the
       beginning of the compartments_text address.

    molecule_text      [molecule_text_length]
    comparmtent_text  [compartment_text_length]
----------
    compartment_list_starts [number_of_reaction_files+1]
      These are indexes to the first element per reaction file 
      for the compartments_list vector.
    compartment_list_starts is only needed as workspace, no need to save it.

    compartment_list   [global_number_of_compartments]
    compartment_list is only needed as workspace, no need to save it.

    Called by: boltzmann
    Calls:     alloc0,
               read_params,
	       open_output_files,
	       vgrng_init,
	       alloc1,
	       echo_params,
	       size_rxns_file,
	       alloc2,
	       parse_reactions_file,
	       echo_reactions_file,
	       sort_compartments,
	       unique_compartments,
	       translate_compartments,
	       sort_molecules,
	       unique_molecules,
	       print_molecules_dictionary,
	       alloc3,
	       read_initial_concentrations,
	       compute_standard_energies,
	       compute_ke,
	       compute_kss,
	       print_rxn_likelihoods_header,
	       print_free_energy_header
  */
  struct super_state_struct sss;
  struct super_state_struct *super_state;
  struct state_struct *boot_state;
  struct state_struct *state;
  struct state_struct *statep;
  struct state_struct local_state;
  struct vgrng_state_struct *vgrng_state;
  struct vgrng_state_struct *vgrng2_state;
  struct molecule_struct mes;
  struct molecule_struct *molecules;
  struct molecule_struct *compartments;
  struct molecule_struct *compartment_sort_ws;
  struct molecule_struct *molecule_sort_ws;
  struct molecule_struct *sorted_compartments;
  struct molecule_struct *unsorted_compartments;
  struct molecule_struct *sorted_molecules;
  struct molecule_struct *unsorted_molecules;
  struct formation_energy_struct *formation_energies;

  double *dg0s;
  double *free_energy;
  double *activities;
  int64_t *molecule_names;
  int64_t *compartment_names;
  char   *rxn_file_list;
  char   *rp;
  char   *global_state_filename;
  char   *molecule_text;
  char   *compartment_text;
  char   *molecule_text_ws;
  char   *compartment_text_ws;
  char    *cstring;
  char    *mstring;
  int64_t *super_statev;
  int64_t *page_fill;
  int64_t *meta_data;
  int64_t *molecule_map;
  int64_t *compartment_map;
  int64_t *state_offsets_sizes;
  int64_t *state_offset_size;
  int64_t *molecule_map_starts;
  /*
  int64_t *compartment_map_indices;
  */
  int64_t *compartment_list_starts;
  int64_t *compartment_list;
  int64_t *global_compartment;
  int64_t *global_molecule;
  int64_t state_offsets_sizes_offset;
  
  int64_t molecule_map_starts_offset;
  int64_t compartment_indices_offset;
  int64_t molecule_map_offset;
  int64_t compartment_map_offset;
  int64_t molecule_names_offset;
  int64_t compartment_names_offset;
  int64_t molecules_text_offset;
  int64_t compartments_text_offset;
  int64_t meta_size;
  int64_t map_size;
  int64_t meta_data_size;
  int64_t tmp_offset;
  int64_t fill_size;
  int64_t cur_size;
  int64_t nw;
  int64_t nr;
  int64_t vgrng_start;
  int64_t i;
  int64_t num_reaction_files;
  int64_t total_length;
  int64_t page_size;
  int64_t page_mask;
  int64_t number_of_pages;
  int64_t align_len;
  int64_t align_mask;
  int64_t global_number_of_molecules;
  int64_t unique_global_molecules;
  int64_t molecule_text_length;
  int64_t global_number_of_compartments;
  int64_t unique_global_compartments;
  int64_t compartment_text_length;
  int64_t rxn_list_line_len;
  int64_t rxn_file_name_len;
  int64_t concs_file_name_len;
  int64_t one_l;
  int64_t nunique_molecules;
  int64_t nunique_compartments;
  int64_t unique_molecules_max;
  int64_t unique_compartments_max;
  int64_t local_state_size;
  int64_t mws_pos;
  int64_t cws_pos;
  int64_t lseek_pos;
  int64_t state_size;
  int64_t offset;
  int64_t csi;
  int64_t msi;
  int64_t compartment_space;
  int64_t molecule_space;
  int64_t pad_size;
  int64_t cmpt_size;
  int64_t mlcl_size;
  int64_t local_state_buff_size_in_pages;
  int64_t ntr;
  int64_t ntw;
  int64_t mmap_file_len;
  int64_t ask_for;
  int64_t mi_base;
  int64_t ci_base;
  int64_t minimum_state_size;
  int64_t maximum_state_size;
  char *local_state_buffer;
  char *rxn_file_name;
  char *concs_file_name;
  char rxn_list_buffer[1024];
  char *rxn_list_buffp;
  char *solvent_string;
  char global_state_filename_buffer[1024];

  int success;
  int vgrng_start_steps;

  int print_output;
  int nzr;

  int j;
  int k;

  int global_state_fd;
  int tmp_state_fd;

  int open_flags;
  int whence;

  int log2_int64_t_size;
  int log2_page_size;

  int ierr;
  int solvent_pos;
  int padi;

  int skip1;
  int skip2;

  FILE *tmp_state_fp;
  FILE *global_state_fp;
  FILE *rxn_list_fp;
  FILE *lfp;
  log2_int64_t_size = 3;
  page_size = (int64_t)4096;
  log2_page_size = 12;
  one_l      = (int64_t)1;
  page_mask  = page_size - one_l;
  maximum_state_size = 0;
  minimum_state_size = 0;
  open_flags = O_RDWR | O_CREAT;
  whence     = SEEK_SET;
  formation_energies = NULL;
  /*
    allocate space for the state struct.
    Allocate space for the reactions line buffer, and the rxn_file keywords.
  */
  success = alloc0(&boot_state);
  if (success) {
    /*
      Read the input parameters file.
    */
    state = boot_state;
    success = read_params(param_file_name,state);
  }
  if (success) {
    print_output = (int)state->print_output;
    if (print_output) {
      success = open_output_files(state);
    }
  }
  if (success) {
    /*
      Echo the input parameters to the log file.
    */
    if (print_output) {
      success = echo_params(state->lfp,state);
    }
  }
  if (success) {
    num_reaction_files = size_rxns_list(state);
    if (num_reaction_files < 1) {
      success = 0;
    }
  }
  /*
    Need to open a temporary file, to store the states.
    Need a page fill buffer of length page_size.
    Need to build a file name for the global state file.
  */
  if (success) {
    rxn_file_list = state->reaction_file;
    tmp_state_fp = fopen("./tmp_global_state.state","w");
    if (tmp_state_fp == NULL) {
      fprintf(stderr,"boltzman_boot: Error, could not open temporary state "
	      "file tmp_global_state.state\n");
      fflush(stderr);
      success = 0;
    }
    tmp_state_fd = fileno(tmp_state_fp);
  }
  if (success) {
    global_state_filename = (char *)&global_state_filename_buffer[0];
    strcpy(global_state_filename,rxn_file_list);
    k = strlen(rxn_file_list);
    i = k;
    for (j=k-1;j>=0;j--) {
      if (global_state_filename[j] == '.') {
	i = j;
	break;
      } 
    }
    if (i > 1017) {
      i = 1017;
    }
    strcpy((char*)&global_state_filename[i],".state");
  }
  if (success) {
    global_state_fp = fopen(global_state_filename,"w");
    if (global_state_fp == NULL) {
      fprintf(stderr,"boltzmann_boot: Error, could not open global state "
	      "file %s\n",global_state_filename);
      fflush(stderr);
      success = 0;
    }
    global_state_fd = fileno(global_state_fp);
  }
  if (success) {
    page_fill = (int64_t *)calloc(one_l,page_size);
    if (page_fill) {
      for (i = 0;i<(page_size>>3);i++) {
        page_fill[i] = 0;
      }
    } else {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld bytes for"
              " page_fill\n",page_size);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    rxn_list_fp       = fopen(rxn_file_list,"r");
    if (rxn_list_fp == NULL) {
      success = 0;
      fprintf(stderr,"boltzmann_boot: Error, could not open "
	      "reaction_file_list, %s\n",state->reaction_file);
      fflush(stderr);
    } 
  }
  if (success) {
    ask_for = (int64_t)(num_reaction_files+1) * (int64_t)sizeof(int64_t);
    compartment_list_starts = (int64_t *)calloc(one_l,ask_for);
    if (compartment_list_starts == NULL) {
      fprintf(stderr,"boltzmann_boot: Error unable to allocate %ld bytes "
	      "for compartment_list_starts\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    meta_data_size = ((sizeof(sss) + sizeof(int64_t)) - 1) >> log2_int64_t_size;    
    state_offsets_sizes_offset = meta_data_size;
    molecule_map_starts_offset    = state_offsets_sizes_offset + 
                                   (num_reaction_files << 1);
    molecule_map_offset           = molecule_map_starts_offset + 
		                     num_reaction_files + one_l;
    
    meta_data_size = molecule_map_offset;
    ask_for = sizeof(int64_t) * meta_data_size;
    meta_data = (int64_t *)calloc(one_l,ask_for);
    if (meta_data == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld bytes for meta_data\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    state_offsets_sizes        = (int64_t *)&meta_data[state_offsets_sizes_offset];
    molecule_map_starts_offset    = state_offsets_sizes_offset + 
                                   (num_reaction_files << 1);
    molecule_map_starts       = (int64_t *)&meta_data[molecule_map_starts_offset];
    molecule_map_starts[0] = (int64_t)0;
    compartment_list_starts[0] = (int64_t)0;

    rxn_list_buffp = (char *)&rxn_list_buffer[0];
    rxn_list_line_len = 1024;
    tmp_offset = (int64_t)0;
    state_offset_size = state_offsets_sizes;
    molecule_text_length = 0;
    compartment_text_length = 0;
    unique_molecules_max = 0;
    unique_compartments_max = 0;
    /*
      Loop over the reaction files.
    */
    for (i=0;i<num_reaction_files;i++) {
      /*rp = fgets(rxn_file_name,rxn_list_line_len,rxn_list_fp);*/
      rp = fgets(rxn_list_buffp,rxn_list_line_len,rxn_list_fp);
      if (rp == NULL) {
	fprintf(stderr,"boltzmann_boot: Error reading %ld'th line of rxns_list file\n",
		i);
	fflush(stderr);
      } else {
	/* 
	   The rest of this loop
	*/
	skip1 = count_ws(rxn_list_buffp);
	rxn_file_name = &rxn_list_buffp[skip1];
	rxn_file_name_len = count_nws(rxn_file_name);
	rxn_file_name[rxn_file_name_len] = '\0';
	skip2 = count_ws(&rxn_file_name[rxn_file_name_len]);
	concs_file_name = &rxn_file_name[rxn_file_name_len+skip2];
	concs_file_name_len = count_nws(concs_file_name);
	if (concs_file_name_len > 0) {
	  concs_file_name[concs_file_name_len] = '\0';
	  strcpy(state->init_conc_file,concs_file_name);
	  success = size_rxns_file(state,rxn_file_name);
	}
	/*
	  At this point in time we can compute how much space the aligned
	  reaction titles, pathway descriptions, compartments and molecules
	  verbage will take.
	  Also we have an upperbound of the number of 
	  molecules, state->number_molecules, and we can allocate space for 
	  the molecules sorting.
	*/
	if (success) {
	  success = alloc2(state);
	}
	/*
	  Read reactions file to count molecules and reactions
	  Then allocate space then read for real.
	*/
	if (success) {
#ifdef DBG_BOLTZMANN_BOOT
	  lfp = state->lfp;
	  if (lfp) {
	    fprintf(lfp,"boltzmann_boot: %s\nboltzmann_boot: "
		    "rxns = %ld, cmpts = %ld, molecules = %ld, "
		    "reaction_file_length = %ld\n",
		    rxn_list_bufferp,
		    state->number_reactions,state->number_compartments,
		    state->number_molecules,
		    state->reaction_file_length);
	    fprintf(lfp,"boltzmann_boot: %s\nboltzmann_boot:molecules_len"
		    " = %ld, reaction_title_len = %ld, "
		    "pathway_len = %ld, compartment_len = %ld\n",
		    rxn_list_bufferp,
		    state->molecules_len,
		    state->rxn_title_len,
		    state->pathway_len,
		    state->compartment_len);
	    fflush(lfp);
	  }
#endif
	}
	if (success) {
	  /*
	    Read and parse the reactions file.
	    NB. parse_reactions_file changes the state->reaction_file
	    to be rxn_file instead of the original reaction file list name.
	  */
	  success = parse_reactions_file(state,rxn_file_name);
	}
	/*
	  First we need to sort the compartments.
	*/
	if (success) {
	  success = sort_compartments(state->unsorted_cmpts,
				      state->sorted_cmpts,
				      state->compartment_text,
				      state->number_compartments);
	}
	/*
	  Then we extract the unique compartments.
	  and fill the compartment_indices vector of the 
	  rxns_matrix structure.
	*/
	if (success) {
	  success = unique_compartments(state);
	}
	/*
	  Now we need to assign the proper compartment numbers to the 
	  unsorted molecules, using the compartment_indices field
	  of the rxns_matrix structure.
	*/
	if (success) {
	  success = translate_compartments(state);
	}
	/*
	  Now we need to sort the molecules, by compartment and name.
	*/
	if (success) {
	  success = sort_molecules(state->unsorted_molecules,
				   state->sorted_molecules,
				   state->molecules_text,
				   state->number_molecules);
	}
	/*
	  Then we extract the unique molecules and set the
	  molecules_indices field of the rxn_matrix struct.
	*/
	if (success) {
	  success = unique_molecules(state);
	}
	/*
	  Print the molecules dictionary and the header lines for 
	  the counts output file.
	*/
	if (success) {
	  if (print_output) {
	    success = print_molecules_dictionary(state);
	  }
	}
	/*
	  Now we need to allocate space for the counts,
	  and read in the intial concentrations.
	*/
	if (success) {
	  success = alloc3(state);
	}
	/*
	  Now in order to enable molecule/compartment lookup for
	  read_initial_concentrations we need to set the compartment 
	  pointers, these are pointers into the list of sorted molecules
	  All the molecules within a compartment are adjacent in the
	  sorted molecules list. This call just sets pointers to the
	  first molecule in each compartment, and one past the last
	  molecule so that molecules in compartment i in the sorted 
	  molecules list are in positions 
	  compartment_ptrs[i]:compartment_ptrs[i+1]-1 inclusive.
	*/
	if (success) {
	  success = set_compartment_ptrs(state);
	}
	/*
	  Read initial concentrations, convert them to counts,
	  and print them to the counts output file.
	*/
	if (success) {
	  success = read_initial_concentrations(state);
	}
	/*
	  Compute the molecules matrix.
	  if (success) {
	  success = form_molecules_matrix(state);
	  }
	*/
	/*
	  Compute the reaction energies of formation if called for.
	*/
	if (success) {
	  if (state->use_pseudoisomers) {
	    success = compute_standard_energies(state,&formation_energies);
	  }
	}
	if (success) {
	  if (print_output) {
	    /*
	      Echo the reactions to the log file.
	    */
	    if (state->rxn_fp) {
	      success = echo_reactions_file(state,state->rxn_fp);
	    }
	  }
	}
	/*
	  Compute the reaction ke's.
	*/
	if (success) {
	  success = compute_ke(state);
	}
	/*
	  Compute the reaction kss's.
	*/
	if (success) {
	  success = compute_kss(state);
	}
	/*
	  Here we won't print the dg0 and ke values as they go to a
	  hardwired output file name and would overwrite eachother.
	if (success) {
	  if (print_output) {
	    success = print_dg0_ke(state);
	  }
	}
	  At this juncture we have echoed the reactions file if requested and
	  need to zero out the coefficients in the reaction matrix that
	  correspond to the solvent molecule (by default H2O) so as not to
	  have it influence the computation of likelihoods, nor change
	  counts (see rxn_likelihood.c and comment in rxn_count_update.c)
	*/
	if (success) {
	  success = zero_solvent_coefficients(state);
	}
	/*
	  If use_activities has not been turned on, set all activities to
	  1.0 so that all reactions are fully active.
	*/
	if (success) {
	  activities = state->activities;
	  if (state->use_activities == 0) {
	    for (j=0;j<state->number_reactions;j++) {
	      activities[j] = 1.0;
	    }
	  }
	}
	if (success) {
	  /*
	    Initialize the random number generators,
	    setting the vgrng_state and vgrng2_state fields of
	    the state structure.
	  */
	  vgrng_state = state->vgrng_state;
	  vgrng_start_steps = 1001;
	  vgrng_start= vgrng_init(vgrng_state,vgrng_start_steps);
	  vgrng2_state = state->vgrng2_state;
	  vgrng_start_steps = 1042;
	  vgrng_start= vgrng_init(vgrng2_state,vgrng_start_steps);
	}
	if (success) {
	  /*
	   *statep = boot_state;
	  */
	  statep = NULL;
	  boot_state->workspace_base = NULL;
	  success = flatten_state(boot_state,&statep);
	  if (success) {
	    success = free_boot_state2(&boot_state);
	  }
	}
	/*
	  Now we need to write the statep block out to a file
	  tracking where we put it for global dictionary determination
	  later.  This stuff needs to be in a routine.
	*/
	molecule_text_length    += statep->molecules_space;
	compartment_text_length += statep->compartment_space;
	cur_size = statep->state_length;
	if (cur_size > maximum_state_size) {
	  maximum_state_size = cur_size;
	}
	if ((minimum_state_size == 0) || (cur_size < minimum_state_size)) {
	  minimum_state_size = cur_size;
	}
	fill_size = (page_size - (cur_size & page_mask)) & page_mask;
	nw = write(tmp_state_fd,statep,cur_size);
	if (nw != cur_size) {
	  fprintf(stderr,"boltzmann_boot: Error writing data for %s state\n",
		  rxn_file_name);
	  fflush(stderr);
	  success = 0;
	}
	if (success) {
	  if (fill_size > 0) {
	    nw = write(tmp_state_fd,(char *)page_fill,fill_size);
	    if (nw != fill_size) {
	      fprintf(stderr,"boltzmann_boot: Error writing data for %s state\n",
		      rxn_file_name);
	      fflush(stderr);
	      success = 0;
	    }
	  }
	}
	if (success) {
	  *state_offset_size = tmp_offset;
	  state_offset_size += 1;
	  *state_offset_size = cur_size;
	  state_offset_size += 1;
	  nunique_molecules = statep->nunique_molecules;
	  nunique_compartments = statep->nunique_compartments;
	  molecule_map_starts[i+1]    = molecule_map_starts[i] + 
	                                   nunique_molecules;
	  compartment_list_starts[i+1] = compartment_list_starts[i] + 
	                                   nunique_compartments;
	  if (nunique_molecules > unique_molecules_max) {
	    unique_molecules_max = nunique_molecules;
	  }
	  if (nunique_compartments > unique_compartments_max) {
	    unique_compartments_max = nunique_compartments;
	  }
	  tmp_offset += cur_size + fill_size;
	}
	
      } /* end else we could open reaction list file */
    } /* end for (i...) */
  }
  if (success) {
    /*
      Need to close and reopen tmp_state_fp to flush the file buffers
      out to the file so that it may be read.
    */
    fclose(tmp_state_fp);
    tmp_state_fp = fopen("./tmp_global_state.state","r");
    if (tmp_state_fp == NULL) {
      fprintf(stderr,"boltzman_boot: Error, could not open temporary state "
	      "file tmp_global_state.state for reading\n");
      fflush(stderr);
      success = 0;
    }
    tmp_state_fd = fileno(tmp_state_fp);
  }
  if (success) {
    global_number_of_compartments = compartment_list_starts[num_reaction_files];
    global_number_of_molecules = molecule_map_starts[num_reaction_files];

    /*
      Now to merge the compartment lists into one global list, and then the 
      molecules lists into one global list.
      To do this we need to sum the compartment_space and molecule_space
      arguments for each workspace, and allocate one large block of storage
      for each of those, copy the strings in with
      reaction_file_number, compartment_number string_offset triples for
      the compartments and 
      reaction_file_number, molecule_number, compartment_number,
      string_offset quadruple for compartments.
      Then we need to sort the compartments building the compartment map
      of length number total_number_compartments indexd by 
      compartment_list_starts, followed by sorting the molecules in
      which we also fill in the global compartment numbers from
      the compartment map.
    	we need to parse the sorted_molecules and sorted_compartment
    	molecule_struct components of the local states, so we
    	need to know their maximum size: unique_molecules_max
    	 and unique_compartments_max.

      So first we need to allocate workspace to store the 
      molecules_text and compartment_text pieces as well as
      double the number of quadruples and triples for each.
    */
    molecule_text_ws = (char*)calloc(one_l,molecule_text_length);
    if (molecule_text_ws == NULL) {
    	fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
    		"bytes for molecule_text_ws\n",molecule_text_length);
    	fflush(stderr);
    	success = 0;
    }
  }
  if (success) {
    compartment_text_ws = (char*)calloc(one_l,compartment_text_length);
    if (compartment_text_ws == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for compartment_text_ws\n",compartment_text_length);
      fflush(stderr);
      success = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(2*sizeof(mes)) * 
      global_number_of_molecules;
    molecule_sort_ws = (struct molecule_struct *)calloc(one_l,ask_for);
    if (molecule_sort_ws == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for molecule_sort_ws\n",ask_for);
      fflush(stderr);
      success  = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(2*sizeof(mes)) * 
      global_number_of_compartments;
    compartment_sort_ws = (struct molecule_struct*)calloc(one_l,ask_for);
    if (compartment_sort_ws == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for compartment_sort_ws\n",ask_for);
      fflush(stderr);
      success  = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(sizeof(int64_t)*global_number_of_compartments);
    compartment_list = (int64_t*)calloc(one_l,ask_for);
    if (compartment_list == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for compartment_list\n",ask_for);
      fflush(stderr);
      success  = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(sizeof(int64_t)*global_number_of_molecules);
    compartment_map = (int64_t*)calloc(one_l,ask_for);
    if (compartment_map == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for compartment_map\n",ask_for);
      fflush(stderr);
      success  = 0;
    }
  }	    
  if (success) {
    ask_for = (int64_t)(sizeof(int64_t)*global_number_of_molecules);
    molecule_map = (int64_t*)calloc(one_l,ask_for);
    if (molecule_map == NULL) {
      fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
	      "bytes for molecule_map\n",ask_for);
      fflush(stderr);
      success  = 0;
    }
  }	    
  /*
    Now loop over the local reaction workspaces building
    compartment_text_ws, and molecule_text_ws strings
    as well as the first half of the sort workspaces.
  */
  if (success) {
    state_offset_size = state_offsets_sizes;
    state_size = (int64_t)sizeof(local_state);
    mws_pos = (int64_t)0;
    cws_pos = (int64_t)0;
    for (i=0;((i<num_reaction_files) && success);i++) {
      /*
	Seek to the position of the metadata for reaction file i.
      */
      offset = *state_offset_size;
      lseek(tmp_state_fd,offset,whence);
      /*
	get the local pointers.
      */
      nr = read(tmp_state_fd,&local_state,state_size);
      if (nr != state_size) {
	fprintf(stderr,"boltzmann_boot: Error i = %ld, could not read "
		"%ld bytes for a local_state\n",i,state_size);
	fflush(stderr);
	/*print some error message.*/
	success =0;
      }
      if (success) {
	/*
	  read in the compartment text 
	*/
	compartment_space = local_state.compartment_space;
	lseek_pos = offset + local_state.compartments_text_offset_in_bytes;
	lseek(tmp_state_fd,lseek_pos,whence);
	nr = read(tmp_state_fd,(char *)&compartment_text_ws[cws_pos],compartment_space);
	if (nr != compartment_space) {
	  fprintf(stderr,"boltzmann_boot: Error i = %ld, could not read "
		  "%ld bytes for compartment text\n",i,compartment_space);
	  fflush(stderr);
	  success = 0;
	}
      }
      if (success) {
	/*
	  read in molecules text.
	*/
	molecule_space    = local_state.molecules_space;
	lseek_pos = offset + local_state.molecules_text_offset_in_bytes;
	lseek(tmp_state_fd,lseek_pos,whence);
	nr = read(tmp_state_fd,(char *)&molecule_text_ws[mws_pos],
		  molecule_space);
	if (nr != molecule_space) {
	  fprintf(stderr,"boltzmann_boot: Error i = %ld, could not read "
		  "%ld bytes for molecules text\n",i,molecule_space);
	  fflush(stderr);
	  success = 0;
	}
      }
      if (success) {
	/*
	  read in the meta data  for the compartments.
	*/
	nunique_compartments = local_state.nunique_compartments;
	compartment_space = nunique_compartments * sizeof(mes);
	lseek_pos = offset + local_state.sorted_compartments_offset_in_bytes;
	lseek(tmp_state_fd,lseek_pos,whence);
	ci_base = compartment_list_starts[i];
	compartments = (struct molecule_struct *)&compartment_sort_ws[ci_base];
	nr = read(tmp_state_fd,compartments,compartment_space);
	if (nr != compartment_space) {
	  fprintf(stderr,"boltzmann_boot: Error i = %ld, could not read "
		  "%ld bytes for sorted compartments meta data\n",i,
		  compartment_space);
	  fflush(stderr);
	  success = 0;
	} else {
	  /*
	    Set up global_compartment string pointers, from the
	    local_compartment string field.
	  */
	  for (j=0;j<nunique_compartments;j++) {
	    compartments->c_index = ci_base + j;
	    compartments->g_index = i;
	    compartments->string += cws_pos;
	    compartments++; /* Caution address arithmetic. */
	  }
	  cws_pos += compartment_space;
	}
      }
      if (success) {
	/*
	  read in the metadata for the molecules;
	*/
	nunique_molecules = local_state.nunique_molecules;
	molecule_space = nunique_molecules * sizeof(mes);
	lseek_pos = offset + local_state.sorted_molecules_offset_in_bytes;
	lseek(tmp_state_fd,lseek_pos,whence);
	mi_base = molecule_map_starts[i];
	molecules = (struct molecule_struct *)&molecule_sort_ws[mi_base];
	nr = read(tmp_state_fd,
		  (struct molecule_struct *)molecules,molecule_space);
	if (nr != molecule_space) {
	  fprintf(stderr,"boltzmann_boot: Error i = %ld, could not read "
		  "%ld bytes for sorted molecules meta data\n",i,
		  molecule_space);
	  fflush(stderr);
	  success = 0;
	} else {
	  /*
	    Set up global_molecule string pointers, from the
	    local_molecule string field.
	  */
	  for (j=0;j<nunique_molecules;j++) {
	    molecules->m_index = mi_base + j;
	    molecules->string  += mws_pos;
	    molecules->g_index = i;
	    molecules++; /* Caution address arithmetic. */
	  }
	  mws_pos += molecule_space;
	}
      }
      state_offset_size += 2;   /* Caution address arithmetic. */
    } /* end for(i...) */
  }
  /*
    at this juncture 
    compartment_sort_ws is a vector of molecule_struct's of length 
    2*global_number_of_compartments, the first half is filled
    with compartments where the c_index field is the global index 
    of the compartment, the g_index field is the reaction file number
    and the string field is an offset into the compartment_text_ws array.
    
    Similarly
    molecule_sort_ws in a vector of molecule_struct's of length
    2*global_number_of_molecules, the first half is filled with
    with molecules where the m_index field is the global index of the
    molecule, the c_index field is the local compartment number
    within the reaction file, the g_index field is the reaction file 
    number, and the string field is an offset into the molecule_text_ws
    array.
  */
  /*
    So now we sort the compartments, within a reaction file they
    are already sorted so its just a matter of merging these lists.
  */
  if (success) {
    unsorted_compartments = (struct molecule_struct *)&compartment_sort_ws[0];
    sorted_compartments   = (struct molecule_struct *)&compartment_sort_ws[global_number_of_compartments];
    success = sort_global_compartments(&unsorted_compartments,
				       &sorted_compartments,
				       compartment_list_starts,
				       compartment_text_ws,
				       num_reaction_files);
  }
  if (success) {
    /*
      Now we need to extract the list of unique global compartments.
    */
    align_len     = state->align_len;
    align_mask     = state->align_mask;
    nzr            = global_number_of_compartments;
    success = unique_compartments_core(nzr,
				       sorted_compartments,
				       compartment_text_ws,
				       compartment_list,
				       &unique_global_compartments,
				       &compartment_text_length,
				       align_len,align_mask);
  }
  if (success) {
    /*
      Replace the local compartment indices with global compartment
      indices in the unsorted_molecules.
    */
    molecules = (struct molecule_struct *)&molecule_sort_ws[0];
    for (i=0;i<global_number_of_molecules;i++) {
      /*
	compartment_list is a list of all unique compartment names within
	each reaction file by reaction file.
	molecules->g_index is the reaction_file number.
	compartment_list_starts[molecules->g_index] points to where
	in the compartment_list is to start, and molecules->c_index gives
	the offset from the start.
      */
      molecules->c_index = compartment_list[compartment_list_starts[molecules->g_index] + molecules->c_index];
      compartment_map[i] = molecules->c_index;
      molecules++; /* Caution address arithmetic */
    }
  }
  /*
    So now we sort the molecules, within a reaction file they
    are already sorted so its just a matter of merging these lists.
  */
  if (success) {
    unsorted_molecules = (struct molecule_struct *)&molecule_sort_ws[0];
    sorted_molecules   = (struct molecule_struct *)&molecule_sort_ws[global_number_of_molecules];
    success = sort_global_molecules(&unsorted_molecules,
				    &sorted_molecules,
				    molecule_map_starts,
				    molecule_text_ws,
				    num_reaction_files);
  }
  if (success) {
    /*
      Now we need to extract the list of unique global molecules.    */
    align_len     = state->align_len;
    align_mask     = state->align_mask;
    nzr            = global_number_of_molecules;
    solvent_string = state->solvent_string;
    success = unique_molecules_core(nzr,
				    sorted_molecules,
				    molecule_text_ws,
				    solvent_string,
				    molecule_map,
				    &unique_global_molecules,
				    &molecule_text_length,
				    &solvent_pos,
				    align_len,align_mask);
  }
  /*
    Now we need to allocate space for the
    compartment_names and the molecule_names character pointer arrays,
    and the condensed molecule_text and compartment text arrays.
  */
  if (success) {
    molecule_text = (char*)calloc(one_l,molecule_text_length);
    if (molecule_text == NULL) {
    	fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
    		"bytes for molecule_text\n",molecule_text_length);
    	fflush(stderr);
    	success = 0;
    }
  }
  if (success) {
    compartment_text = (char*)calloc(one_l,compartment_text_length);
    if (compartment_text == NULL) {
    	fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
    		"bytes for compartment_text\n",compartment_text_length);
    	fflush(stderr);
    	success = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(sizeof(int64_t) * unique_global_molecules);
    molecule_names = (int64_t*)calloc(one_l,ask_for);
    if (molecule_names == NULL) {
    	fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
    		"bytes for molecule_names\n",ask_for);
    	fflush(stderr);
    	success = 0;
    }
  }
  if (success) {
    ask_for = (int64_t)(sizeof(int64_t) * unique_global_compartments);
    compartment_names = (int64_t*)calloc(one_l,ask_for);
    if (compartment_names == NULL) {
    	fprintf(stderr,"boltzmann_boot: Error could not allocate %ld "
    		"bytes for compartment_names\n",ask_for);
    	fflush(stderr);
    	success = 0;
    }
  }
  if (success) {
    cws_pos = (int64_t)0;
    compartments = sorted_compartments;
    for (i=0;i<unique_global_compartments;i++) {
      cstring   = (char*)&compartment_text_ws[compartments->string];
      cmpt_size = (int64_t)strlen(cstring) + one_l;
      pad_size = (align_len  - (cmpt_size & align_mask)) & align_mask;
      strcpy((char*)&(compartment_text[cws_pos]),cstring);
      compartment_names[i] = cws_pos;
      compartments->string = cws_pos;
      cws_pos += (cmpt_size + pad_size);
      compartments++; /* Caution address arithmetic */
    }
    mws_pos = (int64_t)0;
    molecules = sorted_molecules;
    for (i=0;i<unique_global_molecules;i++) {
      mstring = (char*)&molecule_text_ws[molecules->string];
      mlcl_size = (int64_t)strlen(mstring) + one_l;
      pad_size = (align_len  - (mlcl_size & align_mask)) & align_mask;
      strcpy((char*)&(molecule_text[mws_pos]),mstring);
      molecule_names[i] = mws_pos;
      molecules->string = mws_pos;
      mws_pos += (mlcl_size + pad_size);
      molecules++; /* Caution address arithmetic */
    }
    compartment_text_length = cws_pos;
    molecule_text_length = mws_pos;
    molecule_map_offset = meta_data_size;
    compartment_map_offset = meta_data_size + global_number_of_molecules;
    molecule_names_offset  = compartment_map_offset + global_number_of_molecules;
    compartment_names_offset = molecule_names_offset + unique_global_molecules;
    molecules_text_offset    = (compartment_names_offset + unique_global_compartments)<<log2_int64_t_size;
    compartments_text_offset  = molecules_text_offset + mws_pos;
    meta_data_size = (compartments_text_offset + cws_pos);
    fill_size = (page_size - (meta_data_size & page_mask)) & page_mask;
    meta_data_size += fill_size;

    /*
      Fill the meta data array.
    */
    total_length = meta_data_size + tmp_offset;
    super_state = (struct super_state_struct *)&meta_data[0];
    super_state->number_of_reaction_files = num_reaction_files;
    super_state->total_length_in_bytes = total_length;
    super_state->page_size_in_bytes = page_size;
    super_state->number_of_pages = (total_length >> log2_page_size) +
      ((page_size - (total_length & page_mask)) * page_mask);
    super_state->string_align_len = align_len;
    super_state->string_align_mask = align_len - one_l;
    super_state->global_number_of_molecules = global_number_of_molecules;
    super_state->unique_global_molecules;
    super_state->molecule_text_length_in_bytes = mws_pos;
    super_state->global_number_of_compartments = global_number_of_compartments;
    super_state->unique_global_compartments = unique_global_compartments;
    super_state->compartment_text_length_in_bytes = cws_pos;
    super_state->state_offsets_sizes_offset_in_words = state_offsets_sizes_offset;
    super_state->molecule_map_starts_offset_in_words = molecule_map_starts_offset;
    super_state->molecule_map_offset_in_words      = molecule_map_offset;
    super_state->molecule_names_offset_in_words    = molecule_names_offset;
    super_state->compartment_map_offset_in_words   = compartment_map_offset;
    super_state->compartment_names_offset_in_words = compartment_names_offset;
    super_state->molecules_text_offset_in_bytes    = molecules_text_offset;
    super_state->compartments_text_offset_in_bytes = compartments_text_offset;
    super_state->maximum_state_size_in_bytes       = maximum_state_size;
    super_state->minimum_state_size_in_bytes       = minimum_state_size;
    /*
    meta_data[0]  = num_reaction_files;
    meta_data[1]  = meta_data_size + tmp_offset;
    meta_data[2]  = page_size;
    meta_data[3]  = (meta_data[1] >> log2_page_size) + 
      ((page_size  - (meta_data[1] & page_mask)) & page_mask);
    meta_data[4]  = align_len;
    meta_data[5]  = align_len - one_l;
    meta_data[6]  = global_number_of_molecules;
    meta_data[7]  = unique_global_molecules;
    meta_data[8]  = mws_pos;
    meta_data[9]  = global_number_of_compartments;
    meta_data[10] = unique_global_compartments;
    meta_data[11] = cws_pos;
    meta_data[12] = state_offsets_sizes_offset;
    meta_data[13] = molecule_map_starts_offset;
    meta_data[14] = molecule_map_offset;
    meta_data[15] = molecule_names_offset; 
    meta_data[16] = compartment_map_offset;
    meta_data[17] = compartment_names_offset;
    meta_data[18] = molecules_text_offset; 
    meta_data[19] = compartments_text_offset;
    */
    /*
      Increment the state_offset_size offsets by meta_data_size;
    */
    state_offset_size = state_offsets_sizes;
    for (i=0;i<num_reaction_files;i++) {
      *state_offset_size += meta_data_size;
      state_offset_size += 2; /* Caution address arithmetic */
    }
    ntw = molecule_map_offset << log2_int64_t_size;
    nw = write(global_state_fd,(char*)&meta_data[0],ntw);
    if (nw != ntw) {
      success = 0;
      fprintf(stderr,"boltzmann_boot: Error writing meta data to global state file\n");
      fflush(stderr);
    }
  }
  if (success) {
    ntw = global_number_of_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(char*)&molecule_map[0],ntw);
    if (nw!=ntw) {
      success = 0;
      fprintf(stderr,"boltzmann_boot: Error writing molecule_map to global state file\n");
      fflush(stderr);
    }
  }
  if (success) {
    ntw = global_number_of_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(char*)&compartment_map[0],ntw);
    if (nw!=ntw) {
      success = 0;
      fprintf(stderr,"boltzmann_boot: Error writing compartment_map to global state file\n");
      fflush(stderr);
    }
  }  
  if (success) {
    ntw = unique_global_molecules << log2_int64_t_size;
    nw = write(global_state_fd,(int64_t*)&molecule_names[0],ntw);
    if (nw!=ntw) {
      fprintf(stderr,"boltzmann_boot: Error writing molecule_names to global state file\n");
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    ntw = unique_global_compartments << log2_int64_t_size;
    nw = write(global_state_fd,(int64_t*)&compartment_names[0],ntw);
    if (nw!=ntw) {
      fprintf(stderr,"boltzmann_boot: Error writing compartment_names to global state file\n");
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    ntw = mws_pos;
    nw = write(global_state_fd,molecule_text,ntw);
    if (nw!=ntw) {
      fprintf(stderr,"boltzmann_boot: Error writing molecule_text to global state file\n");
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    ntw = cws_pos;
    nw = write(global_state_fd,compartment_text,ntw);
    if (nw!=ntw) {
      fprintf(stderr,"boltzmann_boot: Error writing compartment_text to global state file\n");
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    ntw = fill_size;
    nw = write(global_state_fd,(char*)&page_fill[0],ntw);
    if (nw!=ntw) {
      fprintf(stderr,"boltzmann_boot: Error writing pagefill to global state file\n");
      fflush(stderr);
      success = 0;
    }
  }  
  if (success) {
    lseek_pos = (int64_t)0;
    lseek(tmp_state_fd,lseek_pos,whence);
    /*
      Now copy the local states into the global state file.
      Do this in largest chunks as possible.
    */
    local_state_buff_size_in_pages = 1024;
    ask_for = page_size * local_state_buff_size_in_pages;
    local_state_buffer = (char *)calloc(one_l,ask_for);
    if (local_state_buffer) {
      for (i=0;i<tmp_offset;i+=ask_for) {
	ntr = ask_for;
	if (i + ntr >= tmp_offset) {
	  ntr = tmp_offset - i;
	}
        nr = read(tmp_state_fd,local_state_buffer,ntr);
	if (nr != ntr) {
	  success = 0;
	  fprintf(stderr,"boltzmann_boot: Error reading temporary local state file\n");
	  fflush(stderr);

	  break;
	}
	nw = write(global_state_fd,local_state_buffer,ntr);
	if (nw != ntr) {
	  success = 0;
	  fprintf(stderr,"boltzmann_boot: Error writing global state file\n");
	  fflush(stderr);
	  break;
	}
      }
      fclose(tmp_state_fp);
    }
  }
  if (success) {
    fclose(global_state_fp);
    mmap_file_len   = meta_data[1];
    success = boltzmann_mmap_superstate(global_state_filename,
					mmap_file_len,
					super_statep);
  }
  return(success);
}
