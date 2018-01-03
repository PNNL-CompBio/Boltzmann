/* 
  So we will lay out the state struct a little differently.
  We want to flatten the state structure, and make it self describing.
  At the beginning there will be a set of scalars describing sizes
  and constansts to be used in the code, followed by a partitioning

  Let us make the following partitions
  1, for the size parameters and floating point constants.

  2. for arrays/structs that are both inputs and outputs (4 of these)
     and 2 scalars, dg_forward and entropy for the whole system
     current_concentrations,
     bndry_flux_concs,
     vgrng_state
     vgrng2_state


 
  3. for arrays/structs that are strictly inputs (there are 6 of these
     not counting string pointers).
     reactions,
     reactions_matrix,
     sorted_molecules,
     ke,
     activities,
     sorted_compartments


  4. for work space arrays and structs.
     unsorted_molecules,
     unsorted_compartments.
     future_concentrations,
     rxn_fire,
     no_op_likelihood,
     forward_rxn_likelihood,
     reverse_rxn_likelihood,
     forward_rxn_log_likelihood_ratio,
     reverse_rxn_log_likelihood_ratio,
     rxn_likelihood_ps,
     rxn_view_likelihoods,
     rev_rxn_view_likelihoods,
     dg0s,
     free_energy,


  a partition would have a number of entities as an 8 byte first element,
  a total length for the second 8 byte element followed by 
  offset, length pairs for each of the entities. 
  followed by the entities.

  So the super partition would look like 

    4L 
    total_length_in_8_byte_words_of_the_state_struct
    10L
    total_length_in_8byte_words_of_partition_1_size_params_and_fp_constants

    10 + total_length_in_8byte_words_of_partition_1_size_params_and_fp_constants    
    total_length_in_8byte_words_of_partition_2_2way_data

    10 + total_length_in_8byte_words_of_partition_1_size_params_and_fp_constants    +     total_length_in_8byte_words_of_partition_2_2way_data
    total_length_in_8byte_words_of_partition_3_input_data

    10 + total_length_in_8byte_words_of_partition_1_size_params_and_fp_constants    +     total_length_in_8byte_words_of_partition_2_2way_data
    +     total_length_in_8byte_words_of_partition_3_input_data
    total_length_in_8byte_words_of_partition_4_workspace.


    Partition 1 is a little different in that all the quantities are
    scalars, we will make them all 8 byte quantites.
    First quantity will be num_scalars (includes itself).
*/
/* state_struct.h 
*******************************************************************************
Boltsman

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
#ifndef __STATE_STRUCT__
#define __STATE_STRUCT__ 1
struct state_struct {
  int64_t agent_type;
  int64_t state_length;         	
  int64_t thread_id;
  double  x_coord;
  double  y_coord;
  double  z_coord;
  int64_t num_partitions;       	
  int64_t two_way_data_length;  	
  int64_t two_way_data_offset;  	
  int64_t incoming_data_length; 	
  int64_t incoming_data_offset; 	
  int64_t auxiliary_data_length;
  int64_t auxiliary_data_offset;
  int64_t workspace_length;    	
  int64_t workspace_offset;    	
  int64_t num_scalars;          	
  int64_t number_reactions;     	
  int64_t number_molecules;            
  int64_t unique_molecules;     	
  int64_t use_activities;       	
  int64_t molecules_or_conc;      
  int64_t warmup_steps;           
  int64_t record_steps;           
  int64_t print_output;
  int64_t number_compartments;  	
  int64_t unique_compartments;  	

  int64_t align_len;
  int64_t align_mask;
  int64_t max_filename_len;       
  int64_t max_param_line_len;     
  int64_t num_rxn_file_keywords;  
  int64_t free_energy_format;     
  int64_t rxn_view_freq;          
  int64_t rxn_view_hist_length;    
  int64_t lklhd_view_freq;        
  int64_t conc_view_freq;         
  int64_t reaction_file_length;   
  int64_t rxn_title_space;        
  int64_t pathway_space;          
  int64_t compartment_space;      
  int64_t molecules_space;        
  int64_t num_fixed_concs;      	
  int64_t usage;                  
  int64_t max_molecule_len;     	
  int64_t min_molecule_len;     	
  int64_t max_compartment_len;  	
  int64_t min_compartment_len;  	
  int64_t sum_molecule_len;
  int64_t sum_compartment_len;


  double  ideal_gas_r;
  double  temp_kelvin;
  double  ph;
  double  ionic_strength;
  double  rt;
  double  m_r_rt;
  double  m_rt;
  double  cal_gm_per_joule;
  double  joule_per_cal_gm;
  double  default_initial_conc;
  double  dg_forward;
  double  entropy;
  int64_t *workspace_base;

  /* two way data (modified) */
  double  *dg_forward_p; /* scalar */        /* 8 */             
  double  *entropy_p;    /* scalar */        /* 8 */             
  double  *current_concentrations; /* len = unique_molecules */
  double  *bndry_flux_concs;       /* len = unique_molecules */
  struct  vgrng_state_struct *vgrng_state; /* len = 13 */
  struct  vgrng_state_struct *vgrng2_state;/* len = 13 */ 
  /* 
    2*sizeof(double) + unique_molecules * 2 * sizeof(double) 
    + 2 * sizeof(vgrng_state_struct)
  */
  /*
    Incoming data not modified, depends only on agent type.
    (2*number_reactions + unique_molecules) * sizeof(double)
  */
  double  *dg0s;      /* len = number_reactions */
  double  *ke;        /* len = number_reactions */
  double  *activities; /* len = unique_molecules */
  /* sizeof(rxn_struct) * number of reactions */
  struct rxn_struct *reactions; 
  /* 
     ((number_reactions + 1) * sizeof(int64_t)) +
     number_molecules * (4*sizeof(int64_t))
  */
  struct rxn_matrix_struct *reactions_matrix; 
  /* sizeof(istring_elem_struct) * unique_molecules */
  struct istring_elem_struct *sorted_molecules;   
  /* sizeof(istring_elem_struct) * unique_compartments */
  struct istring_elem_struct *sorted_cmpts;
  /*
    Auxilliary data read in by boltzmann_init, not needed by
    boltzman_run unless printing is enabled.
  */
  char *params_file;      /* max_filename_len */
  char *reaction_file;    /* max_filename_len */
  char *init_conc_file;   /* max_filename_len */
  char *input_dir;        /* max_filename_len */
  char *output_file;      /* max_filename_len */
  char *log_file;         /* max_filename_len */
  char *concs_out_file;   /* max_filename_len */
  char *rxn_lklhd_file;   /* max_filename_len */
  char *free_energy_file; /* max_filename_len */
  char *restart_file;     /* max_filename_len */
  char *rxn_view_file;    /* max_filename_len */
  char *bndry_flux_file;  /* max_filename_len */
  char *output_dir;       /* max_filename_len */
  /* 2*sizeof(double) + 
     ((unique_molecules * 3) + 2*number_reactions) * sizeof(double)+
     2*sizeof(vgrng_struct)  + 
     (sizeof(rxn_struct) * number of reactions) + 
     sizeof(rxn_matrix_struct) +
     ((number_reactions + 1) * sizeof(int64_t)) +
     (number_molecules * ((3*sizeof(int64_t)) + sizeof(char*))) +
     (sizeof(istring_elem_struct) * (unique_molecules+unqiue_compartments))) +
     (num_state_files * max_filename_len) +
     rxn_title_space + 	
     pathway_space + 
     compartment_space + 
     molecules_space +
     
     max(molecules_space+ ,
     (sizeof(istring_elem_struct) * (number_molecules+number_compartments))) +
	     
     ((unqiue_molecules + 6*number_reactions + 1)*sizeof(double)
     + print_output ? 
     ((rxn_view_hist_length *(2*number_reactions+1)) * sizeof(double)
     + ((2*number_reactions+2)*sizeof(int))) : 0))
  */
  char *rxn_title_text;   /* rxn_title_space */
  char *pathway_text;     /* pathway_space   */
  char *compartment_text; /* compartment_space */
  char *molecules_text;   /* molecules_space */
  /*
    Workspace only.
  */
  struct istring_elem_struct *unsorted_molecules;
  struct istring_elem_struct *unsorted_cmpts;
  int64_t *compartment_ptrs;
  int64_t *rxn_file_keyword_lengths;
  char    **rxn_file_keywords; /* 12 */
  char    *rxn_file_keyword_buffer; /* 144 */
  char    *param_buffer; /*  2* max_param_line_len */
  char    *raw_molecules_text; /* molecules_space */
  double  *future_concentrations;  /* unique_molecules */
  double  *free_energy;            /* number_reactions */
  double  *forward_rxn_likelihood; /* number_reactions */
  double  *reverse_rxn_likelihood; /* number_reactions */
  double  *forward_rxn_log_likelihood_ratio; /* number_reactions */
  double  *reverse_rxn_log_likelihood_ratio; /* number_reactions */
  double  *rxn_likelihood_ps;      /* number_reactions + 1 */
  /*
    Workspace only if printing. (debugging);
  */
  double  *no_op_likelihood;       /* rxn_view_hist_length */
  /* rxn_view_hist_length * number_reactions */
  double  *rxn_view_likelihoods;    
  /* rxn_view_hist_length * number_reactions */
  double  *rev_rxn_view_likelihoods; 
  int  *rxn_fire;                  /* (number_reactions * 2) + 2*/
  FILE *rxn_fp;
  FILE *conc_fp;
  FILE *out_fp;
  FILE *concs_out_fp;
  FILE *rxn_lklhd_fp;
  FILE *free_energy_fp;
  FILE *restart_fp;
  FILE *rxn_view_fp;
  FILE *bndry_flux_fp;
  FILE *lfp;
  /*
    Instrumentation.     
  */
  struct timing_struct *timing_data;
}
;
#endif
