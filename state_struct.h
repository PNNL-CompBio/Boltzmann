/* state_struct.h 
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
#ifndef __STATE_STRUCT__
#define __STATE_STRUCT__ 1
/* 
  So we will lay out the state struct a little differently.
  We want to flatten the state structure, and make it self describing.
  At the beginning there will be a set of scalars describing sizes
  and constansts to be used in the code, followed by a partitioning

  Let us make the following partitions
  1, for the size parameters and floating point constants.

  2. for arrays/structs that are both inputs and outputs (4 of these)
     current_concentrations,
     bndry_flux_concs,
     vgrng_state
     vgrng2_state
     and 2 scalars, dg_forward and entropy for the whole system
 
  3. for arrays/structs that are strictly inputs (there are 12 of these
     not counting string pointers).
     reactions,
     reactions_matrix,
     sorted_molecules,
     ke,
     kss,
     kssr,
     activities,
     sorted_compartments

     Each of following 4 arrays of length max_regs_per_rxn * num_reactions
     reg_species,   (int64_t)
     reg_drctn,     (int64_t)
     reg_constant,  (double)
     reg_exponent,  (double)

  4. for work space arrays and structs.
     unsorted_molecules,
     unsorted_cmpts.
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

    This struct is allocated in alloc0.
*/
struct state_struct {
  int64_t agent_type;
  int64_t state_length;         	
  int64_t thread_id;
  int64_t version_no;
  double  x_coord;
  double  y_coord;
  double  z_coord;
  double  epsilon; /* used in choose_rxn */
  int64_t num_partitions;       	
  int64_t two_way_data_length;  	
  int64_t two_way_data_offset;  	
  int64_t incoming_data_length; 	
  int64_t incoming_data_offset; 	
  int64_t auxiliary_data_length;
  int64_t auxiliary_data_offset;
  int64_t workspace_length;    	
  int64_t workspace_offset;    	
  int64_t number_reactions;     	
  int64_t number_molecules;            
  int64_t nunique_molecules;     	
  int64_t use_activities;       	
  int64_t use_deq;
  int64_t no_round_from_deq;
  int64_t adjust_steady_state;
  int64_t molecules_or_conc;      
  int64_t ode_solver_choice;
  int64_t delta_concs_choice;
  int64_t warmup_steps;           
  int64_t record_steps;           
  int64_t print_output;
  int64_t concs_or_counts;
  int64_t number_compartments;  	
  int64_t nunique_compartments;  	
  int64_t align_len;
  int64_t align_mask;
  int64_t max_filename_len;
  int64_t free_energy_format;
  int64_t rxn_view_freq;
  int64_t rxn_view_hist_length;
  int64_t ode_rxn_view_freq;
  int64_t lklhd_view_freq;
  int64_t count_view_freq;
  int64_t fe_view_freq;
  int64_t reaction_file_length;
  int64_t num_fixed_concs;
  int64_t usage;
  int64_t max_molecule_len;
  int64_t min_molecule_len;
  int64_t max_compartment_len;
  int64_t min_compartment_len;
  int64_t sum_molecule_len;
  int64_t sum_compartment_len;
  int64_t num_files;
  int64_t solvent_pos;
  int64_t use_pseudoisomers;
  int64_t use_metropolis;
  int64_t use_regulation;
  int64_t max_regs_per_rxn;
  int64_t rxn_filename_base_length;
  int64_t base_reaction;
  int64_t number_base_reaction_reactants;
  int64_t print_ode_concs;

  /*
    sizes used to self-describe this state vector.
    only needed for parallel version multiple instantiations
    when using boltzmann_boot and boltzmann_load, with flatten_state.
  */
  int64_t keyword_buffer_length;
  int64_t num_rxn_file_keywords;
  int64_t rxn_file_keyword_len;
  int64_t per_keyword_pointer_size;
  int64_t max_param_line_len;
  int64_t reaction_titles_length;
  int64_t pathway_text_length;
  int64_t compartment_text_length;
  int64_t molecule_text_length;
  int64_t regulation_text_length;
  int64_t per_molecule_double_size;
  int64_t per_reaction_double_size;
  int64_t per_compartment_double_size;
  int64_t vgrng_state_size;
  int64_t reaction_data_size;
  int64_t molecule_data_size;
  int64_t compartment_data_size;
  int64_t reactions_matrix_size;
  int64_t molecules_matrix_size;
  int64_t reactions_ptrs_size;
  int64_t molecules_ptrs_size;
  int64_t reactions_matrix_field_size;
  int64_t molecules_matrix_field_size;
  int64_t file_names_size;
  int64_t solvent_string_size;
  /*
    offsets used to self-describe this state vector.
    only needed for parallel version multiple instantiations
    when using boltzmann_boot and boltzmann_load, with flatten_state.
  */
  int64_t state_offset;
  int64_t current_counts_offset;
  int64_t bndry_flux_counts_offset;
  int64_t net_lklhd_bndry_flux_offset;
  int64_t net_likelihood_offset;
  int64_t vgrng_offset;
  int64_t vgrng2_offset;

  int64_t dg0s_offset;
  int64_t ke_offset;
  int64_t rke_offset;
  int64_t kss_offset;
  int64_t kssr_offset;
  int64_t kss_e_val_offset;
  int64_t kss_u_val_offset;
  int64_t molecule_dg0tfs_offset;
  int64_t molecule_probabilities_offset;
  int64_t molecule_chemical_potentials_offset;
  int64_t count_to_conc_offset;
  int64_t conc_to_count_offset;
  int64_t activites_offset;
  int64_t reg_constant_offset;
  int64_t reg_exponent_offset;
  int64_t reg_direction_offset;
  int64_t reg_species_offset;
  
  int64_t reactions_offset;
  int64_t sorted_molecules_offset;
  int64_t sorted_compartments_offset;
  int64_t reactions_matrix_offset;
  int64_t molecules_matrix_offset;
  int64_t reactions_ptrs_offset;
  int64_t molecules_indices_offset;
  int64_t compartment_indices_offset;
  int64_t reactions_coefficients_offset;
  int64_t text_indices_offset;
  int64_t solvent_coefficients_offset;
  int64_t molecules_ptrs_offset;
  int64_t reaction_indices_offset;
  int64_t molecules_coefficients_offset;

  int64_t file_names_offset;
  int64_t solvent_string_offset;
  int64_t reaction_titles_offset;
  int64_t pathway_text_offset;
  int64_t compartment_text_offset;
  int64_t molecules_text_offset;
  int64_t regulation_text_offset;

  int64_t unsorted_molecules_offset;
  int64_t unsorted_compartments_offset;
  int64_t keyword_buff_offset;
  int64_t keyword_lengths_offset;
  int64_t keywords_offset;
  int64_t raw_molecules_text_offset;
  int64_t transpose_workspace_offset;
  int64_t future_counts_offset;
  int64_t free_energy_offset;
  int64_t forward_rxn_likelihood_offset;
  int64_t reverse_rxn_likelihood_offset;
  int64_t forward_rxn_log_likelihood_ratio_offset;
  int64_t reverse_rxn_log_likelihood_ratio_offset;
  int64_t rxn_likelihood_ps_offset;
  int64_t reactant_term_offset;
  int64_t product_term_offset;
  int64_t rxn_q_offset;
  int64_t recip_rxn_q_offset;
  int64_t log_kf_rel_offset;
  int64_t log_kr_rel_offset;
  int64_t ode_counts_offset;
  int64_t ode_concs_offset;
  int64_t ode_forward_lklhds_offset;
  int64_t ode_reverse_lklhds_offset;
  int64_t base_reactant_indicator_offset;
  int64_t base_reactants_offset;
  int64_t rxn_has_flux_offset;
  /*
    Floating point scalars
  */
  double  ideal_gas_r;
  double  temp_kelvin;
  double  avogadro;
  double  recip_avogadro;
  double  min_count;
  double  min_conc;
  double  ph;
  double  ionic_strength;
  double  rt;
  double  m_r_rt;
  double  m_rt;
  double  cals_per_joule;
  double  joules_per_cal;
  double  default_initial_count;
  double  dg_forward;
  double  entropy;
  double  current_concentrations_sum;
  double  default_volume;
  double  recip_default_volume;
  double  conc_units;
  double  ntotal_opt;
  double  ntotal_exp;
  double  ode_t_final;
  /* 
    max_log_g0_sum should be set to maximum argument allowed for
    exp function, for double precision about 704
  */
  double  max_log_g0_sum; 
  double  dg0_scale_factor;
  /*
    Note that flux_scaling is K_f(base_rxn_reaction)*(product of reactant 
    concentrations in base reaction).
  */
  double  kf_base_reaction;
  /*
    Flux scaling, if 0 then use kf_base_raction * product of base reaction 
    reactant species concentratiopns, else use fixed FLUX_SCALING.
  */
  double  flux_scaling;
  /*
    min_molecule_dg0tf is the smallest dg0 of formation from molecules 
    computed in compute_moldeule_dg0tfs and used in 
    compute_moleculear_partition_probability.
  */
  double  min_molecule_dg0tf;
  /*
    The following two fields are not input fields, but are scalars and
    hence included here..
  */
  int64_t *workspace_base;

  /* two way data (modified) */
  double  *current_counts; /* len = unique_molecules */
  double  *bndry_flux_counts;       /* len = unique_molecules */
  double  *net_lklhd_bndry_flux;    /* len = unique_molecules */
  double  *net_likelihood;          /* len = number_reactions */
  struct  vgrng_state_struct *vgrng_state; /* len = 13 */
  struct  vgrng_state_struct *vgrng2_state;/* len = 13 */ 
  /* 
    (2 + (3*unique_molecules)+ number_reactions) * sizeof(double) 
    + 2 * sizeof(vgrng_state_struct)
  */
  /*
    Incoming data not modified, depends only on agent type.
    (9*number_reactions + 7*unique_molecules) * sizeof(double)
    + 2*number_reactions * max_regs_per_rxn * sizeof(double)
    + 2*number_reactions * max_regs_per_rxn * sizeof(int64_t)
        =
    (number_reactions*(5+4*max_regs_per_rxn) + 5*unique_molecules)*sizeof(double)
  */
  double  *dg0s;      /* len = number_reactions  */
  double  *ke;        /* len = number_reactions  */
  double  *rke;       /* len = number_reactions  */
  double  *forward_rc;  /* len = number_reactions */
  double  *reverse_rc;  /* len = number_reactions */
  double  *kss;       /* len = number_reactions  */
  double  *kssr;       /* len = number_reactions  */
  double  *kss_e_val; /* len = unique_molecules  */
  double  *kss_u_val; /* len = unique_molecules  */
  double  *molecule_dg0tfs; /* len = unique_molecules */
  double  *molecule_probabilities; /* len = unique_molecules */
  double  *molecule_chemical_potentials; /* len = unique_molecules */
  double  *count_to_conc; /* len = unique_molecules */
  double  *conc_to_count; /* len = unique_molecules */
  double  *activities;   /* len = number_reactions */
  double  *enzyme_level; /* len = number_reactions */
  double  *reg_constant; /* len = number_reactions * max_regs_per_rxn */
  double  *reg_exponent; /* len = number_reactions * max_regs_per_rxn */
  double  *reg_drctn;    /* len = number_reactions * max_regs_per_rxn */
  int64_t *reg_species;  /* len = number_reactions * max_regs_per_rxn */
  /* 
     sizeof(reaction_struct) * number of reactions. 
     Allocated in alloc2 
  */
  struct reaction_struct *reactions; 
  /* 
    sizeof(molecule_struct) * unique_molecules 
    allocated in alloc2 
  */
  struct molecule_struct *sorted_molecules;  
  /* 
     sizeof(compartment_struct) * unique_compartments 
     allocated in alloc2 
  */ 
  struct compartment_struct *sorted_compartments; 
  /* 
     reactions_matrix size is ((number_reactions + 1) * sizeof(int64_t)) +
     number_molecules * (4*sizeof(int64_t)) 
     allocated in alloc2 
  */
  struct reactions_matrix_struct *reactions_matrix; /* */

  /*
    molecules_matrix is allocated in alloc4.
  */
  struct molecules_matrix_struct *molecules_matrix;
  /*
    Auxilliary data read in by boltzmann_init, not needed by
    boltzman_run unless printing is enabled.
    These strings are allocated in alloc0.
  */
  char *params_file;       /* max_filename_len */
  char *reaction_file;     /* max_filename_len */
  char *init_conc_file;    /* max_filename_len */
  char *input_dir;         /* max_filename_len */
  char *output_file;       /* max_filename_len */
  char *log_file;          /* max_filename_len */
  char *output_dir;        /* max_filename_len */
  char *counts_out_file;   /* max_filename_len */
  char *ode_concs_file;    /* max_filename_len */
  char *net_lklhd_file;    /* max_filename_len */
  char *nl_bndry_flx_file; /* max_filename_len */
  char *rxn_lklhd_file;    /* max_filename_len */
  char *free_energy_file;  /* max_filename_len */
  char *restart_file;      /* max_filename_len */
  char *rxn_view_file;     /* max_filename_len */
  char *bndry_flux_file;   /* max_filename_len */
  char *pseudoisomer_file; /* max_filename_len */
  char *compartment_file;  /* max_filename_len */
  char *sbml_file;         /* max_filename_len */
  char *ms2js_file;        /* max_filename_len */
  char *kg2js_file;        /* max_filename_len */
  char *rxn_echo_file;     /* max_filename_len */
  char *rxn_mat_file;      /* max_filename_len */
  char *dg0ke_file;        /* max_filename_len */
  char *dictionary_file;   /* max_filename_len */
  char *ode_dconcs_file;   /* max_filename_len */
  char *ode_lklhd_file;    /* max_filename_len */
  char *ode_bflux_file;    /* max_filename_len */
  char *concs_out_file;    /* max_filename_len */
  char *solvent_string;    /* Length is 64. Allocated in alloc0 */

  char *rxn_title_text;    /* rxn_title_text_length. Allocated in alloc2  */
  char *pathway_text;      /* pathway_text_length. Allocated in alloc2    */
  char *compartment_text;  /* compartment_text_length Allocated in alloc2 */
  char *molecules_text;    /* molecule_text_length Allocated in alloc2 */
  char *regulation_text;    /* regulation_text_length Allocated in alloc2 */
  /*
    Workspace only.
  */
  struct molecule_struct *unsorted_molecules; /* allocated in alloc2 */
  struct compartment_struct *unsorted_cmpts; /* allocated in alloc2 */
  int64_t *compartment_ptrs;                 /* allocated in alloc3 */
  int64_t *transpose_workspace;  /* used in forming molecules_matrix. */
                                 /* allocated in alloc4 */
  int64_t *rxn_file_keyword_lengths /* allocated in alloc0 */;
  char    **rxn_file_keywords; /* 12, allocated in alloc0 */
  char    *rxn_file_keyword_buffer; /* 144, allocated in alloc0  */
  char    *param_buffer; /*  2* max_param_line_len, allocated in alloc0 */ 
  char    *raw_molecules_text; /* molecule_text_length Allocated in alloc2 */
  /* 
    Allocated in alloc8.
  */
  double  *future_counts;      /* unique_molecules */
  double  *free_energy;            /* number_reactions */
  double  *forward_rxn_likelihood; /* number_reactions */
  double  *reverse_rxn_likelihood; /* number_reactions */
  double  *forward_rxn_log_likelihood_ratio; /* number_reactions */
  double  *reverse_rxn_log_likelihood_ratio; /* number_reactions */
  double  *rxn_likelihood_ps;      /* number_reactions + 1 */

  /* Workspace used by ode routines. Allocated in alloc7 */
  double *reactant_term; /* product of reaction reactant concentrations, length number_reactions */
  double *product_term;  /* product of reaction product concentrations, length number_reactions */
  double *rxn_q; /* Ratio  product_term to reactant_term.  */
  double *recip_rxn_q;/* Ratio of reactant_term to product_term. */
  double *log_kf_rel;
  double *log_kr_rel;
  double *ode_counts; /* counts from concentrations */
  double *ode_concs;  /* concentrations from counts. */
  double *ode_forward_lklhds;
  double *ode_reverse_lklhds;
  int *rxn_has_flux; /* Indicator as to whether a reaction contributes to 
		       species flux  length is number_reactions */
  int  *base_reactants;            /* List of reactant species (by number)
				      in the base reaction */
  int  *base_reactant_indicator;   /* vector of length nunique_molecules
				      elment i 1 for species i in
				      base_reactants list, 0 otherwise. */
  /*
    Workspace only if printing. (debugging);
    allocated in alloc9.
  */
  double  *no_op_likelihood;       /* rxn_view_hist_length */
  /* rxn_view_hist_length * number_reactions */
  double  *rxn_view_likelihoods;    
  /* rxn_view_hist_length * number_reactions */
  double  *rev_rxn_view_likelihoods; 
  int64_t *rxn_fire;                  /* (number_reactions * 2) + 2*/
  int  *rxn_mat_row;               /* (nunique_molecules) */



  FILE *rxn_fp;
  FILE *conc_fp;

  FILE *out_fp;
  FILE *lfp;

  FILE *counts_out_fp;
  FILE *concs_out_fp;

  FILE *ode_concs_fp;
  FILE *net_lklhd_fp;

  FILE *nl_bndry_flx_fp;
  FILE *rxn_lklhd_fp;

  FILE *free_energy_fp;
  FILE *restart_fp;

  FILE *rxn_view_fp;
  FILE *bndry_flux_fp;

  FILE *cmpt_fp;
  FILE *ode_dconcs_fp;

  FILE *ode_lklhd_fp;
  FILE *ode_bflux_fp;
  /*
    Instrumentation.     
  */
  struct timing_struct *timing_data;
}
;
#endif
