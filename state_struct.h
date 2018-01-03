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
struct state_struct {
  struct rxn_struct *reactions; /* sizeof(rxn_struct) * number of reactions */
  struct rxn_matrix_struct *reactions_matrix; 
  /* 
     ((number_reactions + 1) * sizeof(int64_t)) +
     number_molecules * ((3*sizeof(int64_t)) + sizeof(char*))
  */
  struct molecules_matrix_struct *molecules_matrix;
  struct istring_elem_struct *unsorted_molecules; /* 24 * number_molecules */
  struct istring_elem_struct *sorted_molecules;   /* 24 * number_molecules */
  struct istring_elem_struct *unsorted_cmpts;    /* 24 * number_compartments */
  struct istring_elem_struct *sorted_cmpts;      /* 24 * number_compartments */
  struct vgrng_state_struct  *vgrng_state; /* 13 * 8 = 104 bytes */
  struct vgrng_state_struct  *vgrng2_state; /* 13 * 8 = 104 bytes */
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
  char *rxn_buffer;       /* rxn_buff_len */
  char *param_buffer;     /* max_param_line_len */
  char *param_key;        /* max_param_line_len>>1 */
  char *param_value;      /* max_param_line_len>>1 */
  char *rxn_file_keyword_buffer;  /* 144 chars */
  char *rxn_title_text;           /* rxn_title_space */
  char *pathway_text;             /* pathway_space   */
  char *compartment_text;         /* compartment_space */
  char *molecules_text;           /* molecules_space */
  char *molecule_name;            /* max_molecule_len */
  char *compartment_name;         /* max_compartment_len */
  char *raw_molecules_text;       /* molecules_space */
  char **rxn_file_keywords;       /* 12 * sizeof(char*) */
  int64_t *rxn_file_keyword_lengths; /* 12 */
  int64_t *transpose_work;           /* unique_molecules+1 */
  int64_t *compartment_ptrs;         /* unique_compartments + 1 */
  int64_t reaction_file_length;
  int64_t print_output;
  int64_t align_len;
  int64_t align_mask;
  int64_t max_filename_len;
  int64_t max_param_line_len;
  int64_t rxn_buff_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t molecules_len;
  int64_t rxn_title_space;
  int64_t pathway_space;
  int64_t compartment_space;
  int64_t molecules_space;
  /* These are just local variables used in parse_reactions_file.
  int64_t rxn_title_pos;
  int64_t pathway_pos;
  int64_t compartment_pos;
  int64_t molecules_pos;
  int64_t mixed_case_pos;
  */
  int64_t usage;
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
  /*
  double  small_nonzero;
  */
  double  *dg0s;                   /* number_reactions */
  double  *ke;                     /* number_reactions */
  double  *activities;             /* number_reactions */
  double  *current_concentrations; /* unique_molecules */
  double  *bndry_flux_concs;       /* unique_molecules */
  double  *future_concentrations;  /* unique_molecules */
  double  *free_energy;            /* number_reactions */
  double  *forward_rxn_likelihood; /* number_reactions */
  double  *reverse_rxn_likelihood; /* number_reactions */
  double  *rxn_likelihood_ps;      /* number_ractions + 1 */
  double  *forward_rxn_log_likelihood_ratio; /* number_reactions */
  double  *reverse_rxn_log_likelihood_ratio; /* number_reactions */
  double  *rxn_view_likelihoods;    /* rxn_view_hist_length * number_reactions */
  double  *rev_rxn_view_likelihoods; /* rxn_view_hist_length * number_reactions */
  double  *no_op_likelihood;       /* rxn_view_hist_length */
  int  *rxn_fire;                  /* (number_reactions * 2) + 2*/

  int64_t  number_reactions;
  /* free energy_format, 0 for none, 1 for negative log likelihoods,
     2 for KJ/mol, 3 for Kcal/mol.
  */
  int64_t  free_energy_format;
  

  int64_t  number_compartments;
  int64_t  unique_compartments;

  int64_t  number_molecules;
  int64_t  unique_molecules;

  int64_t  max_molecule_len;
  int64_t  min_molecule_len;

  int64_t  max_compartment_len;
  int64_t  min_compartment_len;

  int64_t  molecules_or_conc;
  int64_t  num_rxn_file_keywords;

  int64_t  warmup_steps;
  int64_t  record_steps;

  /*
    reaction view history frequency.
    Default value is 0: don't display.
  */
  int64_t  rxn_view_freq;
  /*
    Likelihood transpose history length:
    rxn_view_hist_lngth = 
       floor((record_steps + rxn_view_freq - 1)/rxn_view_freq);
  */
  int64_t  rxn_view_hist_lngth;

  int64_t  lklhd_view_freq;
  int64_t  conc_view_freq;

  int64_t  num_fixed_concs;
  int64_t  use_activities;

#ifdef TIMING_ON
  struct timing_struct timing_data;
#endif
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
}
;
#endif
