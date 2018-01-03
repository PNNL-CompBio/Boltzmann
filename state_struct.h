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
  struct rxn_struct *reactions;
  struct rxn_matrix_struct *reactions_matrix;
  struct molecules_matrix_struct *molecules_matrix;
  struct istring_elem_struct *unsorted_molecules;
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *unsorted_cmpts;
  struct istring_elem_struct *sorted_cmpts;
  struct vgrng_state_struct  *vgrng_state;
  char *params_file;
  char *reaction_file;
  char *init_conc_file;
  char *input_dir;
  char *output_file;
  char *log_file;
  char *concs_out_file;
  char *rxn_lklhd_file;
  char *free_energy_file;
  char *output_dir;
  char *rxn_buffer;
  char *conc_buffer;
  char *param_buffer;
  char *param_key;
  char *param_value;
  char *rxn_file_keyword_buffer;
  char *rxn_title_text;
  char *pathway_text;
  char *compartment_text;
  char *molecules_text;
  char *molecule_name;
  char *compartment_name;
  char *raw_molecules_text;
  char **rxn_file_keywords;
  int64_t *rxn_file_keyword_lengths;
  int64_t *transpose_work;
  int64_t *compartment_ptrs;
  int64_t reaction_file_length;
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
  int64_t rxn_title_pos;
  int64_t pathway_pos;
  int64_t compartment_pos;
  int64_t molecules_pos;
  int64_t mixed_case_pos;
  int64_t usage;
  double  ideal_gas_r;
  double  temp_kelvin;
  double  rt;
  double  m_r_rt;
  double  m_rt;
  double  cal_gm_per_joule;
  double  joule_per_cal_gm;
  double  default_initial_conc;
  double  small_nonzero;
  double  *current_concentrations;
  double  *future_concentrations;
  double  *dg0s;
  double  *free_energy;
  double  *forward_rxn_likelihood;
  double  *reverse_rxn_likelihood;
  double  *rxn_likelihood_ps;
  double  *ke;
  double  *current_rxn_log_likelihood_ratio;
  double  *future_rxn_log_likelihood_ratio;
  double  *l_thermo;
  int  *cmpts_map;

  int  number_reactions;
  /* free energy_format, 0 for none, 1 for negative log likelihoods,
     2 for KJ/mol, 3 for Kcal/mol.
  */
  int  free_energy_format;
  

  int  number_compartments;
  int  unique_compartments;

  int  number_molecules;
  int  unique_molecules;

  int  max_molecule_len;
  int  min_molecule_len;

  int  max_compartment_len;
  int  min_compartment_len;

  int  molecules_or_conc;
  int  num_rxn_file_keywords;

  int  warmup_steps;
  int  record_steps;

#ifdef TIMING_ON
  struct timing_struct timing_data;
#endif
  FILE *rxn_fp;
  FILE *conc_fp;
  FILE *out_fp;
  FILE *concs_out_fp;
  FILE *rxn_lklhd_fp;
  FILE *free_energy_fp;
  FILE *lfp;
}
;
#endif
