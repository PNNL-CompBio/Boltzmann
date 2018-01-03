/* state_struct.h 
  
sss

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
*/
#ifndef __STATE_STRUCT__
#define __STATE_STRUCT__ 1
struct state_struct {
  char *params_file;
  char *reaction_file;
  char *init_conc_file;
  char *input_dir;
  char *output_file;
  char *log_file;
  char *output_dir;
  char *rxn_buffer;
  char *param_buffer;
  char *param_key;
  char *param_value;
  char *rxn_file_keyword_buffer;
  char **rxn_file_keywords;
  int64_t *rxn_file_keyword_lengths;
  int64_t reaction_file_length;
  int64_t align_len;
  int64_t max_filename_len;
  int64_t max_param_line_len;
  int64_t align_mask;
  int64_t rxn_buff_len;
  int64_t rxn_title_len;
  int64_t pathway_len;
  int64_t compartment_len;
  int64_t species_len;
  double  ideal_gas_r;
  double  temp_kelvin;
  double  cal_gm_per_joule;
  double  joule_per_cal_gm;
  int  number_reactions;
  int  number_species;

  int  number_compartments;
  int  max_species_len;

  int  min_species_len;
  int  molecules_or_conc;

  int  num_rxn_file_keywords;
  int  pad1;
#ifdef TIMING_ON
  struct timing_struct timing_data;
#endif
  FILE *lfp;
}
;
#endif
