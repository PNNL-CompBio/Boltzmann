/* read_params.c
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <float.h>
#include <signal.h>
#include <unistd.h>

#include "boltzmann_structs.h"

#include "read_params.h"
int read_params (char *param_file_name, struct state_struct *state) {
  /*
    Read paramaters for the boltzmann code to determine equilbrium
    concentrations of a set of reactions via Monte Carlo methods.

    Called by: boltzmann_init
    Calls:     fopen, fprintf, fgets, feof, sscanf, strncmp (intrinsic)
  */
  int64_t max_param_line_len;
  char *param_buffer;
  char *key;
  char *value;
  char *rtp;
  int success;
  int sscan_ok;
  FILE *in_fp;
  success = 1;
  if (param_file_name == NULL) {
    fprintf(stdout,"read_params: Warning no param_file_name given, looking for ./boltzmann.in input file.\n");
    fflush(stdout);
    strcpy(state->params_file,"./boltzmann.in");
  } else {
    strcpy(state->params_file,param_file_name);
  }
  in_fp = fopen(state->params_file,"r");
  if (in_fp == NULL) {
    fprintf(stderr,"read_params: Error opening %s.\n",state->params_file);
    fflush(stderr);
    success = 0;
  }
  if (success) {
    /*
      Set defaults.
    */
    strcpy(state->reaction_file,"./rxns.in");
    strcpy(state->init_conc_file,"./concs.in");
    strcpy(state->log_file,"./boltzmann.log");
    strcpy(state->output_file,"./boltzmann.out");
    strcpy(state->concs_out_file,"./concs.out");
    strcpy(state->rxn_lklhd_file,"./rxns.lklhd");
    strcpy(state->free_energy_file,"./rxns.fe");
    strcpy(state->input_dir,"./");
    strcpy(state->output_dir,"./");
    state->align_len        = (int64_t)16;
    state->max_filename_len = 4096;
    state->max_param_line_len = 4096;
    state->align_mask       = state->align_len - 1;
    state->ideal_gas_r      = 0.00198858775;
    state->temp_kelvin      = 298.15;
    state->cal_gm_per_joule = 4.184;
    state->rxn_buff_len     = (int64_t)4194304;
    state->small_nonzero    = 1.e-31;
    state->warmup_steps     = 1000;
    state->record_steps     = 1000;
    state->free_energy_format = 0;
    param_buffer       = state->param_buffer;
    max_param_line_len = state->max_param_line_len;
    key                = state->param_key;
    value              = state->param_value;
    /*
      read in parameters.
    */
    sscan_ok = 0;
    rtp = fgets(param_buffer,max_param_line_len,in_fp);
    if (rtp) {
      sscan_ok = sscanf(param_buffer,"%s %s",key,value);
    }
    while ((!feof(in_fp)) && (sscan_ok == 2)) {
      if (strncmp(key,"RXN_FILE",8) == 0) {
	sscan_ok = sscanf(value,"%s",state->reaction_file);
      } else if (strncmp(key,"INIT_FILE",9) == 0) {
	sscan_ok = sscanf(value,"%s",state->init_conc_file);
      } else if (strncmp(key,"IN_DIR",6) == 0) {
	sscan_ok = sscanf(value,"%s",state->input_dir);
      } else if (strncmp(key,"OUT_DIR",7) == 0) {
	sscan_ok = sscanf(value,"%s",state->output_dir);
      } else if (strncmp(key,"OUT_FILE",8) == 0) {
	sscan_ok = sscanf(value,"%s",state->output_file);
      } else if (strncmp(key,"CONCS_OUT_FILE",13) == 0) {
	sscan_ok = sscanf(value,"%s",state->concs_out_file);
      } else if (strncmp(key,"RXN_LKLHD_FILE",13) == 0) {
	sscan_ok = sscanf(value,"%s",state->rxn_lklhd_file);
      } else if (strncmp(key,"FREE_ENERGY_FILE",15) == 0) {
	sscan_ok = sscanf(value,"%s",state->free_energy_file);
      } else if (strncmp(key,"LOG_FILE",8) == 0) {
	sscan_ok = sscanf(value,"%s",state->log_file);
      } else if (strncmp(key,"ALIGN_LEN",9) == 0) {
	sscan_ok = sscanf(value,"%ld",&(state->align_len));
	if (state->align_len < 0) {
	  state->align_len = 16;
	}
	state->align_mask = state->align_len - 1;
      } else if (strncmp(key,"RXN_BUFF_LEN",9) == 0) {
	sscan_ok = sscanf(value,"%ld",&(state->rxn_buff_len));
      } else if (strncmp(key,"IDEAL_GAS_R",11) == 0) {
	sscan_ok = sscanf(value,"%le",&(state->ideal_gas_r));
      } else if (strncmp(key,"TEMP_KELVIN",11) == 0) {
	sscan_ok = sscanf(value,"%le",&(state->temp_kelvin));
      } else if (strncmp(key,"WARMUP_STEPS",12) == 0) {
	sscan_ok = sscanf(value,"%d",&(state->warmup_steps));
      } else if (strncmp(key,"RECORD_STEPS",12) == 0) {
	sscan_ok = sscanf(value,"%d",&(state->record_steps));
      } else if (strncmp(key,"FREE_ENERGY_FORMAT",12) == 0) {
	if (strncmp(value,"NONE",4) == 0) {
	  state->free_energy_format = 0;
	} else if (strcmp(value,"NEG_LOG_LKLHD") == 0) {
	  state->free_energy_format = 1;
	} else if (strcmp(value,"KJ/MOL") == 0) {
	  state->free_energy_format = 2;
	} else if (strcmp(value,"KCAL/MOL") == 0) {
	  state->free_energy_format = 3;
	} else {
	  sscan_ok = sscanf(value,"%d",&(state->free_energy_format));
	  if (sscan_ok != 1) {
	    fprintf(stderr,"read_params: invalid value for free_energy "
		    "format using 0\n");
	    fflush(stderr);
	    state->free_energy_format = 0;
	  }
	}
      }
      sscan_ok = 0;
      rtp = fgets(param_buffer,max_param_line_len,in_fp);
      if (rtp) {
	sscan_ok = sscanf(param_buffer,"%s %s",key,value);
      }
    }
  }
  return (success);
}
