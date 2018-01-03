/* echo_params.c
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
#include "echo_params.h"
int echo_params (FILE *lfp, struct state_struct *state) {
  /*
    Echo the state paramaters for the boltzmann code to determine equilbrium
    concentrations of a set of reactions via Monte Carlo methods.
    
    Called by: boltzmann main program.
 */
  int success;
  int pad1;
  success = 1;
  if (lfp) {
    fprintf(lfp,"state->params_file    	   = %s\n",state->params_file);
    fprintf(lfp,"state->reaction_file  	   = %s\n",state->reaction_file);
    fprintf(lfp,"state->init_conc_file 	   = %s\n",state->init_conc_file);
    fprintf(lfp,"state->input_dir      	   = %s\n",state->input_dir);
    fprintf(lfp,"state->output_file    	   = %s\n",state->output_file);
    fprintf(lfp,"state->log_file      	   = %s\n",state->log_file);
    fprintf(lfp,"state->output_dir     	   = %s\n",state->output_dir);
    fprintf(lfp,"state->align_len          = %ld\n",state->align_len);
    fprintf(lfp,"state->max_filename_len   = %ld\n",state->max_filename_len);
    fprintf(lfp,"state->max_param_line_len = %ld\n",state->max_param_line_len);
    fprintf(lfp,"state->rxn_buff_len       = %ld\n",state->rxn_buff_len);
    fprintf(lfp,"state->ideal_gas_r        = %le\n",state->ideal_gas_r);
    fprintf(lfp,"state->temp_kelvin        = %le\n",state->temp_kelvin);
  }
  return (success);
}
