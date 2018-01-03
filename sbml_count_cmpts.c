/* sbml_count_cmpts.c
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

int sbml_count_cmpts(struct sbml2bo_struct *sbml_state) {
  /*
    Called by: sbml_to_boltzman and sbml2bo
    Calls:     sbml_find_section, fopen, fprintf, fflush, fclose
  */
  char sbml_buffer_c[2048];
  char *sbml_buffer;
  
  int success;
  int num_cmpts;
  int another;
  int sbml_buffer_len;
    
  FILE *sbml_fp;
  FILE *lfp;
  FILE *error_fp;
  FILE *extra_fp;
  success = 1;
  sbml_buffer_len = 2048;
  sbml_buffer = (char*)&sbml_buffer_c[0];
  lfp = sbml_state->log_fp;
  error_fp = lfp;
  if (error_fp == NULL) {
    error_fp = stderr;
  }
  sbml_fp = fopen(sbml_state->sbml_file,"r");
  if (sbml_fp == NULL) {
    success = 0;
    fprintf(error_fp,
	    "sbml_count_cmpts: Error unable to open sbml_file, %x\n",
	    sbml_state->sbml_file);
    fflush(error_fp);
  }
  if (success) {
    another = 1;
    num_cmpts = 0;
    while (!feof(sbml_fp) && another) {
      another = sbml_find_section(sbml_fp,sbml_buffer,
				  sbml_buffer_len,
				  "<compartment");
      if (another) {
	num_cmpts += 1;
      }
    }
  }
  sbml_state->num_cmpts = num_cmpts;
  fclose(sbml_fp);
  return(success);
}
