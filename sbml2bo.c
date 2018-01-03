/* sbml2bo.c
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

#include "sbml2bo_structs.h"
#include "boltzmann_structs.h"

#include "sbml_alloc0.h"
#include "sbml_set_file_names.h"
#include "sbml_alloc1.h"
#include "sbml_count_cmpts.h"
#include "parse_sbml.h"

int main(int argc, char **argv)
{
  struct sbml2bo_struct *sbml_state;
  int success;
  int base_len;

  int in_len;
  int max_file_name_len;
  char *concs_in_file;
  char *rxns_dat_file;
  char *cmpts_dat_file;
  char *log_file;
  char *tail;
  FILE *lfp;
  FILE *extra_fp;
    
  max_file_name_len = 1024;
  success = 1;
  /*
    Need to allocate space for the sbml_state structure and its
    filename fields.
  */
  success = sbml_alloc0(&sbml_state);
  /*
    Need to read input file name and generate output file names.
  */
  if (argc > 1) {
    strcpy(sbml_state->sbml_file,(char*)argv[1]);
  } else {
    success = 0;
    fprintf(stderr,"sbml_set_file_names: Error, you need to provide as an "
	    "argument the name of the sbml file\n");
    fflush(stderr);
  }
  if (success) {
    success = sbml_set_file_names(sbml_state);
  }
  if (success) {
    lfp = fopen(sbml_state->log_file,"w");
    if (lfp == NULL) {
      fprintf(stderr,"sbml2bo: Error could not open log file for writing.\n");
      success = 0;
    }
    sbml_state->log_fp = lfp;
  }
  if (success) {
    /*
      Count the compartments, setting the num_cmpts field in sbml_state.
    */
    success = sbml_count_cmpts(sbml_state);
  }
  if (success) {
    success = sbml_alloc1(sbml_state);
  }
  if (success) {
    success = parse_sbml(sbml_state);
  }
}
