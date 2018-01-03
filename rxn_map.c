/* rxn_map.c
*******************************************************************************
rxn_mapp

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

/*
#define RXN_MAP_DBG 1
*/
#include "rxn_map_init.h"
#include "rxn_map_run.h"
#include "rxn_map_parse_start_stop_line.h"
int main(int argc, char **argv)
{
  /*
    Calls:
      rxn_map_init
      rxn_map_parse_start_stop_file
      rxn_map_run
  */
  struct state_struct *state;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t *molecules_map_work;
  char *param_file_name;
  int success;
  int padi;

  int mi;
  int mj;

#ifdef LIBUNWIND
#include "luwtb1.h"
#endif
#ifdef LIBUNWIND
#include "luwtb2.h"
#endif
  TIMING_INIT("./timingi.h");
  TIMING_START(TOTAL_TIME);
  TIMING_START(INITIALIZE);
  if (argc > 1) {
    param_file_name = argv[1];
  } else {
    param_file_name = NULL;
  }
  success = rxn_map_init(param_file_name,&state,&molecules_matrix);
  TIMING_STOP(INITIALIZE);
  while (success) {
    success = rxn_map_parse_start_stop_line(state,&mi,&mj);
    if (success) {
      success = rxn_map_run(state,molecules_matrix,mi,mj);
    }
  }
  TIMING_STOP(TOTAL_TIME);
  TIMING_PRINT(stdout);
  fflush(stdout);
  exit(0);
}
