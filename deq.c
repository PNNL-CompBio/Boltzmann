/* deq.c
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

/*
  Initialize the boltzmann system, set up the ode solver and
  solve the ode's to get a steady state restart file printed out.
  Calls:  boltzmann_init, deq_run, print_restart
*/



/*
*/
#define BOLTZMANN_DBG 1
#include "boltzmann_init.h"
#include "deq_run.h"
#include "print_counts.h"
#include "print_restart_file.h"
int main(int argc, char **argv)
{
  struct state_struct *state;
  char *param_file_name;
  int success;
  int j;

  if (argc > 1) {
    param_file_name = argv[1];
  } else {
    param_file_name = NULL;
  }
  success = boltzmann_init(param_file_name,&state);
  if (success) {
    state->print_ode_concs = 1;
    success = deq_run(state);
    j = 1;
    print_counts(state,j);
    print_restart_file(state);
  }
  fflush(stdout);
  exit(0);
}
