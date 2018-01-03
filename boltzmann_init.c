/* boltzmann_init.c
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
#include "alloc0.h"
#include "read_params.h"
#include "boltzmann_init_core.h"
#include "free_boot_state.h"
/*
#define DBG_BOLTZMANN_INIT  
*/
#include "boltzmann_init.h"
int boltzmann_init(char *param_file_name, struct state_struct **statep) {
  /*
    Initialize the reactions and data structures for boltzmann.
    Called by: boltzmann
    Calls:     alloc0,
               read_params,
	       boltzmann_init_core
	       free_boot_state
  */
  struct state_struct *boot_state;
  struct state_struct *state;
  struct state_struct *stateq;
  struct formation_energy_struct *formation_energies;
  struct vgrng_state_struct *vgrng_state;
  struct vgrng_state_struct *vgrng2_state;
  double *dg0s;
  double *free_energy;
  double *activities;
  int64_t vgrng_start;
  int64_t i;

  int success;
  int vgrng_start_steps;

  int print_output;
  int padi;
  
  FILE *rxn_echo_fp;
  FILE *bndry_flux_fp;
  FILE *lfp;
  /*
    allocate space for the state struct.
    Allocate space for the reactions line buffer, and the rxn_file keywords.
  */
  success = alloc0(&boot_state);
  if (success) {
    /*
      Read the input parameters file.
    */
    success = read_params(param_file_name,boot_state);
  }
  success = boltzmann_init_core(boot_state,statep);
  if (success) {
    success = free_boot_state(&boot_state);
  }
  return(success);
}
