/* compute_ke.c
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

#include "compute_ke.h"
int compute_ke(struct state_struct *state) {
  /*
    Fill the dg0s vector from the reactions struct
    and compute the equilibrium constants from the delta G_0's.

    K_eq  = exp(-delta_G_0/RT)
    
    Called by: boltzmann_init
  */
  struct rxn_struct *reactions;
  struct rxn_struct *reaction;
  /*
  struct species_matrix_struct species_matrix;
  */
  double dg0;
  double *dg0s;
  double *ke;
  double m_r_rt;
  double joules_per_cal_gm;
  double ideal_gas_r;
  double temp_kelvin;

  int success;
  int nrxns;
  int i;
  int padi;
  success = 1;
  ideal_gas_r = state->ideal_gas_r;
  temp_kelvin = state->temp_kelvin;
  if (temp_kelvin > 0) {
    m_r_rt = -1.0/(ideal_gas_r * temp_kelvin);
    state->m_r_rt = m_r_rt;
  } else {
    success = 0;
    fprintf (stderr,
	     "compute_ke: Error at temp_kelvin = 0 Ke  = 0, nothing happens.");
    fflush(stderr);
  }
  if (success) {
    nrxns     = state->number_reactions;
    reactions = state->reactions;
    dg0s      = state->dg0s;
    ke        = state->ke;
    joules_per_cal_gm = 1.0/(state->cal_gm_per_joule);
    reaction = reactions;
    for (i=0;i<nrxns;i++) {
      if (reaction->unit_i == 0) {
	/*
	  Delta_g0 units was in calories, convert to joules.
	*/
	dg0 = reaction->delta_g0 * joules_per_cal_gm;
      } else {
	dg0 = reaction->delta_g0;
      }
      dg0s[i] = dg0;
      ke[i] = exp(dg0 * m_r_rt);
      reaction += 1; /* Caution address arithmetic */
    }
  }
  return (success);
}
  
