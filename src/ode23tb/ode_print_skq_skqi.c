/* ode_print_skq_skqi.c
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

#include "ode_print_skq_skqi.h"
void ode_print_skq_skqi(struct state_struct *state, 
			double time, 
			double *skq,
			double *skqi) {
  /* 
    print the scaled KQ and scaled KQ^-1 values per rreaction per time step.
    Prints out the current skq and skqi.
    in a tab delimited row terminated by a newline.

    Called by: boltzmann_monidtor_ode.

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are nrxns,
			                    ode_skq_fp,
                           no fields of state are modified.
    
    time          dsi      time stamp
                  
    skq 	  D*I      vector of scaled KQ to be printed

    skqi 	  D*I      vector of scaled KQ^-1 to be printed

  */
  int nrxns;
  int j;

  FILE *ode_skq_fp;
  nrxns                  = state->number_reactions;
  ode_skq_fp             = state->ode_skq_fp;
  if (ode_skq_fp) {
    fprintf(ode_skq_fp,"%le",time);
    for(j=0;j<nrxns;j++) {
      fprintf(ode_skq_fp,"\t%le\t%le\n",skq[j],skqi[j]);
    }
    fprintf(ode_skq_fp,"\n");
    fflush(ode_skq_fp);
  }
}
