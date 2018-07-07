/* ode_print_kq_kqi.c
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

#include "ode_print_kq_kqi.h"
void ode_print_kq_kqi(struct state_struct *state, double time, double *kq,
		      double *kqi) {
  /* 
    print the gradient forward and reverse parts, kq and kqi, corresponding
    to the forward and reverse likeilhood approximations used by the
    gradient routne.
    Prints out the current kq and kqi.
    in a tab delimited row terminated by a newline.

    Called by: boltzmann_monidtor_ode.

    Arguments:
    
    Name          TMF      Description

    state         G*I      state structure :
                           input fields are nrxns,
			                    ode_kq_fp,
                           no fields of state are modified.
    dconcs 	  D*I      vector of delta concentrations to be printed.

    time          dsi      time stamp
                  
    
  */
  int nrxns;
  int j;

  FILE *ode_kq_fp;
  nrxns                  = state->number_reactions;
  ode_kq_fp              = state->ode_kq_fp;
  if (ode_kq_fp) {
    fprintf(ode_kq_fp,"time = %le\nrxn_no\t   kq   \t   kqi\n",time);
    for(j=0;j<nrxns;j++) {
      fprintf(ode_kq_fp,"%d\t%le\t%le\n",j,kq[j],kqi[j]);
    }
    fflush(ode_kq_fp);
  }
}
