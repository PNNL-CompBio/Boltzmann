/* ode_print_concs.c
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

#include "ode_print_lklhds.h"
void ode_print_lklhds(struct state_struct *state,
		      double t,
		      double *forward_rxn_likelihoods,
		      double *reverse_rxn_likelihoods) {
  int nrxns;
  int j;
  FILE *ode_lklhd_fp;
  FILE *lfp;
  nrxns        = state->number_reactions;
  ode_lklhd_fp = state->ode_lklhd_fp;
  if (ode_lklhd_fp) {
    fprintf(ode_lklhd_fp,"le",time);
    for(j=0;j<nrxns;j++) {
      fprintf(ode_lklhd_fp,"\t%le\t%le",forward_rxn_likelihoods[j],
	      reverse_rxn_likelihoods[j]);
    }
    fprintf(ode_lklhd_fp,"\n");
    fflush(ode_lklhd_fp);
  }
}
