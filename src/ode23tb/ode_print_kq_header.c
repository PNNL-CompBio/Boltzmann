/* ode_print_lklhd_header.c
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
#include "print_rxns_f_and_r_header.h"

#include "ode_print_kq_header.h"
void ode_print_kq_header(struct state_struct *state) {
  /*
    Print the header lines for the ode_kq output file.
    Called by: deq_run .
    Calls    : print_rxns_f_and_r_header, fopen, fprintf, fflush
  */
  struct reaction_struct *reactions;
  char *rxn_title_text;
  char *title;
  int i;
  int nrxns;
  FILE *ode_kq_fp;
  FILE *lfp;
  rxn_title_text = state->rxn_title_text;
  lfp            = state->lfp;
  ode_kq_fp      = fopen(state->ode_kq_file,"w");
  state->ode_kq_fp = ode_kq_fp;
  if (ode_kq_fp) {
    fprintf(ode_kq_fp,"Reactions KQ and KQ^(-1)\n");
    print_rxns_f_and_r_header(state,"Time",ode_kq_fp);
  } else {
    if (lfp) {
      fprintf(lfp,"ode_print_kq_header: Error, unable to open %s for writing\n",
	      state->ode_kq_file);
      fflush(lfp);
    }
  }
  return;
}
