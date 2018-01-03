/* boltzmann_run_sim.c
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

#include "djb_timing_b.h"
#include "boltzmann_structs.h"
/*
#define DBG_BOLTZMANN_RUN_SIM 1
*/

#include "free_energy_changes.h"
#include "vgrng.h"
#include "binary_search_l_u_b.h"
#include "forward_rxn_conc_update.h"
#include "reverse_rxn_conc_update.h"
#include "log_rxn_ratios.h"

#include "boltzmann_run_sim.h"
int boltzmann_run_sim(struct state_struct *state) {
  /*
    Run the boltzmann simulations, to be called after
    boltzmann_init has been called.

    Called by: boltzmann
    Calls:     free_energy_changes,
               vgrng,
	       binary_search_l_u_b,
	       forward_rxn_conc_update,
	       reverse_rxn_conc_update,
	       log_rxn_ratios,
  */
  struct vgrng_state_struct *vgrng_state;
  int64_t choice;
  double dchoice;
  double uni_multiplier;
  double vall;
  double rvall;
  double scaling;
  double dg_forward;
  double *cconcs;
  double *fconcs;
  double *rxn_likelihood_ps;
  double *free_energy;
  double *c_lrr;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *lthermo;
  double m_rt;
  int success;
  int num_state_files;

  int num_rxns;
  int num_rxns_t2;

  int num_rxns_t2_p1;
  int i;

  int n_warmup_steps;
  int n_record_steps;

  int rxn_choice;
  int nu;

  int j;
  int padi;

  FILE *lfp;
  FILE *concs_out_fp;
  FILE *rxn_lklhd_fp;
  success = 1;
  
  n_warmup_steps = state->warmup_steps;
  n_record_steps = state->record_steps;
  vgrng_state    = state->vgrng_state;
  rxn_likelihood_ps         = state->rxn_likelihood_ps;
  free_energy            = state->free_energy;
  c_lrr          = state->current_log_rxn_ratio;
  num_rxns       = state->number_reactions;
  nu             = state->unique_molecules;
  cconcs         = state->current_concentrations;
  fconcs         = state->future_concentrations;
  m_rt           = state->m_rt;
  lthermo        = state->l_thermo;
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  num_rxns_t2    = num_rxns << 1;
  num_rxns_t2_p1 = num_rxns_t2 + 1;
  uni_multiplier = vgrng_state->uni_multiplier;
  concs_out_fp   = state->concs_out_fp;
  rxn_lklhd_fp   = state->rxn_lklhd_fp;
  for (i=0;i<n_warmup_steps;i++) {
    /*
      Comput the reaction likelihoods - forward_rxn_likelihood, 
      and reverse_rxn_likelihood.
    */
    success = free_energy_changes(state);
    /*
      Compute the partial sums of the reaction likelihoods.
    */
    rxn_likelihood_ps[0] = forward_rxn_likelihood[0];
    for (j=1;j<num_rxns;j++) {
      rxn_likelihood_ps[j] = rxn_likelihood_ps[j-1] + forward_rxn_likelihood[j];
    }
    for(j=0;j<num_rxns;j++) {
      rxn_likelihood_ps[num_rxns+j] = rxn_likelihood_ps[num_rxns-1+j] + reverse_rxn_likelihood[j];
    }
    /*
      Doug wonders why 1.0 + here?
    */
    vall = 1.0 + rxn_likelihood_ps[num_rxns+num_rxns-1];
    rxn_likelihood_ps[num_rxns+num_rxns] = vall;
    /*
      Unimultiplier is 1.0/2^31-1
    */
    scaling = vall*uni_multiplier;
    /*
      choice is a pseudo-random integer in [0,2^31-1] (inclusive).
     */
    choice  = vgrng(vgrng_state);
    dchoice = ((double)choice)*scaling;
    /*
      Find index of smallest dg_ps entry that is >= choice.
    */
    rxn_choice = binary_search_l_u_b(rxn_likelihood_ps,dchoice,num_rxns_t2_p1);
    if (rxn_choice < num_rxns) {
      success = forward_rxn_conc_update(rxn_choice,state);
    } else {
      if (rxn_choice < num_rxns_t2) {
	success = reverse_rxn_conc_update(rxn_choice-num_rxns,state);
      } /* else do nothing, no change */
    }
    for (j=0;j<nu;j++) {
      cconcs[j] = fconcs[j];
    }
    /*
      Doug wonders here whether or not we ought to update free_energy based
      on the new concentrations, Absolutelyh we should.
    */
    success = log_rxn_ratios(free_energy,cconcs,c_lrr,c_lrr,state);
    for (j=0;j<num_rxns;j++) {
      free_energy[j] = m_rt * c_lrr[j];
    }
  } /* end for(i...) */
  
  for (i=0;i<n_record_steps;i++) {
    /*
      Comput the reaction likelihoods - forward_rxn_likelihood, and reverse_free_energy.
    */
    success = free_energy_changes(state);
    /*
      Compute the partial sums of the reaction likelihoods.
    */
    rxn_likelihood_ps[0] = forward_rxn_likelihood[0];
    for (j=1;j<num_rxns;j++) {
      rxn_likelihood_ps[j] = rxn_likelihood_ps[j-1] + forward_rxn_likelihood[j];
    }
    for(j=0;j<num_rxns;j++) {
      rxn_likelihood_ps[num_rxns+j] = rxn_likelihood_ps[num_rxns-1+j] + reverse_rxn_likelihood[j];
    }
    /*
      Doug wonders why 1.0 + here?
    */
    vall = 1.0 + rxn_likelihood_ps[num_rxns+num_rxns-1];
    rxn_likelihood_ps[num_rxns+num_rxns] = vall;
    /*
      Unimultiplier is 1.0/2^31-1
    */
    scaling = vall*uni_multiplier;
    /*
      choice is a pseudo-random integer in [0,2^31-1] (inclusive).
     */
    choice  = vgrng(vgrng_state);
    dchoice = ((double)choice)*scaling;
    /*
      Find index of smallest rxn_likelihood_ps entry that is >= choice.
    */
    rxn_choice = binary_search_l_u_b(rxn_likelihood_ps,dchoice,num_rxns_t2_p1);
    if (rxn_choice < num_rxns) {
      success = forward_rxn_conc_update(rxn_choice,state);
    } else {
      if (rxn_choice < num_rxns_t2) {
	success = reverse_rxn_conc_update(rxn_choice-num_rxns,state);
      } /* else do nothing, no change */
    }
    for (j=0;j<nu;j++) {
      cconcs[j] = fconcs[j];
    }
    /*
      Doug wonders here whether or not we ought to update free_energy based
      on the new concentrations. Absolutely we should.
    */
    success = log_rxn_ratios(free_energy,cconcs,lthermo,c_lrr,state);
    dg_forward = 0.0;
    for (j=0;j<num_rxns;j++) {
      free_energy[j] =  m_rt * c_lrr[j];
      dg_forward     += c_lrr[j];
    }
    dg_forward *= m_rt;
    /*
      Now we need to generate the relevant output.
      We want one file for the reaction Likelihoods,
      one file for the reaction concentrations.
      We will include the 
    */
    /* 
      print the concentrations. 
    */
    fprintf(state->concs_out_fp,"%d",i);
    for (j=0;j<num_rxns;j++) {
      fprintf(state->concs_out_fp," %le",cconcs[j]);
    }
    fprintf(state->concs_out_fp,"%\n");
    
    /* 
      print the dg_forward and the reaction likelihoods, in lthermo.
    */
    fprintf(state->rxn_lklhd_fp,"%d %le",i,dg_forward);
    for (j=0;j<num_rxns;j++) {
      fprintf(state->rxn_lklhd_fp," %le",lthermo[j]);
    }
    fprintf(state->rxn_lklhd_fp,"%\n");


  } /* end for(i...) */
  return(success);
}
