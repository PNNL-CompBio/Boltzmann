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
#include <math.h>
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
#include "rxn_log_likelihoods.h"
#include "print_restart_file.h"
#include "print_reactions_view.h"

#include "boltzmann_run_sim.h"
int boltzmann_run_sim(struct state_struct *state) {
  /*
    Run the boltzmann simulations, to be called after
    boltzmann_init has been called.

    Called by: boltzmann_drv
    Calls:     free_energy_changes,
               vgrng,
	       binary_search_l_u_b,
	       forward_rxn_conc_update,
	       reverse_rxn_conc_update,
	       rxn_log_likelihoods,
  */
  struct vgrng_state_struct *vgrng_state;
  struct istring_elem_struct *molecule;
  struct istring_elem_struct *sorted_molecules;

  int64_t choice;
  double dchoice;
  double uni_multiplier;
  double vall;
  double rvall;
  double scaling;
  double dg_forward;
  double sum_likelihood;
  double r_sum_likelihood;
  double entropy;
  double scaled_likelihood;
  double fe;
  double *cconcs;
  double *fconcs;
  double *rxn_likelihood_ps;
  double *free_energy;
  double *c_loglr;
  double *forward_rxn_likelihood;
  double *reverse_rxn_likelihood;
  double *rxn_view_data;
  double *rxn_view_p;
  double *rrxn_view_data;
  double *rrxn_view_p;
  double *bndry_flux_concs;
  double cal_gm_per_joule;
  double delta;
  /*
  double *lthermo;
  */
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
  int forward;

  int view_step;
  int view_pos;

  int view_freq;
  int lthf;

  int lthl;
  int padi;

  FILE *lfp;
  FILE *concs_out_fp;
  FILE *rxn_lklhd_fp;
  FILE *bndry_flux_fp;
  success = 1;
  forward = 1;
  
  n_warmup_steps = state->warmup_steps;
  n_record_steps = state->record_steps;
  vgrng_state    = state->vgrng_state;
  rxn_likelihood_ps         = state->rxn_likelihood_ps;
  free_energy            = state->free_energy;
  c_loglr          = state->current_rxn_log_likelihood_ratio;
  num_rxns       = state->number_reactions;
  nu             = state->unique_molecules;
  cconcs         = state->current_concentrations;
  fconcs         = state->future_concentrations;
  bndry_flux_concs  = state->bndry_flux_concs;

  m_rt           = state->m_rt;
  rxn_view_data  = state->rxn_view_likelihoods;
  rrxn_view_data  = state->rev_rxn_view_likelihoods;
  cal_gm_per_joule = state->cal_gm_per_joule;
  /*
  lthermo        = state->l_thermo;
  */
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  num_rxns_t2    = num_rxns << 1;
  num_rxns_t2_p1 = num_rxns_t2 + 1;
  uni_multiplier = vgrng_state->uni_multiplier;
  concs_out_fp   = state->concs_out_fp;
  rxn_lklhd_fp   = state->rxn_lklhd_fp;
  bndry_flux_fp  = state->bndry_flux_fp;
  sorted_molecules = (struct istring_elem_struct *)state->sorted_molecules;
  lthf           = state->lthf;
  lthl           = state->lthl;
  view_pos       = 0;
  view_step      = 1;
  for (i=0;i<n_warmup_steps;i++) {
    /*
      Compute the reaction likelihoods - forward_rxn_likelihood, 
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
      1.0 is add to the likelihoods to account for the 
      likeilhood that the state does not change.
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
      on the new concentrations, Absolutely we should.
    */
    success = rxn_log_likelihoods(free_energy,cconcs,c_loglr,c_loglr,state,
				  forward);
    dg_forward = 0.0;
    sum_likelihood = 0.0;
    for (j=0;j<num_rxns;j++) {
      dg_forward += c_loglr[j];
      sum_likelihood += forward_rxn_likelihood[j];
    }
    dg_forward *= m_rt;
    entropy = 0.0;
    if (sum_likelihood <= 0.0) {
      fprintf(stderr,"boltzmann_run_sim: Error, nonpositivity sum_likelihood\n");
      fflush(stderr);
      success = 0;
      break;
    }
    if (success) {
      r_sum_likelihood = 1.0/sum_likelihood;
      for (j=0;j<num_rxns;j++) {
	scaled_likelihood = forward_rxn_likelihood[j] * r_sum_likelihood;
	if (scaled_likelihood>0) {
	  entropy -= (scaled_likelihood * log(scaled_likelihood));
	}
      }
    }
  } /* end for(i...) */
  
  if (success) {
    /*
      Reinitialize the boundary fluxes to the current concentrations.
    */
    for (i=0;i<nu;i++) {
      bndry_flux_concs[i] = cconcs[i];
    }
    for (i=0;i<n_record_steps;i++) {
      /*
	Comput the reaction likelihoods - forward_rxn_likelihood, 
	and reverse_free_energy.
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
	rxn_likelihood_ps[num_rxns+j] = rxn_likelihood_ps[num_rxns-1+j] + 
	  reverse_rxn_likelihood[j];
      }
      /*
	1.0 is add to the likelihoods to account for the 
	likeilhood that the state does not change.
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
      rxn_choice = binary_search_l_u_b(rxn_likelihood_ps,dchoice,
				       num_rxns_t2_p1);
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
	success = rxn_log_likelihoods(free_energy,cconcs,lthermo,c_loglr,state);      
      */
      success = rxn_log_likelihoods(free_energy,cconcs,
				    forward_rxn_likelihood,c_loglr,state,
				    forward);
      dg_forward = 0.0;
      sum_likelihood = 0.0;
      for (j=0;j<num_rxns;j++) {
	free_energy[j] =  m_rt * c_loglr[j];
	dg_forward     += c_loglr[j];
	sum_likelihood += forward_rxn_likelihood[j];
      }
      dg_forward *= m_rt;
      entropy = 0.0;
      if (sum_likelihood <= 0.0) {
	fprintf(stderr,
		"boltzmann_run_sim: Error, nonpositivity sum_likelihood\n");
	fflush(stderr);
	success = 0;
	break;
      }
      if (success) {
	r_sum_likelihood = 1.0/sum_likelihood;
	for (j=0;j<num_rxns;j++) {
	  scaled_likelihood = forward_rxn_likelihood[j] * r_sum_likelihood;
	  entropy -= scaled_likelihood * log(scaled_likelihood);
	}
      }
      /*
	Now we need to generate the relevant output.
	We want one file for the reaction Likelihoods,
	one file for the reaction concentrations.
	We will include the 
      */
      /* 
	 print the concentrations. 
      */
      if (concs_out_fp) {
	fprintf(concs_out_fp,"%d",i);
	for (j=0;j<nu;j++) {
	  fprintf(state->concs_out_fp,"\t%le",cconcs[j]);
	}
	fprintf(state->concs_out_fp,"\n");
      }
      /* 
	 print the entropy, dg_forward and the reaction likelihoods, 
      */
      fprintf(state->rxn_lklhd_fp,"%d\t%le\t%le",i,entropy,dg_forward);
      for (j=0;j<num_rxns;j++) {
	fprintf(state->rxn_lklhd_fp,"\t%le",forward_rxn_likelihood[j]);
	fprintf(state->rxn_lklhd_fp,"\t%le",reverse_rxn_likelihood[j]);
      }
      fprintf(state->rxn_lklhd_fp,"\n");
      if (lthf > 0) {
	view_step = view_step - 1;
	/*
	  Save the likelihoods on a per reaction basis.
	*/
	if ((view_step <= 0) || (i == n_record_steps-1)) {
	  rxn_view_p = (double *)&rxn_view_data[view_pos];
	  rrxn_view_p = (double *)&rrxn_view_data[view_pos];
	  for (j = 0; j < num_rxns;j++) {
	    *rxn_view_p = forward_rxn_likelihood[j];
	    *rrxn_view_p = reverse_rxn_likelihood[j];
	    rxn_view_p += lthl; /* Caution address arithmetic here. */
	    rrxn_view_p += lthl; /* Caution address arithmetic here. */
	  }
	  view_step = lthf;
	  view_pos += 1;
	}
      }
      /*
	If user has requested print out free energies as well.
      */
      if (state->free_energy_format > 0) {
	fprintf(state->free_energy_fp,"%d",i);
	if (state->free_energy_format == 1) {
	  for (j=0;j<num_rxns;j++) {
	    fprintf(state->free_energy_fp,"\t%le",
		    -c_loglr[j]);
	  }
	  fprintf(state->free_energy_fp,"\n");
	}
	if (state->free_energy_format == 2) {
	  for (j=0;j<num_rxns;j++) {
	    fe = free_energy[j]*cal_gm_per_joule;
	    fprintf(state->free_energy_fp,"\t%le",fe);
	  }
	  fprintf(state->free_energy_fp,"\n");
	}
	if (state->free_energy_format == 3) {
	  for (j=0;j<num_rxns;j++) {
	    fprintf(state->free_energy_fp,"\t%le",free_energy[j]);
	  }
	  fprintf(state->free_energy_fp,"\n");
	}
      }
    } /* end for(i...) */
    if (state->num_fixed_concs > 0) {
      if (bndry_flux_fp) {
	fprintf(bndry_flux_fp,"  final flux  ");
	molecule = sorted_molecules;
	for (j=0;j<nu;j++) {
	  if (molecule->variable == 0) {
	    fprintf(state->bndry_flux_fp,"\t%le",bndry_flux_concs[j]);
	  }
	  molecule += 1; /* Caution address arithmetic.*/
	}
	fprintf(state->bndry_flux_fp,"\n");
	fprintf(bndry_flux_fp," flux - fixed ");
	molecule = sorted_molecules;
	for (j=0;j<nu;j++) {
	  if (molecule->variable == 0) {
	    delta = bndry_flux_concs[j]-cconcs[j];
	    fprintf(state->bndry_flux_fp,"\t%le",delta);
	  }
	  molecule += 1; /* Caution address arithmetic.*/
	}
	fprintf(state->bndry_flux_fp,"\n");
      }
    }
    if (success) {
      success = print_restart_file(state);
    }
    if (success) {
      if (lthf > 0) {
	success = print_reactions_view(state);
      }
    }
  }
  return(success);
}
