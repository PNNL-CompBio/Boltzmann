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
#include "choose_rxn.h"
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
	       choose_rxn,
	       rxn_log_likelihoods,
	       print_restart_file
	       print_reactions_view
  */
  struct vgrng_state_struct *vgrng_state;
  struct istring_elem_struct *molecule;
  struct istring_elem_struct *sorted_molecules;
  struct istring_elem_struct *cur_cmpts;
  struct istring_elem_struct *cur_cmpt;

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
  double *activities;
  double *no_op_likelihood;
  double cal_gm_per_joule;
  double delta;
  /*
  double *lthermo;
  */
  double m_rt;
  int    *rxn_fire;
  char   *cmpt_string;
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

  int rxn_view_step;
  int rxn_view_pos;

  int rxn_view_freq;
  int rxn_view_hist_lngth;

  int lklhd_view_step;
  int lklhd_view_freq;

  int conc_view_step;
  int conc_view_freq;

  int ci;
  int oi;

  FILE *lfp;
  FILE *concs_out_fp;
  FILE *rxn_lklhd_fp;
  FILE *bndry_flux_fp;
  success = 1;
  forward = 1;
  
  n_warmup_steps    = state->warmup_steps;
  n_record_steps    = state->record_steps;
  vgrng_state       = state->vgrng_state;
  rxn_likelihood_ps = state->rxn_likelihood_ps;
  free_energy       = state->free_energy;
  c_loglr           = state->current_rxn_log_likelihood_ratio;
  num_rxns          = state->number_reactions;
  nu                = state->unique_molecules;
  cconcs            = state->current_concentrations;
  fconcs            = state->future_concentrations;
  bndry_flux_concs  = state->bndry_flux_concs;
  activities        = state->activities;
  m_rt              = state->m_rt;
  rxn_view_data     = state->rxn_view_likelihoods;
  rrxn_view_data    = state->rev_rxn_view_likelihoods;
  cal_gm_per_joule  = state->cal_gm_per_joule;
  rxn_fire          = state->rxn_fire;
  no_op_likelihood  = state->no_op_likelihood;
  /*
  lthermo        = state->l_thermo;
  */
  forward_rxn_likelihood = state->forward_rxn_likelihood;
  reverse_rxn_likelihood = state->reverse_rxn_likelihood;
  num_rxns_t2          	 = num_rxns << 1;
  num_rxns_t2_p1       	 = num_rxns_t2 + 1;
  concs_out_fp         	 = state->concs_out_fp;
  rxn_lklhd_fp         	 = state->rxn_lklhd_fp;
  bndry_flux_fp        	 = state->bndry_flux_fp;
  sorted_molecules     	 = (struct istring_elem_struct *)state->sorted_molecules;
  cur_cmpts              = (struct istring_elem_struct *)state->sorted_cmpts;
  rxn_view_freq        	 = state->rxn_view_freq;
  rxn_view_hist_lngth  	 = state->rxn_view_hist_lngth;
  lklhd_view_freq        = state->lklhd_view_freq;
  conc_view_freq         = state->conc_view_freq;
  rxn_view_pos         	 = 0;
  rxn_view_step        	 = 1;
  conc_view_step         = 1;
  lklhd_view_step 	 = 1;
  for (i=0;i<n_warmup_steps;i++) {
    /*
      Compute the reaction likelihoods - forward_rxn_likelihood, 
      and reverse_rxn_likelihood.
    */
    success = free_energy_changes(state);
    /*
      Compute the partial sums of the reaction likelihoods.
      This call also updates the fconcs vector.
    */
    rxn_choice = choose_rxn(state);
    if (rxn_choice < 0) break;

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
      sum_likelihood += activities[j]*forward_rxn_likelihood[j];
    }
    dg_forward *= m_rt;
    entropy = 0.0;
    if (sum_likelihood <= 0.0) {
      fprintf(stderr,"boltzmann_run_sim: Error, nonpositivity sum_likelihood = %le in warmup loop iteration %d\n",sum_likelihood,i);
      fflush(stderr);
      success = 0;
      break;
    }
    if (success) {
      r_sum_likelihood = 1.0/sum_likelihood;
      for (j=0;j<num_rxns;j++) {
	scaled_likelihood = activities[j]*forward_rxn_likelihood[j] * r_sum_likelihood;
	if (scaled_likelihood>0) {
	  entropy -= (scaled_likelihood * log(scaled_likelihood));
	}
      }
    }
  } /* end for(i...) */
  /* 
    Data collection phase (recording).
  */
  if (success) {
    /*
      Set the rxn_fire counts to 0.
    */
    for (i=0;i<num_rxns_t2;i++) {
      rxn_fire[i] = 0;
    }
    /*
      Initialize the boundary fluxes 0.
    */
    for (i=0;i<nu;i++) {
      bndry_flux_concs[i] = 0.0;
    }
    for (i=0;i<n_record_steps;i++) {
      /*
	Comput the reaction likelihoods - forward_rxn_likelihood, 
	and reverse_free_energy.
      */
      success = free_energy_changes(state);

      rxn_choice = choose_rxn(state);
      if (rxn_choice < 0) break;
      if (rxn_choice <= num_rxns_t2) {
	rxn_fire[rxn_choice] += 1;
      }

      for (j=0;j<nu;j++) {
	cconcs[j] = fconcs[j];
      }
      /*
	Doug wonders here whether or not we ought to update free_energy based
	on the new concentrations. Absolutely we should.
      */
      success = rxn_log_likelihoods(free_energy,cconcs,
				    forward_rxn_likelihood,c_loglr,state,
				    forward);
      dg_forward = 0.0;
      sum_likelihood = 0.0;
      for (j=0;j<num_rxns;j++) {
	free_energy[j] =  m_rt * c_loglr[j];
	dg_forward     += c_loglr[j];
	sum_likelihood += activities[j]*forward_rxn_likelihood[j];
      }
      dg_forward *= m_rt;
      entropy = 0.0;
      if (sum_likelihood <= 0.0) {
	fprintf(stderr,"boltzmann_run_sim: Error, nonpositivity sum_likelihood = %le in recording loop iteration %d\n",sum_likelihood,i);
	fprintf(stderr,
		"boltzmann_run_sim: Error, nonpositivity sum_likelihood\n");
	fflush(stderr);
	success = 0;
	break;
      }
      r_sum_likelihood = 0.0;
      if (success) {
	r_sum_likelihood = 1.0/sum_likelihood;
	for (j=0;j<num_rxns;j++) {
	  scaled_likelihood = activities[j]*forward_rxn_likelihood[j] * r_sum_likelihood;
	  if (scaled_likelihood > 0) {
	    entropy -= scaled_likelihood * log(scaled_likelihood);
	  }
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
	conc_view_step = conc_view_step - 1;
	if ((conc_view_step <= 0) || (i == (n_record_steps-1))) {
	  fprintf(concs_out_fp,"%d",i);
	  for (j=0;j<nu;j++) {
	    fprintf(state->concs_out_fp,"\t%le",cconcs[j]);
	  }
	  fprintf(state->concs_out_fp,"\n");
	  conc_view_step = conc_view_freq;
	}
      }
      /* 
	 print the entropy, dg_forward and the reaction likelihoods, 
      */
      if (state->rxn_lklhd_fp) {
	lklhd_view_step = lklhd_view_step - 1;
	if ((lklhd_view_step <= 0) || (i == (n_record_steps-1))) {
	  fprintf(state->rxn_lklhd_fp,"%d\t%le\t%le",i,entropy,dg_forward);
	  for (j=0;j<num_rxns;j++) {
	    /*
	    fprintf(state->rxn_lklhd_fp,"\t%le",forward_rxn_likelihood[j]);
	    fprintf(state->rxn_lklhd_fp,"\t%le",reverse_rxn_likelihood[j]);
	    */
	    fprintf(state->rxn_lklhd_fp,"\t%le",forward_rxn_likelihood[j]*activities[j]);
	    fprintf(state->rxn_lklhd_fp,"\t%le",reverse_rxn_likelihood[j]*activities[j]);
	  }
	  fprintf(state->rxn_lklhd_fp,"\n");
	  lklhd_view_step = lklhd_view_freq;
	}
      }
      if (rxn_view_freq > 0) {
	rxn_view_step = rxn_view_step - 1;
	/*
	  Save the likelihoods on a per reaction basis.
	*/
	if ((rxn_view_step <= 0) || (i == (n_record_steps-1))) {
	  no_op_likelihood[rxn_view_pos] = r_sum_likelihood;
	  rxn_view_p = (double *)&rxn_view_data[rxn_view_pos];
	  rrxn_view_p = (double *)&rrxn_view_data[rxn_view_pos];
	  for (j = 0; j < num_rxns;j++) {
	    /*
	    *rxn_view_p = forward_rxn_likelihood[j];
	    *rrxn_view_p = reverse_rxn_likelihood[j];
	    */
	    *rxn_view_p  = forward_rxn_likelihood[j]*activities[j];
	    *rrxn_view_p = reverse_rxn_likelihood[j]*activities[j];
	    rxn_view_p   += rxn_view_hist_lngth; /* Caution address arithmetic here. */
	    rrxn_view_p  += rxn_view_hist_lngth; /* Caution address arithmetic here. */
	  }
	  rxn_view_step = rxn_view_freq;
	  rxn_view_pos  += 1;
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
	fprintf(bndry_flux_fp,"final flux\n");
	molecule    = sorted_molecules;
	cmpt_string = NULL;
	oi = -1;
	for (j=0;j<nu;j++) {
	  ci = molecule->c_index;
	  if (ci != oi) {
	    oi = ci;
	    if (ci >= 0) {
	      cur_cmpt = (struct istring_elem_struct *)&(cur_cmpts[ci]);
	      cmpt_string = cur_cmpt->string;
	    } else {
	      cmpt_string = NULL;
	    }
	  }
	  if (molecule->variable == 0) {
	    if (ci >= 0) {
	      fprintf(state->bndry_flux_fp,"%s:%s\t%le\n",molecule->string,
		      cmpt_string,bndry_flux_concs[j]);
	    } else {
	      fprintf(state->bndry_flux_fp,"%s\t%le\n",molecule->string,
		      bndry_flux_concs[j]);
	    }
	  }
	  molecule += 1; /* Caution address arithmetic.*/
	} /* end for (j...) */
	/*
	fprintf(bndry_flux_fp," flux - fixed ");
	molecule = sorted_molecules;
	for (j=0;j<nu;j++) {
	  if (molecule->variable == 0) {
	    delta = bndry_flux_concs[j]-cconcs[j];
	    fprintf(state->bndry_flux_fp,"\t%le",delta);
	  }
	  molecule += 1; 
	}
	fprintf(state->bndry_flux_fp,"\n");
	*/
      } /* end if (bndry_flux_fp) */
    } /* end if (state->num_fixed_concs ...) */
    if (success) {
      success = print_restart_file(state);
    }
    if (success) {
      if (rxn_view_freq > 0) {
	success = print_reactions_view(state);
      }
    }
  }
  return(success);
}
