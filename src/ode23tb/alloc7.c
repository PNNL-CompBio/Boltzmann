/* alloc7.c
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

#include "alloc7.h"

int alloc7(struct state_struct *state) {
  /*
    Allocace space needed by the deq routines.
    Allocate space for a flux vector ( lenth = num_species)
    The Jacobian of the flux vector (length = num_species * num_species);
    the reactant_term vector (length = num_rxns), 
         product_term_vector (length = num_rxns),  
	 recip_rxn_q, p_over_r both of 
	 length num_rxns, 
    and concs of length num_species.
    Also an integer vector, rxn_has_flux of length num_rxns.

    Called by: deq_run
    Calls:     calloc, fprintf, fflush,
    Sets the following fields of state:
      ode_counts,
      ode_concs,
      reactant_term,
      product_term,
      rxn_q,
      recip_rxn_q,
      log_kf_rel,
      log_kr_rel,
      ode_forward_lklhds,
      ode_reverse_lklhds,
      rxn_has_flux,
      base_reactant_indicator
  */
  double *reactant_term;
  double *product_term;
  double *rfc;
  double *deriv_acc;
  double *stable_add_scr;
  double *rxn_q;
  double *recip_rxn_q;
  double *log_kf_rel;
  double *log_kr_rel;
  double *ode_forward_lklhds;
  double *ode_reverse_lklhds;
  double *ode_counts;
  double *ode_concs;
  double *ode_f;
  double *ode_kq;
  double *ode_kqi;
  double *dfdke_dfdmu0_work;
  int64_t ask_for;
  int64_t one_l;
  int64_t usage;
  int64_t run_workspace_bytes;
  int    *rxn_has_flux;
  int    *base_reactants;
  int    *base_reactant_indicator;
  int num_rxns;
  int num_species;
  int success;
  int padi;
  FILE *lfp;
  FILE *efp;

  one_l         = (int64_t)1;
  usage         = state->usage;
  num_species   = (int)state->nunique_molecules;
  num_rxns      = (int)state->number_reactions;
  lfp           = state->lfp;
  success       = 1;
  run_workspace_bytes  = state->run_workspace_bytes;
  /*
    Allocate space for the ode_counts vector to be computed
    from the molecule concentrations and compartment sizes.
  */
  if (success) {
    ask_for = num_species * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_counts = (double *)calloc(one_l,ask_for);
    if (ode_counts == NULL) {
      fprintf(stderr,"alloc7: Error unable to allocate %ld bytes for "
	      "ode_counts\n",ask_for);
      fflush(stderr);
      success = 0;
    } else {
      state->ode_counts = ode_counts;
    }
  }
  if (success) {
    ask_for = num_species * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_concs = (double *)calloc(one_l,ask_for);
    if (ode_concs == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"ode_concs\n",ask_for);
	fflush(lfp);
      }
    } else {
      state->ode_concs = ode_concs;
    }
  }
  /*
    Allocate space for the ode flux vector used in printing.
  */
  if (success) {
    ask_for = num_species * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_f = (double *)calloc(one_l,ask_for);
    if (ode_f == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"ode_f\n",ask_for);
	fflush(lfp);
      }
    } else {
      state->ode_f = ode_f;
    }
  }
  /*
    allocate space for the ode_kq and ode_kqi vectors for printing.
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    ask_for = ask_for + ask_for;
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_kq = (double*)calloc(one_l,ask_for);
    if (ode_kq == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"ode_kq and ode_kqi\n",ask_for);
	fflush(lfp);
      }
    } else {
      ode_kqi = (double*)&ode_kq[num_rxns];
      state->ode_kq = ode_kq;
      state->ode_kqi = ode_kqi;
    }
  }
  /*
    Allocate space for the reactant_terms
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    reactant_term = (double *)calloc(one_l,ask_for);
    if (reactant_term == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
	      "reactant_term\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->reactant_term = reactant_term;
    }
  }
  /*
    Allocate space for the product_term
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    product_term = (double *)calloc(one_l,ask_for);
    if (product_term == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"product_term\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->product_term = product_term;
    }
  }
  /*
    Allocate space for the rxn_q (reaction quotient)
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    rxn_q = (double *)calloc(one_l,ask_for);
    if (rxn_q == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
	      "rxn_q\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->rxn_q = rxn_q;
    }
  }
  /*
    Allocate space for the recip_rxn_q vector
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    recip_rxn_q = (double *)calloc(one_l,ask_for);
    if (recip_rxn_q == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"recip_rxn_q\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->recip_rxn_q = recip_rxn_q;
    }
  }
  /*
    Allocate space for the log_kf_rel vector
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    log_kf_rel = (double *)calloc(one_l,ask_for);
    if (log_kf_rel == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
	      "log_kf_rel\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->log_kf_rel = log_kf_rel;
    }
  }
  /*
    Allocate space for the log_kr_rel vector
  */
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    log_kr_rel = (double *)calloc(one_l,ask_for);
    if (log_kr_rel == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"log_kr_rel\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->log_kr_rel = log_kr_rel;
    }
  }
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_forward_lklhds = (double *)calloc(one_l,ask_for);
    if (ode_forward_lklhds == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"ode_forward_lklhds\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->ode_forward_lklhds = ode_forward_lklhds;
    }
  }
  if (success) {
    ask_for = num_rxns * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    ode_reverse_lklhds = (double *)calloc(one_l,ask_for);
    if (ode_reverse_lklhds == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
		"ode_reverse_lklhds\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->ode_reverse_lklhds = ode_reverse_lklhds;
    }
  }
  /*
    Allocate space for the rxn_has_flux indicator vector.
    This vector is sort of a hack for now for improving mass
    conservation. We need to address mass conservation more carefully.
  */
  if (success) {
    ask_for = (num_rxns + (num_rxns & 1)) * sizeof(int);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    rxn_has_flux = (int*)calloc(one_l,ask_for);
    if (rxn_has_flux == NULL) {
      if (lfp) {
	fprintf(lfp,"alloc7: Error unable to allocate %ld bytes for "
	      "rxn_has_flux\n",ask_for);
	fflush(lfp);
      }
      success = 0;
    } else {
      state->rxn_has_flux = rxn_has_flux;
    }
  }
  if (success) {
    ask_for = (int64_t)(num_species+num_species)*sizeof(int);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    base_reactant_indicator = (int*)calloc(ask_for,one_l);
    if (base_reactant_indicator == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error could not allocate %ld "
		"bytes for int base_reactant scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    base_reactants = (int*)&base_reactant_indicator[num_species];
    state->base_reactant_indicator = base_reactant_indicator;
    state->base_reactants          = base_reactants;
  }
  if (success) {
    ask_for = (int64_t)(6*num_rxns) * sizeof(double);
    usage += ask_for;
    run_workspace_bytes  += ask_for;
    rfc = (double *)calloc(ask_for,one_l);
    if (rfc == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error could not allocate %ld "
		"bytes for int rfc and deriv_acc scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    deriv_acc = &rfc[num_rxns+num_rxns];
    stable_add_scr = &deriv_acc[num_rxns+num_rxns];
    state->rfc = rfc;
    state->deriv_acc = deriv_acc;
    state->stable_add_scr = stable_add_scr;
  }
  if (success) {
    /*
      Perhaps we should allocate this only if print_output = 1
      as the compute_dfdke_dfdmu0 routine only get called when
      print_output = 1.
    */
    ask_for = (int64_t)(3*num_species)*sizeof(double);
    usage += ask_for;
    run_workspace_bytes += ask_for;
    dfdke_dfdmu0_work = (double*)calloc(ask_for,one_l);
    if (dfdke_dfdmu0_work == NULL) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"alloc7: Error could not allocate %ld "
		"bytes for dfdke_dfdmu0_work scratch space.\n",ask_for);
	fflush(lfp);
      }
    }
  }
  if (success) {
    state->dfdke_dfdmu0_work = dfdke_dfdmu0_work;
  }
  state->usage = usage;
  state->run_workspace_bytes  = run_workspace_bytes;

  return(success);
}
