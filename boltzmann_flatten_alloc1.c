#include "boltzmann_structs.h"
#include "alloc2.h"
#include "alloc3.h"
#include "alloc4.h"
#include "alloc7.h"
#include "alloc8.h"
#include "alloc9.h"
#include "boltzmann_flatten_alloc1.h"
int boltzmann_flatten_alloc1(struct state_struct *state) {
  /*
    Allocate one-way, two-way and work fields of state for
    boltzmann_run.
    Called by: boltzmann_flatten_state
    Calls:     alloc2,
               alloc3,
	       alloc4,
	       alloc7,
	       alloc8,
	       alloc9,

    Sets the following fields in state:	       
    reactions,
    reactions_matrix,
      reactions_matrix->rxn_ptrs,
      reactions_matrix->molecules_indices,
      reactions_matrix->compartment_indices,
      reactions_matrix->coefficients,
      reactions_matrix->solvent_coefficients,
      reactions_matrix->text,
      sorted_molecules,
      sorted_compartments,
      activities,
      enzyme_level,
      forward_rc,
      reverse_rc,
      reg_constant,
      reg_exponent,
      reg_drctn,
      reg_species,
      vgrng_state,
      vgnrg2_state

      current_counts         	     
      bndry_flux_counts      	     
      net_lklhd_bndry_flux         
      count_to_conc          	     
      conc_to_count          	     
      net_likelihood               
      dg0s                   	     
      ke                     	     
      rke                          
      kss                    	     
      kssr                   	     
      kss_eval               	     
      kss_uval               	     
      dg0tfs                 	     
      molecule_probabilities 	     
      molecule_chemical_potentials 

      molecules_matrix
      molecules_matrix->molecules_ptrs,
      molecules_matrix->reaction_indices,
      molecules_matrix->coefficients

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

      future_counts,
      free_energy,
      forward_rxn_likelihood,
      reverse_rxn_likelihood,
      forward_rxn_log_likelihood_ratio,
      reverse_rxn_log_likelihood_ratio,
      rxn_likelihood_ps

      no_op_likelihood, 
      rxn_view_likelihoods,
      rev_rxn_view_likelihoods,
      rxn_fire,
      rxn_mat_row 
  */
  int success;
  int setup;
  setup = 0;
  /*
    Allocate the following fields.
    reactions,
    reactions_matrix,
    reactions_matrix->rxn_ptrs,
    reactions_matrix->molecules_indices,
    reactions_matrix->compartment_indices,
    reactions_matrix->coefficients,
    reactions_matrix->recip_coeffs,
    reactions_matrix->solvent_coefficients,
    reactions_matrix->text,
    sorted_molecules,
    sorted_compartments,
    activities,
    enzyme_level,
    forward_rc,
    reverse_rc,
    reg_constant,
    reg_exponent,
    reg_drctn,
    reg_species,
    vgrng_state,
    vgnrg2_state

  */
  success = alloc2(state,setup);
  if (success) {
    /*
  	Allocate the follwing fields:
  	current_counts         	     (nunique_molecules)
  	bndry_flux_counts      	     (nunique_molecules)
  	net_lklhd_bndry_flux         (nunique_molecules)
  	count_to_conc          	     (nunique_molecules)
  	conc_to_count          	     (nunique_molecules)
  	net_likelihood               (number_reactions)
  	dg0s                   	     (number_reactions)  
  	ke                     	     (number_reactions) 
  	rke                          (number_reactions) 
  	kss                    	     (number_reactions)
  	kssr                   	     (number_reactions)
  	kss_eval               	     (nunique_molecules)
  	kss_uval               	     (nunique_molecules)
  	dg0tfs                 	     (nunique_molecules)
  	molecule_probabilities 	     (nunique_molecules)
  	molecule_chemical_potentials (nunique_molecules)
    */
    success = alloc3(state,setup);
  }
  if (success) {
    /* 
  	Allocate space for the 
  	molecules_matrix
  	molecules_matrix->molecules_ptrs,
  	molecules_matrix->reaction_indices,
  	molecules_matrix->coefficients,
	molecules_matrix->recip_coeffs
    */
    success = alloc4(state,setup);
  }
  /*
    We don't call alloc5 or alloc6 as that allocates pseudoisomer space
    only used in setup.

    Allocate space needed by the deq routines.
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
  /*
    Now this routine only gets called if we are loading from a flattened
    state, in which case we need to 
  */
  if (success) {
    if(state->use_deq) {
      success = alloc7(state);
    }
  }
  /*
    Allocate workspace needed in boltzmann_run
    future_counts,
    free_energy,
    forward_rxn_likelihood,
    reverse_rxn_likelihood,
    forward_rxn_log_likelihood_ratio,
    reverse_rxn_log_likelihood_ratio,
    rxn_likelihood_ps
  */
  if (success) {
    success = alloc8(state);
  }
  /*
    if state->print_output)   allocate printing space:
      no_op_likelihood, 
      rxn_view_likelihoods,
      rev_rxn_view_likelihoods,
      rxn_fire,
      rxn_mat_row 
  */
  if (success) {
    if (state->print_output) {
  	success = alloc9(state);
    }
  }
  return(success);
}
