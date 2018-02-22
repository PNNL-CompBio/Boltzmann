#include "boltzmann_structs.h"
#include "boltzmann_cvodes_headers.h"
#include "cvodes_params_struct.h"
#include "boltzmann_free_state.h"
int boltzmann_free_state(struct state_struct *state) {
  /*
    Free the space allocated in boltzmann_init 
    as part of the state struct includeing the state struct itself 
    Called by: boltzmann_test_save_load
    Calls:     free
  */
  struct molecules_matrix_struct *molecules_matrix;
  struct reactions_matrix_struct *reactions_matrix;
  int success;
  int padi;
  success = 1;
  if (state != NULL) {
    /* 
      freeing in reverse order of allocating: alloc9
    */
    if (state->no_op_likelihood != NULL) {
      free((void*)state->no_op_likelihood);
    }
    if (state->rxn_view_likelihoods != NULL) {
      free((void*)state->rxn_view_likelihoods);
    }
    if (state->rev_rxn_view_likelihoods != NULL) {
      free((void*)state->rev_rxn_view_likelihoods);
    }
    if (state->rxn_fire != NULL) {
      free((void*)state->rxn_fire);
    }
    if (state->rxn_mat_row != NULL) {
      free((void*)state->rxn_mat_row);
    }
    /*
      alloc 8
    */
    if (state->future_counts != NULL) {
      free((void*)state->future_counts);
    }
    if (state->free_energy != NULL) {
      free((void*)state->free_energy);
    }
    if (state->forward_rxn_likelihood != NULL) {
      free((void*)state->forward_rxn_likelihood);
    }
    if (state->reverse_rxn_likelihood != NULL) {
      free((void*)state->reverse_rxn_likelihood);
    }
    if (state->forward_rxn_log_likelihood_ratio != NULL) {
      free((void*)state->forward_rxn_log_likelihood_ratio);
    }
    if (state->reverse_rxn_log_likelihood_ratio != NULL) {
      free((void*)state->reverse_rxn_log_likelihood_ratio);
    }
    if (state->rxn_likelihood_ps != NULL) {
      free((void*)state->rxn_likelihood_ps);
    }
    /*
      alloc7
    */
    if (state->ode_counts != NULL) {
      free((void*)state->ode_counts);
    } 
    if (state->ode_concs != NULL) {
      free((void*)state->ode_concs);
    }
    if (state->reactant_term != NULL) {
      free((void*)state->reactant_term);
    }
    if (state->product_term  != NULL) {
      free((void*)state->product_term);
    }
    if (state->rxn_q != NULL) {
      free((void*)state->rxn_q);
    }
    if (state->recip_rxn_q != NULL) {
      free((void*)state->recip_rxn_q);
    }
    if (state->log_kf_rel != NULL) {
      free((void*)state->log_kf_rel);
    }
    if (state->log_kr_rel != NULL) {
      free((void*)state->log_kr_rel);
    }
    if (state->ode_forward_lklhds != NULL) {
      free((void*)state->ode_forward_lklhds);
    }
    if (state->ode_reverse_lklhds != NULL) {
      free((void*)state->ode_reverse_lklhds);
    }
    if (state->rxn_has_flux != NULL) {
      free((void*)state->rxn_has_flux);
    }
    if (state->base_reactant_indicator != NULL) {
      free((void*)state->base_reactant_indicator);
    }
    if (state->rfc != NULL) {
      free((void*)state->rfc);
    }
    if (state->dfdke_dfdmu0_work != NULL) {
      free((void*)state->dfdke_dfdmu0_work);
    }
    /*
      NB ode23tb allocates a large block of worksepace
      that it now frees on termination. ODE23tb really only
      meant to be called in initialization phase.
    */
    /*
      alloc6 has been deprecated.
    */
    /*
      alloc5, allocations done in alloc 5 are freed in
      compute standard_energies. They won't be available here.
    */
    /*
      alloc4
    */
    molecules_matrix = state->molecules_matrix;
    if (molecules_matrix != NULL) {
      if (molecules_matrix->molecules_ptrs != NULL) {
	free((void*)molecules_matrix->molecules_ptrs);
      }
      if (molecules_matrix->reaction_indices != NULL) {
	free((void*)molecules_matrix->reaction_indices);
      }
      if (molecules_matrix->coefficients != NULL) {
	free((void*)molecules_matrix->coefficients);
      }
      if (molecules_matrix->recip_coeffs != NULL) {
	free((void*)molecules_matrix->recip_coeffs);
      }
      free((void*)molecules_matrix);
    }
    if (state->transpose_workspace != NULL) {
      free((void*)state->transpose_workspace);
    }
    /*
      alloc3
    */
    if (state->compartment_ptrs != NULL) {
      free((void*)state->compartment_ptrs);
    }
    if (state->current_counts != NULL) {
      free((void*)state->current_counts);
    }
    if (state->bndry_flux_counts != NULL) {
      free((void*)state->bndry_flux_counts);
    }
    if (state->net_lklhd_bndry_flux != NULL) {
      free((void*)state->net_lklhd_bndry_flux);
    }
    if (state->count_to_conc != NULL) {
      free((void*)state->count_to_conc);
    }
    if (state->conc_to_count != NULL) {
      free((void*)state->conc_to_count);
    }
    if (state->net_likelihood != NULL) {
      free((void*)state->net_likelihood);
    }
    if (state->dg0s != NULL) {
      free((void*)state->dg0s);
    }
    if (state->ke != NULL) {
      free((void*)state->ke);
    }
    if (state->kss != NULL) {
      free((void*)state->kss);
    }
    if (state->kss_e_val != NULL) {
      free((void*)state->kss_e_val);
    }
    if (state->kss_u_val != NULL) {
      free((void*)state->kss_u_val);
    }
    if (state->molecule_dg0tfs != NULL) {
      free((void*)state->molecule_dg0tfs);
    }
    if (state->molecule_probabilities != NULL) {
      free((void*)state->molecule_probabilities);
    }
    if (state->molecule_chemical_potentials != NULL) {
      free((void*)state->molecule_chemical_potentials);
    }
    /*
      alloc2_a:
    */
    if (state->rxn_title_text != NULL) {
      free((void*)state->rxn_title_text);
    }
    /*
      alloc2
    */
    if (state->reactions != NULL) {
      free((void*)state->reactions);
    }
    reactions_matrix = state->reactions_matrix;
    if (reactions_matrix != NULL) {
      if (reactions_matrix->rxn_ptrs != NULL) {
	free((void*)reactions_matrix->rxn_ptrs);
      }
      if (reactions_matrix->molecules_indices != NULL) {
	free((void*)reactions_matrix->molecules_indices);
      }
      if (reactions_matrix->compartment_indices != NULL) {
	free((void*)reactions_matrix->compartment_indices);
      }
      if (reactions_matrix->coefficients != NULL) {
	free((void*)reactions_matrix->coefficients);
      }
      if (reactions_matrix->recip_coeffs != NULL) {
	free((void*)reactions_matrix->recip_coeffs);
      }
      if (reactions_matrix->solvent_coefficients != NULL) {
	free((void*)reactions_matrix->solvent_coefficients);
      }
      if (reactions_matrix->text != NULL) {
	free((void*)reactions_matrix->text);
      }
      free((void*)reactions_matrix);
    }
    if (state->sorted_molecules != NULL) {
      free((void*)state->sorted_molecules);
    }
    if (state->sorted_compartments != NULL) {
      free((void*)state->sorted_compartments);
    }
    if (state->activities != NULL) {
      free((void*)state->activities);
    }
    if (state->enzyme_level != NULL) {
      free((void*)state->enzyme_level);
    }
    if (state->forward_rc != NULL) {
      free((void*)state->forward_rc);
    }
    if (state->reg_constant != NULL) {
      free((void*)state->reg_constant);
    }
    if (state->reg_exponent != NULL) {
      free((void*)state->reg_exponent);
    }
    if (state->reg_drctn != NULL) {
      free((void*)state->reg_drctn);
    }
    if (state->reg_species != NULL) {
      free((void*)state->reg_species);
    }
    if (state->coeff_sum != NULL) {
      free((void*)state->coeff_sum);
    }
    if (state->use_rxn != NULL) {
      free((void*)state->use_rxn);
    }
    /*
      alloc0_a
    */
    if (state->params_file != NULL) {
      free((void*)state->params_file);
    }
    if (state->solvent_string != NULL) {
      free((void*)state->solvent_string);
    }
    /*
      alloc0
    */
    if (state->param_buffer != NULL) {
      free((void*)state->param_buffer);
    }
    if (state->rxn_file_keyword_buffer != NULL) {
      free((void*)state->rxn_file_keyword_buffer);
    }
    if (state->rxn_file_keywords != NULL) {
      free((void*)state->rxn_file_keywords);
    }
    if (state->rxn_file_keyword_lengths != NULL) {
      free((void*)state->rxn_file_keyword_lengths);
    }
    if (state->vgrng_state != NULL) {
      free((void*)state->vgrng_state);
    }
    if (state->vgrng2_state != NULL) {
      free((void*)state->vgrng2_state);
    }
    if (state->cvodes_params != NULL) {
      free((void*)state->cvodes_params);
    }
    free((void*)state);
  }
  return(success);
}
   
    
