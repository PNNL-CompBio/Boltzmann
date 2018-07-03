#include "boltzmann_structs.h"
#include "boltzmann_set_state_ptrs.h"
void boltzmann_set_state_ptrs(struct state_struct *state,
			      int *length_changed) {
  /*
    Set the state pointer fields for a flattened (contiguous storage) state.
    If the computed sizes and offsets meta_size, two_way_data_offset, 
    incoming_data_offset, auxiliary_data_offset, workspace_offset, and
    workspace_length don't match up, the 1, 2, 4, 8, 16, or 32 bits
    correspondingly in *length_changed are set.
    Called by: flatten_state.

    Calls:     
  */
  struct state_struct ss;
  struct vgrng_state_struct vss;
  struct molecule_struct ies;
  struct compartment_struct ces;
  struct reaction_struct rs;
  struct reactions_matrix_struct rms;
  struct reactions_matrix_struct *reactions_matrix;
  struct molecules_matrix_struct mms;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t print_output;
  int64_t maX_regs_per_rxn;
  int64_t number_reactions;
  int64_t number_molecules;
  int64_t number_compartments;
  int64_t unique_compartments;
  int64_t unique_molecules;
  int64_t max_filename_len;
  int64_t num_files;
  int64_t meta_size;
  int64_t data_pad;
  int64_t data_len;
  int64_t align_len;
  int64_t align_mask;

  int64_t per_molecule_double_size;
  int64_t per_molecule_data_pad;
  int64_t padded_per_molecule_size;
  int64_t per_molecule_int_pad;
  int64_t per_molecule_int_size;
  int64_t padded_per_molecule_int_size;
  int64_t per_reaction_double_size;
  int64_t per_reaction_data_pad;
  int64_t padded_per_reaction_size;
  int64_t per_reg_size;
  int64_t per_reg_data_pad;
  int64_t padded_per_reg_size;

  int64_t per_compartment_double_size;
  int64_t vgrng_state_size;
  int64_t reaction_data_size;
  int64_t molecule_data_size;
  int64_t compartment_data_size;
  int64_t reactions_matrix_size;
  int64_t molecules_matrix_size;
  int64_t reactions_ptrs_size;
  int64_t molecules_ptrs_size;
  int64_t reactions_matrix_field_size;
  int64_t molecules_matrix_field_size;
  int64_t file_names_size;
  int64_t solvent_string_size;
  int64_t reaction_titles_length;
  int64_t pathway_text_length;
  int64_t compartment_text_length;
  int64_t molecule_text_length;
  int64_t regulation_text_lengtjh;
  int64_t keyword_buffer_length;
  int64_t num_keywords;
  int64_t per_keyword_pointer_size;
  int64_t max_param_line_len;

  int64_t state_offset;
  int64_t two_way_data_offset;
  int64_t incoming_data_offset;
  int64_t auxiliary_data_offset;
  int64_t workspace_offset;
  int64_t workspace_length;
  int64_t current_counts_offset;
  int64_t bndry_flux_counts_offset;
  int64_t net_lklhd_bndry_flux_offset;
  int64_t net_likelihood_offset;
  int64_t vgrng_offset;
  int64_t vgrng2_offset;

  int64_t dg0s_offset;
  int64_t ke_offset;
  int64_t rke_offset;
  int64_t kss_offset;
  int64_t kssr_offset;
  int64_t kss_e_val_offset;
  int64_t kss_u_val_offset;
  int64_t molecule_dg0tfs_offset;
  /*
  int64_t molecule_probabilities_offset;
  int64_t molecule_chemical_potentials_offset;
  */
  int64_t count_to_conc_offset;
  int64_t conc_to_count_offset;
  int64_t activites_offset;
  int64_t reg_constant_offset;
  int64_t reg_exponent_offset;
  int64_t reg_direction_offset;
  int64_t reg_species_offset;
  int64_t coeff_sum_offset;
  int64_t use_rxn_offset;
  int64_t dg0tfs_set_offset;
  
  int64_t reactions_offset;
  int64_t sorted_molecules_offset;
  int64_t sorted_compartments_offset;
  int64_t reactions_matrix_offset;
  int64_t molecules_matrix_offset;
  int64_t reactions_ptrs_offset;
  int64_t molecules_indices_offset;
  int64_t compartment_indices_offset;
  int64_t reactions_coefficients_offset;
  int64_t text_indices_offset;
  int64_t solvent_coefficients_offset;
  int64_t molecules_ptrs_offset;
  int64_t reaction_indices_offset;
  int64_t molecules_coefficients_offset;

  int64_t file_names_offset;
  int64_t solvent_string_offset;
  int64_t reaction_titles_offset;
  int64_t pathway_text_offset;
  int64_t compartment_text_offset;
  int64_t molecules_text_offset;
  int64_t regulation_text_offset;

  int64_t unsorted_molecules_offset;
  int64_t unsorted_compartments_offset;
  int64_t keyword_buff_offset;
  int64_t keyword_lengths_offset;
  int64_t keywords_offset;

  int64_t rxn_has_flux_length;
  int64_t workspace_end;


  *length_changed      = 0;
  print_output         = state->print_output;
  max_regs_per_rxn     = state->max_regs_per_rxn;
  number_reactions     = state->number_reactions;
  number_molecules     = state->number_molecules;
  number_compartments  = state->number_compartments;
  unique_compartments  = state->nunique_compartments;
  unique_molecules     = state->nunique_molecules;
  max_filename_len     = state->max_filename_len;
  num_files            = state->num_files;
  align_len            = state->align_len;
  align_mask           = state->align_mask;

  meta_size    = (int64_t)sizeof(ss);
  if (meta_size != state->meta_size) {
    *length_changed += 1;
  }
  state->meta_size = meta_size;

  vgrng_state_size = (int64_t) sizeof(vss);
  state->vgrng_state_size = vgrng_state_size;

  reaction_data_size    = number_reactions * sizeof(rs);
  molecule_data_size    = unique_molecules  * sizeof(ies);
  compartment_data_size = unique_compartments * sizeof(ces);

  reactions_matrix_size = sizeof(rms);
  molecules_matrix_size = sizeof(mms);

  reactions_ptrs_size   = (number_reactions + 1) * sizeof(int64_t);

  molecules_ptrs_size   = (unique_molecules + 1) * sizeof(int64_t);

  reactions_matrix_field_size = number_molecules * sizeof(int64_t);
  
  molecules_matrix_field_size = number_molecules * sizeof(int64_t);

  solvent_string_size = state->solvent_string_size;

  reaction_titles_length  = state->reaction_titles_length;
  pathway_text_length     = state->pathway_text_length;
  compartment_text_length = state->compartment_text_length;
  molecule_text_length    = state->molecule_text_length;
  regulation_text_length  = state->regulation_text_length;

  keyword_buffer_length   = state->keyword_buffer_length; /*256*/
  num_keywords            = state->max_num_rxn_file_keywords;    /*32 */
  per_keyword_pointer_size = num_keywords * sizeof(int64_t);
  max_param_line_len      = state->max_param_line_len;

  per_molecule_double_size = unique_molecules * sizeof(double);
  state->per_molecule_double_size = per_molecule_double_size;


  

  per_molecule_data_pad = (align_len -
		   (per_molecule_double_size & align_mask)) & align_mask;

  padded_per_molecule_size = per_molecule_double_size + per_molecule_data_pad;


  per_molecule_int_size = (unique_molecules + (unique_molecules & 1))*sizeof(int);
  state->per_molecule_int_size = per_molecule_int_size;
  per_molecule_int_pad = (align_len - (per_molecule_int_size & align_mask)) & align_mask;
  padded_per_molecule_int_size = per_molecule_int_size + per_molecule_int_pad;

  per_reaction_double_size = number_reactions * sizeof(double);
  state->per_reaction_double_size = per_reaction_double_size;

  per_reaction_data_pad = (align_len -
		   (per_reaction_double_size & align_mask)) & align_mask;

  padded_per_reaction_size = per_reaction_double_size + per_reaction_data_pad;


  per_reg_size = number_reactions * max_regs_per_rxn * sizeof(double);
  per_reg_data_pad = (align_len = (per_reg_size & align_mask)) & align_mask;
  padded_per_reg_size = per_reg_size + per_reg_data_pad;
    
  per_compartment_double_size = number_compartments * sizeof(double);


  state_offset = (int64_t)0;
  state->state_offset = state_offset;

  data_pad     = (align_len - (align_mask & meta_size)) & align_mask;

  two_way_data_offset = state_offset + meta_size + data_pad;
  if (two_way_data_offset != state->two_way_data_offset) {
    *length_changed += 2;
  }
  state->two_way_data_offset = two_way_data_offset;

  current_counts_offset = two_way_data_offset;
  state->current_counts_offset = current_counts_offset;
  state->current_counts = (double *)((void *)state + current_counts_offset);

  bndry_flux_counts_offset = current_counts_offset + padded_per_molecule_size;
  state->bndry_flux_counts_offset = bndry_flux_counts_offset;
  state->bndry_flux_counts = 
    (double *)((void *)state + bndry_flux_counts_offset);

  net_lklhd_bndry_flux_offset = bndry_flux_counts_offset + 
    padded_per_molecule_size;
  state->net_lklhd_bndry_flux_offset = net_lklhd_bndry_flux_offset;
  state->net_lklhd_bndry_flux = 
    (double *)((void *)state + net_lklhd_bndry_flux_offset);

  net_likelihood_offset = net_lklhd_bndr_flux_offset + padded_per_molecule_size;
  state->net_likelihood_offset = net_likelihood_offset;
  state->net_likelihood = (douible *)((void*)state + net_likelihood_offset);


  vgrng_state_offset = net_likelihood_offset + padded_per_reaction_size;
  state->vgrng_state_offset = vgrng_state_offset;
  state->vgnrg_state        = (struct vgrng_state_struct *)((void*)state + vgrng_state_offset);
  data_pad = (align_len - (vgrng_state_size & align_mask)) & align_mask;

  vgrng2_state_offset = vgrng_state_offset + vgrng_state_size + data_pad;
  state->vgrng2_state_offset = vgrng2_state_offset;
  state->vgrng2_state        = (struct vgrng_state_struct *)((void*)state + vgrng2_state_offset);

  /*
    Incoming fields.
  */

  incoming_data_offset = vgrng2_state_offset + vgrng_state_size + data_pad;
  if (incoming_data_offset != state->incoming_data_offset) {
    *length_changed += 4;
  }
  state->incoming_data_offset = incoming_data_offset;
  
  dg0s_offset = incoming_data_offset;
  state->dg0s_offset = dg0s_offset;
  state->dg0s  = (double *)((void*)state + dg0s_offset);
  
  ke_offset        = dg0s_offset + padded_per_reaction_size;
  state->ke_offset = ke_offset;
  state->ke        = (double *)((void*)state + ke_offset);
  
  rke_offset        = ke_offset + padded_per_reaction_size;
  state->rke_offset = rke_offset;
  state->rke        = (double *)((void*)state + rke_offset);
  
  kss_offset        = rke_offset + padded_per_reaction_size;
  state->kss_offset = kss_offset;
  state->kss        = (double *)((void*)state + kss_offset);

  kssr_offset        = kss_offset + padded_per_reaction_size;
  state->kssr_offset = kssr_offset;
  state->kssr        = (double *)((void*)state + kssr_offset);

  kss_e_val_offset        = kssr_offset + padded_per_reaction_size;
  state->kss_e_val_offset = kss_e_val_offset;
  state->kss_e_val        = (double *)((void*)state + kss_e_val_offset);

  kss_u_val_offset        = kss_e_val_offset + padded_per_molecule_size;
  state->kss_u_val_offset = kss_u_val_offset;
  state->kss_u_val        = (double *)((void*)state + kss_u_val_offset);

  molecule_dg0tfs_offset        = kss_u_val_offset + padded_per_molecule_size;
  state->molecule_dg0tfs_offset = molecule_dg0tfs_offset;
  state->molecule_dg0tfs        = 
    (double *)((void*)state + molecule_dg0tfs_offset);
  /*
  molecule_probabilities_offset = molecule_dg0tfs_offset + 
    padded_per_molecule_size;
  state->molecule_probabilities_offset = molecule_probabilities_offset;
  state->molecule_probabilities = 
    (double *)((void*)state + molecule_probabilities_offset);

  molecule_chemical_potentials_offset = molecule_probabilities_offset +
    padded_per_molecule_size;
  state->molecule_chemical_potentials_offset = 
    molecule_chemical_potentials_offset;
  state->molecule_chemical_potentials = 
    (double *)((void*)state + molecule_chemical_potentials_offset);

  count_to_conc_offset = molecule_chemical_potentials_offset + 
    padded_per_molecule_size;
  */
  count_to_conc_offset = molecule_dg0tfs_offset + 
    padded_per_molecule_size;

  state->count_to_conc_offset = count_to_conc_offset;
  state->count_to_conc  = (double *)((void*)state + count_to_conc_offset);
  
  conc_to_count_offset        = count_to_conc_offset + padded_per_molecule_size;
  state->conc_to_count_offset = conc_to_count_offset;
  state->conc_to_count        = (double*)((void*)state + conc_to_count_offset);

  activities_offset         = conc_to_count_offset + padded_per_molecule_size;
  state->activities_offset  = activities_offset;
  state->activities         = (double*)((void*)state + activities_offset);

  
  reg_constant_offset        = activities_offset + padded_per_reaction_size;
  state->reg_constant_offset = reg_constant_offset;
  state->reg_constant        = (double*)((void*)state + reg_constant_offset);

  reg_exponent_offset        = reg_constant_offset + padded_per_reg_size;
  state->reg_exponent_offset = reg_exponent_offset;
  state->reg_exponent        = (double*)((void*)state + reg_exponent_offset);

  reg_direction_offset        = reg_exponent_offset + padded_per_reg_size;
  state->reg_direction_offset = reg_direction_offset;
  state->reg_direction        = (double*)((void*)state + reg_direction_offset);

  reg_species_offset          = reg_direction_offset + padded_per_reg_size;
  state->reg_species_offset   = reg_species_offset;
  state->reg_species          = (double*)((void*)state + reg_species_offset);

  coeff_sum_offset            = reg_species_offset + padded_per_reg_size;
  state->coeff_sum_offset     = coeff_sum_offset;
  state->coeff_sum            = (double*)((void*)state + coeff_sum_offset);
  
  
  use_rxn_offset              = coeff_sum_offset + padded_per_molecule_size;
  state->use_rxn_offset       = use_rxn_offset;
  state->use_rxn              = (int*)((void*)state + use_rxn_offset);

  dg0tfs_set_offset           = use_rxn_offset + padded_per_molecule_int_size;
  state->dg0tfs_set_offset    = dg0tfs_set_offset;
  state->dg0tfs               = (int*)((void*)state + dg0tfs_set_offset);
    
  reactions_offset            = dg0tfs_set_offset + padded_per_molecule_int_size;
  state->reactions_offset     = reactions_offset;
  state->reactions            = (struct reaction_struct *) ((void*)state + reactions_offset);
  

  data_pad = (align_len - (reaction_data_size & align_mask)) & align_mask;
  sorted_molecules_offset = reactions_offset + reaction_data_size + data_pad;
  state->sorted_molecules_offset = sorted_molecules_offset;
  state->sorted_molecules = (struct molecule_struct *)((void*)state + sorted_molecules_offset);

  data_pad = (align_len - (molecule_data_size & align_mask)) & align_mask;
  sorted_compartments_offset = sort_molecules_offset + molecule_data_size + data_pad;
  state->sorted_compartments_offset = sorted_compartments_offset;
  state->sorted_compartments = (struct compartment_struct *)((void*)state + sorted_compartments_offset);
  
  data_pad = (align_len - (compartment_data_size & align_mask)) & align_mask;
  reactions_matrix_offset = sorted_compartments_offset + compartment_data_size + data_pad;
  state->reactions_matrix_offset = reactions_matrix_offset;
  state->reactions_matrix = (struct reactions_matrix_struct *)((void*)state+reactions_matrix_offset);

  data_pad = (align_len - (reactions_matrix_size & align_mask)) & align_mask;
  molecules_matrix_offset = reactions_matrix_offset + reactions_matrix_size  + data_pad;
  state->molecules_matrix_offset = molecules_matrix_offset;
  state->molecules_matrix        = (struct molecules_matrix_struct *)((void*state) + molecules_matrix_offset);
  
  /* 
    Fields of the reaction matrix.
  */
  reactions_matrix = state->reactions_matrix;

  data_pad = (align_len - (molecules_matrix_size & align_mask)) & align_mask;
  reactions_ptrs_offset = molecules_matrix_offset + molecules_matrix_size + data_pad;
  state->reactions_ptrs_offset = reactions_ptrs_offset;
  reactions_matrix->reactions_ptrs = (int64_t *)((void*)state + reactions_ptrs_offset);


  data_pad = (align_len - (reactions_ptrs_size & align_mask)) & align_mask;
  molecules_indices_offset = reactions_ptrs_offset + reactions_ptrs_size + data_pad;
  state->molecules_indices_offset     = molecules_indices_offset;
  reactions_matrix->molecules_indices = (int64_t*)((void*)state + molecules_indicies_offset);
  

  data_pad = (align_len - (reactions_matrix_field_size & align_mask)) & align_mask;
  compartment_indices_offset = molecules_indices_offset + reactions_matrix_field_size + data_pad;
  state->compartment_indices_offset = compartment_indices_offset;
  reactions_matrix->compartment_indices = (int64_t *)((void*)state + compartment_indices_offset);

  reactions_coefficients_offset = compartment_indices_offset + reactions_matrix_field_size + data_pad;
  state->reactions_coefficients_offset = reactions_coefficients_offset;
  reactions_matrix->coefficients = (double*)((void*)state + reactions_coefficients_offset);
  

  text_indices_offset = reactions_coefficients_offset + reactions_matrix_field_size + data_pad;
  state->text_indices_offset = text_indiceds_offset;
  reactions_matrix->text = (int64_t *)((void*)state + text_indices_offset);

  solvent_coefficients_offset = 
    reactions_text_offset + reactions_matrix_field_size + data_pad;
  state->solvent_coefficients_offset = solvent_coefficients_offset;
  reactions_matrix->solvent_coefficients = (double*)((void*)state + solvent_coefficients_offset);
  /*
    Molecules_matrix_fields.
  */
  molecules_matrix      = state->molecules_matrix;
  molecules_ptrs_offset = solvent_coefficients_offset + padded_per_reaction_size;
  state->molecules_ptrs_offset = molecules_ptrs_offset;
  molecules_matrix->molecules_ptrs = (int64_t*)((void*)state + molecules_ptrs_offset);
  
  data_pad = (align_len - (molecules_ptrs_size & align_mask)) & align_mask;
  reaction_indices_offset = molecules_ptrs_offset + molecules_ptrs_size + data_pad;
  state->reaction_indices_offset = reaction_indices_offset;
  molecules_matrix->reaction_indices = (int64_t *)((void*)state + reaction_indices_offset);
  
  data_pad = (align_len - (molecules_matrix_field_size & align_mask)) & align_mask;
  molecules_coefficients_offset = reaction_indices_offset + moleucles_matrix_field_size + data_pad;
  state->molecules_coefficients_offset = molecules_coefficients_offset;
  molecules_matrix->coefficients = (double*)((void*)state + molecules_coefficients_offset);

  auxiliary_data_offset = molecules_coefficients_offset + molecules_matrix_field_size + data_pad;
  if (auxiliary_data_offset != state->auxiliary_data_offset) {
    *length_changed += 8;
  }
  state->auxiliary_data_offset = auxiliary_data_offset;

  file_names_offset = auxiliary_data_offset;
  
  state->file_names_offset = file_names_offset;
  state->param_file_name   = (char *)((void*) state+ file_names_offset);
  boltzmann_set_filename_ptrs(state);
						   
  file_names_size = num_files * max_filename_len;
  data_pad        = (align_len - (file_names_size & align_mask)) & align_mask;
  
  solvent_string_offset = file_names_offset + file_names_size + data_pad;

  state->solvent_string_offset = solvent_string_offset;
  state->solvent_string = (char *)((void*)state+solvent_string_offset);

  data_pad        = (align_len - (solvent_string_size & alignmask)) & align_mask;
  
  reaction_titles_offset = solvent_string_offset + solvent_string_size + data_pad;
  state->reaction_titles_offset = reaction_titles_offset;
  state->reaction_titles = (char *)((void*)state+reaction_titles_offset);

  data_pad = (aligh_len - (reaction_titles_length & align_mask)) & align_mask;
  pathway_text_offset  = reaction_titles_offset + reaction_titles_length + data_pad;
  state->pathway_text_offset = pathway_test_offset;
  state->pathway_text = (char*)((void*)state+pathway_text_offset);

  data_pad = (align_len - (pathway_text_length & align_mask)) & align_mask;
  compartment_text_offset = pathway_text_offset + pathway_text_length + data_pad;
  state->compartment_text_offset = compartment_text_offset;
  state->compartment_text = (char*)((void*)state+compartment_text_offset);

  data_pad = (align_len - (compartment_text_length & align_mask)) & align_mask;
  molecule_text_offset = compartment_text_offset + compartment_text_length + data_pad;
  state->molecule_text_offset = molecule_text_offset;
  state->molecule_text = (char*)((void*)state+molecule_text_offset);
  
  data_pad = (align_len - (molecule_text_length & align_mask)) & align_mask;
  regulation_text_offset = molecule_text_offset + molecule_text_length + data_pad;
  state->regulation_text_offset = regulation_text_offset;
  state->reglation_text  = (char*)((void*)state+regulation_text_offset);

  data_pad = (align_len - (regulation_text_length & align_mask)) & align_mask;

  workspace_offset = regulation_text_offset + regulation_text_length + data_pad;
  if (workspace_offset != state->workspace_offset) {
    *length_changed += 16;
  }

  state->workspace_offset = workspace_offset;

  /*
    these fields are only used in setup never, running, so wew
    don't need to preserver those.
  */
  /*
  unsorted_molecules_offset = workspace_offset;
  state->unsorted_molecules_offset = unsorted_molecules_offset;
  state->unsorted_molecules = (struct molecule_struct *)((void*)state+ unsorted_molecules_offset);
  
  data_pad = (align_len - (molecule_data_size & align_mask)) & align_mask;
  unsorted_compartments_offset = unsorted_molecules_offset + molecule_data_size + data_pad;
  state->unsorted_compartments_offset = unsorted_compartments_offset;
  state->unsorted_compartments = (struct compartment_struct *)((void*)state+unsorted_compartments_offset);

  data_pad = (align_len - (compartment_data_size & align_mask)) & align_mask;
  compartment_ptrs_offset = unsorted_compartments_offset + compartment_data_size + data_pad;
  state->compartment_ptrs_offset = compartment_ptrs_offset;
  state->compartment_ptrs =  (int64_t *)((void*)state + compartment_ptrs_offset);
  data_pad = (align_len - (per_compartment_double_size & align_mask)) & align_mask;
  keyword_buff_offset = compartment_ptrs_offset + per_compartment_double_size + data_pad;
  state->keyword_buff_offset = rxn_file_keyword_buff_offset;
  state->keyword_buff = (char*)((void*)state+keyword_buff_offset);


  data_pad = (align_len - (keyword_buffer_length & align_mask)) & align_mask;
  keyword_lengths_offset = keyword_buff_offset + keyword_buffer_length + data_pad;
  state->keyword_lengths_offset = keyword_lengths_offset;
  state->keyword_lengths = (int64_t*)((void*)state+keyword_lengths_offset);


  data_pad = (align_len - (per_keyword_pointer_size & align_mask)) & align_mask;
  
  keywords_offset  = keyword_lengths_offset + per_keyword_pointer_size + data_pad;
  state->keywords_offset = keywords_offset;
  state->keywords = (char **)((void*)state + keywords_offset);

  param_buffer_offset = keywords_offset + per_keyword_pointer_size + data_pad;
  state->param_buffer_offset = param_buffer_offset;
  state->param_buffer = (char **)((void*)state+param_buffer_offset);
  
  data_pad = (align_len - ((max_param_line_length*2) & max_align)) & max_align;
  
  raw_molecules_text_offset = param_buffer_offset + (max_param_line_length*2) +
    data_pad;
  state->raw_molecules_text_offset = raw_molecules_text_offset;
  data_pad = (align_len - (molecules_text_length & align_mask)) & align_mask;
  
  transpose_workspace_offset = raw_molecules_text_offset + molecules_text_length + data_pad;
  state->transpose_workspace_offset = transpose_workspace_offset;
  state->transpose_workspace_offset = (int64_t*)((void*)state + transpose_workspace_offset);
  
  data_pad = (align_len - (molecules_ptrs_size &align_mask)) & align_mask;
  future_counts_offset = transpose_workspace_offset + molecules_ptrs_size + darta_pad;
  */
  future_counts_offset = workspace_offset;
  state->future_counts_offset = future_counts_offset;
  state->future_counts = (double *)((void*)state + future_counts_offset);

  free_energy_offset = future_counts_offset + padded_per_molecule_size;
  state->free_energy_offset = free_energy_offset;
  state->free_energy = (double*)((void*)state + free_energy_offset);

  forward_rxn_likelihood_offset = free_energy_offset + padded_per_reaction_size;
  state->forward_rxn_likelihood_offset = forward_rxn_likelihood_offset;
  state->forward_rxn_likelihood = (double*)((void*)state+ forwrard_rxn_likelihood_offset);
  
  reverse_rxn_likelihood_offset = forward_rxn_likelihood_offset + padded_per_reaction_size;
  state->reverse_rxn_likelihood_offset = reverse_rxn_likelihood_offset;
  state->reverse_rxn_likelihood = (double*)((void*)state+ reverse_rxn_likelihood_offset);

  forward_rxn_log_likelihood_ratio_offset = reverse_rxn_likelihood_offset + padded_per_reaction_size;
  state->forward_rxn_log_likelihood_ratio_offset = 
    forward_rxn_log_likelihood_ratio_offset;
  state->forward_rxn_log_likelihood_ratio = 
    (double*)((void*)state + forward_rxn_log_likelihood_ratio_offset);

  reverse_rxn_log_likelihood_ratio_offset = 
    forward_rxn_log_likelihood_ratio_offset + padded_per_reaction_size;
  state->reverse_rxn_log_likelihood_ratio_offset = 
    reverse_rxn_log_likelihood_ratio_offset;
  state->reverse_rxn_log_likelihood_ratio = 
    (double*)((void*)state + reverse_rxn_log_likelihood_ratio_offset);

  rxn_likelihood_ps_offset = 
    reverse_rxn_log_likelihood_ratio_offset + padded_per_reaction_size;
  state->rxn_likelihood_ps_offset = rxn_likelihood_ps_offset;
  state->rxn_likelihood_ps = 
    (double*)((void*)state + rxn_likelihood_ps_offset);

  data_pad = (align_len - (reactions_ptrs_size & align_mask)) & align_mask;
  reactant_term_offset = rxn_liklihood_ps_offset + 
                         reactions_ptrs_size + data_pad;
  state->reactant_term_offset = reactant_term_offset;
  state->reactant_term = (double*)((void*)state+ reactant_term_offset);

  product_term_offset  = reactant_term_offset + padded_per_reaction_size;
  state->product_term_offset = product_term_offset;
  state->product_tem   = (double*)((void*)state + product_term_offset);
  
  rxn_q_offset = product_term_offset + padded_per_reaction_size;
  state->rxn_q_offset = rxn_q_offset;
  state->rxn_q = (double*)((void*)state + rxn_q_offset);

  recip_rxn_q_offset = rxn_q_offset + padded_per_reaction_size;
  state->recip_rxn_q_offset = recip_rxn_q_offset;
  state->recip_rxn_q = (double*)((void*)state + recip_rxn_q_offset);

  log_kf_rel_offset = recip_rxn_q_offset + padded_per_reaction_size;
  state->log_kf_rel_offset = log_kf_rel_offset;
  state->log_kf_rel = (double*)((void*)state + log_kf_rel_offset);

  log_kr_rel_offset = log_kf_rel_offset + padded_per_reaction_size;
  state->log_kr_rel_offset = log_kr_rel_offset;
  state->log_kr_rel = (double*)((void*)state + log_kr_rel_offset);

  ode_counts_offset = log_kr_rel_oofset + padded_per_reaction_size;
  state->ode_counts_offset = ode_counts_offset;
  state->ode_counts = (double*)((void*)state + ode_counts_offset);

  ode_concs_offset  = ode_counts_offset + padded_per_molecule_size;
  state->ode_concs_offset = ode_concs_offset;
  state->ode_concs = (double*)((void*)state + ode_concs_offset);

  ode_forward_lklhds_offset = ode_concs_offset + padded_per_molecule_size;
  state->ode_forward_lklhds_offset = ode_forward_lklhds_offset;
  state->ode_forward_lklhds = (double*)((void*)state + ode_forward_lklhds_offset);

  ode_reverse_lklhds_offset = ode_forward_lklhds_offset + padded_per_reaction_size;
  state->ode_reverse_lklhds_offset = ode_reverse_lklhds_offset;
  state->ode_reverse_lklhds = (double*)((void*)state + ode_reverse_lklhds_offset);

  /*
    Here we assuem that the two integer vectors each of length
    unique_molecules can fit in the padded_per_molecule_size number of bytes.
  */
  base_reactant_indicator_offset = ode_reverse_lklhds_offset + padded_per_reaction_size;
  state->base_reactant_indicator_offset = base_reactant_indicator;
  state->base_reactant_indicator = (int*)((void*)state + base_reactant_indicator_offset);
  base_reactants_offset = base_reactant_indicator_offset + (unique_molecules * sizeof(int));
  state->base_reactants_offset = base_reactants_offset;
  state->base_reactants = (int*)((void*)state + base_reactants_offset);

  rxn_has_flux_offset = base_reactant_indicator_offset + padded_per_molecule_size;
  state->rxn_has_flux_offset = rxn_has_flux_offset;
  state->rxn_has_flux = (int*)((void*)state + rxn_has_flux_offset);

  rxn_has_flux_length = (int64_t)(number_reactions*sizeof(int));
  data_pad = (align_len - (rxn_has_flux_length &align_mask)) & align_mask;
  workspace_end = rxn_has_flux_offset + rxn_has_flux_length + data_pad;
  workspace_length = workspace_end - workspace_offset;
  if (workspace_length > state->workspace_length) {
    *length_changed += 32;
  }
  /* 
    We don't reserve the printing workspace, as that is unlikely to happen
    in the restart/parallel case. We let boltzmann_run create the printing
    space for 
    no_op_likelihood (length rxn_view_hist_length),
    rxn_view_likelihoods (length rxn_view_hist_length* number_reactions),
    rev_rxn_view_likelihoods (length rxn_view_hist_length* number_reactions),
    rxn_fire (length (number_reactions+1)*2),
    rxn_mat_row (length unique_molecules);
    These fields are now set in alloc8 called by boltzmann_run.
  */
}
