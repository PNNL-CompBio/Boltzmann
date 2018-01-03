/* flatten_state.c
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

#include "boltzmann_set_filename_ptrs.h"
#include "flatten_state.h"

int flatten_state(struct state_struct *boot_state, 
		  struct state_struct **new_state_p) {
  /*
    If new_state_p is null,
      allocate space for a contiguous piece of memory to store state,
      and copy the pieces of the boot_state into the new_state and return
      a pointer to the new_state. Also set the self referential pointers 
      in the new_state appropriately.
    else 
      just set the self referential pointers in the new_state appropriately.

      
    Called by run_init and boltzmann_run

  */
  struct state_struct *new_state;
  struct state_struct ss;
  struct vgrng_state_struct vss;
  struct molecule_struct ies;
  struct compartment_struct ces;
  struct rxn_struct rs;
  struct rxn_matrix_struct rms;
  struct rxn_matrix_struct *reactions_matrix;
  struct molecules_matrix_struct mms;
  struct molecules_matrix_struct *molecules_matrix;
  int64_t *new_state_l;
  int64_t *workspace_base_l;
  int64_t *workspace_base;
  int64_t ask_for;
  int64_t one_l;
  int64_t align_len;
  int64_t align_mask;
  int64_t word_len;
  int64_t log2_word_len;
  int64_t state_offset;
  int64_t state_length;
  int64_t meta_size;
  int64_t meta_pad;
  int64_t two_way_data_offset;
  int64_t two_way_data_length;
  int64_t incoming_data_offset;
  int64_t incoming_data_length;
  int64_t auxiliary_data_offset;
  int64_t auxiliary_data_length;
  int64_t auxiliary_end;
  int64_t workspace_offset;
  int64_t workspace_length;
  int64_t print_output;
  int64_t unique_molecules;
  int64_t number_reactions;
  int64_t number_molecules;
  int64_t number_compartments;
  int64_t unique_compartments;
  int64_t max_filename_len;
  int64_t num_files;
  int64_t solvent_pos;
  int64_t dg_forward_offset;
  int64_t dg_forward_size;
  int64_t entropy_offset;
  int64_t entropy_size;
  int64_t current_counts_offset;
  int64_t current_counts_size;
  int64_t bndry_flux_counts_offset;
  int64_t bndry_flux_counts_size;
  int64_t net_lklhd_bndry_flux_offset;
  int64_t net_lklhd_bndry_flux_size;
  int64_t vgrng_offset;
  int64_t vgrng_size;
  int64_t vgrng_pad;
  int64_t vgrng2_offset;
  int64_t net_likelihood_offset;
  int64_t net_likelihood_size;
  int64_t dg0s_offset;
  int64_t dg0s_size;
  int64_t ke_offset;
  int64_t ke_size;
  int64_t rke_offset;
  int64_t rke_size;
  int64_t kss_size;
  int64_t kss_offset;
  int64_t kss_e_val_size;
  int64_t kss_e_val_offset;
  int64_t kss_u_val_size;
  int64_t kss_u_val_offset;
  int64_t molecule_dg0tfs_offset;
  int64_t molecule_dg0tfs_size;
  int64_t molecule_probabilities_offset;
  int64_t molecule_probabilities_size;
  int64_t molecule_chemical_potentials_offset;
  int64_t molecule_chemical_potentials_size;
  int64_t count_to_conc_offset;
  int64_t count_to_conc_size;
  int64_t conc_to_count_offset;
  int64_t conc_to_count_size;
  
  int64_t activities_offset;
  int64_t activities_size;

  int64_t reg_constant_offset;
  int64_t reg_constant_size;
  int64_t reg_exponent_offset;
  int64_t reg_exponent_size;
  int64_t reg_species_offset;
  int64_t reg_species_size;
  int64_t reg_drctn_offset;
  int64_t reg_drctn_size;

  int64_t reactions_offset;
  int64_t reactions_size;
  int64_t reactions_pad;
  int64_t reactions_matrix_offset;
  int64_t reactions_matrix_size;
  int64_t reactions_matrix_pad;
  int64_t rxn_ptrs_offset;
  int64_t rxn_ptrs_size;
  int64_t molecules_offset;
  int64_t molecules_size;
  int64_t pathway_offset;
  int64_t pathway_size;
  int64_t compartment_offset;
  int64_t compartment_size;
  int64_t molecules_indices_offset;
  int64_t molecules_indices_size;
  int64_t compartment_indices_offset;
  int64_t compartment_indices_size;
  int64_t coefficients_offset;
  int64_t coefficients_size;
  int64_t solvent_coefficients_offset;
  int64_t solvent_coefficients_size;
  int64_t text_offset;
  int64_t text_size;
  int64_t sorted_molecules_offset;
  int64_t sorted_molecules_size;
  int64_t sorted_molecules_pad;
  int64_t sorted_cmpts_offset;
  int64_t sorted_cmpts_size;
  int64_t sorted_cmpts_pad;
  int64_t file_names_offset;
  int64_t file_names_size;
  int64_t file_names_pad;
  int64_t rxn_title_offset;
  int64_t rxn_title_size;
  int64_t future_counts_offset;
  int64_t future_counts_size;
  int64_t free_energy_offset;
  int64_t free_energy_size;
  int64_t forward_rxn_offset;
  int64_t forward_rxn_size;
  int64_t reverse_rxn_offset;
  int64_t reverse_rxn_size;
  int64_t forward_rxn_log_offset;
  int64_t forward_rxn_log_size;
  int64_t reverse_rxn_log_offset;
  int64_t reverse_rxn_log_size;
  int64_t rxn_likelihood_ps_offset;
  int64_t rxn_likelihood_ps_size;
  int64_t workspace_end;
  int64_t rxn_view_hist_length;
  int64_t no_op_likelihood_offset;
  int64_t no_op_likelihood_size;
  int64_t rxn_view_offset;
  int64_t rxn_view_size;
  int64_t rev_rxn_view_offset;
  int64_t rev_rxn_view_size;
  int64_t rxn_fire_offset;
  int64_t rxn_fire_size;
  int64_t rxn_mat_offset;
  int64_t rxn_mat_size;
  int64_t move_size;
  int64_t sorted_molecules_offset_in_bytes;
  int64_t sorted_compartments_offset_in_bytes;
  int64_t molecules_text_offset_in_bytes;
  int64_t compartments_text_offset_in_bytes;
  int64_t max_regs_per_rxn;
  /*
    Adding stuff for molecules matrix.
  sorted_molecules_offset = text_offset + text_size;
  */
  int64_t molecules_matrix_offset;
  int64_t molecules_matrix_size;
  int64_t molecules_matrix_pad;
  int64_t molecules_ptrs_offset;
  int64_t molecules_ptrs_size;
  int64_t rxn_indices_offset;
  int64_t rxn_indices_size;
  int64_t mlcls_coeffs_offset;
  int64_t mlcls_coeffs_size;
  int64_t transpose_workspace_offset;
  int64_t transpose_workspace_size;
  int success;
  int load_from_boot;
  one_l   = (int64_t)1;
  success     = 1;
  word_len    = 8;
  log2_word_len = 3;
  align_len  = word_len + word_len;
  align_mask = align_len - 1;
  state_offset = (int64_t)0;
  meta_size    = (int64_t)sizeof(ss);
  meta_pad     = (align_len - (align_mask & meta_size)) & align_mask;
  two_way_data_offset = state_offset + ((meta_size + meta_pad) >> log2_word_len);
  print_output         = boot_state->print_output;
  max_regs_per_rxn     = boot_state->max_regs_per_rxn;
  number_reactions     = boot_state->number_reactions;
  number_molecules     = boot_state->number_molecules;
  number_compartments  = boot_state->number_compartments;
  unique_compartments  = boot_state->nunique_compartments;
  unique_molecules     = boot_state->nunique_molecules;
  max_filename_len     = boot_state->max_filename_len;
  num_files            = boot_state->num_files;

  dg_forward_offset  = two_way_data_offset;
  dg_forward_size    = 1;
  entropy_offset     = dg_forward_offset + dg_forward_size;
  entropy_size       = 1;
  current_counts_offset    = entropy_offset + entropy_size;
  current_counts_size      = unique_molecules + (unique_molecules & 1);
  bndry_flux_counts_offset = current_counts_offset + 
    current_counts_size;
  bndry_flux_counts_size   = current_counts_size;
  net_lklhd_bndry_flux_offset = bndry_flux_counts_offset + bndry_flux_counts_size;
  net_lklhd_bndry_flux_size   = current_counts_size;
  number_reactions            = boot_state->number_reactions;
  net_likelihood_offset       = net_lklhd_bndry_flux_offset + net_lklhd_bndry_flux_size;
  net_likelihood_size         = number_reactions + (number_reactions & 1);
  vgrng_offset                = net_likelihood_offset + net_likelihood_size;
  vgrng_size                  = (int64_t)sizeof(vss);
  vgrng_pad = ((align_len - (vgrng_size & align_mask)) & align_mask);
  vgrng_size = (vgrng_size + vgrng_pad) >> log2_word_len;
  vgrng2_offset = vgrng_offset + vgrng_size;
  two_way_data_length  = vgrng2_offset + vgrng_size - two_way_data_offset;
  two_way_data_length += two_way_data_length & 1;

  incoming_data_offset = two_way_data_offset + two_way_data_length;
  solvent_pos          = boot_state->solvent_pos;
  dg0s_offset          = incoming_data_offset;
  dg0s_size            = number_reactions + (number_reactions & 1);
  ke_offset            = dg0s_offset + dg0s_size;
  ke_size              = dg0s_size;
  rke_offset           = ke_offset + ke_size;
  rke_size             = ke_size;
  kss_offset           = rke_offset + rke_size;
  kss_size             = ke_size + ke_size;
  kss_e_val_size       = unique_molecules + (unique_molecules & 1);
  kss_e_val_offset     = kss_offset + kss_size;
  kss_u_val_size       = kss_e_val_size;
  kss_u_val_offset     = kss_e_val_offset + kss_e_val_size;
  
  molecule_dg0tfs_offset = kss_u_val_offset + kss_u_val_size;
  molecule_dg0tfs_size = unique_molecules;
  molecule_probabilities_offset =  molecule_dg0tfs_offset + molecule_dg0tfs_size;
  molecule_probabilities_size   = unique_molecules;
  molecule_chemical_potentials_offset = molecule_probabilities_offset +
    molecule_probabilities_size;
  molecule_chemical_potentials_size   = unique_molecules;
  
  count_to_conc_offset = molecule_chemical_potentials_offset +
    molecule_chemical_potentials_size;
  count_to_conc_size  = unique_molecules;

  conc_to_count_offset = count_to_conc_offset + count_to_conc_size;
  conc_to_count_size   = unique_molecules;

  activities_offset    = conc_to_count_offset + conc_to_count_size;
  activities_size      = number_reactions;

  reg_constant_size    = number_reactions * max_regs_per_rxn;
  reg_constant_offset  = activities_offset + activities_size;

  reg_exponent_size    = reg_constant_size;
  reg_exponent_offset  = reg_constant_offset + reg_constant_size;

  reg_drctn_size       = reg_constant_size;
  reg_drctn_offset     = reg_exponent_offset + reg_exponent_size;

  reg_species_size     = reg_constant_size;
  reg_species_offset   = reg_drctn_offset + reg_drctn_size;

  reactions_offset     = reg_species_offset + reg_species_size;
  reactions_size       = (int64_t)sizeof(rs) * number_reactions; 
  reactions_pad        = (align_len - (reactions_size & align_mask)) & align_mask;
  reactions_size       = (reactions_size + reactions_pad) >> log2_word_len;
  reactions_matrix_offset = reactions_offset + reactions_size;
  reactions_matrix_size   = sizeof(rms);
  reactions_matrix_pad    = (align_len - (reactions_matrix_size & align_mask)) & align_mask;
  reactions_matrix_size   = (reactions_matrix_size + reactions_matrix_pad) >> log2_word_len;
  rxn_ptrs_offset        = reactions_matrix_offset + reactions_matrix_size;
  rxn_ptrs_size          = number_reactions + 1;
  rxn_ptrs_size          = rxn_ptrs_size + (rxn_ptrs_size&1);
  molecules_indices_offset  = rxn_ptrs_offset + rxn_ptrs_size;
  molecules_indices_size    = number_molecules + (number_molecules&1);
  compartment_indices_offset = molecules_indices_offset + 
                               molecules_indices_size;
  compartment_indices_size = molecules_indices_size;
  coefficients_offset      = compartment_indices_offset + compartment_indices_size;
  coefficients_size      = molecules_indices_size;
  solvent_coefficients_offset = coefficients_offset + coefficients_size;
  solvent_coefficients_size   = number_reactions + (number_reactions & 1);
  text_offset            = solvent_coefficients_offset + solvent_coefficients_size;
  text_size              = molecules_indices_size;

  molecules_matrix_offset = text_offset + text_size;
  molecules_matrix_size   = sizeof(mms);
  molecules_matrix_pad    = (align_len - (molecules_matrix_size & align_mask)) & align_mask;
  molecules_matrix_size   = (molecules_matrix_size + molecules_matrix_pad) >> log2_word_len;
  molecules_ptrs_offset   = molecules_matrix_offset + molecules_matrix_size;
  molecules_ptrs_size     = unique_molecules + 1;
  molecules_ptrs_size     = molecules_ptrs_size + (molecules_ptrs_size & 1);
  rxn_indices_offset      = molecules_ptrs_offset + molecules_ptrs_size;
  rxn_indices_size        = number_molecules + (number_molecules&1);
  mlcls_coeffs_offset     = rxn_indices_offset + rxn_indices_size;
  mlcls_coeffs_size       = rxn_indices_size;

  sorted_molecules_offset = mlcls_coeffs_offset + mlcls_coeffs_size;
  sorted_molecules_offset_in_bytes = sorted_molecules_offset << log2_word_len;
  sorted_molecules_size   = sizeof(ies) * number_molecules;
  sorted_molecules_pad    = (align_len - (sorted_molecules_size & align_mask)) & align_mask;
  sorted_molecules_size   = (sorted_molecules_size + sorted_molecules_pad) >> log2_word_len;
  sorted_cmpts_offset     = sorted_molecules_offset + sorted_molecules_size;
  sorted_compartments_offset_in_bytes = sorted_cmpts_offset << log2_word_len;
  sorted_cmpts_size       = sizeof(ces) * unique_compartments;
  sorted_cmpts_pad        = (align_len - (sorted_cmpts_size & align_mask)) & align_mask;
  sorted_cmpts_size       = (sorted_cmpts_size + sorted_cmpts_pad) >> log2_word_len;
  auxiliary_data_offset   = sorted_cmpts_offset + sorted_cmpts_size;
  incoming_data_length    = auxiliary_data_offset - incoming_data_offset;  
  
  file_names_offset       = auxiliary_data_offset;
  file_names_size         = num_files * max_filename_len;
  file_names_pad          = (align_len - (file_names_size & align_mask)) & align_mask;
  file_names_size         = (file_names_size + file_names_pad) >> log2_word_len;
  rxn_title_offset        = file_names_offset + file_names_size;
  rxn_title_size          = boot_state->rxn_title_text_length >> log2_word_len;
  pathway_offset          = rxn_title_offset + rxn_title_size;
  pathway_size            = boot_state->pathway_text_length >> log2_word_len;
  compartment_offset      = pathway_offset + pathway_size;
  compartments_text_offset_in_bytes = compartment_offset << log2_word_len;
  compartment_size        = boot_state->compartment_text_length >> log2_word_len;
  molecules_offset        = compartment_offset + compartment_size;
  molecules_text_offset_in_bytes = molecules_offset << log2_word_len;
  molecules_size          = boot_state->molecule_text_length >> log2_word_len;

  auxiliary_end           = molecules_offset + molecules_size;
  auxiliary_data_length   = auxiliary_end - auxiliary_data_offset;
  /*
    Overlap workspace with auxiliary_data as the auxiliary_data is
    not needed for running the simulations.
  */
  /*
  if (print_output) {
    workspace_offset       = auxiliary_end;
  } else {
    workspace_offset       = auxiliary_data_offset;
  }
  */
  /*
    Here we are just computing the size of the workspace.
  */
  workspace_base          = boot_state->workspace_base;
  workspace_offset        = (int64_t)0;
  future_counts_offset     = workspace_offset;
  future_counts_size       = unique_molecules + (unique_molecules&1);
  free_energy_offset      = future_counts_offset + future_counts_size;
  free_energy_size        = number_reactions + (number_reactions &1);
  forward_rxn_offset      = free_energy_offset + free_energy_size;
  forward_rxn_size        = free_energy_size;
  reverse_rxn_offset      = forward_rxn_offset + forward_rxn_size;
  reverse_rxn_size        = free_energy_size;
  forward_rxn_log_offset  = reverse_rxn_offset + reverse_rxn_size;
  forward_rxn_log_size    = free_energy_size;
  reverse_rxn_log_offset  = forward_rxn_log_offset + forward_rxn_log_size;
  reverse_rxn_log_size    = free_energy_size;
  rxn_likelihood_ps_offset = reverse_rxn_log_offset + reverse_rxn_log_size;
  rxn_likelihood_ps_size  = (free_energy_size + 1) << 1;
  transpose_workspace_offset = rxn_likelihood_ps_offset + rxn_likelihood_ps_size;
  transpose_workspace_size = unique_molecules+1;
  workspace_end =  transpose_workspace_offset + transpose_workspace_size;
  if (print_output) {
    rxn_view_hist_length     = boot_state->rxn_view_hist_length;
    no_op_likelihood_offset = workspace_end;
    no_op_likelihood_size   = rxn_view_hist_length;
    rxn_view_offset         = no_op_likelihood_offset + no_op_likelihood_size;
    rxn_view_size           = rxn_view_hist_length * number_reactions;
    rxn_view_size           += (rxn_view_size & 1);
    rev_rxn_view_offset     = rxn_view_offset + rxn_view_size;
    rev_rxn_view_size       = rxn_view_size;
    rxn_fire_offset         = rev_rxn_view_offset + rev_rxn_view_size;
    rxn_fire_size           = number_reactions+1;
    rxn_fire_size           += (rxn_fire_size&1);
    rxn_mat_offset          = rxn_fire_offset + rxn_fire_size;
    rxn_mat_size            = (unique_molecules + (unique_molecules & 1)) >> 1;
    workspace_end           = rxn_mat_offset + rxn_mat_size;
  }
  workspace_length = workspace_end - workspace_offset;
  /*
  state_length = workspace_end;
  if (auxiliary_end > workspace_end) {
    state_length = auxiliary_end;
  }
  */
  state_length = auxiliary_end;
  if (*new_state_p == NULL) {
    ask_for = state_length * sizeof(int64_t);
    new_state_l = (int64_t*)calloc(one_l,ask_for);
    if (new_state_l) {
      new_state = (struct state_struct *)new_state_l;
      *new_state_p = new_state;
      /*
	Copy boot_state.
      */
      memcpy(new_state,boot_state,meta_size);
      new_state->state_length          = state_length<<log2_word_len;
      new_state->two_way_data_offset   = two_way_data_offset;
      new_state->two_way_data_length   = two_way_data_length;
      new_state->incoming_data_offset  = incoming_data_offset;
      new_state->incoming_data_length  = incoming_data_length;
      new_state->auxiliary_data_offset = auxiliary_data_offset;
      new_state->auxiliary_data_length = auxiliary_data_length;
      new_state->workspace_base        = workspace_base;
      new_state->workspace_offset      = workspace_offset;
      new_state->workspace_length      = workspace_length;
      new_state->sorted_molecules_offset_in_bytes = sorted_molecules_offset_in_bytes;
      new_state->sorted_compartments_offset_in_bytes = sorted_compartments_offset_in_bytes;
      new_state->molecules_text_offset_in_bytes = molecules_text_offset_in_bytes;
      new_state->compartments_text_offset_in_bytes = compartments_text_offset_in_bytes;
      load_from_boot = 1;
    } else {
      fprintf(stderr,"flatten_state: Error unable to allocate %lld bytes for new_state\n",ask_for);
      fflush(stderr);
      success = 0;
    }
  } else {
    new_state = boot_state;
    new_state_l = (int64_t *)new_state;
    load_from_boot = 0;
  }
  if (success  && new_state) {
    /*
      Allocate space for a workspace if on input it was null.
    */
    if (workspace_base == NULL) {
      ask_for = workspace_length * sizeof(int64_t);
      workspace_base_l = (int64_t*)calloc(one_l,ask_for);
      if (workspace_base_l) {
	new_state->workspace_base = workspace_base_l;
      } else {
	new_state->workspace_base = 0;
	success = 0;
	fprintf(stderr,"flatten_state: could not allocate %lld bytes "
		"for workspace\n",ask_for);
	fflush(stderr);
      }
    } else {
      workspace_base_l = workspace_base;
    }
  }
  if (success  && new_state) {
    /*
      Set pointers.
    */
    new_state->dg_forward_p = (double*)&new_state_l[dg_forward_offset];
    new_state->entropy_p = (double*)&new_state_l[entropy_offset];
    new_state->current_counts = (double*)&new_state_l[current_counts_offset];
    new_state->bndry_flux_counts = (double *)&new_state_l[bndry_flux_counts_offset];
    new_state->net_lklhd_bndry_flux = (double*)&new_state_l[net_lklhd_bndry_flux_offset];
    new_state->net_likelihood = (double*)&new_state_l[net_likelihood_offset]; 
    new_state->vgrng_state = (struct vgrng_state_struct *)&new_state_l[vgrng_offset];
    new_state->vgrng2_state = (struct vgrng_state_struct *)&new_state_l[vgrng2_offset];
    new_state->dg0s = (double*)&new_state_l[dg0s_offset];
    new_state->ke   = (double*)&new_state_l[ke_offset];
    new_state->rke  = (double*)&new_state_l[rke_offset];
    new_state->kss  = (double*)&new_state_l[kss_offset];
    new_state->kssr = (double*)&new_state_l[kss_offset+ke_size];
    new_state->kss_e_val  = (double*)&new_state_l[kss_e_val_offset];
    new_state->kss_u_val  = (double*)&new_state_l[kss_u_val_offset];
    new_state->molecule_dg0tfs   = (double*)&new_state_l[molecule_dg0tfs_offset];
    new_state->molecule_probabilities   = (double*)&new_state_l[molecule_probabilities_offset];
    new_state->molecule_chemical_potentials   = (double*)&new_state_l[molecule_chemical_potentials_offset];
    new_state->count_to_conc = (double*)&new_state_l[count_to_conc_offset];
    new_state->conc_to_count = (double*)&new_state_l[conc_to_count_offset];

    new_state->activities    = (double*)&new_state_l[activities_offset];
    new_state->reg_constant  = (double*)&new_state_l[reg_constant_offset];
    new_state->reg_exponent  = (double*)&new_state_l[reg_exponent_offset];
    new_state->reg_drctn     = (double*)&new_state_l[reg_drctn_offset];
    new_state->reg_species   = (int64_t*)&new_state_l[reg_species_offset];


    new_state->reactions   = (struct rxn_struct *)&new_state_l[reactions_offset];
    new_state->reactions_matrix = (struct rxn_matrix_struct*)&new_state_l[reactions_matrix_offset];
    reactions_matrix = new_state->reactions_matrix;
    reactions_matrix->rxn_ptrs = (int64_t *)&new_state_l[rxn_ptrs_offset];
    reactions_matrix->molecules_indices = (int64_t *)&new_state_l[molecules_indices_offset];
    reactions_matrix->compartment_indices = (int64_t *)&new_state_l[compartment_indices_offset];
    reactions_matrix->coefficients = (int64_t *)&new_state_l[coefficients_offset];
    reactions_matrix->solvent_coefficients = (int64_t *)&new_state_l[solvent_coefficients_offset];
    reactions_matrix->text = (int64_t*)&new_state_l[text_offset];
    new_state->molecules_matrix = (struct molecules_matrix_struct*)&new_state_l[molecules_matrix_offset];
    molecules_matrix = new_state->molecules_matrix;
    molecules_matrix->molecules_ptrs = (int64_t*)&new_state_l[molecules_ptrs_offset];
    molecules_matrix->rxn_indices    = (int64_t*)&new_state_l[rxn_indices_offset];
    molecules_matrix->coefficients   = (int64_t*)&new_state_l[mlcls_coeffs_offset];
    new_state->sorted_molecules = (struct molecule_struct*)&new_state_l[sorted_molecules_offset];
    new_state->sorted_cmpts = (struct compartment_struct*)&new_state_l[sorted_cmpts_offset];
    new_state->params_file  = (char*)&new_state_l[file_names_offset];
    new_state->max_filename_len = max_filename_len;
    boltzmann_set_filename_ptrs(new_state);
    new_state->rxn_title_text = (char *)&new_state_l[rxn_title_offset];
    new_state->pathway_text   = (char *)&new_state_l[pathway_offset];
    new_state->compartment_text   = (char *)&new_state_l[compartment_offset];
    new_state->molecules_text   = (char *)&new_state_l[molecules_offset];
    new_state->future_counts = (double *)&workspace_base_l[future_counts_offset];
    new_state->free_energy = (double *)&workspace_base_l[free_energy_offset];
    new_state->forward_rxn_likelihood = (double *)&workspace_base_l[forward_rxn_offset];
    new_state->reverse_rxn_likelihood = (double *)&workspace_base_l[reverse_rxn_offset];
    new_state->forward_rxn_log_likelihood_ratio = (double *)&workspace_base_l[forward_rxn_log_offset];
    new_state->reverse_rxn_log_likelihood_ratio = (double *)&workspace_base_l[reverse_rxn_log_offset];
    new_state->rxn_likelihood_ps = (double*)&workspace_base_l[rxn_likelihood_ps_offset];
    if (print_output) {
      new_state->no_op_likelihood = (double *)&workspace_base_l[no_op_likelihood_offset];
      new_state->rxn_view_likelihoods = (double *)&workspace_base_l[rxn_view_offset];
      new_state->rev_rxn_view_likelihoods = (double *)&workspace_base_l[rev_rxn_view_offset];
      new_state->rxn_fire = (int *)&workspace_base_l[rxn_fire_offset];
      new_state->rxn_mat_row  = (int *)&workspace_base_l[rxn_mat_offset];
    }
  }
  if (success) {
    if (load_from_boot) {
      move_size = unique_molecules * sizeof(double);
      memcpy(new_state->current_counts,boot_state->current_counts,move_size);
      move_size = (int64_t)sizeof(vss);
      memcpy(new_state->vgrng_state,boot_state->vgrng_state,move_size);
      memcpy(new_state->vgrng2_state,boot_state->vgrng2_state,move_size);
      move_size = number_reactions*sizeof(double);
      memcpy(new_state->dg0s,boot_state->dg0s,move_size);
      memcpy(new_state->ke,boot_state->ke,move_size);
      memcpy(new_state->rke,boot_state->rke,move_size);
      move_size = kss_size * sizeof(double);
      memcpy(new_state->kss,boot_state->kss,move_size);
      move_size = unique_molecules * sizeof(double);
      memcpy(new_state->kss_e_val,boot_state->kss_e_val,move_size);
      memcpy(new_state->kss_u_val,boot_state->kss_u_val,move_size);
      memcpy(new_state->molecule_dg0tfs,boot_state->molecule_dg0tfs,move_size);
      memcpy(new_state->molecule_probabilities,
	     boot_state->molecule_probabilities,move_size);
      memcpy(new_state->molecule_chemical_potentials,
	     boot_state->molecule_chemical_potentials,move_size);
      memcpy(new_state->count_to_conc,boot_state->count_to_conc,move_size);
      memcpy(new_state->conc_to_count,boot_state->conc_to_count,move_size);
      move_size = number_reactions * sizeof(double);
      memcpy(new_state->activities,boot_state->activities,move_size);
      move_size = number_reactions * max_regs_per_rxn * sizeof(double);
      memcpy(new_state->reg_constant,boot_state->reg_constant,move_size);
      memcpy(new_state->reg_exponent,boot_state->reg_exponent,move_size);
      memcpy(new_state->reg_drctn,boot_state->reg_drctn,move_size);
      memcpy(new_state->reg_species,boot_state->reg_species,move_size);

      move_size = (int64_t)sizeof(rs) * number_reactions;
      memcpy(new_state->reactions,boot_state->reactions,move_size);
      move_size = (number_reactions + 1) * sizeof(int64_t);      
      memcpy(new_state->reactions_matrix->rxn_ptrs,
	     boot_state->reactions_matrix->rxn_ptrs,move_size);
      move_size = number_molecules * sizeof(int64_t);
      memcpy(new_state->reactions_matrix->molecules_indices,
	     boot_state->reactions_matrix->molecules_indices,move_size);
      memcpy(new_state->reactions_matrix->compartment_indices,
	     boot_state->reactions_matrix->compartment_indices,move_size);
      memcpy(new_state->reactions_matrix->coefficients,
	     boot_state->reactions_matrix->coefficients,move_size);
      memcpy(new_state->reactions_matrix->text,
	     boot_state->reactions_matrix->text,move_size);

      move_size = (unique_molecules + 1) * sizeof(int64_t);
      memcpy(new_state->molecules_matrix->molecules_ptrs,
	      boot_state->molecules_matrix->molecules_ptrs,move_size);
      move_size = number_molecules * sizeof(int64_t);
      memcpy(new_state->molecules_matrix->rxn_indices,
	     boot_state->molecules_matrix->rxn_indices,move_size);
      memcpy(new_state->molecules_matrix->coefficients,
             boot_state->molecules_matrix->coefficients,move_size);

      move_size = number_reactions * sizeof(int64_t);
      memcpy(new_state->reactions_matrix->solvent_coefficients,
	     boot_state->reactions_matrix->solvent_coefficients,move_size);
      move_size = unique_molecules * sizeof(ies);
      memcpy(new_state->sorted_molecules,
	     boot_state->sorted_molecules,move_size);
      move_size = unique_compartments * sizeof(ces);
      memcpy(new_state->sorted_cmpts,
	     boot_state->sorted_cmpts,move_size);
      move_size = num_files * max_filename_len;
      memcpy(new_state->params_file,
	     boot_state->params_file,move_size);
      move_size  = (int64_t)64;
      memcpy(new_state->solvent_string,boot_state->solvent_string,move_size);
      move_size = boot_state->rxn_title_text_length;
      memcpy(new_state->rxn_title_text,
	     boot_state->rxn_title_text,move_size);
      move_size = boot_state->pathway_text_length;
      memcpy(new_state->pathway_text,
	     boot_state->pathway_text,move_size);
      move_size = boot_state->compartment_text_length;
      memcpy(new_state->compartment_text,
	     boot_state->compartment_text,move_size);
      move_size = boot_state->molecule_text_length;
      memcpy(new_state->molecules_text,
	     boot_state->molecules_text,move_size);
    }
  }
  return(success);
}

    




