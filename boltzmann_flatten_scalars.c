#include "boltzmann_structs.h"
#include "boltzmann_flatten_scalars.h"
int boltzmann_flatten_scalars(struct state_struct *state,
			      void *flattened, int direction, int *word_pos_p, 
			      FILE *lfp) {
  /* 
    Transfer scalar variabless, to or from flattened
    Called by: boltzmann_flatten_state
    
    Right now this section has length 96 words, 
    2 meta data words, 2 meta data words, 52 int64_t's 
    2 meta_data words, 38 doubles,
 
    Word_pos (*word_pos_p) is assumed to point to the last set/read 
    word position in flattened.
  */
  int64_t *lflattened;
  double  *dflattened;
  int success;
  int word_pos;
  success  = 1;
  word_pos = *word_pos_p;
  lflattened = (int64_t *)flattened;
  dflattened = (double  *)flattened;
  word_pos += 1; /* 32 */
  if (direction == 0) {
    lflattened[word_pos] = (int64_t)96;
  }
  word_pos += 1; /* 33 */
  if (direction == 0) {
    lflattened[word_pos] = (int64_t)2;
  }
  word_pos += 1; /* 34 */
  if (direction == 0) {
    lflattened[word_pos] = (int64_t)54;
  }
  word_pos += 1; /* 35 */
  if (direction == 0) {
    lflattened[word_pos] = (int64_t)-2;
  }
  word_pos += 1; /* 36 */
  if (direction == 0) {
    lflattened[word_pos] = state->align_len;
  } else {
    state->align_len = lflattened[word_pos];
  }
  word_pos += 1; /* 37 */
  if (direction == 0) {
    lflattened[word_pos] = state->align_mask;
  } else {
    state->align_mask = lflattened[word_pos];
  }
  word_pos += 1; /* 38 */
  if (direction == 0) {
    lflattened[word_pos] = state->number_reactions;
  } else {
    state->number_reactions = lflattened[word_pos];
  }
  word_pos += 1; /* 39 */
  if (direction == 0) {
    lflattened[word_pos] = state->number_molecules;
  } else {
    state->number_molecules = lflattened[word_pos];
  }
  word_pos += 1; /* 40 */
  if (direction == 0) {
    lflattened[word_pos] = state->nunique_molecules;
  } else {
    state->nunique_molecules = lflattened[word_pos];
  }
  word_pos += 1; /* 41 */
  if (direction == 0) {
    lflattened[word_pos] = state->max_molecule_len;
  } else {
    state->max_molecule_len = lflattened[word_pos];
  }
  word_pos += 1; /* 42 */
  if (direction == 0) {
    lflattened[word_pos] = state->min_molecule_len;
  } else {
    state->min_molecule_len = lflattened[word_pos];
  }
  word_pos += 1; /* 43 */
  if (direction == 0) {
    lflattened[word_pos] = state->sum_molecule_len;
  } else {
    state->sum_molecule_len = lflattened[word_pos];
  }
  word_pos += 1; /* 44 */
  if (direction == 0) {
    lflattened[word_pos] = state->number_compartments;
  } else {
    state->number_compartments = lflattened[word_pos];
  }
  word_pos += 1; /* 45 */
  if (direction == 0) {
    lflattened[word_pos] = state->nunique_compartments;
  } else {
    state->nunique_compartments = lflattened[word_pos];
  }
  word_pos += 1; /* 46 */
  if (direction == 0) {
    lflattened[word_pos] = state->max_compartment_len;
  } else {
    state->max_compartment_len = lflattened[word_pos];
  }
  word_pos += 1; /* 47 */
  if (direction == 0) {
    lflattened[word_pos] = state->min_compartment_len;
  } else {
    state->min_compartment_len = lflattened[word_pos];
  }
  word_pos += 1; /* 48 */
  if (direction == 0) {
    lflattened[word_pos] = state->sum_compartment_len;
  } else {
    state->sum_compartment_len = lflattened[word_pos];
  }
  word_pos += 1; /* 49 */
  if (direction == 0) {
    lflattened[word_pos] = state->solvent_pos;
  } else {
    state->solvent_pos = lflattened[word_pos];
  }
  word_pos += 1; /* 50 */
  if (direction == 0) {
    lflattened[word_pos] = state->num_fixed_concs;
  } else {
    state->num_fixed_concs = lflattened[word_pos];
  }
  word_pos += 1; /* 51 */
  if (direction == 0) {
    lflattened[word_pos] = state->num_files;
  } else {
    state->num_files = lflattened[word_pos];
  }
  word_pos += 1; /* 52 */
  if (direction == 0) {
    lflattened[word_pos] = state->max_filename_len;
  } else {
    state->max_filename_len = lflattened[word_pos];
  }
  word_pos += 1; /* 53 */
  if (direction == 0) {
    lflattened[word_pos] = state->reaction_file_length;
  } else {
    state->reaction_file_length = lflattened[word_pos];
  }
  word_pos += 1; /* 54 */
  if (direction == 0) {
    lflattened[word_pos] = state->warmup_steps;
  } else {
    state->warmup_steps = lflattened[word_pos];
  }
  word_pos += 1; /* 55 */
  if (direction == 0) {
    lflattened[word_pos] = state->record_steps;
  } else {
    state->record_steps = lflattened[word_pos];
  }
  word_pos += 1; /* 56 */
  if (direction == 0) {
    lflattened[word_pos] = state->ode_solver_choice;
  } else {
    state->ode_solver_choice = lflattened[word_pos];
  }
  word_pos += 1; /* 57 */
  if (direction == 0) {
    lflattened[word_pos] = state->delta_concs_choice;
  } else {
    state->delta_concs_choice = lflattened[word_pos];
  }
  word_pos += 1; /* 58 */
  if (direction == 0) {
    lflattened[word_pos] = state->base_reaction;
  } else {
    state->base_reaction = lflattened[word_pos];
  }
  word_pos += 1; /* 59 */
  if (direction == 0) {
    lflattened[word_pos] = state->number_base_reaction_reactants;
  } else {
    state->number_base_reaction_reactants = lflattened[word_pos];
  }
  word_pos += 1; /* 60 */
  if (direction == 0) {
    lflattened[word_pos] = state->adjust_steady_state;
  } else {
    state->adjust_steady_state = lflattened[word_pos];
  }
  word_pos += 1; /* 61 */
  if (direction == 0) {
    lflattened[word_pos] = state->molecules_or_conc;
  } else {
    state->molecules_or_conc = lflattened[word_pos];
  }
  word_pos += 1; /* 62 */
  if (direction == 0) {
    lflattened[word_pos] = state->concs_or_counts;
  } else {
    state->concs_or_counts = lflattened[word_pos];
  }
  word_pos += 1; /* 63 */
  if (direction == 0) {
    lflattened[word_pos] = state->free_energy_format;
  } else {
    state->free_energy_format = lflattened[word_pos];
  }
  word_pos += 1; /* 64 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_activities;
  } else {
    state->use_activities = lflattened[word_pos];
  }
  word_pos += 1; /* 65 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_deq;
  } else {
    state->use_deq = lflattened[word_pos];
  }
  word_pos += 1; /* 66 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_pseudoisomers;
  } else {
    state->use_pseudoisomers = lflattened[word_pos];
  }
  word_pos += 1; /* 67 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_metropolis;
  } else {
    state->use_metropolis = lflattened[word_pos];
  }
  word_pos += 1; /* 68 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_regulation;
  } else {
    state->use_regulation = lflattened[word_pos];
  }
  word_pos += 1; /* 69 */
  if (direction == 0) {
    lflattened[word_pos] = state->max_regs_per_rxn;
  } else {
    state->max_regs_per_rxn = lflattened[word_pos];
  }
  word_pos += 1; /* 70 */
  if (direction == 0) {
    lflattened[word_pos] = state->no_round_from_deq;
  } else {
    state->no_round_from_deq = lflattened[word_pos];
  }
  word_pos += 1; /* 71 */
  if (direction == 0) {
    lflattened[word_pos] = state->adjust_steady_state;
  } else {
    state->adjust_steady_state = lflattened[word_pos];
  }
  word_pos += 1; /* 72 */
  if (direction == 0) {
    lflattened[word_pos] = state->print_output;
  } else {
    state->print_output = lflattened[word_pos];
  }
  word_pos += 1; /* 73 */
  if (direction == 0) {
    lflattened[word_pos] = state->rxn_view_freq;
  } else {
    state->rxn_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 74 */
  if (direction == 0) {
    lflattened[word_pos] = state->rxn_view_hist_length;
  } else {
    state->rxn_view_hist_length = lflattened[word_pos];
  }
  word_pos += 1; /* 75 */
  if (direction == 0) {
    lflattened[word_pos] = state->ode_rxn_view_freq;
  } else {
    state->ode_rxn_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 76 */
  if (direction == 0) {
    lflattened[word_pos] = state->lklhd_view_freq;
  } else {
    state->lklhd_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 77 */
  if (direction == 0) {
    lflattened[word_pos] = state->lklhd_view_freq;
  } else {
    state->lklhd_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 78 */
  if (direction == 0) {
    lflattened[word_pos] = state->count_view_freq;
  } else {
    state->count_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 79 */
  if (direction == 0) {
    lflattened[word_pos] = state->fe_view_freq;
  } else {
    state->fe_view_freq = lflattened[word_pos];
  }
  word_pos += 1; /* 80 */
  if (direction == 0) {
    lflattened[word_pos] = state->print_ode_concs;
  } else {
    state->print_ode_concs = lflattened[word_pos];
  }
  word_pos += 1; /* 81 */
  if (direction == 0) {
    lflattened[word_pos] = state->usage;
  } else {
    state->usage = lflattened[word_pos];
  }
  word_pos += 1; /* 82 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_dgzero;
  } else {
    state->use_dgzero = lflattened[word_pos];
  }
  word_pos += 1; /* 83 */
  if (direction == 0) {
    lflattened[word_pos] = state->use_bulk_water;
  } else {
    state->use_bulk_water = lflattened[word_pos];
  }
  /*
    Leave a litle extra space, round up to 52 data words 
  */
  word_pos += 4;
  if (direction == 0) {
    lflattened[word_pos] = 40;
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = -1;
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ideal_gas_r;
  } else {
    state->ideal_gas_r = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->temp_kelvin;
  } else {
    state->temp_kelvin = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->avogadro;
  } else {
    state->avogadro = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->recip_avogadro;
  } else {
    state->recip_avogadro = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->min_count;
  } else {
    state->min_count = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->min_conc;
  } else {
    state->min_conc = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ph;
  } else {
    state->ph = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ionic_strength;
  } else {
    state->ionic_strength = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->rt;
  } else {
    state->rt = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->m_r_rt;
  } else {
    state->m_r_rt = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->m_rt;
  } else {
    state->m_rt = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->cals_per_joule;
  } else {
    state->cals_per_joule = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->default_initial_count;
  } else {
    state->default_initial_count = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->dg_forward;
  } else {
    state->dg_forward = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->entropy;
  } else {
    state->entropy = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->current_concentrations_sum;
  } else {
    state->current_concentrations_sum = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->default_volume;
  } else {
    state->default_volume = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->recip_default_volume;
  } else {
    state->recip_default_volume = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->conc_units;
  } else {
    state->conc_units = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->conc_units;
  } else {
    state->conc_units = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ntotal_opt;
  } else {
    state->ntotal_opt = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ntotal_exp;
  } else {
    state->ntotal_exp = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->ode_t_final;
  } else {
    state->ode_t_final = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->max_log_g0_sum;
  } else {
    state->max_log_g0_sum = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->dg0_scale_factor;
  } else {
    state->dg0_scale_factor = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->kf_base_reaction;
  } else {
    state->kf_base_reaction = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->flux_scaling;
  } else {
    state->flux_scaling = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->min_molecule_dg0tf;
  } else {
    state->min_molecule_dg0tf = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->epsilon;
  } else {
    state->epsilon = dflattened[word_pos];
  }
  /*
    Leave some extra space.
  */
  word_pos += 9;
  *word_pos_p  = word_pos;
  return(success);
}
