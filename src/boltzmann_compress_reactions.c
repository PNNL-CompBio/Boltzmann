#include "boltzmann_structs.h"
#include "boltzmann_compress_reactions.h"
int boltzmann_compress_reactions(struct state_struct *state) {
  /*
    Remove reactions that have a use_rxn value of 0.
    Called by: parse_reactions_file

    Reads the following_fields in state:
      number_reactions,
      reg_constant,
      reg_exponent,
      reg_species,
      reg_drctn,
      coeff_sum,
      use_rxn,
      reactions,
      reactions_matrix, and subfields
      pathway_text,
      compartment_text,
      molecules_text,
      raw_molecules_text,
      regulation_text
      activities,
      enzyme_level

    Sets the following fields in state:
      number_reactions,
      reg_constant,
      reg_exponent,
      reg_species,
      reg_drctn,
      coeff_sum,
      reactions,
      reactions_matrix, and subfields
      pathway_text,
      activities,
      enzyme_level
    
  */
  struct reaction_struct *reactions;
  struct reaction_struct *reaction;
  struct reaction_struct *next_reaction;
  struct reactions_matrix_struct *rxns_matrix;
  double *activities;
  double *enzyme_level;
  double *reg_constant;
  double *reg_exponent;
  double *reg_drctn;
  double *forward_rc;
  double *reverse_rc;
  double *recip_coeffs;
  double  *coefficients;
  double  *coeff_sum;
  int64_t *rxn_ptrs;
  int64_t *molecules_indices;
  int64_t *reg_species;
  int  *use_rxn;

  int success;
  int rxns;

  int mol_pos;
  int reg_pos;

  int next_mol_pos;
  int next_reg_pos;

  int i;
  int j;

  int max_regs_per_rxn;
  int num_reactions;

  int num_species;
  int padi;


  FILE *lfp;
  FILE *efp;

  success = 1;
  lfp          	   = state->lfp;
  reg_constant 	   = state->reg_constant;
  reg_exponent 	   = state->reg_exponent;
  reg_species  	   = state->reg_species;
  reg_drctn    	   = state->reg_drctn;
  coeff_sum        = state->coeff_sum;
  use_rxn          = state->use_rxn;
  max_regs_per_rxn = (int)state->max_regs_per_rxn;
  num_reactions    = state->number_reactions;

  activities                  = state->activities;
  enzyme_level                = state->enzyme_level;
  reactions                   = state->reactions;
  rxns_matrix                 = state->reactions_matrix;
  forward_rc                  = state->forward_rc;
  reverse_rc                  = state->reverse_rc;
  rxn_ptrs                    = rxns_matrix->rxn_ptrs;
  molecules_indices           = rxns_matrix->molecules_indices;
  coefficients                = rxns_matrix->coefficients;
  recip_coeffs                = rxns_matrix->recip_coeffs;
  reaction                    = reactions;
  next_reaction               = reaction;

  rxns = 0;
  mol_pos = 0;
  next_mol_pos = 0;
  reg_pos = 0;
  next_reg_pos = 0;
  for (i=0;i<num_reactions;i++) {
    next_mol_pos = rxn_ptrs[i];
    if (use_rxn[i] == 1) {
      rxn_ptrs[rxns]     = mol_pos;
      activities[rxns]   = activities[i];
      enzyme_level[rxns] = enzyme_level[i];
      coeff_sum[rxns]    = coeff_sum[i];
      forward_rc[rxns]   = forward_rc[i];
      reverse_rc[rxns]   = reverse_rc[rxns];
      reaction->self_id = rxns;
      num_species = rxn_ptrs[i+1] - rxn_ptrs[i];
      /*
	Shift the molecule_indices, coefficients, and recip_coeffs field.
      */
      for (j=0;j<num_species;j++) {
	molecules_indices[mol_pos+j] = molecules_indices[next_mol_pos+j];
	coefficients[mol_pos+j] = coefficients[next_mol_pos+j];
	recip_coeffs[mol_pos+j] = recip_coeffs[next_mol_pos+j];
      }
      /*
	Move fields of reaction struct.
      */
      reaction->title = next_reaction->title;
      reaction->pathway = next_reaction->pathway;
      reaction->lcompartment = next_reaction->lcompartment;
      reaction->rcompartment = next_reaction->rcompartment;
      reaction->delta_g0     = next_reaction->delta_g0;
      reaction->unit_v       = next_reaction->unit_v;
      reaction->k_epsilon    = next_reaction->k_epsilon;
      reaction->activity     = next_reaction->activity;
      reaction->enzyme_level = next_reaction->enzyme_level;
      reaction->forward_rc   = next_reaction->forward_rc;
      reaction->reverse_rc   = next_reaction->reverse_rc;
      
      reaction->ph           = next_reaction->ph;
      reaction->temp_kelvin  = next_reaction->temp_kelvin;
      reaction->ionic_strength = next_reaction->ionic_strength;
      reaction->num_reactants  = next_reaction->num_reactants;
      reaction->num_products   = next_reaction->num_products;
      reaction->self_id        = rxns;
      reaction->unit_i         = next_reaction->unit_i;
      reaction->left_compartment  = next_reaction->left_compartment;
      reaction->right_compartment = next_reaction->right_compartment;
      reaction->deltag0_computed  = next_reaction->deltag0_computed;
      reaction->num_regulators    = next_reaction->num_regulators;
      reaction->coefficient_sum   = next_reaction->coefficient_sum;
      /*
	Move regulation fields.
      */
      for (j=0;j<max_regs_per_rxn;j++) {
	reg_constant[reg_pos + j] = reg_constant[next_reg_pos + j];
	reg_exponent[reg_pos + j] = reg_exponent[next_reg_pos + j];
	reg_drctn[reg_pos + j]    = reg_drctn[next_reg_pos + j];
	reg_species[reg_pos + j]  = reg_species[next_reg_pos + j];
      }
      /*
	Advance mol_pos and reg_base and set rxn_ptrs field, and reaction
      */
      mol_pos  += num_species;
      reg_pos  += max_regs_per_rxn;
      reaction += 1; /* Caution address arithmetic here. */
      rxns     += 1;
    } /* end if rxn_use[i] */
    /*
      advance next_reaction and next_reg_base;
    */
    next_reaction += 1; /* Caution address arithmetic here */
    next_reg_pos  += max_regs_per_rxn;
  }
  /*
    Set the last rxn_ptrs field.
  */
  rxn_ptrs[rxns] = mol_pos;
  /*
    reset the number_reactions field in state.
  */
  state->number_reactions = rxns;
  return(success);
}
