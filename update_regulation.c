#include "boltzmann_structs.h"
#include "update_regulation.h"
int update_regulation (struct state_struct *state, int rxn) {
  /*
    Compute the activity for a regulated reaction.
    Called by: update_regulations, metropolis
    Calls:     log, exp, fprintf, fflush
  */
  struct molecule_struct *sorted_molecules;  
  struct molecule_struct *molecule;  
  struct compartment_struct *sorted_cmpts;
  struct compartment_struct *compartment;
  double multiplier;
  double *reg_constant;
  double *reg_exponent;
  double *reg_drctn;
  double *current_counts;
  double *activities;
  double direction;
  double ndirection;
  double constant;
  double exponent;
  double activity;
  double conc;
  double exp_arg;
  double c_to_exp;
  double denom;
  double numer;
  int64_t *reg_species;
  int64_t max_regs_per_rxn;
  int64_t reg_base;
  int64_t loop_lim;
  int64_t species;
  int64_t count;
  int64_t i;
  int64_t one_l;


  int success;
  int c_index;
  FILE *lfp;
  FILE *efp;
  success             = 1;
  one_l               = (int64_t)1;
  lfp                 = state->lfp;
  max_regs_per_rxn    = state->max_regs_per_rxn;
  current_counts      = state->current_counts;
  reg_base            = max_regs_per_rxn * rxn;
  reg_constant        = state->reg_constant;
  reg_exponent        = state->reg_exponent;
  reg_species         = state->reg_species;
  reg_drctn           = state->reg_drctn;
  sorted_molecules    = state->sorted_molecules;
  sorted_cmpts        = state->sorted_compartments;
  activities          = state->activities;
  
  activity = 1.0;
  loop_lim = reg_base+max_regs_per_rxn;
  for (i=reg_base;((i<loop_lim) && success);i++) {
    species = reg_species[i];
    if (species < 0) break;
    direction = reg_drctn[i];
    ndirection = 1.0 - direction;
    count = current_counts[species];
    if (count > 0.0) {
      molecule = (struct molecule_struct *)&sorted_molecules[rxn];
      c_index   = molecule->c_index;
      if (c_index >= 0) {
	compartment = (struct compartment_struct *)&sorted_cmpts[c_index];
	multiplier  = compartment->count_to_conc;
	conc        = count * multiplier;
      } else {
	conc = count;
      }
      if (conc <=0) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"update_regulation non positive multiplier for compartment %d, reaction was %d\n",c_index,rxn);
	  fflush(lfp);
	}
      }
      if (success) {
	exponent = reg_exponent[i];
	constant = reg_constant[i];
	exp_arg  = exponent * log(conc);
	c_to_exp = exp(exp_arg);
	numer = (ndirection * constant) + (direction * c_to_exp);
	denom = constant + c_to_exp;
	if (denom != 0.0) {
	  activity *= (numer/denom);
	}
      }
    } else { 
      /*
	No regulator so regulator does not happen.
	If regulator is a positive regulator, activity is 0,
	if its a negative regulator activity is 1.
      */
      activity *= ndirection;
    }
  } /* end for (i...) */
  if (activity > 0.5) {
    activity = 1.0;
  } else {
    activity = 0.0;
  }
  activities[rxn] = activity;
  return (success);
}
