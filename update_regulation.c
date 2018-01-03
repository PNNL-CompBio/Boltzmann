#include "boltzmann_structs.h"
#include "update_regulation.h"
int update_regulation (struct state_struct *state, int rxn,
		       double *counts_or_concs, int count_or_conc) {
  /*
    Compute the activity for a regulated reaction.
    Called by: update_regulations, metropolis
    Calls:     log, exp, fprintf, fflush
    Recently add a  counts_or_concs vector and a count_or_conc indicator,
    which is 1 if counts_or_concs is a vector of counts, and 0, if
    it is a vector of concentrations.
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
  double *enzyme_level;
  double direction;
  double ndirection;
  double constant;
  double power;
  double activity;
  double conc;
  double count;
  double exp_arg;
  double conc_to_power;
  double constant_to_power;
  double denom;
  double numer;
  int64_t *reg_species;
  int64_t max_regs_per_rxn;
  int64_t reg_base;
  int64_t loop_lim;
  int64_t species;
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
  enzyme_level        = state->enzyme_level;

  /*
  activity = 1.0;
  */
  activity = enzyme_level[rxn];
  loop_lim = reg_base+max_regs_per_rxn;
  for (i=reg_base;((i<loop_lim) && success);i++) {
    species = reg_species[i];
    if (species < 0) break;
    direction = reg_drctn[i];
    ndirection = 1.0 - direction;
    count = counts_or_concs[species];
    if (count_or_conc == 1) {
      molecule = (struct molecule_struct *)&sorted_molecules[species];    
      c_index   = molecule->c_index;
      if (c_index >= 0) {
	compartment = (struct compartment_struct *)&sorted_cmpts[c_index];
	multiplier  = compartment->count_to_conc;
	conc        = count * multiplier;
      } else {
	conc = count;
      }
    } else {
      conc = count;
    }
    if (conc < 0) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"update_regulation non positive conc for reaction %d, regulation %lld, species %lld\n",rxn,i,species);
	fflush(lfp);
      }
    }
    if (success) {
      power = reg_exponent[i];
      constant = reg_constant[i];
      if (constant <=0) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"update_regulation non positive regulation constant for reaction %d, regulation %lld\n",rxn,i);
	}
	fflush(lfp);
      }
    }
    if (success) {
      if (conc > 0) {
	/*
	  we want to use conc^exponent and constant^exponent as the pieces.
	  Here instead of using the pow function which is very inefficient
	  we us a^x  = exp(x*log(a))
	*/
	exp_arg  = power * log(conc);
	conc_to_power = exp(exp_arg);
	exp_arg  = power * log(constant);
	constant_to_power = exp(exp_arg);
	numer = (ndirection * constant_to_power) + (direction * conc_to_power);
	denom = constant_to_power + conc_to_power;
	if (denom != 0.0) {
	  activity *= (numer/denom);
	}
      } else { 
	/*
	  No regulator concentration so regulation does not happen.
	  If regulator is a positive regulator, activity is 0,
	  if its a negative regulator activity is 1.
	*/
	activity *= ndirection;
      }
    }
  } /* end for (i...) */
  /*
  if (activity > 0.5) {
    activity = 1.0;
  } else {
    activity = 0.0;
  }
  */
  activities[rxn] = activity;
  return (success);
}
