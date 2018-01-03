#include "boltzmann_structs.h"
#include "boltzmann_flatten_oneway_vecs.h"
int boltzmann_flatten_oneway_vecs(struct state_struct *state,
				  void *fstate, 
				  int direction, 
				  int *word_pos_p,
				  FILE *lfp) {
  /*
    Flatten the oneway input vector fields
    Called by: boltzmann_flatten_oneway_data
    Calls:     fprintf, fflush

  */
  double  *dfstate;
  int64_t *lfstate;
  int64_t number_reactions;
  int64_t nunique_molecules;
  int64_t rxn_len;
  int64_t mlcl_len;
  int64_t cs_len;

  void *dg0s;
  void *ke;
  void *rke;
  void *forward_rc;
  void *reverse_rc;
  void *kss;
  void *kssr;
  void *activities;
  void *enzyme_level;
  void *reg_constant;
  void *reg_exponent;
  void *reg_drctn;
  void *reg_species;
  void *kss_e_val;
  void *kss_u_val;
  void *molecule_dg0tfs;
  void *molecule_probabilities;
  void *molecule_chemical_potentials;
  void *count_to_conc;
  void *conc_to_count;
  void *coeff_sum;

  int success;
  int word_pos;

  int oneway_vecs_len;
  int cs_words;
  
  success = 1;
  nunique_molecules = state->nunique_molecules;
  number_reactions  = state->number_reactions;

  word_pos = *word_pos_p;
  lfstate  = (int64_t *)fstate;
  dfstate  = (double *)fstate;

  rxn_len           = number_reactions * sizeof(double);
  mlcl_len          = nunique_molecules * sizeof(double);
  cs_words          = (number_reactions + (number_reactions & 1))/2;
  cs_len            = cs_words * sizeof(double);

  oneway_vecs_len = (int64_t)2 + (13 * number_reactions) + 
    (7*nunique_molecules) + cs_words;
  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = oneway_vecs_len;
  } else {
    if (lfstate[word_pos] != oneway_vecs_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_oneway_vecs: Error oneway_vecs lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],oneway_vecs_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    }
    
    word_pos += 1;
    dg0s = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(dg0s,state->dg0s,rxn_len);
    } else {
      memcpy(state->dg0s,dg0s,rxn_len);
    }
    word_pos += number_reactions;

    ke = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(ke,state->ke,rxn_len);
    } else {
      memcpy(state->ke,ke,rxn_len);
    }
    word_pos += number_reactions;

    rke = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(rke,state->rke,rxn_len);
    } else {
      memcpy(state->rke,rke,rxn_len);
    }
    word_pos += number_reactions;

    forward_rc = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(forward_rc,state->forward_rc,rxn_len);
    } else {
      memcpy(state->forward_rc,forward_rc,rxn_len);
    }
    word_pos += number_reactions;

    reverse_rc = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(reverse_rc,state->reverse_rc,rxn_len);
    } else {
      memcpy(state->reverse_rc,reverse_rc,rxn_len);
    }
    word_pos += number_reactions;

    kss = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(kss,state->kss,rxn_len);
    } else {
      memcpy(state->kss,kss,rxn_len);
    }
    word_pos += number_reactions;

    kssr = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(kssr,state->kssr,rxn_len);
    } else {
      memcpy(state->kssr,kssr,rxn_len);
    }
    word_pos += number_reactions;

    activities = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(activities,state->activities,rxn_len);
    } else {
      memcpy(state->activities,activities,rxn_len);
    }
    word_pos += number_reactions;

    enzyme_level = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(enzyme_level,state->enzyme_level,rxn_len);
    } else {
      memcpy(state->enzyme_level,enzyme_level,rxn_len);
    }
    word_pos += number_reactions;

    reg_constant = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(reg_constant,state->reg_constant,rxn_len);
    } else {
      memcpy(state->reg_constant,reg_constant,rxn_len);
    }
    word_pos += number_reactions;

    reg_exponent = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(reg_exponent,state->reg_exponent,rxn_len);
    } else {
      memcpy(state->reg_exponent,reg_exponent,rxn_len);
    }
    word_pos += number_reactions;

    reg_drctn = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(reg_drctn,state->reg_drctn,rxn_len);
    } else {
      memcpy(state->reg_drctn,reg_drctn,rxn_len);
    }
    word_pos += number_reactions;

    reg_species = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(reg_species,state->reg_species,rxn_len);
    } else {
      memcpy(state->reg_species,reg_species,rxn_len);
    }
    word_pos += number_reactions;

    kss_e_val = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(kss_e_val,state->kss_e_val,mlcl_len);
    } else {
      memcpy(state->kss_e_val,kss_e_val,mlcl_len);
    }
    word_pos += nunique_molecules;

    kss_u_val = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(kss_u_val,state->kss_u_val,mlcl_len);
    } else {
      memcpy(state->kss_u_val,kss_u_val,mlcl_len);
    }
    word_pos += nunique_molecules;

    molecule_dg0tfs = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(molecule_dg0tfs,state->molecule_dg0tfs,mlcl_len);
    } else {
      memcpy(state->molecule_dg0tfs,molecule_dg0tfs,mlcl_len);
    }
    word_pos += nunique_molecules;

    molecule_probabilities = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(molecule_probabilities,state->molecule_probabilities,mlcl_len);
    } else {
      memcpy(state->molecule_probabilities,molecule_probabilities,mlcl_len);
    }
    word_pos += nunique_molecules;

    molecule_chemical_potentials = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(molecule_chemical_potentials,
	     state->molecule_chemical_potentials,mlcl_len);
    } else {
      memcpy(state->molecule_chemical_potentials,
	     molecule_chemical_potentials,mlcl_len);
    }
    word_pos += nunique_molecules;

    count_to_conc = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(count_to_conc,state->count_to_conc,mlcl_len);
    } else {
      memcpy(state->count_to_conc,count_to_conc,mlcl_len);
    }
    word_pos += nunique_molecules;

    conc_to_count = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(conc_to_count,state->conc_to_count,mlcl_len);
    } else {
      memcpy(state->conc_to_count,conc_to_count,mlcl_len);
    }
    /*
      NB we need to leave word_pos pointing at the last word set.
    */
    word_pos += nunique_molecules;

    coeff_sum = (void*)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(coeff_sum,state->coeff_sum,cs_len);
    } else {
      memcpy(state->coeff_sum,coeff_sum,cs_len);
    }
    word_pos += cs_words - 1;
  } /* end if (success) */
  *word_pos_p = word_pos;
  return(success);
}
