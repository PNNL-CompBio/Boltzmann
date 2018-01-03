#include "boltzmann_structs.h"
#include "boltzmann_flatten_reactions.h"
int boltzmann_flatten_reactions(struct state_struct *state,
				void *fstate, 
				int direction, 
				int *word_pos_p,
				FILE *lfp) {
  /*
    Flatten the reactions_struct.
    Called by: boltzmann_flatten_oneway_data
  */
  struct reaction_struct *reaction;
  int64_t *lfstate;
  double  *dfstate;
  int64_t number_reactions;
  int *ltoi;

  int success;
  int word_pos;

  int reactions_len;
  int words_per_rxn;

  int i;
  int padi;
  
  success = 1;
  number_reactions = state->number_reactions;
  reaction = state->reactions;
  word_pos = *word_pos_p;

  lfstate = (int64_t *)fstate;
  dfstate = (double  *)fstate;

  words_per_rxn = 18;
  reactions_len = 2 + (words_per_rxn * number_reactions);

  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = (int64_t)reactions_len;
  } else {
    if (lfstate[word_pos] != (int64_t)reactions_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_reactions: Error reactions struct lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],reactions_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    } 
    for (i=0; i<number_reactions; i++) {
      word_pos += 1;
      if (direction == 0) {
        lfstate[word_pos] = reaction->title;
      } else {
        reaction->title = lfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        lfstate[word_pos] = reaction->pathway;
      } else {
        reaction->pathway = lfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        lfstate[word_pos] = reaction->lcompartment;
      } else {
        reaction->lcompartment = lfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        lfstate[word_pos] = reaction->rcompartment;
      } else {
        reaction->rcompartment = lfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->delta_g0;
      } else {
        reaction->delta_g0 = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->unit_v;
      } else {
        reaction->unit_v = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->k_epsilon;
      } else {
        reaction->k_epsilon = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->activity;
      } else {
        reaction->activity = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->enzyme_level;
      } else {
        reaction->enzyme_level = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->forward_rc;
      } else {
        reaction->forward_rc = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->reverse_rc;
      } else {
        reaction->reverse_rc = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->ph;
      } else {
        reaction->ph = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->temp_kelvin;
      } else {
        reaction->temp_kelvin = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
        dfstate[word_pos] = reaction->ionic_strength;
      } else {
        reaction->ionic_strength = dfstate[word_pos];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
        ltoi[0] = reaction->num_reactants;
        ltoi[1] = reaction->num_products;
      } else {
        reaction->num_reactants = ltoi[0];
        reaction->num_products  = ltoi[1];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
        ltoi[0] = reaction->self_id;
        ltoi[1] = reaction->unit_i;
      } else {
        reaction->self_id = ltoi[0];
        reaction->unit_i  = ltoi[1];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
        ltoi[0] = reaction->left_compartment;
        ltoi[1] = reaction->right_compartment;
      } else {
        reaction->left_compartment = ltoi[0];
        reaction->right_compartment = ltoi[1];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
        ltoi[0] = reaction->deltag0_computed;
        ltoi[1] = reaction->num_regulators;
      } else {
        reaction->deltag0_computed = ltoi[0];
        reaction->num_regulators = ltoi[1];
      }
      reaction += 1; /* Caution address arithmetic here */
    } /* end for(i... )*/
  }
  *word_pos_p = word_pos;
  return(success);
}
 


