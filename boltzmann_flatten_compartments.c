#include "boltzmann_structs.h"
#include "boltzmann_flatten_compartments.h"
int boltzmann_flatten_compartments(struct state_struct *state,
				   void *fstate, 
				   int direction, 
				   int *word_pos_p,
				   FILE *lfp) {
  /*
    Flatten the compartment_structs.
    Called by: boltzmann_flatten_oneway_data
    Calls:     fprintf, fflush
  */
  struct compartment_struct *compartment;
  int64_t *lfstate;
  double  *dfstate;
  int64_t  nunique_compartments;
  int *ltoi;
  int success;
  int word_pos;

  int compartments_len;
  int words_per_compartment;

  int i;
  int padi;
  
  success = 1;
  nunique_compartments = state->nunique_compartments;
  compartment          = state->sorted_compartments;;
  word_pos = *word_pos_p;

  lfstate = (int64_t *)fstate;
  dfstate = (double *)fstate;

  words_per_compartment = 8;
  compartments_len = (int64_t)2 + (nunique_compartments * words_per_compartment);

  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = compartments_len;
  } else {
    if (lfstate[word_pos] != compartments_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_compartments: Error compartment struct lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],compartments_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    } 
    for (i=0;i<nunique_compartments;i++) {
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->volume;
      } else {
	compartment->volume = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->recip_volume;
      } else {
	compartment->recip_volume = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->ntotal_exp;
      } else {
	compartment->ntotal_exp = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->ntotal_opt;
      } else {
	compartment->ntotal_opt = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->conc_to_count;
      } else {
	compartment->conc_to_count = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	dfstate[word_pos] = compartment->count_to_conc;
      } else {
	compartment->count_to_conc = dfstate[word_pos];
      }
      word_pos += 1;
      if (direction == 0) {
	lfstate[word_pos] = compartment->string;
      } else {
	compartment->string = lfstate[word_pos];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
	ltoi[0] = compartment->c_index;
	ltoi[1] = compartment->g_index;
      } else {
	compartment->c_index = ltoi[0];
	compartment->g_index = ltoi[1];
      }
      compartment += 1; /* Caution address arithmetic here */
    } /* end for(i...) */
  } /* end if(success) */
  *word_pos_p = word_pos;
  return(success);
}
