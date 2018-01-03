#include "boltzmann_structs.h"
#include "boltzmann_flatten_external.h"
int boltzmann_flatten_external(struct state_struct *state,
			       void *flattened, 
			       int direction, 
			       int *word_pos_p, 
			       FILE *lfp) {
  /* 
    Transfer external constants, to or from flattened
    Called by: boltzmann_flatten_state
    While there are only 7 external constants right now, we 
    allocate space for 10 for ease of expansion.
    Word_pos is assumed to point to the last set/read word position 
    in flattened.
  */
  int64_t *lflattened;
  double  *dflattened;
  int success;
  int word_pos;
  success  = 1;
  word_pos = *word_pos_p;
  lflattened = (int64_t *)flattened;
  dflattened = (double  *)flattened;
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = 12; /* 10 words + 2 for meta data*/
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = -5; /* packed fields */
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = state->version_no;
  } else {
    state->version_no = lflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = state->agent_type;
  } else {
    state->agent_type = lflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = state->mpi_rank;
  } else {
    state->mpi_rank = lflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = state->thread_id;
  } else {
    state->thread_id = lflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = state->current_counts_offset;
  } else {
    state->current_counts_offset = lflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->x_coord;
  } else {
    state->x_coord = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->y_coord;
  } else {
    state->y_coord = dflattened[word_pos];
  }
  word_pos += 1;
  if (direction == 0) {
    dflattened[word_pos] = state->z_coord;
  } else {
    state->z_coord = dflattened[word_pos];
  }
  /*
    Skip 2 words reserved for future external variables.
  */
  word_pos += 2;
  *word_pos_p = word_pos;
  return(success);
}

