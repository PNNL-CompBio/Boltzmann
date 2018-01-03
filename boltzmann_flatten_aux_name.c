#include "boltzmann_structs.h"
#include "boltzmann_flatten_aux_name.h"
int boltzmann_flatten_aux_name(struct state_struct *state,
			       void *flattened, int direction, int *word_pos_p, 
			       FILE *lfp) {
  /* 
    Transfer auxiliary file name, to or from flattened

    While there are only 7 external constants right now, we 
    allocate space for 10 for ease of expansion.
    Word_pos is assumed to point to the last set/read word position 
    in flattened.

    Called by: boltzmann_flatten
    Calls:     strcpy
  */
  int64_t *lflattened;
  char *aux_filename;
  int success;
  int word_pos;
  success  = 1;
  word_pos = *word_pos_p;
  lflattened = (int64_t *)flattened;
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = 18;
  }
  word_pos += 1;
  if (direction == 0) {
    lflattened[word_pos] = -4;
  }
  word_pos += 1;
  if (direction == 0) {
    aux_filename = (char *) &lflattened[word_pos];
    if (strlen(state->aux_data_file) < 127) {
      strcpy(aux_filename,state->aux_data_file);
    } else {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_aux_name, aux file name is to long\n");
	fflush(lfp);
      }
    }
  }
  word_pos += 15;
  *word_pos_p = word_pos;
  return(success);
}
  
