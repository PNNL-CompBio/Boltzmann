#include "boltzmann_structs.h"
#include "compute_flattened_state_size.h"
#include "alloc0.h"
#include "boltzmann_flatten_alloc.h"
int boltzmann_flatten_alloc (struct state_struct **state_p,
			     void **flattened_state_p,
			     int direction, int *word_pos_p, FILE *lfp) {
  /*
    Depending on direction allocate flattened_state 
    (if direction == 0 and *flattened_state == NULL) 
    or allocate and load state (if direction == 1)
    Called by: boltzmann_flatten_state
    Calls:     compute_flattened_state_size,
               calloc,
	       alloc0,
	       fprintf,
	       fflush
  */
  struct state_struct *state;
  void *flattened_state;
  int64_t *lfstate;
  int64_t word_size;
  int64_t flattened_state_size;
  int64_t one_l;
  int64_t ask_for;
  int success;
  int setup;
  success = 1;
  one_l   = (int64_t)1;
  if (direction == 0) {
    /*
      Saving a state.
      Allocating flattened_state
      Loading fields from state into flattened_state.
    */
    state = *state_p;
    flattened_state_size = compute_flattened_state_size(state);
    if (*flattened_state_p == NULL) {
      word_size = (int64_t)8;
      ask_for = flattened_state_size * word_size;
      flattened_state = (void*)calloc(one_l,ask_for);
      if (flattened_state == NULL) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_flatten_alloc: could not allocate %ld bytes for flattened_state\n", ask_for);
	  fflush(lfp);
	}
      } else {
	lfstate = (int64_t*)flattened_state;
	lfstate[0] = flattened_state_size;
	/*
	  Basic groupings: externals, auxiliary-filename, scalars, oneway, twoway, 
	*/
	lfstate[1] = 5;
	*flattened_state_p = flattened_state;
      }
    } else {
      /*
	We are filling flattened state, but it has been already allocated.
      */
      flattened_state = *flattened_state_p;
      lfstate   = (int64_t*)flattened_state;
      lfstate[0] = flattened_state_size;
      lfstate[1] = 5;
    }
  } else {
    /*
      Loading a state from flatten_state.
      Allocate state.
      Loading fields from flattened_state into state.
    */
    setup = 0;
    success = alloc0(state_p,setup);
  }
  *word_pos_p = 1;
  return(success);
}
