#include "boltzmann_structs.h"
#include "boltzmann_flatten_vgrng_state.h"
int boltzmann_flatten_vgrng_state(int64_t *fvgrng_state,
				  struct vgrng_state_struct *vgrng_state,
				  int direction,
				  int *vgrng_state_size) {
  /*
    Flatten a random number generator state (vgrng_state).
    Called by: boltzmann_flatten_twoway_data
  */
  double *uni_multiplier_p;
  int *ltoi;
  int success;
  int padi;

  success = 1;

  uni_multiplier_p = (double*)&fvgrng_state[11];
  ltoi             = (int*)&fvgrng_state[12];
  if (direction == 0) {
    fvgrng_state[0] = vgrng_state->fib_history[0];
  } else {
    vgrng_state->fib_history[0] = fvgrng_state[0];
  }
  if (direction == 0) {
    fvgrng_state[1] = vgrng_state->fib_history[1];
  } else {
    vgrng_state->fib_history[1] = fvgrng_state[1];
  }
  if (direction == 0) {  
    fvgrng_state[2] = vgrng_state->lcg_history;
  } else {
    vgrng_state->lcg_history = fvgrng_state[2];
  }
  if (direction == 0) {  
    fvgrng_state[3] = vgrng_state->fib_constant;
  } else {
    vgrng_state->fib_constant = fvgrng_state[3];
  }
  if (direction == 0) {  
    fvgrng_state[4] = vgrng_state->lcg_constant;
  } else {
    vgrng_state->lcg_constant = fvgrng_state[4];
  }
  if (direction == 0) {  
    fvgrng_state[5] = vgrng_state->lcg_multiplier;
  } else {
    vgrng_state->lcg_multiplier = fvgrng_state[5];
  }
  if (direction == 0) {  
    fvgrng_state[6] = vgrng_state->mask;
  } else {
    vgrng_state->mask = fvgrng_state[6];
  }
  if (direction == 0) {  
    fvgrng_state[7] = vgrng_state->fib_seed[0];
  } else {
    vgrng_state->fib_seed[0] = fvgrng_state[7];
  }
  if (direction == 0) {  
    fvgrng_state[8] = vgrng_state->fib_seed[1];
  } else {
    vgrng_state->fib_seed[1] = fvgrng_state[8];
  }
  if (direction == 0) {  
    fvgrng_state[9] = vgrng_state->lcg_seed;
  } else {
    vgrng_state->lcg_seed = fvgrng_state[9];
  }
  if (direction == 0) {  
    fvgrng_state[10] = vgrng_state->xor_lcg_fib;
  } else {
    vgrng_state->xor_lcg_fib = fvgrng_state[10];
  }
  if (direction == 0) {  
    *uni_multiplier_p = vgrng_state->uni_multiplier;
  } else {
    vgrng_state->uni_multiplier = *uni_multiplier_p;
  }
  if (direction == 0) {
    ltoi[0] = vgrng_state->fib_cur_ptr;
  } else {
    vgrng_state->fib_cur_ptr = ltoi[0];
  }
  if (direction == 0) {
    /*
      pad the flattened vgrng_state to 16 elements.
    */
    fvgrng_state[13] = (int64_t)0;
    fvgrng_state[14] = (int64_t)0;
    fvgrng_state[15] = (int64_t)0;
  }
  *vgrng_state_size = 16;
  return(success);
}
