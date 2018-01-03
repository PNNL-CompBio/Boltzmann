#include "boltzmann_structs.h"
#include "boltzmann_flatten_vgrng_state.h"
#include "boltzmann_flatten_twoway_data.h"
int boltzmann_flatten_twoway_data(struct state_struct *state,
				 void *fstate, int direction,
				 int *word_pos_p, FILE *lfp) {
  /*
    Flatten the incoming only structs and vectors.
    Called by: boltzmann_flatten_state
    Calls:     boltzmann_flatten_vgrng_state,
               memcpy, fprintf, fflush
  */
  struct cvodes_params_struct cps;
  struct cvodes_params_struct *cvodes_params;
  int64_t *fvgrng_state;
  int64_t *lfstate;
  double  *dfstate;
  double  *current_counts;
  double  *bndry_flux_counts;
  double  *net_lklhd_bndry_flux;
  double  *net_likelihood;
  int nunique_molecules;
  int word_pos;

  int vgrng_state_size;
  int vec_len;

  int two_way_len;
  int success;

  int number_reactions;
  int cvodes_params_size;

  int cvodes_words;
  int padi;

  success = 1;
  nunique_molecules = state->nunique_molecules;
  number_reactions  = state->number_reactions;
  lfstate = (int64_t *)fstate;
  dfstate = (double  *)fstate;
  word_pos = *word_pos_p;
  vgrng_state_size = 16;
  cvodes_params_size = sizeof(cps);
  two_way_len = 2 + 3*nunique_molecules + number_reactions + 
    2*vgrng_state_size + cvodes_params_size;
  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = two_way_len;
  } else {
    if (lfstate[word_pos] != two_way_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_twoway_data: Error two_way data lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],two_way_len);
      }
      
    }
  }
  if (success) {
    word_pos += 1;
    lfstate[word_pos] = -6;
    
    word_pos += 1;

    fvgrng_state = (int64_t *)&lfstate[word_pos];
    boltzmann_flatten_vgrng_state(fvgrng_state,state->vgrng_state,direction,&vgrng_state_size);
    word_pos += vgrng_state_size;

    fvgrng_state = (int64_t *)&lfstate[word_pos];
    boltzmann_flatten_vgrng_state(fvgrng_state,state->vgrng2_state,direction,&vgrng_state_size);
    word_pos += vgrng_state_size;

    cvodes_params   = (void *)&lfstate[word_pos];
    cvodes_words    = (cvodes_params_size + 7) >> 3;
    if (direction == 0) {
      memcpy(cvodes_params, (void*)state->cvodes_params,cvodes_params_size);
    } else {
      memcpy((void*)state->cvodes_params,cvodes_params,cvodes_params_size);
    }
    word_pos += cvodes_words;
    current_counts = (double*)&dfstate[word_pos];
    if (direction == 0) {
      /*
	lfstate[8] is currently hardwired to be the position that
	points the the first word in current_counts when stored
	in the flattened_state array.
      */
      lfstate[8] = word_pos;
      vec_len     = nunique_molecules * sizeof(double);
      state->current_counts_offset = word_pos;
      memcpy(current_counts,state->current_counts,vec_len);
    } else {
      if (lfstate[8] != word_pos) {
	success = 0;
	if (lfp) {
	  fprintf(lfp,"boltzmann_flatten_twoway_data: Error current_counts_offset, = %ld,  does not match word_pos = %d\n",
		  lfstate[8],word_pos);
	  fflush(lfp);
	}
      } else {
	vec_len     = nunique_molecules * sizeof(double);
	memcpy(state->current_counts,current_counts,vec_len);
      }
    }
  }
  if (success) {
    word_pos += nunique_molecules;
    bndry_flux_counts = (double*)&dfstate[word_pos];
    vec_len     = nunique_molecules * sizeof(double);
    if (direction == 0) {
      memcpy(bndry_flux_counts,state->bndry_flux_counts,vec_len);
    } else {
      memcpy(state->bndry_flux_counts,bndry_flux_counts,vec_len);
    }
    word_pos += nunique_molecules;
    net_lklhd_bndry_flux = (double*)&dfstate[word_pos];
    vec_len     = nunique_molecules * sizeof(double);
    if (direction == 0) {
      memcpy(net_lklhd_bndry_flux,state->net_lklhd_bndry_flux,vec_len);
    } else {
      memcpy(state->net_lklhd_bndry_flux,net_lklhd_bndry_flux,vec_len);
    }
    word_pos += nunique_molecules;
    net_likelihood = (double*)&dfstate[word_pos];
    vec_len     = number_reactions * sizeof(double);
    if (direction == 0) {
      memcpy(net_likelihood,state->net_likelihood,vec_len);
    } else {
      memcpy(state->net_likelihood,net_likelihood,vec_len);
    }
    word_pos += number_reactions;


  } /* end if (success) */
  /*
    Set word_pos to be the last word set.
  */
  word_pos = word_pos - 1;
  *word_pos_p = word_pos;
  return(success);
}
