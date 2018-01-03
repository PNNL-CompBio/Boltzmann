#include "boltzmann_structs.h"
#include "boltzmann_flatten_mmatrix.h"
int boltzmann_flatten_mmatrix(struct state_struct *state,
			      void *fstate, 
			      int direction, 
			      int *word_pos_p,
			      FILE *lfp) {
  /*
    Flatten the molecules_matrix struct
    Called by: boltzmann_flatten_oneway_data
    Calls:     fprintf, fflush
  */
  struct molecules_matrix_struct *mmatrix;
  int64_t *lfstate;
  int64_t nunique_molecules;
  int64_t molecules_ptrs_len;
  int64_t nz_len;
  void *molecules_ptrs;
  void *reaction_indices;
  void *coefficients;
  void *recip_coeffs;

  int success;
  int word_pos;

  int mmatrix_len;
  int nz;

  success = 1;
  nunique_molecules = state->nunique_molecules;
  mmatrix          = state->molecules_matrix;
  nz               = state->number_molecules;

  word_pos = *word_pos_p;

  lfstate = (int64_t *)fstate;

  mmatrix_len = (int64_t)3 + nunique_molecules + (3*nz);

  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = mmatrix_len;
  } else {
    if (lfstate[word_pos] != mmatrix_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_mmatrix: Error mmatrix struct lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],mmatrix_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    }
    molecules_ptrs_len = (nunique_molecules + 1) * sizeof(int64_t);
    nz_len       = nz * sizeof(double);
    
    word_pos += 1;
    molecules_ptrs = (void*)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(molecules_ptrs,mmatrix->molecules_ptrs,molecules_ptrs_len);
    } else {
      memcpy(mmatrix->molecules_ptrs,molecules_ptrs,molecules_ptrs_len);
    }
    word_pos += nunique_molecules + 1;
    reaction_indices = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(reaction_indices,mmatrix->reaction_indices,nz_len);
    } else {
      memcpy(mmatrix->reaction_indices,reaction_indices,nz_len);
    }
    word_pos += nz;
    coefficients = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(coefficients,mmatrix->coefficients,nz_len);
    } else {
      memcpy(mmatrix->coefficients,coefficients,nz_len);
    }
    word_pos += nz;
    recip_coeffs = (void*)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(recip_coeffs,mmatrix->recip_coeffs,nz_len);
    } else {
      memcpy(mmatrix->recip_coeffs,recip_coeffs,nz_len);
    }
    /*
      NB we need to leave word_pos pointing at the last word set.
    */
    word_pos += (nz - 1);
  } /* end if (success) */
  *word_pos_p = word_pos;
  return(success);
}
