#include "boltzmann_structs.h"
#include "boltzmann_flatten_rmatrix.h"
int boltzmann_flatten_rmatrix(struct state_struct *state,
			      void *fstate, 
			      int direction, 
			      int *word_pos_p,
			      FILE *lfp) {
  /*
    Flatten the reactions_matrix struct
    Called by: boltzmann_flatten_oneway_data
    Calls:     fprintf, fflush

    We don't need to save the text field of the reaction matrix as
    that is only used in setup.
  */
  struct reactions_matrix_struct *rmatrix;
  double  *dfstate;
  int64_t *lfstate;
  int64_t  number_reactions;
  int64_t rxn_ptrs_len;
  int64_t solventc_len;
  int64_t nz_len;
  void *rxn_ptrs;
  void *molecules_indices;
  void *compartment_indices;
  void *coefficients;
  void *recip_coeffs;  
  void *solvent_coefficients;

  int success;
  int word_pos;

  int rmatrix_len;
  int nz;

  success = 1;
  number_reactions = state->number_reactions;
  rmatrix          = state->reactions_matrix;
  nz               = state->number_molecules;

  word_pos = *word_pos_p;

  lfstate = (int64_t *)fstate;
  dfstate = (double  *)fstate;

  rmatrix_len = (int64_t)3 + (3*number_reactions)  + (4*nz);

  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = rmatrix_len;
  } else {
    if (lfstate[word_pos] != rmatrix_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_rmatrix: Error rmatrix struct lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],rmatrix_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    }
    rxn_ptrs_len = (number_reactions + 1) * sizeof(int64_t);
    nz_len       = nz * sizeof(double);
    solventc_len = 2* number_reactions *sizeof(int64_t);
    
    word_pos += 1;
    rxn_ptrs = (void*)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(rxn_ptrs,rmatrix->rxn_ptrs,rxn_ptrs_len);
    } else {
      memcpy(rmatrix->rxn_ptrs,rxn_ptrs,rxn_ptrs_len);
    }
    word_pos += number_reactions + 1;
    molecules_indices = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(molecules_indices,rmatrix->molecules_indices,nz_len);
    } else {
      memcpy(rmatrix->molecules_indices,molecules_indices,nz_len);
    }
    word_pos += nz;
    compartment_indices = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(compartment_indices,rmatrix->compartment_indices,nz_len);
    } else {
      memcpy(rmatrix->compartment_indices,compartment_indices,nz_len);
    }
    word_pos += nz;
    coefficients = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(coefficients,rmatrix->coefficients,nz_len);
    } else {
      memcpy(rmatrix->coefficients,coefficients,nz_len);
    }
    word_pos += nz;
    recip_coeffs = (void *)&dfstate[word_pos];
    if (direction == 0) {
      memcpy(recip_coeffs,rmatrix->recip_coeffs,nz_len);
    } else {
      memcpy(rmatrix->recip_coeffs,recip_coeffs,nz_len);
    }
    word_pos += nz;
    solvent_coefficients = (void *)&lfstate[word_pos];
    if (direction == 0) {
      memcpy(solvent_coefficients,rmatrix->solvent_coefficients,solventc_len);
    } else {
      memcpy(rmatrix->solvent_coefficients,solvent_coefficients,solventc_len);
    }
    /*
      NB we need to leave word_pos pointing at the last word set.
    */
    word_pos += (number_reactions + number_reactions - 1);
  } /* end if (success) */
  *word_pos_p = word_pos;
  return(success);
}
