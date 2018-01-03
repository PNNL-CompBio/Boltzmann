#include "boltzmann_structs.h"
#include "boltzmann_flatten_reactions.h"
#include "boltzmann_flatten_molecules.h"
#include "boltzmann_flatten_compartments.h"
#include "boltzmann_flatten_rmatrix.h"
#include "boltzmann_flatten_mmatrix.h"
#include "boltzmann_flatten_oneway_vecs.h"

#include "boltzmann_flatten_oneway_data.h"
int boltzmann_flatten_oneway_data(struct state_struct *state,
				 void *fstate, int direction,
				 int *word_pos_p, FILE *lfp) {
  /*
    Flatten the incoming only structs and vectors.
    Called by: boltzmann_flatten_state
    Calls:     boltzmann_flatten_reactions,
               boltzmann_flatten_molecules,
               boltzmann_flatten_compartments,
               boltzmann_flatten_rmatrix,
               boltzmann_flatten_mmatrix,
               boltzmann_flatten_oneway_vecs
  */
  int64_t *lfstate;
  double  *dfstate;
  int word_pos;
  int oneway_len;

  int oneway_len_pos;
  int reactions_len_pos;

  int molecules_len_pos;
  int compartments_len_pos;

  int rmatrix_len_pos;
  int mmatrix_len_pos;

  int oneway_vec_len_pos;
  int success;

  lfstate = (int64_t *)fstate;
  dfstate = (double  *)fstate;
  word_pos = *word_pos_p;

  word_pos += 1;
  oneway_len_pos = word_pos;
  
  oneway_len  = 2; /* meta data - total length */
  /*
    6 component: reactions, molecules, compartments, reactions_matrix,
    molecules_matrix, and vectors:
  */
  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = 6;
  }
  reactions_len_pos = word_pos + 1;
  success = boltzmann_flatten_reactions(state,fstate,
					direction,&word_pos,lfp);

  if (success) {
    oneway_len += lfstate[reactions_len_pos];
    molecules_len_pos = word_pos + 1;
    success = boltzmann_flatten_molecules(state,fstate,direction,&word_pos,
					  lfp);
  }
  if (success) {
    oneway_len += lfstate[molecules_len_pos];
    compartments_len_pos = word_pos + 1;
    success = boltzmann_flatten_compartments(state,fstate,
					     direction,&word_pos,lfp);
  }
  if (success) {
    oneway_len += lfstate[compartments_len_pos];
    rmatrix_len_pos = word_pos + 1;
    success = boltzmann_flatten_rmatrix(state,fstate,direction,
					&word_pos,lfp);
  }
  if (success) {
    oneway_len += lfstate[rmatrix_len_pos];
    mmatrix_len_pos = word_pos + 1;
    success = boltzmann_flatten_mmatrix(state,fstate,direction,
					&word_pos,lfp);
  }
  if (success) {
    oneway_len += lfstate[mmatrix_len_pos];
    oneway_vec_len_pos = word_pos + 1;
    success = boltzmann_flatten_oneway_vecs(state,fstate,direction,
					    &word_pos,lfp);
  }
  if (success) {
    oneway_len += lfstate[oneway_vec_len_pos];
    if (direction == 0) {
      lfstate[oneway_len_pos] = oneway_len;
    }
  }
  *word_pos_p = word_pos;
  return(success);
}

