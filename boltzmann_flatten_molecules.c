#include "boltzmann_structs.h"
#include "boltzmann_flatten_molecules.h"
int boltzmann_flatten_molecules(struct state_struct *state,
				void *fstate, 
				int direction, 
				int *word_pos_p,
				FILE *lfp) {
  /*
    Flatten the molecule_structs.
    Called by: boltzmann_flatten_oneway_data
    Calls:     fprintf, fflush
  */
  struct molecule_struct *molecule;
  int64_t *lfstate;
  int64_t  nunique_molecules;
  int *ltoi;
  int success;
  int word_pos;

  int molecules_len;
  int words_per_molecule;

  int i;
  int padi;
  
  success = 1;
  nunique_molecules = state->nunique_molecules;
  molecule = state->sorted_molecules;;
  word_pos = *word_pos_p;

  lfstate = (int64_t *)fstate;

  words_per_molecule = 4;
  molecules_len = (int64_t)2 + (nunique_molecules * words_per_molecule);

  word_pos += 1;
  if (direction == 0) {
    lfstate[word_pos] = molecules_len;
  } else {
    if (lfstate[word_pos] != molecules_len) {
      success = 0;
      if (lfp) {
	fprintf(lfp,"boltzmann_flatten_molecules: Error molecules struct lengths do not match, input value was %ld, expected value was %d\n",
		lfstate[word_pos],molecules_len);
      }
    }
  }
  if (success) {
    word_pos += 1;
    if (direction == 0) {
      lfstate[word_pos] = -5; /* packed structure. */
    } 
    for (i=0;i<nunique_molecules;i++) {
      word_pos += 1;
      if (direction == 0) {
	lfstate[word_pos] = molecule->string;
      } else {
	molecule->string = lfstate[word_pos];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
	ltoi[0] = molecule->m_index;
	ltoi[1] = molecule->c_index;
      } else {
	molecule->m_index = ltoi[0];
	molecule->c_index = ltoi[1];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
	ltoi[0] = molecule->variable;
	ltoi[1] = molecule->g_index;
      } else {
	molecule->variable = ltoi[0];
	molecule->g_index  = ltoi[1];
      }
      word_pos += 1;
      ltoi = (int*)&lfstate[word_pos];
      if (direction == 0) {
	ltoi[0] = molecule->compute_init_conc;
	ltoi[1] = molecule->solvent;
      } else {
	molecule->compute_init_conc = ltoi[0];
	molecule->solvent           = ltoi[1];
      }
      molecule += 1; /* Caution address arithmetic here */
    } /* end for(i...) */
  } /* end if(success) */
  *word_pos_p = word_pos;
  return(success);
}
