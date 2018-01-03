#ifndef _BOLTZMANN_FLATTEN_COMPARTMENTS_H_
#define _BOLTZMANN_FLATTEN_COMPARTMENTS_H_ 1
extern int boltzmann_flatten_compartments(struct state_struct *state,
					  void *fstate, 
					  int direction, 
					  int *word_pos_p,
					  FILE *lfp);
#endif
