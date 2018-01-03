#ifndef _BOLTZMANN_FLATTEN_RMATRIX_H_
#define _BOLTZMANN_FLATTEN_RMATRIX_H_ 1
extern int boltzmann_flatten_rmatrix(struct state_struct *state,
				     void *fstate, 
				     int direction, 
				     int *word_pos_p,
				     FILE *lfp);
#endif
