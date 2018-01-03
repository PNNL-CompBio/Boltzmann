#ifndef _BOLTZMANN_FLATTEN_SCALARS_H_
#define _BOLTZMANN_FLATTEN_SCALARS_H_ 1
extern int boltzmann_flatten_scalars(struct state_struct *state,
				     void *flattened, 
				     int direction, 
				     int *word_pos_p, 
				     FILE *lfp);
#endif
