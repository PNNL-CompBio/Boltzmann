#ifndef _BOLTZMANN_FLATTEN_EXTERNAL_H_
#define _BOLTZMANN_FLATTEN_EXTERNAL_H_ 1
extern int boltzmann_flatten_external(struct state_struct *state,
				      void *flattened, 
				      int direction, 
				      int *word_pos_p, 
				      FILE *lfp);
#endif
