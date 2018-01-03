#ifndef _BOLTZMANN_FLATTEN_ALLOC_H_
#define _BOLTZMANN_FLATTEN_ALLOC_H_ 1
extern int boltzmann_flatten_alloc (struct state_struct **state_p,
				    void **flattened_state_p,
				    int direction, 
				    int *word_pos_p, 
				    FILE *lfp);
#endif
