#ifndef _FLATTEN_ONEWAY_VECS_H_
#define _FLATTEN_ONEWAY_VECS_H_ 1
extern int boltzmann_flatten_oneway_vecs(struct state_struct *state,
					 void *fstate, 
					 int direction, 
					 int *word_pos_p,
					 FILE *lfp);
#endif
