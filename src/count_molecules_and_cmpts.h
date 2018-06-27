#ifndef _COUNT_MOLECULES_AND_CMPTS_H_
#define _COUNT_MOLECULES_AND_CMPTS_H_ 1
extern void count_molecules_and_cmpts(char* molecules_line, 
				      int *num_molecules_p, 
				      int *num_compartments_p, 
				      int64_t *molecules_len_p, 
				      int64_t *compartment_len_p,
				      FILE *lfp);
#endif
