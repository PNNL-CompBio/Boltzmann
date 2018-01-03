#ifndef _FILL_FLUX_PIECES_H_
#define _FILL_FLUX_PIECES_H_ 1
extern int fill_flux_pieces(struct state_struct *state, 
			    struct molecules_matrix_struct *molecules_matrix,
			    int base_reaction,
			    int fill_jacobi);
#endif
