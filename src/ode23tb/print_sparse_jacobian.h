#ifndef _PRINT_SPARSE_JACOBIAN_H_
#define _PRINT_SPARSE_JACOBIAN_H_ 1
extern void print_sparse_jacobian(int ny, double *dfdy_a, 
				  int *dfdy_ia, int *dfdy_ja, 
				  FILE *fp);
#endif
