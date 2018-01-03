#ifndef _BOLTZMANN_DENSE_TO_SPARSE_H_
#define _BOLTZMANN_DENSE_TO_SPARSE_H_ 1
extern void boltzmann_dense_to_sparse(int ny,
				      double *dfdy, double *dfdy_a,
				      int *dfdy_ia, int *dfdy_ja);
#endif
