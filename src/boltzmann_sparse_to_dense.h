#ifndef _BOLTZMANN_SPARSE_TO_DENSE_H_
#define _BOLTZMANN_SPARSE_TO_DENSE_H_ 1
extern void boltzmann_sparse_to_dense(int ny,
				      double *dfdy_a,
				      int *dfdy_ia,
				      int *dfdy_ja,
				      double *dfdy);
#endif
