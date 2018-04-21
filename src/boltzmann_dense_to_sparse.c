#include "boltzmann_structs.h"
#include "boltzmann_dense_to_sparse.h"
void boltzmann_dense_to_sparse(int ny,
			       double *dfdy, double *dfdy_a,
			       int *dfdy_ia, int *dfdy_ja) {
  /*
    Convert a dense matrix in column-major order to a sparse matrix in 
    crs format.
    Called by: approximate_jacobian.
  */
  double *rowi;
  double *colj;
  int i;
  int j;
  int k;
  int padi;
  k = 0;
  rowi = dfdy;
  for (i=0;i<ny;i++) {
    colj = rowi;
    dfdy_ia[i] = k;
    for (j=0;j<ny;j++) {
      if (i == j) {
	/*
	  Always have a diagonal element even if its 0.
	*/
	dfdy_a[k] = *colj;
	dfdy_ja[k] = j;
	k+= 1;
      } else {
	if (*colj != 0.0) {
	  dfdy_a[k] = *colj;
	  dfdy_ja[k] = j;
	  k += 1;
	}
      }
      colj += ny; /* caution address arithmetic. */
    } /* end for (j...) */
    rowi += 1; /* caution address arithmetic. */
  } /* end for (i...) */
  dfdy_ia[ny] = k;
}



