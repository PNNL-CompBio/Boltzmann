#include "boltzmann_structs.h"
#include "vec_set_constant.h"
#include "boltzmann_sparse_to_dense.h"
void boltzmann_sparse_to_dense(int ny,
			       double *dfdy_a,
			       int *dfdy_ia,
			       int *dfdy_ja,
			       double *dfdy) {
  /*
    Convert a sparse matrix in crs format to a dense column-major ordered 
    matrix.
    Called by: approximate_jacobian.
  */
  double zero;
  double *colj;
  double *rowi;
  int ny2;
  int i;
  int j;
  int k;
  ny2 = ny * ny;
  zero = 0.0;
  vec_set_constant(ny2, dfdy, zero);
  rowi = dfdy;
  for (i=0;i<ny;i++) {
    for (k=dfdy_ia[i];k<dfdy_ia[i+1];k++) {
      j = dfdy_ja[k];
      colj = rowi + (j * ny); /* caution address arithmetic here */
      *colj = dfdy_a[k];
    }
    rowi += 1; /* Caution address arithmetic here. */
  }
}
