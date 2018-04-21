#include "boltzmann_structs.h"
#include "boltzmann_sparse_mvp.h"
int  boltzmann_sparse_mvp(int n, double *a, int *ia, int *ja, double *v, 
			  double *av) {
  /*
    Form matirx vector product av = A*v where A is stored in compressed
    row sparse matrix format represented by a,ia,ja.

    Called by: boltzmann_cvodes_jtimes
    Calls:
  */
  double avi;
  int i;
  int j;
  int k;
  int success;
  success = 1;
  for (i=0;i<n;i++) {
    avi = 0.0;
    for (k=ia[i];k<ia[i+1];k++) {
      j = ja[k];
      avi += a[k] * v[j];
    }
    av[i] = avi;
  }
  return(success);
}
