#include "boltzmann_structs.h"
#include "print_dense_jacobian.h"
void print_sparse_jacobian(int ny, double *dfdy_a, int *dfdy_ia, int *dfdy_ja, 
                           FILE *fp) {
  /*
    Print the entries of a full jacobian matrix.
    The jacobian is expected to be in column major order,
    but we want to print in row major order.
  */
  int i;
  int j;
  if (fp) {
    for (i=0;i<ny;i++) {
      fprintf(fp,"row %d\t ja |",i);
      for (j=dfdy_ia[i];j<dfdy_ia[i+1];j++) {
	fprintf(fp,"\t%d",dfdy_ja[j]);
      }
      fprintf(fp,"\n");
      fprintf(fp,"row %d\t  a |",i);
      for (j=dfdy_ia[i];j<dfdy_ia[i+1];j++) {
	fprintf(fp,"\t%le",dfdy_a[j]);
      }
      fprintf(fp,"\n");
    }
  }
}
