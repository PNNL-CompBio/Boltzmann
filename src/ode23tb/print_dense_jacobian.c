#include "boltzmann_structs.h"
#include "print_dense_jacobian.h"
void print_dense_jacobian(int ny, double *dfdy, FILE *fp) {
  /*
    Print the entries of a full jacobian matrix.
    The jacobian is expected to be in column major order,
    but we want to print in row major order.
  */
  double *di;
  double *dj;
  int i;
  int j;
  if (fp) {
    fprintf(fp,"Dense Jacobian: ny = %d\n",ny);
    di = dfdy;
    for (i=0;i<ny;i++) {
      dj = di;
      fprintf(fp,"%d\t|",i);
      for (j=0;j<ny;j++) {
	fprintf(fp,"\t%le",*dj);
	dj += ny; /* Caution address arithmetic */
      }
      fprintf(fp,"\n");
      di += 1; /* Caution address arithmetic */
    }
  }
}
