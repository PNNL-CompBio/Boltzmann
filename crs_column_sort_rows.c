#include "boltzmann_structs.h"
#include "crs_column_sort_rows.h"
void crs_column_sort_rows(int m,
			  int n,
			  double *a, 
			  int *ia, 
			  int *ja, 
			  double *at, 
			  int *iat,
			  int *jat) {
  /*
    Column sort the rows of a compressed row storage format sparse matrix
    of an m-rows by n-column matrix.
    Use the double transpose method as it is linear in the number
    of nonzero lements.
    Called by: lr8_approximate_jacobian
  */
  /*
    Count the elements per column
  */
  double val;

  int i;
  int j;

  int nz;
  int tpos;

  int this_row;
  int next_row;

  int col;
  int row;

  nz = ia[m];
  for (j=0;j<n;j++) {
    iat[j] = 0;
  }
  for (j=0;j<nz;j++) {
    iat[ja[j]] += 1;
  }
  /*
    Form partial sums of iat in place so
    that iat are the row pointers for at.
  */
  this_row = iat[0];
  iat[0]   = 0;
  for (j=1;j<=n;j++ ) {
    next_row = iat[j];
    iat[j]   = iat[j-1] + this_row;
    this_row = next_row;
  }
  /*
    Now form at and jat
  */
  for (i=0;i<m;i++) {
    for (j=ia[i];j<ia[i+1];j++) {
      col = ja[j];
      val = a[j];
      tpos = iat[col];
      at[tpos] = val;
      jat[tpos] = i;
      iat[col] += 1;
    }
  }
  /*
    Now we need to shift right iat as row pointer for row i is now in iat[i-1],
    and we want it to be in iat[i]
  */
  for (i=n-1;i>0;i--) {
    iat[i] = iat[i-1];
  }
  iat[0] = 0;
  /*
    Now we transpose back
  */
  for (j=0;j<n;j++) {
    for (i=iat[j];i<iat[j+1];i++) {
      row = jat[i];
      val = at[i];
      tpos = ia[row];
      a[tpos] = val;
      ja[tpos] = j;
      ia[row] += 1;
    }
  }
  /*
    Now we need to shift right ia as row pointer for row i is now in ia[i-1],
    and we want it to be in ia[i]
  */
  for (i=n-1;i>0;i--) {
    ia[i] = ia[i-1];
  }
  ia[0] = 0;
}
  
