#include "system_includes.h"
/*
  dlaswp.c addapted from lapack dlaswp.f
*/
#include "dlaswp.h"
void dlaswp_( int *n_p,
	     double *a,
	     int *lda_p,
	     int *k1_p,
	     int *k2_p,
	     int *ipiv,
	     int *incx_p) {
  /*
    DLASWP performs a series of row interchanges on the matrix A.
    One row interchange is initiated for each of rows K1 through K2 of A.
     A is double array, dimension (LDA,N)
       On entry, the matrix of column dimension N to which the row
       interchanges will be applied.
       On exit, the permuted matrix.
     LDA is the leading dimension of the array A.
     K1 is the first element of IPIV for which a row interchange will be done.
     K2 is the last element of IPIV for which a row interchange will
     IPIV is integer array, dimension (K2*abs(INCX))
        The vector of pivot indices.  Only the elements in positions
        K1 through K2 of IPIV are accessed.
        IPIV(K) = L implies rows K and L are to be interchanged.
     INCX is the increment between successive values of IPIV.  If INCX
         is negative, the pivots are applied in reverse order.
  */
  double temp;
  int64_t lda32;
  int64_t jlda;
  int64_t ai;
  int64_t aip;
  int n;
  int lda;

  int k1;
  int k2;

  int incx;
  int ix0;

  int ix;
  int ip;

  int i1;
  int i2;

  int inc;
  int i;

  int j;
  int k;

  int n32;
  int padi;

  n    = *n_p;
  lda  = *lda_p;
  k1   = *k1_p;
  k2   = *k2_p;
  incx = *incx_p;
  /*
    Interchange row I with row IPIV(I) for each of rows K1 through K2.
  */
  if( incx > 0 ) {
    ix0 = k1;
    i1 = k1;
    i2 = k2;
    inc = 1;
  } else {
    if (incx < 0) {
      /*
      ix0 = 1 + ( 1-k2 )*incx
	The above line is is what is in the dlaswap.f, but if what you wanted
	was to do the pivots IPIV(k1), IPIV(k1+inc), IPIV(k1+2*inc), ..
	IPIV(k1 + (k2-k1)*inc)) in reverse order the formula below
	is correct.
      */
      ix0 = k1 + (k1-k2)*incx;
      i1 = k2;
      i2 = k1;
      inc = -1;
    }
  }
  if (incx != 0) {
    /*
    n32 = ( n / 32 )*32;
    */
    n32 = n >> 5;
    n32 = n32 << 5;
    lda32 = lda << 5;
    jlda  = 0;
    if (n32 > 0) {
      for (j=0;j<n32;j+=32) {
	ix = ix0;
	for (i=i1;i != (i2 + inc); i += inc) {
	  ip = ipiv[ix] - 1;
	  if (ip != i)  {
	    /*
	      a[i,j] is at a[i + j*lda]
	    */
	    ai = i + jlda;
	    aip = ip + jlda;
	    for (k=j;k<j+32;k++) {
	      temp = a[ai];
	      a[ai] = a[aip];
	      a[aip] = temp;
	      ai += lda;
	      aip += lda;
	    } /* end for (k...) */
	  }  /* end if (ip != i) */
	  ix += incx;
	} /* end for (i...) */
	jlda =jlda + lda32;
      } /* end for (j...) */
    } /* end if n32 > 0 */
    /*
      Clean up loop.
    */
    if (n32 != n) {
      ix = ix0;
      jlda = n32 * lda;
      for (i=i1;i != (i2+inc);i += inc) {
	ip = ipiv[ix] - 1;
	if (ip != i) {
	  ai = i + jlda;
	  aip = ip + jlda;
	  for (k=n32;k<n;k++) {
	    temp = a[ai];
	    a[ai] = a[aip];
	    a[aip] = temp;
	    ai += lda;
	    aip += lda;
	  } /* end for (k ...) */
	} /* end if (ip != i) */
	ix += incx;
      } /* end for (i...) */
    } /* end if (n32 != n) */
  } /* end if (incx != 0) */
  return;
} /* end of dlaswp */
