#include "system_includes.h"
#include "blas.h"
#include "dgemm.h"
#include "dlaswp.h"
#include "dtrsm.h"
#include "dgetrf2.h"
void dgetrf2(int *m_p, 
	     int *n_p, 
	     double *a, 
	     int *lda_p, 
	     int *ipiv, 
	     int *info_p) {
  /*
    Adapted from dgetrf2.f from lapack source.
  */
  /*
   
    DGETRF2 computes an LU factorization of a general M-by-N matrix A
    using partial pivoting with row interchanges.
   
    The factorization has the form
       A = P * L * U
    where P is a permutation matrix, L is lower triangular with unit
    diagonal elements (lower trapezoidal if m > n), and U is upper
    triangular (upper trapezoidal if m < n).
   
    This is the recursive version of the algorithm. It divides
    the matrix into four submatrices:
               
           [  A11 | A12  ]  where A11 is n1 by n1 and A22 is n2 by n2
       A = [ -----|----- ]  with n1 = min(m,n)/2
           [  A21 | A22  ]       n2 = n-n1
               
                                          [ A11 ]
    The subroutine calls itself to factor [ --- ],
                                          [ A12 ]
                    [ A12 ]
    do the swaps on [ --- ], solve A12, update A22,
                    [ A22 ]
   
    then calls itself to factor A22 and do the swaps on A21.
  */
  double one;
  double mone;
  double zero;
  double temp;
  double sfmin;
  double *dbl_ptr;
  int64_t sfmin_hex;
  int64_t a12_start;
  int64_t a21_start;
  int64_t a22_start;
  int m;
  int n;

  int lda;
  int info;

  int iinfo;
  int i;

  int n1;
  int n2;

  int inc1;
  int mm1;

  int izero;
  int n1m1;

  int mmn1;
  int min_m_n;

  int min_m_n_m1;
  int padi;

  char l_char;
  char n_char;
  char u_char;
  char padc;
  m   = *m_p;
  n   = *n_p;
  lda = *lda_p;
  mm1 = m-1;
  sfmin_hex       = 0x0004000000000001L;
  dbl_ptr = (double*)&sfmin_hex;
  sfmin   = *dbl_ptr;
  one     = 1.0;
  zero    = 0.0;
  inc1    = 1;
  izero   = 0;
  l_char  = 'L';
  n_char  = 'N';
  u_char  = 'U';
  info = 0;
  if (m < 0) {
    info = -1;
  } else {
    if (n < 0) {
      info = -2;
    } else {
      if ((lda < 1) || (lda < m)) {
	info = -4;
      }
    }
  }
  if (info != 0) {
    fprintf(stderr,"dgetrf2: error: info = %d\n",info);
    fflush(stderr);
  } else {
    if ((m > 0) && (n > 0)) {
      if (m == 1) {
	/*
	  Use unblocked code for one row case
	  Just need to handle IPIV and INFO
	*/
	ipiv[0] = 0;
	if (a[0] == 0) {
	  info = 1;
	} 
      } else {
	if (n == 1) {
	  /*
	    Use unblocked code for one column case
	  */
	  i = idamax(&m,a,&inc1);
	  ipiv[0] = i;
	  if (a[i] != zero) {
	    /*
	      apply the interchange.
	    */
	    if (i != 0) {
	      temp = a[0];
	      a[0] = a[i];
	      a[i] = temp;
	    }
	    /*
	      Compute elements 1:m-1 of the column
	    */
	    temp = fabs(a[0]);
	    if (temp > sfmin) {
	      temp = one/temp;
	      dscal(&mm1,&temp,&a[1],&inc1);
	    } else {
	      for (i=1;i<m;i++) {
		a[i] = a[i]/temp;
	      }
	    }
	  } else {
	    /*
	      Zero pivot.
	    */
	    info = 1;
	  }
	} else {
	  /*
	    m > 1 and n > 1
	    Use recursive mode.
	  */
	  min_m_n = m;
	  if (n < m) {
	    min_m_n = n;
	  }
	  n1 = min_m_n >> 1;
	  n2 = n - n1;
	  /*
	    Recursive call.
	    Factor [ A11 ]
	           [ --- ]
	           [ A21 ]
	  */
	  dgetrf2(&m,&n1,a,&lda,ipiv,&iinfo);
	  if ((info == 0) && (iinfo > 0)) {
	    info = iinfo;
	  }
	  /*
  	                        [ A12 ]
          Apply interchanges to [ --- ]
                                [ A22 ]
	  */
	  n1m1 = n1 - 1;
	  a12_start = n1 * lda;
	  a22_start = a12_start + n1;
	  a21_start = n1;
	  dlaswp(&n2,&a[a12_start], &lda, &izero, &n1m1, ipiv, &inc1 );
	  /*
	    Solve A12
	  */
	  dtrsm(&l_char,&l_char, &n_char, &u_char,&n1,&n2,&one,a,&lda,
		&a[a12_start],&lda);
	  /*
	    Update A22
	  */
	  mmn1 = m - n1;
	  mone = -1.0;
	  dgemm(&n_char, &n_char, &mmn1, &n2, &n1, &mone, &a[a21_start],
		&lda, &a[a12_start], &lda, &one, &a[a22_start], &lda);
	  /*
	    Factor A22
	    Recursive Call.
	  */
	  dgetrf2(&mmn1, &n2, &a[a22_start], &lda, &ipiv[n1],&iinfo);
	  /*
	    Adjust INFO and the pivot indices
	  */
	  if ((info == 0) && (iinfo > 0)) {
	    info = iinfo + n1;
	  }
	  for (i=n1;i<min_m_n;i++) {
	    ipiv[i] = ipiv[i] + n1;
	  }
	  /*
	    Apply interchanges to A21
	  */
	  min_m_n_m1 = min_m_n - 1;
	  dlaswp(&n1, a, &lda, &n1, &min_m_n_m1, ipiv,&inc1);
	} /* end else (n > 1 and m > 1) */
      } /* end else m> 1 */
    } /* end if (m > 0 & n > 0) */
  } /* end else (info == 0) */
  *info_p = info;
}
