#include "system_includes.h"
#include "blas.h"
#include "lsame.h"

#include "dtrsm.h"
void dtrsm_(char *side_p, char *uplo_p, char *transa_p, char *diag_p,
	   int *m_p, int *n_p, double *alpha_p, double *a,
	   int *lda_p, double *b, int *ldb_p) {
  /*
    Adapted from LAPACK dtrsm.f
    DTRSM  solves one of the matrix equations
   
      op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
   
      where alpha is a scalar, X and B are m by n matrices, A is a unit, or
      non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
   
      op( A ) = A   or   op( A ) = A**T.
   
      The matrix X is overwritten on B.
      .. Scalar Arguments ..
        double alpha
        int lda,ldb,m,n
`       char diag,side,transa,uplo
      .. Array Arguments ..
       double  a(lda,*),b(ldb,*)

      Arguments:
      ==========
      SIDE is char 
        On entry, SIDE specifies whether op( A ) appears on the left
        or right of X as follows:
              SIDE = 'L' or 'l'   op( A )*X = alpha*B.

              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.

      UPLO is char
        On entry, UPLO specifies whether the matrix A is an upper or
        lower triangular matrix as follows:

              UPLO = 'U' or 'u'   A is an upper triangular matrix.
              UPLO = 'L' or 'l'   A is a lower triangular matrix.

      TRANSA is char
        On entry, TRANSA specifies the form of op( A ) to be used in
        the matrix multiplication as follows:

              TRANSA = 'N' or 'n'   op( A ) = A.
              TRANSA = 'T' or 't'   op( A ) = A**T.
              TRANSA = 'C' or 'c'   op( A ) = A**T.

      DIAG is char
        On entry, DIAG specifies whether or not A is unit triangular
        as follows:
              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
              DIAG = 'N' or 'n'   A is not assumed to be unit
                                  triangular.
      M is int
        On entry, M specifies the number of rows of B. M must be at
        least zero.

      N is int
        On entry, N specifies the number of columns of B.  N must be
        at least zero.

      alpha is double
        On entry,  ALPHA specifies the scalar  alpha. When  alpha is
        zero then  A is not referenced and  B need not be set before
        entry.

      A is double array of DIMENSION ( LDA, k ),
        where k is m when SIDE = 'L' or 'l'  
        and k is n when SIDE = 'R' or 'r'.
        Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
        upper triangular part of the array  A must contain the upper
        triangular matrix  and the strictly lower triangular part of
        A is not referenced.
        Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
        lower triangular part of the array  A must contain the lower
        triangular matrix  and the strictly upper triangular part of
        A is not referenced.
        Note that when  DIAG = 'U' or 'u',  the diagonal elements of
        A  are not referenced either,  but are assumed to be  unity.

        LDA is int
          On entry, LDA specifies the first dimension of A as declared
          in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
          LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
          then LDA must be at least max( 1, n ).

        B is double array of DIMENSION ( LDB, n ).
          Before entry,  the leading  m by n part of the array  B must
          contain  the  right-hand  side  matrix  B,  and  on exit  is
          overwritten by the solution matrix  X.

        LDB is INTEGER
          On entry, LDB specifies the first dimension of B as declared
          in  the  calling  (sub)  program.   LDB  must  be  at  least
          max( 1, m ).

        Called by dgetrf, dgetrf2, dgetrs
	Calls lsame, dscal, daxpy, ddot
  */
  double temp;
  double alpha;
  double one;
  double zero;
  double *acoli;
  double *acolj;
  double *acolk;
  double *bcolj;
  double *bcolk;
  int64_t acoli_pos;
  int64_t acolj_pos;
  int64_t acolk_pos;
  int64_t bcolj_pos;
  int64_t bcolk_pos;
  int64_t mm1_lda;
  int  m;
  int  n;

  int  lda;
  int  ldb;

  int  i;
  int  j;

  int  k;
  int  info;

  int  nrowa;
  int  lside;

  int  rside;
  int  nounit;

  int  upper;
  int  lower;

  int  trans;
  int  conjt;

  int  notrans;
  int  udiag;

  int  notudiag;
  int  inc1;

  int  min_lda;
  int  min_ldb;

  int  mmkm1;
  int  mmim1;

  char side;
  char uplo;
  char transa;
  char diag;
  char c_char;
  char l_char;
  char n_char;
  char r_char;
  char t_char;
  char u_char;
  char padc1;
  char padc2;
  char padc3;
  char padc4;
  char padc5;
  char padc6;


  alpha  = *alpha_p;
  m      = *m_p;
  n      = *n_p;
  lda    = *lda_p;
  ldb    = *ldb_p;
  side   = *side_p;
  uplo   = *uplo_p;
  transa = *transa_p;
  diag   = *diag_p;
  one    = 1.0;
  zero   = 0.0;
  inc1   = 1;
  c_char = 'C';
  l_char = 'L';
  n_char = 'N';
  r_char = 'R';
  t_char = 'T';
  u_char = 'U';
  /*
    Test the input parameters.
  */
  lside = lsame_(&side,&l_char);
  if (lside) {
    nrowa = m;
  } else {
    nrowa = n;
  }
  nounit   = lsame_(&diag,&n_char);
  upper    = lsame_(&uplo,&u_char);
  lower    = lsame_(&uplo,&l_char);
  notrans  = lsame_(&transa,&n_char);
  trans    = lsame_(&transa,&t_char);
  conjt    = lsame_(&transa,&c_char);
  notudiag = lsame_(&diag,&n_char);
  udiag    = lsame_(&diag,&u_char);
  info = 0;
  rside = lsame_(&side,&r_char);
  if ((lside == 0) && (rside == 0)) {
    info = 1;
  } else {
    if ((upper == 0) && (lower == 0)) {
      info = 2;
    } else {
      if ((notrans == 0) && (trans == 0) && (conjt == 0)) {
	info = 3;
      } else {
	if ((udiag == 0) && (notudiag == 0)) {
	  info = 4;
	} else {
	  if (m < 0) {
	    info = 5;
	  } else {
	    if (n < 0) {
	      info = 6;
	    } else {
	      min_lda = 1;
	      if (nrowa > min_lda) {
		min_lda = nrowa;
	      }
	      if (lda < min_lda) {
		info = 9;
	      } else {
		min_ldb = 1;
		if (m > min_ldb) {
		  min_ldb = m;
		}
		if (ldb < min_ldb) {
		  info = 11;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  if (info == 0) {
    mm1_lda = (m-1) * lda;
    /*
      When alpha is zero, all solutions are zero.
    */
    if (alpha == zero) {
      bcolj_pos = 0;
      for (j=0;j<n;j++) {
	bcolj = &b[bcolj_pos];
	for (i=0;i<m;i++) {
	  bcolj[i] = 0;
	}
	bcolj_pos += ldb;
      }
    } else {
      if (lside) {
	/*
	  Left side.
	*/
	if (notrans) {
	  /*
	    Form  B := alpha*inv( A )*B.
	  */
	  if (upper) {
	    bcolj_pos = 0;
	    for (j=0;j<n;j++) {
	      bcolj = &b[bcolj_pos];
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolj,&inc1);
	      }
	      acolk_pos = mm1_lda;
	      for (k=m-1;k>=0;k-=1) {
		acolk = &a[acolk_pos];
		temp = bcolj[k];
		if (temp != zero) {
		  if (nounit) {
		    temp = temp/acolk[k];
		    bcolj[k] = temp;
		  }
		  temp = - temp;
		  daxpy_(&k,&temp,acolk,&inc1,bcolj,&inc1);
		} /* end if (bcolk != zero) */
		acolk_pos = acolk_pos - lda;
	      } /* end for (k...) */
	      bcolj_pos = bcolj_pos + ldb;
	    } /* end for (j ...) */
	    /* end if (upper) */
	  } else {
	    /* lower */
	    bcolj_pos = 0;
	    for (j=0;j<n;j++) {
	      bcolj = &b[bcolj_pos];
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolj,&inc1);
	      }
	      acolk_pos = 0;
	      mmkm1 = m - 1;
	      for (k=0;k<m;k++) {
		acolk = &a[acolk_pos];
		temp = bcolj[k];
		if (temp != zero) {
		  if (nounit) {
		    temp = temp/acolk[k];
		    bcolj[k] = temp;
		  }
		  temp = -temp;
		  daxpy_(&mmkm1,&temp,&acolk[k+1],&inc1,&bcolj[k+1],&inc1);
		}
		mmkm1 = mmkm1 - 1;
		acolk_pos = acolk_pos + lda;
	      } /* end for (k...) */
	      bcolj_pos = bcolj_pos + ldb;
	    } /* end for (j ...) */
	  } /* end else lower */
	  /* end if (notrans) */
	} else {
	  /*
	    Form  B := alpha*inv( A**T )*B.
	  */
	  if (upper) {
	    bcolj_pos = 0;
	    for (j=0;j<n;j++) {
	      bcolj = &b[bcolj_pos];
	      acoli_pos = 0;
	      for (i=0;i<m;i++) {
		acoli = &a[acoli_pos];
		temp = ddot_(&i,acoli,&inc1,bcolj,&inc1);
		temp = (alpha * bcolj[i]) - temp;
		if (nounit) {
		  temp = temp / acoli[i];
		}
		bcolj[i] = temp;
		acoli_pos = acoli_pos + lda;
	      } /* end for (i...) */
	      bcolj_pos = bcolj_pos + ldb;
	    } /* end for (j...) */
	  } else {
	    /*
	      Lower
	    */
	    bcolj_pos = 0;
	    for (j=0;j<n;j++) {
	      bcolj = &b[bcolj_pos];
	      acoli_pos = mm1_lda;
	      mmim1 = 0;
	      for (i=m-1;i>=0;i-=1) {
		acoli = &a[acoli_pos];
		temp = ddot_(&mmim1,&acoli[i+1],&inc1,&bcolj[i+1],&inc1);
		temp = (alpha * bcolj[i]) - temp;
		if (nounit) {
		  temp = temp/acoli[i];
		}
		bcolj[i] = temp;
		acoli_pos = acoli_pos - lda;
		mmim1 += 1;
	      } /* end for (i...) */
	      bcolj_pos = bcolj_pos + ldb;
	    } /* end for (j...) */
	  } /* end else lower */
	} /* end else transpose */
      } else {
	/*
	  rside.
	*/
	if (notrans) {
	  /*
	    Form B = alpha*B*inv(A)
	  */
	  if (upper) {
	    bcolj_pos = 0;
	    acolj_pos = 0;
	    for (j=0;j<n;j++) {
	      bcolj = &b[bcolj_pos];
	      acolj = &a[acolj_pos];
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolj,&inc1);
	      }
	      bcolk_pos = 0;
	      for (k=0;k<j;k++) {
		bcolk = &b[bcolk_pos];
		temp = - acolj[k];
		if (temp != zero) {
		  daxpy_(&m,&temp,bcolk,&inc1,bcolj,&inc1);
		}
		bcolk_pos = bcolk_pos + ldb;
	      } /* end for (k...) */
	      if (nounit) {
		temp = one/acolj[j];
		dscal_(&m,&temp,bcolj,&inc1);
	      }
	      bcolj_pos = bcolj_pos + ldb;
	      acolj_pos = acolj_pos + lda;
	    }
	    /*
	      end if (upper).
	    */
	  } else {
	    /*
	      lower
	    */
	    acolj_pos = (n-1)*lda;
	    bcolj_pos = (n-1)*ldb;
	    for (j=n-1;j>=0;j--) {
	      acolj = &a[acolj_pos];
	      bcolj = &b[bcolj_pos];
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolj,&inc1);
	      }
	      bcolk_pos = bcolj_pos + ldb;
	      for(k=j+1;k<n;k++) {
		temp = - acolj[k];
		bcolk = &b[bcolk_pos];
		if (temp != zero) {
		  daxpy_(&m,&temp,bcolk,&inc1,bcolj,&inc1);
		}
		bcolk_pos = bcolk_pos + ldb;
	      } /* end for (k...) */
	      if (nounit) {
		temp = one/acolj[j];
		dscal_(&m,&temp,bcolj,&inc1);
	      }
	      bcolj_pos = bcolj_pos - ldb;
	      acolj_pos = acolj_pos - lda;
	    } /* end for (j...) */
	  } /* end else lower */
	  /*
	    end if (notrans)
          */
	} else {
	  /*
	    Transpose.
	    Form  B := alpha*B*inv( A**T ).
	  */
	  if (upper) {
	    acolk_pos = (n-1) * lda;
	    bcolk_pos = (n-1) * ldb;
	    for (k=n-1;k>=0;k--) {
	      acolk = &a[acolk_pos];
	      bcolk = &b[bcolk_pos];
	      if (nounit) {
		temp = one/acolk[k];
		dscal_(&m,&temp,bcolk,&inc1);
	      }
	      bcolj_pos = 0;
	      for (j = 0;j<k;j++) {
		bcolj = &b[bcolj_pos];
		temp = -acolk[j];
		if (temp != zero) {
		  daxpy_(&m,&temp,bcolk,&inc1,bcolj,&inc1);
		}
		bcolj_pos = bcolj_pos + ldb;
	      } /* end for (j...) */
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolk,&inc1);
	      }
	      acolk_pos = acolk_pos - lda;
	      bcolk_pos = bcolk_pos - ldb;
	    } /* end for (k ...) */
	    /* end if (upper) */
	  } else {
	    /*
	      Lower
	    */
	    acolk_pos = 0;
	    bcolk_pos = 0;
	    for (k=0;k<n;k++) {
	      acolk = &a[acolk_pos];
	      bcolk = &b[bcolk_pos];
	      if (nounit) {
		temp = one/acolk[k];
		dscal_(&m,&temp,bcolk,&inc1);
	      }
	      bcolj_pos = bcolk_pos + ldb;
	      for (j=k+1;j<n;j++) {
		bcolj = &b[bcolj_pos];
		temp = -acolk[j];
		if (temp != zero) {
		  daxpy_(&m,&temp,bcolk,&inc1,bcolj,&inc1);
		}
		bcolj_pos = bcolj_pos + ldb;
	      } /* end for (j...) */
	      if (alpha != one) {
		dscal_(&m,&alpha,bcolk,&inc1);
	      }
	      acolk_pos = acolk_pos + lda;
	      bcolk_pos = bcolk_pos + ldb;
	    } /* end for (k...) */
	  } /* end else lower */
	} /* end else transpose */
      } /* end else rside */
    } /* end else alpha != 0 */
    /* end if (info == 0) */
  } else {
    fprintf(stderr,"dtrsm: Error info = %d\n",info);
    fflush(stderr);
  }
} /* end dtrsm */
