/* dtbsv.c
*> \brief \b DTBSV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,K,LDA,N
*       CHARACTER DIAG,TRANS,UPLO
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DTBSV  solves one of the systems of equations
*>
*>    A*x = b,   or   A**T*x = b,
*>
*> where b and x are n element vectors and A is an n by n unit, or
*> non-unit, upper or lower triangular band matrix, with ( k + 1 )
*> diagonals.
*>
*> No test for singularity or near-singularity is included in this
*> routine. Such tests must be performed before calling this routine.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] UPLO
*> \verbatim
*>          UPLO is CHARACTER*1
*>           On entry, UPLO specifies whether the matrix is an upper or
*>           lower triangular matrix as follows:
*>
*>              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*>
*>              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*> \endverbatim
*>
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the equations to be solved as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   A*x = b.
*>
*>              TRANS = 'T' or 't'   A**T*x = b.
*>
*>              TRANS = 'C' or 'c'   A**T*x = b.
*> \endverbatim
*>
*> \param[in] DIAG
*> \verbatim
*>          DIAG is CHARACTER*1
*>           On entry, DIAG specifies whether or not A is unit
*>           triangular as follows:
*>
*>              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*>
*>              DIAG = 'N' or 'n'   A is not assumed to be unit
*>                                  triangular.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the order of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>           On entry with UPLO = 'U' or 'u', K specifies the number of
*>           super-diagonals of the matrix A.
*>           On entry with UPLO = 'L' or 'l', K specifies the number of
*>           sub-diagonals of the matrix A.
*>           K must satisfy  0 .le. K.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*>           by n part of the array A must contain the upper triangular
*>           band part of the matrix of coefficients, supplied column by
*>           column, with the leading diagonal of the matrix in row
*>           ( k + 1 ) of the array, the first super-diagonal starting at
*>           position 2 in row k, and so on. The top left k by k triangle
*>           of the array A is not referenced.
*>           The following program segment will transfer an upper
*>           triangular band matrix from conventional full matrix storage
*>           to band storage:
*>
*>                 DO 20, J = 1, N
*>                    M = K + 1 - J
*>                    DO 10, I = MAX( 1, J - K ), J
*>                       A( M + I, J ) = matrix( I, J )
*>              10    CONTINUE
*>              20 CONTINUE
*>
*>           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*>           by n part of the array A must contain the lower triangular
*>           band part of the matrix of coefficients, supplied column by
*>           column, with the leading diagonal of the matrix in row 1 of
*>           the array, the first sub-diagonal starting at position 1 in
*>           row 2, and so on. The bottom right k by k triangle of the
*>           array A is not referenced.
*>           The following program segment will transfer a lower
*>           triangular band matrix from conventional full matrix storage
*>           to band storage:
*>
*>                 DO 20, J = 1, N
*>                    M = 1 - J
*>                    DO 10, I = J, MIN( N, J + K )
*>                       A( M + I, J ) = matrix( I, J )
*>              10    CONTINUE
*>              20 CONTINUE
*>
*>           Note that when DIAG = 'U' or 'u' the elements of the array A
*>           corresponding to the diagonal elements of the matrix are not
*>           referenced, but are assumed to be unity.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           ( k + 1 ).
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is DOUBLE PRECISION array of dimension at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ).
*>           Before entry, the incremented array X must contain the n
*>           element right-hand side vector b. On exit, X is overwritten
*>           with the solution vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee 
*> \author Univ. of California Berkeley 
*> \author Univ. of Colorado Denver 
*> \author NAG Ltd. 
*
*> \date November 2011
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DTBSV(UPLO,TRANS,DIAG,N,K,A,LDA,X,INCX)
*/
#include "system_includes.h"
#include "lsame.h"
#include "dtbsv.h"
void dtbsv_(char *uplop, char *transp, char *diagp, int *np, int *kp,
	    double *ap, int *ldap, double *xp, int *incxp) {
  /*
*
*  -- Reference BLAS level2 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,K,LDA,N
      CHARACTER DIAG,TRANS,UPLO
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*)
*     ..
*
*  =====================================================================
*
*/
  double *a;
  double *x;
  double zero;
  double temp;

  int    n;
  int    k;

  int    lda;
  int    incx;

  int    i;
  int    info;

  int    ix;
  int    j;

  int    jx;
  int    kplus1;

  int    kx;
  int    l;

  int    nounit;
  int    max_1_jmk;

  int    j_lda;
  int    min_n_jpk;

  char   uplo;
  char   trans;
  char   diag;
  char   uchar;
  char   lchar;
  char   nchar;
  char   tchar;
  char   cchar;
  
  /*  
  
*     .. Parameters ..
      DOUBLE PRECISION ZERO
      PARAMETER (ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,J,JX,KPLUS1,KX,L
      LOGICAL NOUNIT
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX,MIN
*     ..
  */
  uplo  = *uplop;
  trans = *transp;
  diag  = *diagp;
  n     = *np;
  k     = *kp;
  lda   = *ldap;
  incx  = *incxp;
  /*
    Caution address arithmetic here to account for starting index 
    differencess between c and fortran.
  */
  a     = ap - 1 - lda;
  x     = xp - 1;
  zero  = 0.0;
  uchar = 'U';
  lchar = 'L';
  tchar = 'T';
  nchar = 'N';
  cchar = 'C';
  /*
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(UPLO,'U') .AND. .NOT.LSAME(UPLO,'L')) THEN
          INFO = 1
      ELSE IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +         .NOT.LSAME(TRANS,'C')) THEN
          INFO = 2
      ELSE IF (.NOT.LSAME(DIAG,'U') .AND. .NOT.LSAME(DIAG,'N')) THEN
          INFO = 3
      ELSE IF (N.LT.0) THEN
          INFO = 4
      ELSE IF (K.LT.0) THEN
          INFO = 5
      ELSE IF (LDA.LT. (K+1)) THEN
          INFO = 7
      ELSE IF (INCX.EQ.0) THEN
          INFO = 9
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DTBSV ',INFO)
          RETURN
      END IF
  */
  info = 0;
  if ((lsame_(&uplo,&uchar) || lsame_(&uplo,&lchar)) == 0) {
    info = 1;
  } else if ((lsame_(&trans,&nchar) || lsame_(&trans,&tchar) || 
	      lsame_(&trans,&cchar)) == 0) {
    info = 2;
  } else if ((lsame_(&diag,&uchar) || lsame_(&diag,&nchar)) == 0) {
    info = 3;
  } else if (n < 0) {
    info = 4;
  } else if (k < 0) {
    info = 5;
  } else if (lda < k+1) {
    info = 7;
  } else if (incx == 0) {
    info = 9;
  }
  if (info != 0) {
    fprintf(stderr,"dtbsv_: Error in argument %d\n",info);
    fflush(stderr);
  } else {
    /*    
*
*     Quick return if possible.
*
      IF (N.EQ.0) RETURN
    */
    if (n > 0) {
      /*
      NOUNIT = LSAME(DIAG,'N')
      */
      nounit = lsame_(&diag,&nchar);
      /*
	Set up the start point in X if the increment is not unity. This
        will be  ( N - 1 )*INCX  too small for descending loops.
      IF (INCX.LE.0) THEN
          KX = 1 - (N-1)*INCX
      ELSE IF (INCX.NE.1) THEN
          KX = 1
      END IF
      */
      if (incx <= 0) {
	kx = 1 - ((n-1)*incx);
      } else if (incx != 1) {
	kx = 1;
      }
      /*
*
*     Start the operations. In this version the elements of A are
*     accessed by sequentially with one pass through A.
*
      IF (LSAME(TRANS,'N')) THEN
      */
      if (lsame_(&trans,&nchar)) {
	/*
*
*        Form  x := inv( A )*x.
*
          IF (LSAME(UPLO,'U')) THEN
	*/
	if (lsame_(&uplo,&uchar)) {
	  /*
              KPLUS1 = K + 1
              IF (INCX.EQ.1) THEN
	  */
	  kplus1 = k + 1;
	  if (incx == 1) {
	    /*
	    DO 20 J = N,1,-1
                      IF (X(J).NE.ZERO) THEN
                          L = KPLUS1 - J
                          IF (NOUNIT) X(J) = X(J)/A(KPLUS1,J)
                          TEMP = X(J)
                          DO 10 I = J - 1,MAX(1,J-K),-1
                              X(I) = X(I) - TEMP*A(L+I,J)
   10                     CONTINUE
                      END IF
   20             CONTINUE
	    */
	    j_lda = n * lda;
	    for (j=n;j>=1;j--) {
	      if (x[j] != zero) {
		l = kplus1 - j;
		if (nounit) {
		  x[j] = x[j]/a[j_lda + kplus1];
		}
		temp = x[j];
		max_1_jmk = j - k;
		if (1 > max_1_jmk) {
		  max_1_jmk = 1;
		}
		for (i=j-1;i>=max_1_jmk;i--) {
		  x[i] = x[i] - (temp * a[l+i+j_lda]);
		}
	      }
	      j_lda -= lda;
	    }
	    /*	      
              ELSE
                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 40 J = N,1,-1
                      KX = KX - INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = KPLUS1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(KPLUS1,J)
                          TEMP = X(JX)
                          DO 30 I = J - 1,MAX(1,J-K),-1
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX - INCX
   30                     CONTINUE
                      END IF
                      JX = JX - INCX
   40             CONTINUE
              END IF
	    */
	  } else {
	    /* 
	      incx != 1
	    */
	    kx = kx + (n-1) * incx;
	    jx = kx;
	    j_lda = n*lda;
	    for (j = n;j>=1;j--) {
	      kx = kx - incx;
	      if (x[jx] != zero) {
		ix = kx;
		l = kplus1 - j;
		if (nounit) {
		  x[jx] = x[jx]/a[j_lda + kplus1];
		}
		temp = x[jx];
		max_1_jmk = j - k;
		if (1 > max_1_jmk) {
		  max_1_jmk = 1;
		}
		for (i=j-1;i>=max_1_jmk;i--) {
		  x[ix] = x[ix] - (temp * a[j_lda + l+i]);
		  ix = ix - incx;
		} /* end for (i... ) */
	      } /* end if (x[jx] != zero) */
	      j_lda -= lda;
	    } /* end for j */
	    /*
          ELSE
              IF (INCX.EQ.1) THEN
                  DO 60 J = 1,N
                      IF (X(J).NE.ZERO) THEN
                          L = 1 - J
                          IF (NOUNIT) X(J) = X(J)/A(1,J)
                          TEMP = X(J)
                          DO 50 I = J + 1,MIN(N,J+K)
                              X(I) = X(I) - TEMP*A(L+I,J)
   50                     CONTINUE
                      END IF
   60             CONTINUE
	    */
	  }
	} else {
	  /* 
	    uplo != U
	  */
	  if (incx == 1) {
	    j_lda = lda;
	    for (j=1;j<=n;j++) {
	      if (x[j] != zero) {
		l = 1-j;
		if (nounit) {
		  x[j] = x[j]/a[j_lda + 1];
		}
		temp = x[j];
		min_n_jpk = j + k;
		if (n < min_n_jpk) {
		  min_n_jpk = n;
		}
		for (i=j+1;i<=min_n_jpk;i++) {
		  x[i] = x[i] - (temp * a[j_lda+l+i]);
		}
	      }
	      j += lda;
	    } /* end for (j...) */
	    /*
              ELSE
                  JX = KX
                  DO 80 J = 1,N
                      KX = KX + INCX
                      IF (X(JX).NE.ZERO) THEN
                          IX = KX
                          L = 1 - J
                          IF (NOUNIT) X(JX) = X(JX)/A(1,J)
                          TEMP = X(JX)
                          DO 70 I = J + 1,MIN(N,J+K)
                              X(IX) = X(IX) - TEMP*A(L+I,J)
                              IX = IX + INCX
   70                     CONTINUE
                      END IF
                      JX = JX + INCX
   80             CONTINUE
              END IF
          END IF
	    */
	  } else {
	    /*
	      incx != 1;
	    */
	    jx = kx;
	    j_lda = lda;
	    for (j=1;j<=n;j++) {
	      kx = kx + incx;
	      if (x[jx] != zero) {
		ix = kx;
		l = 1-j;
		if (nounit) {
		  x[jx] = x[jx]/a[j_lda+1];
		}
		temp = x[jx];
		min_n_jpk = j + k;
		if (n < min_n_jpk) {
		  min_n_jpk = n;
		}
		for (i=j+1;i<=min_n_jpk;i++) {
		  x[ix] = x[ix] - (temp * a[j_lda + l + i]);
		  ix = ix + incx;
		}
	      }
	      jx = jx + incx;
	      j_lda += lda;
	    } /* end for (j...) */
	  }
	}
	/*
      ELSE
*
*        Form  x := inv( A**T)*x.
*
        */
      } else {
	/*
          IF (LSAME(UPLO,'U')) THEN
              KPLUS1 = K + 1
	*/
	if (lsame_(&uplo,&uchar)) {
	  kplus1 = k+1;
	  /*	    
	    IF (INCX.EQ.1) THEN
	      DO 100 J = 1,N
                      TEMP = X(J)
                      L = KPLUS1 - J
                      DO 90 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(I)
   90                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(J) = TEMP
  100          CONTINUE
	  */
	  if (incx == 1) {
	    j_lda = lda;
	    for (j=1;j<=n;j++) {
	      temp = x[j];
	      l = kplus1 - j;
	      max_1_jmk = j-k;
	      if (1 > max_1_jmk) {
		max_1_jmk = 1;
	      }
	      for (i=max_1_jmk;i<=j-1;i++) {
		temp = temp - a[j_lda+l+i]*x[i];
	      } /* end for (i...) */
	      if (nounit) {
		temp = temp/a[j_lda+kplus1];
	      }
	      x[j] = temp;
	      j_lda += lda;
	    } /* end for (j...) */
	      /*
              ELSE
	      */
	  } else {
	    /*
	      incx != 1
                  JX = KX
                  DO 120 J = 1,N
                      TEMP = X(JX)
                      IX = KX
                      L = KPLUS1 - J
                      DO 110 I = MAX(1,J-K),J - 1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX + INCX
  110                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(KPLUS1,J)
                      X(JX) = TEMP
                      JX = JX + INCX
                      IF (J.GT.K) KX = KX + INCX
  120             CONTINUE
              END IF
	    */
	    jx = kx;
	    j_lda = lda;
	    for (j=1;j<=n;j++) {
	      temp = x[jx];
	      ix = ix;
	      l = kplus1 - j;
	      max_1_jmk = j-k;
	      if (1 > max_1_jmk) {
		max_1_jmk = 1;
	      }
	      for (i=max_1_jmk;i<=j-1;i++) {
		temp = temp - (a[j_lda+l+i]*x[ix]);
		ix = ix + incx;
	      } /* end for (i...) */
	      if (nounit) {
		temp = temp/a[j_lda + kplus1];
	      }
	      x[jx] = temp;
	      jx = jx + incx;
	      if (j>k) {
		kx = kx + incx;
	      }
	      j_lda += lda;
	    } /* end for j */
	  }
	    /*
          ELSE
	    */
	} else {
	  /*
	    uplo != U
	  */
	  /*
              IF (INCX.EQ.1) THEN
                  DO 140 J = N,1,-1
                      TEMP = X(J)
                      L = 1 - J
                      DO 130 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(I)
  130                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(J) = TEMP
  140             CONTINUE
	  */
	  if (incx == 1) {
	    j_lda = n * lda;
	    for (j=n;j>=1;j--) {
	      temp = x[j];
	      l = 1 - j;
	      min_n_jpk = j+k;
	      if (n < min_n_jpk) {
		min_n_jpk = n;
	      }
	      for (i=min_n_jpk;i>=j+1;i--) {
		temp = temp - (a[j_lda+l+i]*x[i]);
	      } /* end for (i...) */
	      if (nounit) {
		temp = temp/a[j_lda+1];
	      }
	      x[j] = temp;
	      j_lda -= lda;
	    } /* end for j */
	    /*
              ELSE
	    */
	  } else {
	    /*
    	      incx != 1;

                  KX = KX + (N-1)*INCX
                  JX = KX
                  DO 160 J = N,1,-1
                      TEMP = X(JX)
                      IX = KX
                      L = 1 - J
                      DO 150 I = MIN(N,J+K),J + 1,-1
                          TEMP = TEMP - A(L+I,J)*X(IX)
                          IX = IX - INCX
  150                 CONTINUE
                      IF (NOUNIT) TEMP = TEMP/A(1,J)
                      X(JX) = TEMP
                      JX = JX - INCX
                      IF ((N-J).GE.K) KX = KX - INCX
  160             CONTINUE
              END IF
	    */
	    kx = kx + (n-1)*incx;
	    jx = kx;
	    j_lda = n*lda;
	    for (j=n;j>=1;j--) {
	      temp = x[jx];
	      ix = kx;
	      l = 1-j;
	      min_n_jpk = j+k;
	      if (n < min_n_jpk) {
		min_n_jpk = n;
	      }
	      for (i=min_n_jpk;i>=j+1;i--) {
		temp = temp - (a[j_lda+l+i]*x[ix]);
		ix = ix - incx;
	      }
	      if (nounit) {
		temp = temp/a[j_lda+1];
	      }
	      x[jx] = temp;
	      jx = jx -incx;
	      if ((n-j) >= k) {
		kx = kx - incx;
	      }
	      j_lda -= lda;
	    } /* end for (j...) */
	      /*
          END IF
      END IF
*
      RETURN
*
*     End of DTBSV .
*
      END
	      */
	  } /* end else incx != 1 */
	} /* end else uplo != u */
      } /* end else solving the transpose problem. */
    } /* end if not quick return */
  } /* end else valid parameters */
}

    
