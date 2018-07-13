/*
  dgemv.c
  Translated from LAPACK dgemv.f
*/
/*
*> \brief \b DGEMV
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION ALPHA,BETA
*       INTEGER INCX,INCY,LDA,M,N
*       CHARACTER TRANS
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DGEMV  performs one of the matrix-vector operations
*>
*>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
*>
*> where alpha and beta are scalars, x and y are vectors and A is an
*> m by n matrix.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] TRANS
*> \verbatim
*>          TRANS is CHARACTER*1
*>           On entry, TRANS specifies the operation to be performed as
*>           follows:
*>
*>              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*>
*>              TRANS = 'T' or 't'   y := alpha*A**T*x + beta*y.
*>
*>              TRANS = 'C' or 'c'   y := alpha*A**T*x + beta*y.
*> \endverbatim
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>           On entry, M specifies the number of rows of the matrix A.
*>           M must be at least zero.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>           On entry, N specifies the number of columns of the matrix A.
*>           N must be at least zero.
*> \endverbatim
*>
*> \param[in] ALPHA
*> \verbatim
*>          ALPHA is DOUBLE PRECISION.
*>           On entry, ALPHA specifies the scalar alpha.
*> \endverbatim
*>
*> \param[in] A
*> \verbatim
*>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*>           Before entry, the leading m by n part of the array A must
*>           contain the matrix of coefficients.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>           On entry, LDA specifies the first dimension of A as declared
*>           in the calling (sub) program. LDA must be at least
*>           max( 1, m ).
*> \endverbatim
*>
*> \param[in] X
*> \verbatim
*>          X is DOUBLE PRECISION array of DIMENSION at least
*>           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*>           Before entry, the incremented array X must contain the
*>           vector x.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>           On entry, INCX specifies the increment for the elements of
*>           X. INCX must not be zero.
*> \endverbatim
*>
*> \param[in] BETA
*> \verbatim
*>          BETA is DOUBLE PRECISION.
*>           On entry, BETA specifies the scalar beta. When BETA is
*>           supplied as zero then Y need not be set on input.
*> \endverbatim
*>
*> \param[in,out] Y
*> \verbatim
*>          Y is DOUBLE PRECISION array of DIMENSION at least
*>           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*>           and at least
*>           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*>           Before entry with BETA non-zero, the incremented array Y
*>           must contain the vector y. On exit, Y is overwritten by the
*>           updated vector y.
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>           On entry, INCY specifies the increment for the elements of
*>           Y. INCY must not be zero.
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
*> \date November 2015
*
*> \ingroup double_blas_level2
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  Level 2 Blas routine.
*>  The vector and matrix arguments are not referenced when N = 0, or M = 0
*>
*>  -- Written on 22-October-1986.
*>     Jack Dongarra, Argonne National Lab.
*>     Jeremy Du Croz, Nag Central Office.
*>     Sven Hammarling, Nag Central Office.
*>     Richard Hanson, Sandia National Labs.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)
*
*  -- Reference BLAS level2 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION ALPHA,BETA
      INTEGER INCX,INCY,LDA,M,N
      CHARACTER TRANS
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION A(LDA,*),X(*),Y(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION TEMP
      INTEGER I,INFO,IX,IY,J,JX,JY,KX,KY,LENX,LENY
*     ..
*     .. External Functions ..
      LOGICAL LSAME
      EXTERNAL LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MAX
*     ..
*
*     Test the input parameters.
*
      INFO = 0
      IF (.NOT.LSAME(TRANS,'N') .AND. .NOT.LSAME(TRANS,'T') .AND.
     +    .NOT.LSAME(TRANS,'C')) THEN
          INFO = 1
      ELSE IF (M.LT.0) THEN
          INFO = 2
      ELSE IF (N.LT.0) THEN
          INFO = 3
      ELSE IF (LDA.LT.MAX(1,M)) THEN
          INFO = 6
      ELSE IF (INCX.EQ.0) THEN
          INFO = 8
      ELSE IF (INCY.EQ.0) THEN
          INFO = 11
      END IF
      IF (INFO.NE.0) THEN
          CALL XERBLA('DGEMV ',INFO)
          RETURN
      END IF
*
*     Quick return if possible.
*
      IF ((M.EQ.0) .OR. (N.EQ.0) .OR.
     +    ((ALPHA.EQ.ZERO).AND. (BETA.EQ.ONE))) RETURN
*
*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
*     up the start points in  X  and  Y.
*
      IF (LSAME(TRANS,'N')) THEN
          LENX = N
          LENY = M
      ELSE
          LENX = M
          LENY = N
      END IF
      IF (INCX.GT.0) THEN
          KX = 1
      ELSE
          KX = 1 - (LENX-1)*INCX
      END IF
      IF (INCY.GT.0) THEN
          KY = 1
      ELSE
          KY = 1 - (LENY-1)*INCY
      END IF
*
*     Start the operations. In this version the elements of A are
*     accessed sequentially with one pass through A.
*
*     First form  y := beta*y.
*
      IF (BETA.NE.ONE) THEN
          IF (INCY.EQ.1) THEN
              IF (BETA.EQ.ZERO) THEN
                  DO 10 I = 1,LENY
                      Y(I) = ZERO
   10             CONTINUE
              ELSE
                  DO 20 I = 1,LENY
                      Y(I) = BETA*Y(I)
   20             CONTINUE
              END IF
          ELSE
              IY = KY
              IF (BETA.EQ.ZERO) THEN
                  DO 30 I = 1,LENY
                      Y(IY) = ZERO
                      IY = IY + INCY
   30             CONTINUE
              ELSE
                  DO 40 I = 1,LENY
                      Y(IY) = BETA*Y(IY)
                      IY = IY + INCY
   40             CONTINUE
              END IF
          END IF
      END IF
      IF (ALPHA.EQ.ZERO) RETURN
      IF (LSAME(TRANS,'N')) THEN
*
*        Form  y := alpha*A*x + y.
*
          JX = KX
          IF (INCY.EQ.1) THEN
              DO 60 J = 1,N
                  TEMP = ALPHA*X(JX)
                  DO 50 I = 1,M
                      Y(I) = Y(I) + TEMP*A(I,J)
   50             CONTINUE
                  JX = JX + INCX
   60         CONTINUE
          ELSE
              DO 80 J = 1,N
                  TEMP = ALPHA*X(JX)
                  IY = KY
                  DO 70 I = 1,M
                      Y(IY) = Y(IY) + TEMP*A(I,J)
                      IY = IY + INCY
   70             CONTINUE
                  JX = JX + INCX
   80         CONTINUE
          END IF
      ELSE
*
*        Form  y := alpha*A**T*x + y.
*
          JY = KY
          IF (INCX.EQ.1) THEN
              DO 100 J = 1,N
                  TEMP = ZERO
                  DO 90 I = 1,M
                      TEMP = TEMP + A(I,J)*X(I)
   90             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  100         CONTINUE
          ELSE
              DO 120 J = 1,N
                  TEMP = ZERO
                  IX = KX
                  DO 110 I = 1,M
                      TEMP = TEMP + A(I,J)*X(IX)
                      IX = IX + INCX
  110             CONTINUE
                  Y(JY) = Y(JY) + ALPHA*TEMP
                  JY = JY + INCY
  120         CONTINUE
          END IF
      END IF
*
      RETURN
*
*     End of DGEMV .
*
      END
*/
#include "system_includes.h"
#include "blas.h"
void dgemv_(char *trans_p,
	    int *m_p,
	    int *n_p,
	    double *alpha_p,
	    double *a_p,
	    int *lda_p,
	    double *x_p,
	    int *incx_p,
	    double *beta_p,
	    double *y_p,
	    int *incy_p,
	    int len_trans) {
  /*
    Called by:
    Calls:     lsame_
  */
  double *a;
  double *x;
  double *y;
  double alpha;
  double beta;

  double d_one;
  double d_zero;
  double temp;

  int m;
  int n;
  
  int lda;
  int incx;

  int incy;
  int i;

  int info;
  int ix;

  int iy;
  int j;

  int jx;
  int jy;

  int kx;
  int ky;

  int lenx;
  int leny;

  int valid_trans;
  int i_one;

  int max_1_m;
  int no_trans;

  int jlda;

  char trans;
  char cchar;
  char tchar;
  char nchar;
  /*
    Unpack arguments
  */
  cchar = 'C';
  nchar = 'N';
  tchar = 'T';
  trans = *trans_p;
  m     = *m_p;
  n     = *n_p;
  lda   = *lda_p;
  incx  = *incx_p;
  incy  = *incy_p;
  alpha = *alpha_p;
  beta  = *beta_p;
  a     = a_p - 1 - lda; /* Caution address arithmetic */
  x     = x_p - 1;       /* Caution address arithmetic */
  y     = y_p - 1;       /* Caution address arithmetic */

  i_one = 1;
  d_zero = 0.0;
  d_one  = 1.0;
  /*
    Check for valid arguments.
  */
  info = 0;
  
  if (incy == 0) {
    info = 11;
  }
  if (incx == 0) {
    info = 8;
  }
  max_1_m = 1;
  if (m > max_1_m ) {
    max_1_m = m;
  }
  if (lda < max_1_m) {
    info = 6;
  }
  if (n < 0) {
    info = 3;
  }
  if (m < 0) {
    info = 2;
  } 
  valid_trans = 0;
  no_trans    = lsame_(&trans,&nchar,i_one,i_one);
  valid_trans = no_trans + lsame_(&trans,&cchar,i_one,i_one) +
    lsame_(&trans,&tchar,i_one,i_one);
  if (valid_trans == 0) {
    info = 1;
  } 
  if (info != 0) {
    fprintf(stderr,"dgemv_ : Error argument %d is invalid \n",info);
  } else {
    /*
      Check for trivial input and quic return
    */
    if ((m > 0) && (n > 0) && ((alpha != d_zero) || (beta != d_one))) {
      /*
        Set  LENX  and  LENY, the lengths of the vectors x and y, and set
        up the start points in  X  and  Y.
      */
      if (no_trans) {
	lenx = n;
	leny = m;
      } else {
	lenx = m;
	leny = n;
      }
      if (incx > 0) {
	kx = 1;
      } else {
	kx = 1 - (lenx-1)*incx;
      }
      if (incy > 0) {
	ky = 1;
      } else {
	ky = 1 - (leny-1)*incy;
      }
      /*

        Start the operations. In this version the elements of A are
        accessed sequentially with one pass through A.

        First form  y := beta*y.

      */
      if (beta != d_one) {
	if (incy == 1) {
	  if (beta == d_zero) {
	    for (i=1;i<=leny;i++) {
	      y[i] = d_zero;
	    }
	  } else {
	    for (i=1;i<=leny;i++) {
	      y[i] = beta * y[i];
	    }
	  }
	} else {
	  iy = ky;
	  if (beta == d_zero) {
	    for (i=1;i<=leny;i++) {
	      y[iy] = d_zero;
	      iy += incy;
	    }
	  } else {
	    for (i=1;i<=leny;i++) {
	      y[iy] = beta *y[iy];
	      iy += incy;
	    }
	  }
	}
      }
      if (alpha != d_zero) {
	if (no_trans) {
	  /*
	    Form  y := alpha*A*x + y.
	  */
          jx = kx;
          if (incy == 1) {
	    jlda = lda;
	    for (j=1;j<=n;j++) {
	      temp = alpha * x[jx];
	      for (i=1;i<=m;i++) {
		y[i] = y[i] + temp * a[i+jlda];
	      }
	      jlda += lda;
	      jx += incx;
	    }
	  } else {
	    jlda = lda;
	    for (j=1;j<=n;j++) {
	      temp = alpha * x[jx];
	      iy = ky;
	      for (i=1;i<=m;i++) {
		y[iy] += temp * a[i+jlda];
		iy    += incy;
	      }
	      jlda += lda;
	      jx += incx;
	    }
	  }
	} else {
	  /*
	    Form  y <= alpha*A^T*x + y.
	  */
          jy = ky;
	  if (incx == 1) {
	    jlda = lda;
	    for (j=1;j<=n;j++) {
	      temp = d_zero;
	      for (i=1;i<=m;i++) {
		temp += a[i+jlda] * x[i];
	      }
	      y[jy] += alpha * temp;
	      jy += incy;
	      jlda += lda;
	    }
	  } else {
	    jlda = lda;
	    for (j=1;j<=n;j++) {
	      temp = d_zero;
	      ix = kx;
	      for (i=1;i<=m;i++) {
		temp = temp + a[i+jlda] * x[ix];
		ix += incx;
	      }
	      y[jy] += alpha * temp;
	      jy += incy;
	      jlda += lda;
	    }
	  } /* end else incx != 1 */
	} /* end else a transposse */
      } /* end else alpha != 0. */
    } /* end if nontrivial input*/
  } /* end else valid input */
  return;
}
