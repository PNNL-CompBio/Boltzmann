#include "boltzmann_structs.h"
#include "dger.h"
/*
    1 *> \brief \b DGER
    2 *
    3 *  =========== DOCUMENTATION ===========
    4 *
    5 * Online html documentation available at 
    6 *            http://www.netlib.org/lapack/explore-html/ 
    7 *
    8 *  Definition:
    9 *  ===========
   10 *
   11 *       SUBROUTINE DGER(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
   12 * 
   13 *       .. Scalar Arguments ..
   14 *       DOUBLE PRECISION ALPHA
   15 *       INTEGER INCX,INCY,LDA,M,N
   16 *       ..
   17 *       .. Array Arguments ..
   18 *       DOUBLE PRECISION A(LDA,*),X(*),Y(*)
   19 *       ..
   20 *  
   21 *
   22 *> \par Purpose:
   23 *  =============
   24 *>
   25 *> \verbatim
   26 *>
   27 *> DGER   performs the rank 1 operation
   28 *>
   29 *>    A := alpha*x*y**T + A,
   30 *>
   31 *> where alpha is a scalar, x is an m element vector, y is an n element
   32 *> vector and A is an m by n matrix.
   33 *> \endverbatim
   34 *
   35 *  Arguments:
   36 *  ==========
   37 *
   38 *> \param[in] M
   39 *> \verbatim
   40 *>          M is INTEGER
   41 *>           On entry, M specifies the number of rows of the matrix A.
   42 *>           M must be at least zero.
   43 *> \endverbatim
   44 *>
   45 *> \param[in] N
   46 *> \verbatim
   47 *>          N is INTEGER
   48 *>           On entry, N specifies the number of columns of the matrix A.
   49 *>           N must be at least zero.
   50 *> \endverbatim
   51 *>
   52 *> \param[in] ALPHA
   53 *> \verbatim
   54 *>          ALPHA is DOUBLE PRECISION.
   55 *>           On entry, ALPHA specifies the scalar alpha.
   56 *> \endverbatim
   57 *>
   58 *> \param[in] X
   59 *> \verbatim
   60 *>          X is DOUBLE PRECISION array of dimension at least
   61 *>           ( 1 + ( m - 1 )*abs( INCX ) ).
   62 *>           Before entry, the incremented array X must contain the m
   63 *>           element vector x.
   64 *> \endverbatim
   65 *>
   66 *> \param[in] INCX
   67 *> \verbatim
   68 *>          INCX is INTEGER
   69 *>           On entry, INCX specifies the increment for the elements of
   70 *>           X. INCX must not be zero.
   71 *> \endverbatim
   72 *>
   73 *> \param[in] Y
   74 *> \verbatim
   75 *>          Y is DOUBLE PRECISION array of dimension at least
   76 *>           ( 1 + ( n - 1 )*abs( INCY ) ).
   77 *>           Before entry, the incremented array Y must contain the n
   78 *>           element vector y.
   79 *> \endverbatim
   80 *>
   81 *> \param[in] INCY
   82 *> \verbatim
   83 *>          INCY is INTEGER
   84 *>           On entry, INCY specifies the increment for the elements of
   85 *>           Y. INCY must not be zero.
   86 *> \endverbatim
   87 *>
   88 *> \param[in,out] A
   89 *> \verbatim
   90 *>          A is DOUBLE PRECISION array of DIMENSION ( LDA, n ).
   91 *>           Before entry, the leading m by n part of the array A must
   92 *>           contain the matrix of coefficients. On exit, A is
   93 *>           overwritten by the updated matrix.
   94 *> \endverbatim
   95 *>
   96 *> \param[in] LDA
   97 *> \verbatim
   98 *>          LDA is INTEGER
   99 *>           On entry, LDA specifies the first dimension of A as declared
  100 *>           in the calling (sub) program. LDA must be at least
  101 *>           max( 1, m ).
  102 *> \endverbatim
  103 *
  104 *  Authors:
  105 *  ========
  106 *
  107 *> \author Univ. of Tennessee 
  108 *> \author Univ. of California Berkeley 
  109 *> \author Univ. of Colorado Denver 
  110 *> \author NAG Ltd. 
  111 *
  112 *> \date November 2011
  113 *
  114 *> \ingroup double_blas_level2
  115 *
  116 *> \par Further Details:
  117 *  =====================
  118 *>
  119 *> \verbatim
  120 *>
  121 *>  Level 2 Blas routine.
  122 *>
  123 *>  -- Written on 22-October-1986.
  124 *>     Jack Dongarra, Argonne National Lab.
  125 *>     Jeremy Du Croz, Nag Central Office.
  126 *>     Sven Hammarling, Nag Central Office.
  127 *>     Richard Hanson, Sandia National Labs.
  128 *> \endverbatim
  129 *>
  130 *  =====================================================================
  131       SUBROUTINE dger(M,N,ALPHA,X,INCX,Y,INCY,A,LDA)
  132 *
  133 *  -- Reference BLAS level2 routine (version 3.4.0) --
  134 *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
  135 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  136 *     November 2011
  137 *
  138 *     .. Scalar Arguments ..
  139       DOUBLE PRECISION ALPHA
  140       INTEGER INCX,INCY,LDA,M,N
  141 *     ..
  142 *     .. Array Arguments ..
  143       DOUBLE PRECISION A(lda,*),X(*),Y(*)
  144 *     ..
  145 *
  146 *  =====================================================================
  147 *
  148 *     .. Parameters ..
  149       DOUBLE PRECISION ZERO
  150       parameter(zero=0.0d+0)
  151 *     ..
  152 *     .. Local Scalars ..
  153       DOUBLE PRECISION TEMP
  154       INTEGER I,INFO,IX,J,JY,KX
  155 *     ..
  156 *     .. External Subroutines ..
  157       EXTERNAL xerbla
  158 *     ..
  159 *     .. Intrinsic Functions ..
  160       INTRINSIC max
  161 *     ..
  */
void dger_(int *mp, int *np, double *alphap, double *xp, int *incxp, 
	   double *yp, int *incyp, double *ap, int *ldap) {
  double *x;
  double *y;
  double *a;
  double alpha;
  double zero;
  double temp;

  int m;
  int n;

  int incx;
  int incy;

  int lda;
  int i;

  int info;
  int ix;
  
  int j;
  int jy;

  int kx;
  int j_lda;

  int lda_inc;

  m    	= *mp;
  n    	= *np;
  incx 	= *incxp;
  incy 	= *incyp;
  lda  	= *ldap;
  alpha = *alphap;
  zero  = 0.0;

  /*
    Caution address arithmetic here to handle fortran array indices
    start at 1 instead of zero.
  */
  a = ap - (lda + 1);
  x = xp - 1;
  y = yp - 1;
  /*
  162 *
  163 *     Test the input parameters.
  164 *
  165       info = 0
  166       IF (m.LT.0) THEN
  167           info = 1
  168       ELSE IF (n.LT.0) THEN
  169           info = 2
  170       ELSE IF (incx.EQ.0) THEN
  171           info = 5
  172       ELSE IF (incy.EQ.0) THEN
  173           info = 7
  174       ELSE IF (lda.LT.max(1,m)) THEN
  175           info = 9
  176       END IF
  177       IF (info.NE.0) THEN
  178           CALL xerbla('DGER  ',info)
  179           RETURN
  180       END IF
  */
  info = 0;
  if (m < 0) {
    info = 1;
  } else if (n < 0) {
    info = 2;
  } else if (incx == 0) {
    info = 5;
  } else if (incy == 0) {
    info = 7;
  } else if ((lda < 1) || (lda < m)) {
    info = 9;
  } 
  if (info != 0) {
    fprintf(stderr,"Error: from dger, info = %d\n",info);
    fflush(stderr);
  } else {
    /*  
  181 *
  182 *     Quick return if possible.
  183 *
  184       IF ((m.EQ.0) .OR. (n.EQ.0) .OR. (alpha.EQ.zero)) RETURN
  185 *
    */
    if ( (m != 0) && (n != 0) && (alpha != zero)) {
      /*
  186 *     Start the operations. In this version the elements of A are
  187 *     accessed sequentially with one pass through A.
  188 *
  189       IF (incy.GT.0) THEN
  190           jy = 1
  191       ELSE
  192           jy = 1 - (n-1)*incy
  193       END IF
      */
      if (incy > 0) {
	jy = 1;
	lda_inc = lda;
      } else {
	jy = 1 - (n-1) * incy;
	lda_inc = -lda;
      }
	/*
  194       IF (incx.EQ.1) THEN
  195           DO 20 j = 1,n
  196               IF (y(jy).NE.zero) THEN
  197                   temp = alpha*y(jy)
  198                   DO 10 i = 1,m
  199                       a(i,j) = a(i,j) + x(i)*temp
  200    10             CONTINUE
  201               END IF
  202               jy = jy + incy
  203    20     CONTINUE
	*/
      if (incx == 1) {
	j_lda = lda;
	for (j = 1;j<=n;j++) {
	  if (y[jy] != zero) {
	    temp = alpha *y[jy];
	    for (i=1;i<=m;i++) {
	      a[i + j_lda] += x[i]*temp;
	    }
	  }
	  jy += incy;
	  j_lda += lda;
	} /* end for (j...) */
      } else {
	/*
  204       ELSE
  205           IF (incx.GT.0) THEN
  206               kx = 1
  207           ELSE
  208               kx = 1 - (m-1)*incx
  209           END IF
  210           DO 40 j = 1,n
  211               IF (y(jy).NE.zero) THEN
  212                   temp = alpha*y(jy)
  213                   ix = kx
  214                   DO 30 i = 1,m
  215                       a(i,j) = a(i,j) + x(ix)*temp
  216                       ix = ix + incx
  217    30             CONTINUE
  218               END IF
  219               jy = jy + incy
  220    40     CONTINUE
  221       END IF
  222 *
	*/
	if (incx > 0) {
	  kx = 1;
	} else {
	  kx = 1 - incx * (m-1);
	}
	j_lda = lda;
	for (j=1;j<=n;j++) {\
	  if (y[jy] != 0) {
	    temp = alpha * y[jy];
	    ix = kx;
	    for (i =1;i<=m;i++) {
	      a[i+j_lda] += x[ix]*temp;
	      ix += incx;
	    }
	  }
	  jy += incy;
	  j_lda += lda;
	}
      } /* end else non-unit incx */
    } /* end else not a quick exit */
  } /* end else valid arguments */
  /*
  223       RETURN
  224 *
  225 *     End of DGER  .
  226 *
  227       END
  */
}
