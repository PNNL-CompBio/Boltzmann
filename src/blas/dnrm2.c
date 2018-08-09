/*
  Translation to c of LAPACK dnrm2.f
*/
/*
*> \brief \b DNRM2
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION X(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DNRM2 returns the euclidean norm of a vector via the function
*> name, so that
*>
*>    DNRM2 := sqrt( x'*x )
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
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  -- This version written on 25-October-1982.
*>     Modified on 14-October-1993 to inline the call to DLASSQ.
*>     Sven Hammarling, Nag Ltd.
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION DNRM2(N,X,INCX)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION X(*)
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION ONE,ZERO
      PARAMETER (ONE=1.0D+0,ZERO=0.0D+0)
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION ABSXI,NORM,SCALE,SSQ
      INTEGER IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS,SQRT
*     ..
      IF (N.LT.1 .OR. INCX.LT.1) THEN
          NORM = ZERO
      ELSE IF (N.EQ.1) THEN
          NORM = ABS(X(1))
      ELSE
          SCALE = ZERO
          SSQ = ONE
*        The following loop is equivalent to this call to the LAPACK
*        auxiliary routine:
*        CALL DLASSQ( N, X, INCX, SCALE, SSQ )
*
          DO 10 IX = 1,1 + (N-1)*INCX,INCX
              IF (X(IX).NE.ZERO) THEN
                  ABSXI = ABS(X(IX))
                  IF (SCALE.LT.ABSXI) THEN
                      SSQ = ONE + SSQ* (SCALE/ABSXI)**2
                      SCALE = ABSXI
                  ELSE
                      SSQ = SSQ + (ABSXI/SCALE)**2
                  END IF
              END IF
   10     CONTINUE
          NORM = SCALE*SQRT(SSQ)
      END IF
*
      DNRM2 = NORM
      RETURN
*
*     End of DNRM2.
*
      END
*/
#include "system_includes.h"
#include "blas.h"
double dnrm2_(int *n_p,double *x_p,int *incx_p) {
  double *x;
  double d_one;
  double d_zero;
  double absxi;
  double norm;
  double scale;
  double ssq;
  double ratio;
  int    n;
  int    incx;

  int    ix;
  int    ilim;

  d_zero = 0.0;
  d_one  = 1.0;
  n      = *n_p;
  incx   = *incx_p;
  x      = x_p - 1; /* Caution address arithmetic */

  norm = d_zero;
  if ((n > 0) && (incx > 0)) {
    if (n == 1) {
      norm = x[1];
      if (norm < d_zero) {
	norm = d_zero - norm;
      }
    } else {
      scale = d_zero;
      ssq   = d_one;
      ilim = 1 + (n-1)*incx;
      for (ix=1;ix<=ilim;ix+=incx) {
	absxi = x[ix];
	if (absxi != d_zero) {
	  if (absxi < d_zero) {
	    absxi = d_zero - absxi;
	  }
	  if (scale < absxi) {
	    ratio = scale/absxi;
	    ssq   = d_one + ((ssq*ratio)*ratio);
	    scale = absxi;
	  } else {
	    ratio = absxi/scale;
	    ssq += ratio * ratio;
	  }
	}
      }
      norm = scale * sqrt(ssq);
    }
  }
  return(norm);
}
