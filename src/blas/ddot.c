/*
  Translation to c of LAPACK ddot.f
*/
/*
*> \brief \b DDOT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*),DY(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DDOT forms the dot product of two vectors.
*>    uses unrolled loops for increments equal to one.
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
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DTEMP
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      DDOT = 0.0d0
      DTEMP = 0.0d0
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DTEMP = DTEMP + DX(I)*DY(I)
            END DO
            IF (N.LT.5) THEN
               DDOT=DTEMP
            RETURN
            END IF
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
          DTEMP = DTEMP + DX(I)*DY(I) + DX(I+1)*DY(I+1) +
     $            DX(I+2)*DY(I+2) + DX(I+3)*DY(I+3) + DX(I+4)*DY(I+4)
         END DO
      ELSE
*
*        code for unequal increments or equal increments
*          not equal to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = DTEMP + DX(IX)*DY(IY)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      DDOT = DTEMP
      RETURN
      END
*/
#include "system_includes.h"
#include "blas.h"
double ddot_(int *n_p,double *dx_p,int *incx_p,double *dy_p,int *incy_p) {
  double *dx;
  double *dy;
  double dtemp;

  int    n;
  int    incx;
  int    incy;

  int    m;
  int    nmm;

  int    nmmp1;
  int    ix;

  int    iy;
  int    i;

  dtemp = 0.0;
  n     = *n_p;
  incx  = *incx_p;
  incy  = *incy_p;
  dx    = dx_p - 1; /* Caution address arithmetic here */
  dy    = dy_p - 1; /* Caution address arithmetic here */

  if (n > 0) {
    if ((incx == 1) && (incy == 1)) {
      /*
	At least one non-unit stride.
      */
      m = n & 7;
      nmm = n-m;
      nmmp1 = nmm + 1;
      for (i=1;i<=nmm;i+=8) {
	dtemp += (dx[i] * dy[i]) + (dx[i+1] * dy[i+1]) +
		 (dx[i+2] * dy[i+2]) + (dx[i+3] * dy[i+3]) +
		 (dx[i+4] * dy[i+4]) + (dx[i+5] * dy[i+5]) +
	         (dx[i+6] * dy[i+6]) + (dx[i+7] * dy[i+7]);
      }          
      if (m > 0) {
	for (i=nmmp1;i<=n;i++) {
	  dtemp += dx[i] * dy[i];
	}
      }
    } else {
      /*
	Code for any non-unit increment
      */
      ix = 1;
      iy = 1;
      if (incx < 0) {
	ix = 1 + (1-n) * incx;
      }
      if (incy < 0) {
	iy = 1 + (1-n) * incy;
      }
      for (i=1;i<=n;i++) {
	dtemp += dx[ix] * dy[iy];
	ix += incx;
	iy += incy;
      }
    }
  }
  return(dtemp);
}
