/*
  Translation to c of LAPACK dcopy.f
*/
/*
*> \brief \b DCOPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
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
*>    DCOPY copies a vector, x, to a vector, y.
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
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
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
      INTEGER I,IX,IY,M,MP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,7)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DX(I)
            END DO
            IF (N.LT.7) RETURN
         END IF   
         MP1 = M + 1
         DO I = MP1,N,7
            DY(I) = DX(I)
            DY(I+1) = DX(I+1)
            DY(I+2) = DX(I+2)
            DY(I+3) = DX(I+3)
            DY(I+4) = DX(I+4)
            DY(I+5) = DX(I+5)
            DY(I+6) = DX(I+6)
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
            DY(IY) = DX(IX)
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*/
#include "system_includes.h"
#include "blas.h"
void dcopy_(int *n_p,
	    double *dx_p,
	    int *incx_p,
	    double *dy_p,
	    int *incy_p) {
  double *dx;
  double *dy;
  size_t size;
  int    n;
  int    incx;
  int    incy;
  int    i;
  int    ix;
  int    iy;
  int    nmm;
  int    nmmp1;
  

  n    = *n_p;
  incx = *incx_p;
  incy = *incy_p;
  dx = dx_p - 1;
  dy = dy_p - 1;

  if (n > 0) {
    if ((incx == 1) && (incy == 1)) {
      size = ((size_t)(n))*sizeof(double);
      /*
	If both increments are one use memcpy.
      */
      memcpy(dy_p,dx_p,size);
    } else {
      /*
	Both increments are not 1.
      */
      ix = 1;
      iy = 1;
      if (incx < 0) {
	ix = (1-n)*incx + 1;
      }
      if (incy < 0) {
	iy = (1-n)*incy + 1;
      }
      for (i=1;i<=n;i++) {
	dy[iy] = dx[ix];
	ix += incx;
	iy += incy;
      }
    }
  }
  return;
}
