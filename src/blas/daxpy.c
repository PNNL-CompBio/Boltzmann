/*
  Translation to c of LAPACK daxpy.f
*/
/*
*> \brief \b DAXPY
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
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
*>    DAXPY constant times a vector plus a vector.
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
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
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
      IF (DA.EQ.0.0d0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*        code for both increments equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,4)
         IF (M.NE.0) THEN
            DO I = 1,M
               DY(I) = DY(I) + DA*DX(I)
            END DO
         END IF
         IF (N.LT.4) RETURN
         MP1 = M + 1
         DO I = MP1,N,4
            DY(I) = DY(I) + DA*DX(I)
            DY(I+1) = DY(I+1) + DA*DX(I+1)
            DY(I+2) = DY(I+2) + DA*DX(I+2)
            DY(I+3) = DY(I+3) + DA*DX(I+3)
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
          DY(IY) = DY(IY) + DA*DX(IX)
          IX = IX + INCX
          IY = IY + INCY
         END DO
      END IF
      RETURN
      END
*/
void daxpy_(int *n_p,
	    double *da_p,
	    double *dx_p,
	    int *incx_p,
	    double *dy_p,
	    int *incy_p) {
  double d_zero;
  double da;
  double *dx;
  double *dy;

  int n;
  int incx;

  int incy;
  int i;

  int ix;
  int iy;

  int m;
  int nmm;

  int nmmp1;
  int padi;

  n = *n_p;
  incx = *incx_p;
  incy = *incy_p;
  da   = *da_p;
  dx   = dx_p - 1; /* Caution address arithmetic */
  dy   = dy_p - 1; /* Caution address arithmetic */

  d_zero = 0.0;

  if ((n > 0) && (da != d_zero)) {
    if ((incx == 1) && (incy == 1)) {
      m = n & 7;
      nmm = n - m;
      nmmp1 = nmm + 1;
      for (i=1;i<=nmm;i+=8) {
	dy[i]   += (da*dx[i]);
	dy[i+1] += (da*dx[i+1]);
	dy[i+2] += (da*dx[i+2]);
	dy[i+3] += (da*dx[i+3]);
	dy[i+4] += (da*dx[i+4]);
	dy[i+5] += (da*dx[i+5]);
	dy[i+6] += (da*dx[i+6]);
	dy[i+7] += (da*dx[i+7]);
      }
      /*
	Clean up lop
      */
      if (m > 0) {
	for (i=nmmp1;i<=n;i++) {
	  dy[i]   += (da*dx[i]);
	}
      }
    } else {
      /*
        code for unequal increments or equal increments
        not equal to 1
      */
      ix = 1;
      iy = 1;
      if (incx < 0) {
	ix = 1 + (1-n)*incx;
      }
      if (incy < 0) {
	iy = 1 + (1-n)*incy;
      }
      for (i=1;i<=n;i++) {
	dy[iy]   += (da*dx[ix]);
	ix += incx;
	iy += incy;
      }
    }
  }
  return;
}
