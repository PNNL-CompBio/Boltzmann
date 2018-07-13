/*
  Translation to c of LAPACK dscal.f
*/
/*

*> \brief \b DSCAL
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       SUBROUTINE DSCAL(N,DA,DX,INCX)
* 
*       .. Scalar Arguments ..
*       DOUBLE PRECISION DA
*       INTEGER INCX,N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION DX(*)
*       ..
*  
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DSCAL scales a vector by a constant.
*>    uses unrolled loops for increment equal to one.
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
*>     modified 3/93 to return if incx .le. 0.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DSCAL(N,DA,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.4.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2011
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION DA
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER I,M,MP1,NINCX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC MOD
*     ..
      IF (N.LE.0 .OR. INCX.LE.0) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
*
*        clean-up loop
*
         M = MOD(N,5)
         IF (M.NE.0) THEN
            DO I = 1,M
               DX(I) = DA*DX(I)
            END DO
            IF (N.LT.5) RETURN
         END IF
         MP1 = M + 1
         DO I = MP1,N,5
            DX(I) = DA*DX(I)
            DX(I+1) = DA*DX(I+1)
            DX(I+2) = DA*DX(I+2)
            DX(I+3) = DA*DX(I+3)
            DX(I+4) = DA*DX(I+4)
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         NINCX = N*INCX
         DO I = 1,NINCX,INCX
            DX(I) = DA*DX(I)
         END DO
      END IF
      RETURN
      END
*/
#include "system_includes.h"
#include "blas.h"
void dscal_(int *n_p, double *da_p, double *dx_p, int *incx_p) {
  double da;
  double *dx;
  int n;
  int incx;
  int i;
  int m;
  int nmm;
  int nmmp1;
  int nincx;
  
  n = *n_p;
  incx = *incx_p;
  dx   = dx_p - 1; /* caution address arithmetic */
  da   = *da_p;

  if ((n > 0) && (incx > 0)) {
    if (incx == 1) {
      m     = n & 7;
      nmm   = n-m;
      nmmp1 = nmm+1;
      for (i=1; i<=nmm; i+=8) {
	dx[i]   = da * dx[i];
	dx[i+1] = da * dx[i+1];
	dx[i+2] = da * dx[i+2];
	dx[i+3] = da * dx[i+3];
	dx[i+4] = da * dx[i+4];
	dx[i+5] = da * dx[i+5];
	dx[i+6] = da * dx[i+6];
	dx[i+7] = da * dx[i+7];
      }
      /*
	cleanup loop.
      */
      if (m > 0) {
	for (i=nmmp1;i<=n;i++) {
	  dx[i] = da * dx[i];
	}
      }
    } else {
      /*
	non unit incx;
      */
      nincx = n*incx;
      for (i=1;i<=nincx;i+=incx) {
	dx[i] = da * dx[i];
      }
    }
  } /* end if nontrivial input.*/
  return;
}
