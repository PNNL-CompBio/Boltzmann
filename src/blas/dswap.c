/* dswap.c
    1 *> \brief \b DSWAP
    2 *
    3 *  =========== DOCUMENTATION ===========
    4 *
    5 * Online html documentation available at 
    6 *            http://www.netlib.org/lapack/explore-html/ 
    7 *
    8 *  Definition:
    9 *  ===========
   10 *
   11 *       SUBROUTINE DSWAP(N,DX,INCX,DY,INCY)
   12 * 
   13 *       .. Scalar Arguments ..
   14 *       INTEGER INCX,INCY,N
   15 *       ..
   16 *       .. Array Arguments ..
   17 *       DOUBLE PRECISION DX(*),DY(*)
   18 *       ..
   19 *  
   20 *
   21 *> \par Purpose:
   22 *  =============
   23 *>
   24 *> \verbatim
   25 *>
   26 *>    interchanges two vectors.
   27 *>    uses unrolled loops for increments equal one.
   28 *> \endverbatim
   29 *
   30 *  Authors:
   31 *  ========
   32 *
   33 *> \author Univ. of Tennessee 
   34 *> \author Univ. of California Berkeley 
   35 *> \author Univ. of Colorado Denver 
   36 *> \author NAG Ltd. 
   37 *
   38 *> \date November 2011
   39 *
   40 *> \ingroup double_blas_level1
   41 *
   42 *> \par Further Details:
   43 *  =====================
   44 *>
   45 *> \verbatim
   46 *>
   47 *>     jack dongarra, linpack, 3/11/78.
   48 *>     modified 12/3/93, array(1) declarations changed to array(*)
   49 *> \endverbatim
   50 *>
   51 *  =====================================================================
   52       SUBROUTINE dswap(N,DX,INCX,DY,INCY)
*/
#include "system_includes.h"
#include "dswap.h"
void dswap_(int *np, double *dxp, int *incxp, double *dyp, int *incyp) {
  double *dx;
  double *dy;
  double dtemp;
  int    n;
  int    incx;
  int    incy;
  int    i;
  int    ix;
  int    iy;
  int    m;
  int    mp1;
  /*
   53 *
   54 *  -- Reference BLAS level1 routine (version 3.4.0) --
   55 *  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
   56 *  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
   57 *     November 2011
   58 *
   59 *     .. Scalar Arguments ..
   60       INTEGER INCX,INCY,N
   61 *     ..
   62 *     .. Array Arguments ..
   63       DOUBLE PRECISION DX(*),DY(*)
   64 *     ..
   65 *
   66 *  =====================================================================
   67 *
   68 *     .. Local Scalars ..
   69       DOUBLE PRECISION DTEMP
   70       INTEGER I,IX,IY,M,MP1
   71 *     ..
   72 *     .. Intrinsic Functions ..
   73       INTRINSIC mod
   74 *     ..
  */
  n = *np;
  incx = *incxp;
  incy = *incxp;
  dx   = dxp - 1;
  dy   = dyp - 1;
  /*
   75       IF (n.LE.0) RETURN
  */
  if (n > 0) {
    /*
   76       IF (incx.EQ.1 .AND. incy.EQ.1) THEN
   */
    if ((incx == 1) && (incy == 1)) {
      /*
   77 *
   78 *       code for both increments equal to 1
   79 *
   80 *
   81 *       clean-up loop
   82 *
      */
      /*
   83          m = mod(n,3)
      */
      m = n/3;
      m = n - (m*3);
      /*
   84          IF (m.NE.0) THEN
   85             DO i = 1,m
   86                dtemp = dx(i)
   87                dx(i) = dy(i)
   88                dy(i) = dtemp
   89             END DO
   90             IF (n.LT.3) RETURN
   91          END IF
      */
      if (m > 0) {
	for (i=1;i<=m;i++) {
	  dtemp = dx[i];
	  dx[i] = dy[i];
	  dy[i] = dtemp;
	}
      }
      if (n >= 3) {
	/*
   92          mp1 = m + 1
   93          DO i = mp1,n,3
   94             dtemp = dx(i)
   95             dx(i) = dy(i)
   96             dy(i) = dtemp
   97             dtemp = dx(i+1)
   98             dx(i+1) = dy(i+1)
   99             dy(i+1) = dtemp
  100             dtemp = dx(i+2)
  101             dx(i+2) = dy(i+2)
  102             dy(i+2) = dtemp
  103          END DO
	*/
	mp1 = m + 1;
	for (i=mp1;i<=n;i+=3) {
	  dtemp   = dx[i];
	  dx[i]   = dy[i];
	  dy[i]   = dtemp;
	  dtemp   = dx[i+1];
	  dx[i+1] = dy[i+1];
	  dy[i+1] = dtemp;
	  dtemp   = dx[i+2];
	  dx[i+2] = dy[i+2];
	  dy[i+2] = dtemp;
	} /* end for i */
      } /* end if (n >= 3) */
	/*
  104       ELSE
	*/
    } else {
      /* 
	 incx != 1 | incy != 1
  105 *
  106 *       code for unequal increments or equal increments not equal
  107 *         to 1
  108 *
  109          ix = 1
  110          iy = 1
  111          IF (incx.LT.0) ix = (-n+1)*incx + 1
  112          IF (incy.LT.0) iy = (-n+1)*incy + 1
  113          DO i = 1,n
  114             dtemp = dx(ix)
  115             dx(ix) = dy(iy)
  116             dy(iy) = dtemp
  117             ix = ix + incx
  118             iy = iy + incy
  119          END DO
  
  120       END IF
  121       RETURN
  122       END
    BLAS
    SRC
    dswap.f
    Generated on Sun Jun 19 2016 20:52:08 for LAPACK by doxygen 1.8.10
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
	dtemp = dx[ix];
	dx[ix] = dy[iy];
	dy[iy] = dtemp;
	ix += incx;
	iy += incy;
      } /* end for (i...) */
    } /* end else incx | incy != 1 */
  } /* end if (n > 0) */
}
