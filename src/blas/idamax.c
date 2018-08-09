/*
  Translation to c of LAPACK idamax.f
*/
/*
*> \brief \b IDAMAX
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at 
*            http://www.netlib.org/lapack/explore-html/ 
*
*  Definition:
*  ===========
*
*       INTEGER FUNCTION IDAMAX(N,DX,INCX)
* 
*       .. Scalar Arguments ..
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
*>    IDAMAX finds the index of the first element having maximum absolute value.
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
*> \ingroup aux_blas
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
      INTEGER FUNCTION IDAMAX(N,DX,INCX)
*
*  -- Reference BLAS level1 routine (version 3.6.0) --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2015
*
*     .. Scalar Arguments ..
      INTEGER INCX,N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION DX(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DMAX
      INTEGER I,IX
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DABS
*     ..
      IDAMAX = 0
      IF (N.LT.1 .OR. INCX.LE.0) RETURN
      IDAMAX = 1
      IF (N.EQ.1) RETURN
      IF (INCX.EQ.1) THEN
*
*        code for increment equal to 1
*
         DMAX = DABS(DX(1))
         DO I = 2,N
            IF (DABS(DX(I)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(I))
            END IF
         END DO
      ELSE
*
*        code for increment not equal to 1
*
         IX = 1
         DMAX = DABS(DX(1))
         IX = IX + INCX
         DO I = 2,N
            IF (DABS(DX(IX)).GT.DMAX) THEN
               IDAMAX = I
               DMAX = DABS(DX(IX))
            END IF
            IX = IX + INCX
         END DO
      END IF
      RETURN
      END
*/
int idamax_(int *n_p,double *dx_p,int *incx_p) {
  double *dx;
  double dmax;
  double absxi;
  double d_zero;

  int    n;
  int    incx;

  int    i;
  int    iresult;

  int    ix;
  int    padi;
  
  iresult = 0;
  d_zero  = 0.0;
  n       = *n_p;
  incx    = *incx_p;
  dx      = dx_p - 1; /* Caution address arithmetic */
  if ((n > 0) && incx > 0) {
    iresult = 1;
    if (n > 1) {
      dmax = dx[1];
      if (dmax < d_zero) {
	dmax = d_zero - dmax;
      }
      if (incx == 1) {
	for (i=2;i<=n;i++) {
	  absxi = dx[i];
	  if (absxi < d_zero) {
	    absxi = d_zero - absxi;
	  }
	  if (absxi > dmax) {
	    iresult = i;
	    dmax    = absxi;
	  }
	}
      } else {
	ix = 1;
	ix += incx;
	for (i=2;i<=n;i++) {
	  absxi = dx[ix];
	  if (absxi < d_zero) {
	    absxi = d_zero - absxi;
	  }
	  if (absxi > dmax) {
	    iresult = i;
	    dmax    = absxi;
	  }
	  ix += incx;
	}
      }
    }
  }
  return(iresult);
}
