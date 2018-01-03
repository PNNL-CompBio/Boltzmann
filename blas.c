#include "boltzmann_structs.h"

#include "blas.h"
void dscal(int *nx, double *a_p, double *x, int *inc) {
  /*
    Dscal routine that ignores inc, assuming its 1.
    x <- a * x
    Called by: ode23tb
  */
  double a;
  int    i;
  int    padi;
  a = *a_p;
  for (i=0;i<(*nx);i++) {
    x[i] = a * x[i];
  }
}

void dcopy(int *nx, double *x, int *incx, double *y, int *incy) {
  /*
    dcopy routine that ignores incx and incy assuming them both to be 1.
    y <- x
    Called by: ode_num_jac
  */
  int i;
  int padi;
  for (i=0;i<(*nx);i++) {
    y[i] = x[i];
  }
}

double dnrm2(int *nx, double *x, int *inc) {
  /*
    Dnrm2 routine that ignores inc, assuming its 1.
    This still needs to be rewritten for stability and 
    scaling.
    Called by: ode23tb
    Calls:     sqrt
  */
  double sum;
  double twonorm;
  int i;
  int padi;
  sum = 0.0;
  for (i=0;i<(*nx);i++) {
    sum = sum + (x[i]*x[i]);
  }
  twonorm = sqrt(sum);
  return(twonorm);
}


void dgemv(int *trans, int *m_p, int *n_p, double *alpha_p, double *a, 
	   int *lda_p, double *x, int *incx_p, double *beta_p, 
	   double *y, int *incy_p) {
   /*
     Simplified version of dgemv that assumes trans = 0,
     and the scalars alpha and beta are both 1.0 and 
     incx and incy are both 1.
     a is an m x n  matrix stored contiguously,
     x is a n x 1 vector, and y is an  m x 1 vector
     We compute y = y + Ax using daxpy operations.
     Called by: ode23tb
     Calls:     daxpy
   */
  double *a_col;
  double xv;
  int m;
  int n;
  int lda;
  int j;
  int incx;
  int incy;
  a_col = a;
  m   = *m_p;
  n   = *n_p;
  incx = 1;
  incy = 1;
  /*
    if beta were not 1, we would want:
    for (i=0;i<m;i++) {
       y[i] = beta * y[i]
    }       
  */
  for (j=0;j<n;j++) {
    /*
      if alpha wer not wone we would want 
      xv = alpha * xj
    */
    xv = x[j];
    daxpy(&m,&xv,a_col,&incx,y,&incy);
    a_col += m; /* Caution address arithmetic here */
    /* in general we would replace the m above with lda. */
  }
}

void daxpy(int *n_p, double *alpha_p, double *x, int *incx_p, 
	   double *y, int *incy_p) {
  /*
    Simplified daxpy routine where we assume incx and incy are both 1.
    Called by dgemv.
  */
  double alpha;
  int i;
  int n;
  n     = *n_p;
  alpha = *alpha_p;
  for (i=0;i<n;i++) {
    y[i] = y[i] + alpha * x[i];
  }
}
